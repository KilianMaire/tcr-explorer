"""Multiple sequence alignment engine.

Two backends:
  - clustalo, if the binary is present on PATH (checked via `clustalo_available`,
    which is monkeypatchable so tests can force the fallback path).
  - center_star, a pure-biopython progressive center-star merge used whenever
    clustalo is absent or fails.

The center-star merge is a real progressive alignment, not a naive pad: each
sequence is aligned pairwise to a running profile (seeded by the center
sequence), and whenever a new alignment inserts a gap into a position of the
profile that already has aligned rows, that gap is propagated into every row
already in the profile. This keeps column coordinates consistent across all
rows, which a simple right-pad-to-max-length approach does not guarantee.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from collections import Counter

from Bio.Align import PairwiseAligner, substitution_matrices

from .cdr_enricher import _translate
from .dossier_models import AlignedRecord, AlignRequest, DossierWarning, MSAResult, Provenance
from .germline_sets import resolve_sequences


def clustalo_available() -> bool:
    return shutil.which("clustalo") is not None


_BLOSUM62 = substitution_matrices.load("BLOSUM62")
_BLOSUM62_ALPHABET = set(_BLOSUM62.alphabet)


def _sanitize_aa(seq: str) -> str:
    """Map any residue absent from the BLOSUM62 alphabet to a neutral in-alphabet
    residue so PairwiseAligner never raises on an unknown letter. The mouse
    frame-0 translations here happen to stay in-alphabet (the alphabet includes
    the 20 aa plus B, Z, X and the stop '*'), but genuinely foreign characters
    (e.g. U, O, J, lowercase, whitespace) would otherwise crash the engine."""
    return "".join(c if c in _BLOSUM62_ALPHABET else ("X" if "X" in _BLOSUM62_ALPHABET else "A") for c in seq)


def _aligner(seq_type: str) -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    if seq_type == "aa":
        a.substitution_matrix = _BLOSUM62
        a.open_gap_score, a.extend_gap_score = -10, -0.5
    else:
        a.match_score, a.mismatch_score = 2, -1
        a.open_gap_score, a.extend_gap_score = -5, -0.5
    return a


def _parse_pairwise(center_gapped: str, seq_gapped: str):
    """Decompose one pairwise alignment of the center (ungapped) against a
    sequence into center-coordinate form.

    Returns:
      - col_at: list of length = number of center residues; col_at[p] is the
        character of the sequence aligned to center residue p (a residue or '-').
      - ins: list of length = number of center residues + 1; ins[slot] is the
        run of sequence residues inserted relative to the center at that slot
        (slot 0 = before the first center residue, slot p = after center
        residue p-1).
    """
    col_at: list[str] = []
    ins: list[list[str]] = [[]]
    for c, s in zip(center_gapped, seq_gapped):
        if c == "-":
            # sequence residue inserted where the center has no residue
            ins[-1].append(s)
        else:
            col_at.append(s)
            ins.append([])
    return col_at, ["".join(x) for x in ins]


def center_star_align(seqs, seq_type):
    """Progressive center-star MSA. seqs: list[(name, seq)]. Returns
    list[(name, gapped)], all gapped strings the same length and in a shared
    column coordinate system.

    Every sequence is aligned pairwise to the (ungapped) center, then all the
    pairwise alignments are merged onto a shared center coordinate: each center
    residue owns a fixed column, and insertions relative to the center are
    left-justified into a per-slot block sized to the widest insertion at that
    slot. Only ungapped input sequences are ever handed to PairwiseAligner (for
    aa, additionally sanitized to the BLOSUM62 alphabet), so the engine never
    trips over a gap character or a foreign residue. Deterministic.
    """
    if len(seqs) < 2:
        return [(n, s) for n, s in seqs]

    aligner = _aligner(seq_type)
    is_aa = seq_type == "aa"
    n = len(seqs)

    # Alignment-safe copies (aa: sanitized to the matrix alphabet). Names and
    # residue counts are preserved, so column geometry maps back 1:1 to the
    # original characters, which is what we emit.
    raw = [s for _, s in seqs]
    safe = [_sanitize_aa(s) if is_aa else s for s in raw]

    # Pick the center = sequence with the max summed pairwise score against
    # all others (the classic center-star choice), scored on the safe copies.
    scores = [0.0] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                scores[i] += aligner.score(safe[i], safe[j])
    c = max(range(n), key=lambda i: scores[i])
    center_safe, center_raw = safe[c], raw[c]
    center_len = len(center_raw)

    # Align every sequence (including the center against itself) to the center
    # and decompose into center-coordinate form on the ORIGINAL characters.
    parsed = []  # (name, col_at, ins)
    max_ins = [0] * (center_len + 1)
    for (name, _), s_safe, s_raw in zip(seqs, safe, raw):
        al = aligner.align(center_safe, s_safe)[0]
        idx = al.indices  # 2 x L: -1 marks a gap in that row
        center_gapped = "".join(center_raw[i] if i >= 0 else "-" for i in idx[0])
        seq_gapped = "".join(s_raw[j] if j >= 0 else "-" for j in idx[1])
        col_at, ins = _parse_pairwise(center_gapped, seq_gapped)
        parsed.append((name, col_at, ins))
        for slot in range(center_len + 1):
            if len(ins[slot]) > max_ins[slot]:
                max_ins[slot] = len(ins[slot])

    out = []
    for name, col_at, ins in parsed:
        row: list[str] = []
        for slot in range(center_len + 1):
            blk = ins[slot]
            row.append(blk + "-" * (max_ins[slot] - len(blk)))
            if slot < center_len:
                row.append(col_at[slot])
        out.append((name, "".join(row)))
    return out


def _consensus(aligned):
    if not aligned:
        return ""
    L = len(aligned[0][1])
    cols = []
    for i in range(L):
        col = [a[1][i] for a in aligned if i < len(a[1])]
        counts = Counter(c for c in col if c != "-")
        if not counts:
            cols.append("-")
        else:
            top, n = counts.most_common(1)[0]
            cols.append(top if n > len(col) / 2 else "X")
    return "".join(cols)


def _mean_identity(aligned):
    if len(aligned) < 2:
        return 100.0 if aligned else 0.0
    seqs = [a[1] for a in aligned]
    L = len(seqs[0])
    pairs = 0
    total = 0.0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            same = cols = 0
            for k in range(L):
                a, b = seqs[i][k], seqs[j][k]
                if a == "-" and b == "-":
                    continue
                cols += 1
                if a == b:
                    same += 1
            total += (100.0 * same / cols) if cols else 0.0
            pairs += 1
    return round(total / pairs, 2) if pairs else 0.0


def _pick_frame(nt_list):
    """Reading frame (0/1/2) minimizing total stop/ambiguous codons across the set."""
    best_f, best_pen = 0, None
    for f in (0, 1, 2):
        pen = sum(_translate(nt[f:]).count("*") + _translate(nt[f:]).count("?") for nt in nt_list)
        if best_pen is None or pen < best_pen:
            best_f, best_pen = f, pen
    return best_f


def _translate_in_frame(nt, frame):
    """Return (aa, coding_nt) where coding_nt is the in-frame slice, len == 3*len(aa)."""
    usable = (len(nt) - frame) // 3 * 3
    coding = nt[frame:frame + usable]
    return _translate(coding), coding


def _expand_to_nt(gapped_aa, coding_nt):
    """Thread the observed codons under a gapped aa row; gap columns become '---'."""
    out, pos = [], 0
    for ch in gapped_aa:
        if ch == "-":
            out.append("---")
        else:
            out.append(coding_nt[pos:pos + 3])
            pos += 3
    return "".join(out)


def _conservation(aligned):
    """Per-column fraction that is the most common residue (gaps count against)."""
    if not aligned:
        return []
    L = len(aligned[0][1])
    cons = []
    for i in range(L):
        col = [a[1][i] for a in aligned]
        non_gap = [c for c in col if c != "-"]
        cons.append(round(Counter(non_gap).most_common(1)[0][1] / len(col), 3) if non_gap else 0.0)
    return cons


def _align_codon_aware(request, prov, warnings):
    raw = request.model_copy(update={"translate": False})
    seqs, more = resolve_sequences(raw)
    warnings = warnings + more
    if len(seqs) < 2:
        warnings.append(DossierWarning(code="too_few_sequences", block="alignment",
                                       message="at least two sequences are required to align"))
        return MSAResult(engine="none", seq_type=request.seq_type, n_sequences=len(seqs),
                         alignment_length=0, records=[], view="aa_nt", provenance=prov, warnings=warnings)
    # Per-sequence reading frame: each row threads its own codons under its own aa,
    # and the aa alignment supplies the cross-row column correspondence. A single
    # set-wide frame would mistranslate germline segments that sit in different
    # phases (e.g. mouse TRBJ), breaking the aa alignment.
    triples = [(n, *_translate_in_frame(nt, _pick_frame([nt]))) for n, nt in seqs]
    try:
        aligned_aa = center_star_align([(n, aa) for n, aa, _ in triples], "aa")
    except Exception:
        warnings.append(DossierWarning(code="alignment_failed", block="alignment",
                                       message="alignment engine failed"))
        return MSAResult(engine="none", seq_type=request.seq_type, n_sequences=len(seqs),
                         alignment_length=0, records=[], view="aa_nt", provenance=prov, warnings=warnings)
    L = max(len(s) for _, s in aligned_aa)
    aligned_aa = [(n, s + "-" * (L - len(s))) for n, s in aligned_aa]
    coding_by = {n: c for n, _, c in triples}
    records = [AlignedRecord(name=n, aligned=gaa, aligned_aa=gaa,
                             aligned_nt=_expand_to_nt(gaa, coding_by[n])) for n, gaa in aligned_aa]
    return MSAResult(engine="center_star", seq_type=request.seq_type, n_sequences=len(records),
                     alignment_length=L, records=records, consensus=_consensus(aligned_aa),
                     mean_pct_identity=_mean_identity(aligned_aa), conservation=_conservation(aligned_aa),
                     view="aa_nt", provenance=prov, warnings=warnings)


def _run_clustalo(seqs):
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "in.fa")
        with open(inp, "w") as fh:
            for name, s in seqs:
                fh.write(f">{name}\n{s}\n")
        out = subprocess.run(
            ["clustalo", "-i", inp, "--outfmt=fasta"],
            capture_output=True, text=True, timeout=60,
        )
        if out.returncode != 0:
            return None
        recs, name, buf = [], None, []
        for line in out.stdout.splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        return recs


def align(request: AlignRequest) -> MSAResult:
    prov = []
    if request.chain or request.genes:
        prov.append(Provenance(block="alignment", source="cdr_enricher",
                                confidence="high", kind="germline_lookup"))
    # Codon-aware AA+NT view: only when translating nucleotide input.
    if request.translate and request.seq_type != "aa":
        return _align_codon_aware(request, prov, [])
    seqs, warnings = resolve_sequences(request)
    if len(seqs) < 2:
        warnings.append(DossierWarning(code="too_few_sequences", block="alignment",
                                        message="at least two sequences are required to align"))
        return MSAResult(engine="none", seq_type=request.seq_type, n_sequences=len(seqs),
                          alignment_length=0, records=[], provenance=prov, warnings=warnings)
    engine = "center_star"
    aligned = None
    if clustalo_available():
        try:
            aligned = _run_clustalo(seqs)
        except (subprocess.SubprocessError, OSError):
            aligned = None
        if aligned:
            engine = "clustalo"
        else:
            # clustalo was present but failed or timed out; degrade to the bundled
            # engine and record that the authoritative backend was attempted.
            warnings.append(DossierWarning(code="alignment_failed", block="alignment",
                                           message="clustalo failed; used the bundled center-star aligner"))
    if not aligned:
        try:
            aligned = center_star_align(seqs, request.seq_type)
        except Exception:
            warnings.append(DossierWarning(code="alignment_failed", block="alignment",
                                            message="alignment engine failed"))
            return MSAResult(engine="none", seq_type=request.seq_type, n_sequences=len(seqs),
                              alignment_length=0, records=[], provenance=prov, warnings=warnings)
        engine = "center_star"
    L = max(len(s) for _, s in aligned)
    aligned = [(n, s + "-" * (L - len(s))) for n, s in aligned]
    view = "aa" if request.seq_type == "aa" else "nt"
    return MSAResult(engine=engine, seq_type=request.seq_type, n_sequences=len(aligned),
                      alignment_length=L, records=[AlignedRecord(name=n, aligned=s) for n, s in aligned],
                      consensus=_consensus(aligned), mean_pct_identity=_mean_identity(aligned),
                      conservation=_conservation(aligned), view=view,
                      provenance=prov, warnings=warnings)


def to_fasta(result: MSAResult) -> str:
    return "\n".join(f">{r.name}\n{r.aligned}" for r in result.records)
