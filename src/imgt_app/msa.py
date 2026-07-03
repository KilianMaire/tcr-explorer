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

from .dossier_models import AlignedRecord, AlignRequest, DossierWarning, MSAResult, Provenance
from .germline_sets import resolve_sequences


def clustalo_available() -> bool:
    return shutil.which("clustalo") is not None


def _aligner(seq_type: str) -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    if seq_type == "aa":
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -10, -0.5
    else:
        a.match_score, a.mismatch_score = 2, -1
        a.open_gap_score, a.extend_gap_score = -5, -0.5
    return a


def _merge_into_profile(profile: list[list[str]], new_profile_gapped: str, new_seq_gapped: str) -> list[str]:
    """Merge one pairwise alignment (profile-row-vs-new-seq) into a growing
    profile.

    `profile` is a list of rows, each row a list of single characters (one
    per existing profile column), all rows the same length. `new_profile_gapped`
    is the profile's first row as it came out of this pairwise alignment
    (with possibly-new gaps inserted); `new_seq_gapped` is the new sequence
    aligned against it, same length as `new_profile_gapped`.

    Walks the old profile's first row in lockstep with `new_profile_gapped`:
    wherever the pairwise alignment introduced a gap that isn't in the old
    profile (an insertion relative to the profile), a gap column is spliced
    into every existing row at that position. Returns the new sequence's row,
    now expressed in the (possibly widened) profile's column coordinates.
    """
    if not profile:
        # First sequence becomes the profile verbatim.
        return list(new_seq_gapped)

    old_ref = profile[0]
    new_row: list[str] = []
    # widened rows we're building, one per existing profile row, plus the new one
    widened = [[] for _ in profile]
    old_i = 0
    for np_char, ns_char in zip(new_profile_gapped, new_seq_gapped):
        if np_char == "-":
            # A gap introduced by this pairwise step. Is it at an old-profile
            # gap position (already there) or a genuinely new insertion?
            if old_i < len(old_ref) and old_ref[old_i] == "-":
                # matches an existing profile gap column; consume it
                for r, row in enumerate(profile):
                    widened[r].append(row[old_i])
                old_i += 1
            else:
                # brand-new gap column: insert '-' into every existing row,
                # do NOT consume an old_ref position
                for r in range(len(profile)):
                    widened[r].append("-")
            new_row.append(ns_char)
        else:
            # consumes exactly one old_ref residue/gap position
            for r, row in enumerate(profile):
                widened[r].append(row[old_i])
            old_i += 1
            new_row.append(ns_char)
    # any trailing old-profile columns not consumed (shouldn't happen with
    # global alignment, but guard anyway)
    while old_i < len(old_ref):
        for r, row in enumerate(profile):
            widened[r].append(row[old_i])
        new_row.append("-")
        old_i += 1

    for r in range(len(profile)):
        profile[r] = widened[r]
    profile.append(new_row)
    return new_row


def center_star_align(seqs, seq_type):
    """Progressive center-star MSA. seqs: list[(name, seq)]. Returns
    list[(name, gapped)], all gapped strings the same length and in a shared
    column coordinate system."""
    if len(seqs) < 2:
        return [(n, s) for n, s in seqs]

    aligner = _aligner(seq_type)
    n = len(seqs)

    # Pick the center = sequence with the max summed pairwise score against
    # all others (the classic center-star choice).
    scores = [0.0] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                scores[i] += aligner.score(seqs[i][1], seqs[j][1])
    c = max(range(n), key=lambda i: scores[i])

    order = [c] + [i for i in range(n) if i != c]

    profile: list[list[str]] = []
    names: list[str] = []
    for idx in order:
        name, s = seqs[idx]
        if not profile:
            profile.append(list(s))
            names.append(name)
            continue
        ref = "".join(profile[0])
        al = aligner.align(ref, s)[0]
        ref_gapped, seq_gapped = str(al[0]), str(al[1])
        _merge_into_profile(profile, ref_gapped, seq_gapped)
        names.append(name)

    rows = ["".join(row) for row in profile]
    return list(zip(names, rows))


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
    seqs, warnings = resolve_sequences(request)
    prov = []
    if request.chain or request.genes:
        prov.append(Provenance(block="alignment", source="cdr_enricher",
                                confidence="high", kind="germline_lookup"))
    if len(seqs) < 2:
        warnings.append(DossierWarning(code="too_few_sequences", block="alignment",
                                        message="at least two sequences are required to align"))
        return MSAResult(engine="none", seq_type=request.seq_type, n_sequences=len(seqs),
                          alignment_length=0, records=[], provenance=prov, warnings=warnings)
    engine = "center_star"
    aligned = None
    if clustalo_available():
        aligned = _run_clustalo(seqs)
        if aligned:
            engine = "clustalo"
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
    return MSAResult(engine=engine, seq_type=request.seq_type, n_sequences=len(aligned),
                      alignment_length=L, records=[AlignedRecord(name=n, aligned=s) for n, s in aligned],
                      consensus=_consensus(aligned), mean_pct_identity=_mean_identity(aligned),
                      provenance=prov, warnings=warnings)


def to_fasta(result: MSAResult) -> str:
    return "\n".join(f">{r.name}\n{r.aligned}" for r in result.records)
