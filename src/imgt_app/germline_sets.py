"""
Germline segment resolution for multiple-sequence alignment requests.

Reuses stitchr's germline FASTA (via cdr_enricher) to resolve V/J/C segment
sets or explicit gene names. stitchr's data source carries no D-segment
records (the last pipe field never marks ~DIVERSITY), so a germline D
request always resolves to an empty set plus a `segment_unavailable`
warning -- this is a deliberate honesty behaviour, not a bug.
"""
from __future__ import annotations

from .cdr_enricher import _stitchr_data_dir, _translate, _SPECIES_STITCHR
from .dossier_models import AlignRequest, DossierWarning

_SEG_TOKEN = {"V": "~VARIABLE", "J": "~JOINING", "C": "~CONSTANT"}  # D absent from stitchr


def load_segment_map(chain: str, species: str, segment: str) -> dict[str, str]:
    """Return {gene_name: nt} for a segment from the stitchr chain FASTA. Empty if
    the segment token is absent (e.g. D) or data missing. Prefers *01 allele."""
    token = _SEG_TOKEN.get(segment.upper())
    if token is None:
        return {}
    data_dir = _stitchr_data_dir()
    if data_dir is None:
        return {}
    fa = data_dir / species.upper() / f"{chain.upper()}.fasta"
    if not fa.exists():
        return {}

    best: dict[str, tuple[str, str]] = {}   # gene -> (best_allele, seq)
    header, parts_seq = "", []

    def _commit() -> None:
        if not header:
            return
        parts = header.split("|")
        if len(parts) < 2:
            return
        allele = parts[1].strip()
        seg_field = parts[-1].strip().upper()
        if token not in seg_field:
            return
        gene = allele.split("*")[0].upper()
        seq = "".join(parts_seq).upper().strip()
        if not seq:
            return
        existing = best.get(gene)
        if existing is None or "*01" in allele:
            best[gene] = (allele, seq)

    with fa.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                _commit()
                header, parts_seq = line[1:], []
            else:
                parts_seq.append(line)
    _commit()

    return {gene: seq for gene, (_, seq) in best.items()}


def _best_frame_translate(nt: str) -> str:
    """Translate in the reading frame that yields the fewest stop/ambiguous codons.

    Germline segments (notably J) are not in frame 0, so a naive frame-0 translation
    is biologically meaningless. Trying all three frames and keeping the cleanest one
    gives a sensible protein for mutual alignment without fabricating a frame.
    """
    best, best_penalty = _translate(nt), None
    for frame in (0, 1, 2):
        aa = _translate(nt[frame:])
        penalty = aa.count("*") + aa.count("?")
        if best_penalty is None or penalty < best_penalty:
            best, best_penalty = aa, penalty
    return best


def _maybe_translate(seqs: list[tuple[str, str]], translate: bool) -> list[tuple[str, str]]:
    if not translate:
        return seqs
    return [(n, _best_frame_translate(s)) for n, s in seqs]


def resolve_sequences(request: AlignRequest) -> tuple[list[tuple[str, str]], list[DossierWarning]]:
    """Resolve an AlignRequest into a list of (name, sequence) pairs plus warnings.

    Three mutually-exclusive input modes, checked in order:
      1. request.sequences  -- provided sequences, passed through verbatim.
      2. request.genes      -- explicit gene/allele names resolved against germline.
      3. request.chain + request.segment -- a full germline segment set (V/J/C only).
    """
    warns: list[DossierWarning] = []

    if request.sequences:
        return [(s.name, s.seq) for s in request.sequences], warns

    if request.genes:
        seqs: list[tuple[str, str]] = []
        for name in request.genes:
            base = name.split("*")[0].upper()
            chain = base[:3]
            seg = base[3:4]  # V/D/J/C letter
            sp = _SPECIES_STITCHR.get(request.species.lower(), "HUMAN")
            m = load_segment_map(chain, sp, seg)
            hit = m.get(base)
            if hit:
                seqs.append((base, hit))
            else:
                warns.append(DossierWarning(
                    code="ambiguous_gene", block="alignment",
                    message=f"gene {name} not resolvable from germline",
                ))
        return _maybe_translate(seqs, request.translate), warns

    if request.chain and request.segment:
        sp = _SPECIES_STITCHR.get(request.species.lower(), "HUMAN")
        m = load_segment_map(request.chain, sp, request.segment)
        if not m:
            warns.append(DossierWarning(
                code="segment_unavailable", block="alignment",
                message=(
                    f"segment {request.segment} not available in the germline source "
                    f"for {request.chain}/{request.species}"
                ),
            ))
            return [], warns
        return _maybe_translate(sorted(m.items()), request.translate), warns

    return [], warns
