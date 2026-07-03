from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

from Bio.Align import PairwiseAligner

from .germline_db import Allele, germline_alleles
from .reconstructor import (
    _j_frame_and_fw,
    _translate,
    _v_cys_cut,
    reconstruct_tcr,
)
from .records import infer_vj_from_cdr3

_NT = set("ACGTN")
_CHAINS = ("TRA", "TRB", "TRG", "TRD")


def detect_alphabet(seq: str) -> str:
    s = (seq or "").strip().upper()
    return "nucleotide" if s and set(s) <= _NT else "amino_acid"


@lru_cache(maxsize=1)
def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -1
    a.open_gap_score = -6
    a.extend_gap_score = -1
    return a


@dataclass
class AlleleCall:
    alleles: list[str]
    identity: float
    aligned_span: tuple[int, int]


def _score_identity_span(query: str, ref: str):
    """Return (score, identity, (qstart, qend)) for the best local alignment."""
    if not query or not ref:
        return None
    aln = _aligner().align(ref, query)[0]
    tb, qb = aln.aligned  # target (ref) blocks, query blocks
    matches = cols = 0
    for (ts, te), (qs, qe) in zip(tb, qb):
        for i in range(te - ts):
            cols += 1
            matches += ref[ts + i] == query[qs + i]
    identity = matches / cols if cols else 0.0
    qstart = int(qb[0][0])
    qend = int(qb[-1][1])
    return float(aln.score), identity, (qstart, qend)


def best_alleles(query: str, alleles: list[Allele], level: str):
    q = (query or "").strip().upper()
    scored = []
    for al in alleles:
        ref = al.nt if level == "nt" else al.aa
        r = _score_identity_span(q, ref)
        if r is not None:
            scored.append((al.name, r[0], r[1], r[2]))
    if not scored:
        return None
    top = max(s[1] for s in scored)
    winners = [s for s in scored if s[1] == top]
    identity = max(s[2] for s in winners)
    span = winners[0][3]
    return AlleleCall(alleles=[w[0] for w in winners], identity=identity, aligned_span=span)


def detect_chain(seq: str, species: str, level: str):
    best_chain = None
    best_score = float("-inf")
    q = (seq or "").strip().upper()
    for chain in _CHAINS:
        for segment in ("V", "J"):
            alleles = germline_alleles(species, chain, segment)
            for al in alleles:
                ref = al.nt if level == "nt" else al.aa
                r = _score_identity_span(q, ref)
                if r and r[0] > best_score:
                    best_score = r[0]
                    best_chain = chain
    return best_chain


# IMGT amino acid region ranges, 1-based inclusive. Converted below to 0-based
# half-open reference coordinates (and scaled by 3 for the nucleotide level).
_REGIONS_AA = {
    "FR1": (0, 26),
    "CDR1": (26, 38),
    "FR2": (38, 55),
    "CDR2": (55, 65),
    "FR3": (65, 104),
}
_V_FRAMEWORK_MIN = 30  # germline residues upstream of Cys104 for a determinable V


@dataclass
class Assignment:
    input_kind: str
    species: str
    chain: str | None
    v_call: dict | None = None
    j_call: dict | None = None
    d_call: dict | None = None
    constant_call: dict | None = None
    regions: dict = field(default_factory=dict)
    cdr3_aa: str | None = None
    v_determinable: bool = True
    v_reason: str | None = None
    v_db_inference: list | None = None
    reconstruction: dict | None = None
    warnings: list = field(default_factory=list)


def _to_call_dict(call: AlleleCall | None) -> dict | None:
    if call is None:
        return None
    return {
        "alleles": list(call.alleles),
        "identity": float(call.identity),
        "aligned_span": [int(call.aligned_span[0]), int(call.aligned_span[1])],
    }


def _winner(call: AlleleCall | None, alleles: list[Allele]) -> Allele | None:
    if call is None or not call.alleles:
        return None
    name = call.alleles[0]
    return next((a for a in alleles if a.name == name), None)


def _map_target_to_query(aln, target_pos: int):
    """Map a target (reference) position to the aligned query position, or None
    when the target position falls outside every aligned block."""
    tblocks, qblocks = aln.aligned
    for (ts, te), (qs, qe) in zip(tblocks, qblocks):
        if ts <= target_pos < te:
            return qs + (target_pos - ts)
    return None


def _region_identity(aln, ref: str, query: str, start: int, end: int):
    """Identity over a reference region span, matches divided by region length.
    Query positions that do not align to the region count as mismatches.

    Returns None when no query position covers the region at all (a partial
    input missing that region), so the caller omits it rather than reporting a
    misleading 0.0 that would read as fully divergent."""
    end = min(end, len(ref))
    if start >= end:
        return None
    matches = 0
    covered = 0
    for t in range(start, end):
        q = _map_target_to_query(aln, t)
        if q is not None and q < len(query):
            covered += 1
            if ref[t] == query[q]:
                matches += 1
    if covered == 0:
        return None
    return matches / (end - start)


def assign(seq: str, species: str | None = None, chain: str | None = None, want_d: bool = False) -> Assignment:
    """Turn a raw TCR sequence into a full allele assignment.

    Honesty rules enforced here: a bare CDR3 (no framework upstream of Cys104)
    never yields an allele-level V call, it refuses and attaches a database
    frequency inference instead; every call carries identity and span; a D call
    is always flagged low confidence; a segment absent from the germline is a
    null call, not a guess.
    """
    seq = (seq or "").strip().upper()
    kind = detect_alphabet(seq)
    level = "nt" if kind == "nucleotide" else "aa"
    species = species or "human"
    chain = chain or detect_chain(seq, species, level)

    if chain is None:
        return Assignment(
            input_kind=kind,
            species=species,
            chain=None,
            warnings=["no germline alignment"],
        )

    warnings: list[str] = []
    v_alleles = germline_alleles(species, chain, "V")
    j_alleles = germline_alleles(species, chain, "J")
    v_best = best_alleles(seq, v_alleles, level)
    j_best = best_alleles(seq, j_alleles, level)
    win_v = _winner(v_best, v_alleles)
    win_j = _winner(j_best, j_alleles)

    # Query reading frame for translation. For an amino acid input the query is
    # already residues; for a nucleotide input the frame comes from where the V
    # alignment starts on the query.
    if level == "aa":
        query_aa = seq
    else:
        frame = 0
        if win_v:
            aln_frame = _aligner().align(win_v.nt, seq)[0]
            target_blocks, query_blocks = aln_frame.aligned
            tstart = int(target_blocks[0][0])
            qstart = int(query_blocks[0][0])
            # The query reading frame is where the V allele's own frame 0 falls
            # on the query. A mid-V fragment aligns at target start tstart > 0,
            # so the frame is (qstart - tstart) % 3, not qstart % 3 (the two
            # agree only when tstart is a multiple of 3).
            frame = (qstart - tstart) % 3
        query_aa = _translate(seq[frame:])

    # V refusal rule and per region identity share one alignment of the query
    # to the winning V allele at the input's own level.
    v_determinable = True
    v_reason: str | None = None
    regions: dict[str, float] = {}
    if win_v:
        ref_v = win_v.nt if level == "nt" else win_v.aa
        aln_v = _aligner().align(ref_v, seq)[0]
        cys_ref = _v_cys_cut(win_v.nt) if level == "nt" else _v_cys_cut(win_v.nt) // 3
        covered_upstream = 0
        for (ts, te), (_qs, _qe) in zip(*aln_v.aligned):
            hi = min(te, cys_ref)
            if hi > ts:
                covered_upstream += hi - ts
        upstream_residues = covered_upstream / 3 if level == "nt" else covered_upstream
        if upstream_residues < _V_FRAMEWORK_MIN:
            v_determinable = False
            v_reason = "no framework in input; V allele not determinable from a CDR3 alone"

        if v_determinable:
            scale = 3 if level == "nt" else 1
            for name, (a, b) in _REGIONS_AA.items():
                ident = _region_identity(aln_v, ref_v, seq, a * scale, b * scale)
                if ident is not None:
                    regions[name] = ident

    # CDR3 extraction. Both anchors are mapped through amino acid alignments so
    # a nucleotide and an amino acid input share one code path.
    cdr3_aa: str | None = None
    if win_v and win_j and query_aa:
        cys_ref_aa = _v_cys_cut(win_v.nt) // 3
        aln_v_aa = _aligner().align(win_v.aa, query_aa)[0]
        cys_q = _map_target_to_query(aln_v_aa, cys_ref_aa)
        frame_j, fw_nt = _j_frame_and_fw(win_j.nt)
        fw_q = None
        if fw_nt >= 0:
            fw_ref_aa = (fw_nt - frame_j) // 3
            aln_j_aa = _aligner().align(win_j.aa, query_aa)[0]
            fw_q = _map_target_to_query(aln_j_aa, fw_ref_aa)
        if cys_q is not None and fw_q is not None and fw_q >= cys_q:
            cdr3_aa = query_aa[cys_q : fw_q + 1]
        else:
            warnings.append("CDR3 anchor could not be mapped")

    # V call: refused inputs carry no allele call, only a database inference.
    v_db_inference: list | None = None
    if v_determinable:
        v_call = _to_call_dict(v_best)
    else:
        v_call = None
        v_db_inference = infer_vj_from_cdr3(cdr3_aa or "", species)

    j_call = _to_call_dict(j_best)

    # Constant identification: align the 3' remainder past the J span.
    constant_call = None
    if j_best:
        j_qend = int(j_best.aligned_span[1])
        remainder = seq[j_qend:]
        # Only attempt a constant call when the 3' remainder is long enough to
        # plausibly be a constant region. A short remainder (e.g. a V region
        # only input, which carries no constant) would otherwise produce a
        # spurious confident call, so we refuse it. Scale by the input level.
        min_remainder = 60 if level == "nt" else 20
        if len(remainder) >= min_remainder:
            c_best = best_alleles(remainder, germline_alleles(species, chain, "C"), level)
            constant_call = _to_call_dict(c_best)

    # D call: only for TRB nucleotide input, always low confidence. When the
    # caller asked for D and none is available for this chain and species, say
    # why rather than returning a silent null.
    d_call = None
    if want_d:
        d_alleles = germline_alleles(species, chain, "D")
        if not d_alleles:
            if chain == "TRB":
                warnings.append(f"D germline not vendored for {species}")
            else:
                warnings.append(f"no D segment for {chain}")
        elif chain == "TRB" and level == "nt" and v_best and j_best:
            seg = seq[int(v_best.aligned_span[1]) : int(j_best.aligned_span[0])]
            d_best = best_alleles(seg, d_alleles, "nt") if seg else None
            if d_best:
                d_call = _to_call_dict(d_best)
                d_call["low_confidence"] = True
                warnings.append("D call is low confidence")

    # Reconstruction from the assigned parts.
    reconstruction = None
    if v_determinable and cdr3_aa and v_best and j_best:
        reconstruction = reconstruct_tcr(v_best.alleles[0], j_best.alleles[0], cdr3_aa, species)

    if v_determinable and v_best and v_best.identity < 0.90:
        warnings.append("low V identity")
    if j_best and j_best.identity < 0.90:
        warnings.append("low J identity")
    if constant_call and constant_call["identity"] < 0.90:
        warnings.append("low constant identity")
    if v_determinable and v_best and len(v_best.alleles) > 1:
        warnings.append(f"ambiguous V allele ({len(v_best.alleles)} co-optimal)")

    return Assignment(
        input_kind=kind,
        species=species,
        chain=chain,
        v_call=v_call,
        j_call=j_call,
        d_call=d_call,
        constant_call=constant_call,
        regions=regions,
        cdr3_aa=cdr3_aa,
        v_determinable=v_determinable,
        v_reason=v_reason,
        v_db_inference=v_db_inference,
        reconstruction=reconstruction,
        warnings=warnings,
    )
