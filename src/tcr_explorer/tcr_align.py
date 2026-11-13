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
_CHAIN_SYNONYMS = {
    "alpha": "TRA", "beta": "TRB", "gamma": "TRG", "delta": "TRD",
    "tra": "TRA", "trb": "TRB", "trg": "TRG", "trd": "TRD",
}


def _normalize_chain(chain: str) -> str | None:
    """Map a caller-supplied chain to a canonical TRA/TRB/TRG/TRD, or None if it
    is not recognized. Accepts the IMGT codes and the alpha/beta/gamma/delta
    names, case-insensitively. Prevents a stray value like "beta" from silently
    yielding an all-null assignment (germline_alleles would return nothing)."""
    key = (chain or "").strip().lower()
    if key in _CHAIN_SYNONYMS:
        return _CHAIN_SYNONYMS[key]
    up = (chain or "").strip().upper()
    return up if up in _CHAINS else None


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
_CDR3_MAX_AA = 30  # longer than any plausible Cys104 to Phe/Trp118 junction


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


def _cdr3_start(win_v: Allele, v_alleles: list[Allele]) -> int | None:
    """Nucleotide offset in the V allele of the codon that opens the CDR3 (IMGT
    position 104), or None when it cannot be placed.

    A few alleles have lost the conserved Cys (TRBV7-3*03 carries an Arg), and
    _v_cys_cut then falls back to the last codon of the V-REGION, several
    residues too far. Position 104 sits at the same offset in every allele of a
    gene, so such an allele takes it from one that kept the Cys."""
    gene = win_v.name.split("*")[0]
    for al in [win_v] + [a for a in v_alleles if a.name.split("*")[0] == gene]:
        cut = _v_cys_cut(al.nt)
        if _translate(al.nt[cut : cut + 3]) == "C":
            return cut
    return None


def _confident_j(call: AlleleCall | None, min_len: int) -> bool:
    if call is None:
        return False
    start, end = call.aligned_span
    return call.identity >= 0.90 and end - start >= min_len


def _joined_j(seq: str, j_alleles: list[Allele], level: str, v_end: int) -> AlleleCall | None:
    """J call for the J segment joined to the V, or None when no confident J
    lies downstream of the V alignment end.

    A transcript can carry a second J besides the recombined one (read-through
    into a downstream J of the locus, or a J upstream of the V). That J is
    untrimmed germline, so it outscores the recombined, 5' trimmed J and the
    best score over the whole input names the wrong gene. The search is
    therefore limited to the sequence after the V, and a confident J hit lying
    between the V and the current best hit takes over until none is left."""
    min_len = 25 if level == "nt" else 8
    call = best_alleles(seq[v_end:], j_alleles, level)
    # The first hit after the V may be short: an input can stop inside its J.
    if not _confident_j(call, 15 if level == "nt" else 5):
        return None
    while call.aligned_span[0] >= min_len:
        nearer = best_alleles(seq[v_end : v_end + call.aligned_span[0]], j_alleles, level)
        if not _confident_j(nearer, min_len):
            break
        call = nearer
    start, end = call.aligned_span
    return AlleleCall(call.alleles, call.identity, (start + v_end, end + v_end))


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

    chain_warning: str | None = None
    if chain:
        norm = _normalize_chain(chain)
        if norm is None:
            chain_warning = (
                f"unrecognized chain '{chain}'; expected one of "
                f"{', '.join(_CHAINS)} (or alpha/beta/gamma/delta). Auto-detecting instead."
            )
            chain = None
        else:
            chain = norm
    chain = chain or detect_chain(seq, species, level)

    if chain is None:
        warns = ["no germline alignment"]
        if chain_warning:
            warns.insert(0, chain_warning)
        return Assignment(
            input_kind=kind,
            species=species,
            chain=None,
            warnings=warns,
        )

    warnings: list[str] = []
    if chain_warning:
        warnings.append(chain_warning)
    v_alleles = germline_alleles(species, chain, "V")
    j_alleles = germline_alleles(species, chain, "J")
    v_best = best_alleles(seq, v_alleles, level)
    win_v = _winner(v_best, v_alleles)

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

    # J call: the J joined to the V when the input carries a determinable V,
    # otherwise (bare CDR3, J only input) the best J over the whole input.
    j_best = None
    if win_v and v_determinable:
        j_best = _joined_j(seq, j_alleles, level, int(v_best.aligned_span[1]))
    if j_best is None:
        j_best = best_alleles(seq, j_alleles, level)
    win_j = _winner(j_best, j_alleles)

    # CDR3 extraction. Both anchors are mapped at the input's own level, through
    # the V alignment and through an alignment of the J allele to its own span,
    # so a nucleotide input needs no single reading frame for the whole query
    # (an indel upstream in the V would shift the frame before the junction).
    cdr3_aa: str | None = None
    if win_v and win_j:
        step = 3 if level == "nt" else 1
        cys_nt = _cdr3_start(win_v, v_alleles)
        cys_q = None
        if cys_nt is not None:
            cys_q = _map_target_to_query(aln_v, cys_nt if level == "nt" else cys_nt // 3)
        frame_j, fw_nt = _j_frame_and_fw(win_j.nt)
        fw_q = None
        if fw_nt >= 0:
            j_start, j_end = j_best.aligned_span
            ref_j = win_j.nt if level == "nt" else win_j.aa
            aln_j = _aligner().align(ref_j, seq[j_start:j_end])[0]
            fw_q = _map_target_to_query(aln_j, fw_nt if level == "nt" else (fw_nt - frame_j) // 3)
            if fw_q is not None:
                fw_q += j_start
        if cys_q is not None and fw_q is not None and (fw_q - cys_q) % step:
            # A one base indel inside the J, at or ahead of its anchor (typically
            # a homopolymer miscall), leaves the junction readable in the V
            # frame: take the in-frame neighbour codon when it still encodes the
            # J anchor residue, and say so. A shift that arises in the junction
            # itself leaves the J alignment ungapped and is not rescued.
            j_indel = any(ts <= fw_nt for ts, _te in aln_j.aligned[0][1:])
            near = fw_q + (1 if (fw_q + 1 - cys_q) % step == 0 else -1)
            anchor = _translate(win_j.nt[fw_nt : fw_nt + 3])
            if j_indel and _translate(seq[near : near + 3]) == anchor:
                fw_q = near
                warnings.append("frameshift in the J at the CDR3 anchor")
        if cys_q is None or fw_q is None or fw_q < cys_q:
            warnings.append("CDR3 anchor could not be mapped")
        elif (fw_q - cys_q) % step:
            warnings.append("V and J are out of frame; no CDR3")
        else:
            cdr3_aa = seq[cys_q : fw_q + step]
            if level == "nt":
                cdr3_aa = _translate(cdr3_aa)
            if len(cdr3_aa) > _CDR3_MAX_AA:
                warnings.append(f"implausible CDR3 length ({len(cdr3_aa)} residues)")

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
