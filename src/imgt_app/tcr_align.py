from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from Bio.Align import PairwiseAligner

from .germline_db import Allele, germline_alleles

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
