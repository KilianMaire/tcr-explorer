"""
K-mer germline aligner: raw-sequence annotation backend.

Annotates a raw TCR sequence (nt or protein) against the vendored stitchr
germline V/J FASTA using shared-k-mer scoring. This is the always-available
fallback used when IgBLAST is absent.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
from .cdr_enricher import _cached_v_map, _translate, _SPECIES_STITCHR
from .reconstructor import _cached_j_map

_K = 7
_CHAINS = {"TRA": "alpha", "TRB": "beta", "TRG": "gamma", "TRD": "delta"}
_HAS_D = {"TRB", "TRD"}


@dataclass
class KmerAnnotation:
    v_call: Optional[str] = None
    j_call: Optional[str] = None
    d_call: Optional[str] = None
    v_score: Optional[float] = None
    j_score: Optional[float] = None
    chain: str = "unknown"
    warnings: list[tuple[str, str]] = field(default_factory=list)


def _kmers(s: str, k: int = _K) -> set[str]:
    s = s.upper()
    return {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else set()


def _best(query_kmers: set[str], gene_map: dict[str, str]) -> tuple[Optional[str], float]:
    best_gene, best = None, 0.0
    for gene, nt in gene_map.items():
        gk = _kmers(nt)
        if not gk:
            continue
        shared = len(query_kmers & gk)
        score = 100.0 * shared / max(len(gk), 1)
        if score > best:
            best_gene, best = gene, score
    return best_gene, round(best, 2)


def annotate_sequence(seq: str, species: str, is_protein: bool) -> KmerAnnotation:
    sp = _SPECIES_STITCHR.get(species.lower(), "HUMAN")
    ann = KmerAnnotation()
    if is_protein:
        ann.warnings.append(("aa_annotation_limited",
            "protein input yields V-region annotation only; no D and no junction nt"))
    qk = _kmers(seq)
    # Try each chain's V map; the chain with the best V hit wins.
    best_overall = (None, 0.0, None)  # (v_gene, score, chain_key)
    for chain_key, chain_name in _CHAINS.items():
        vmap = _cached_v_map(chain_key, sp)
        g, sc = _best(qk, vmap)
        if sc > best_overall[1]:
            best_overall = (g, sc, chain_key)
    v_gene, v_score, chain_key = best_overall
    if chain_key is None:
        return ann
    ann.v_call, ann.v_score, ann.chain = v_gene, v_score, _CHAINS[chain_key]
    jmap = _cached_j_map(chain_key, sp)
    ann.j_call, ann.j_score = _best(qk, jmap)
    if chain_key in _HAS_D:
        ann.d_call = None
        ann.warnings.append(("d_segment_unresolved",
            "k-mer backend does not resolve the D segment; use IgBLAST for D"))
    return ann
