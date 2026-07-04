"""Vendored TCR beta D-REGION alleles (human), which the stitchr germline set
does not carry. Source: IMGT/GENE-DB TRBD. D assignment is inherently low
confidence (the segment is 12 to 16 nt) and callers must flag it as such. Mouse
D is not vendored in this version, so mouse D calls are null with a reason.
"""
from __future__ import annotations

_HUMAN_TRBD = {
    "TRBD1*01": "GGGACAGGGGGC",
    "TRBD2*01": "GGGACTAGCGGGGGGG",
    "TRBD2*02": "GGGACTAGCGGGAGGG",
}


def d_alleles(species: str) -> dict[str, str]:
    """Return {allele_name: nt} for D-REGION alleles, or {} when not vendored."""
    return dict(_HUMAN_TRBD) if (species or "").lower() == "human" else {}
