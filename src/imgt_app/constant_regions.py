"""Canonical membrane-bound TCR constant-region proteins.

The constant region is the invariant portion of a TCR chain that follows the
variable domain (after FR4). It is curated reference data per species and chain,
not derived per query, because the genomic germline does not cleanly carry the
membrane-bound protein isoform. The mouse sequences below are validated byte
exact against the paired-chain ground-truth fixture and are consistent with
IMGT TRBC and TRAC. Human constants are not yet vendored (full_chain is mouse
only for now).
"""
from __future__ import annotations

from typing import Optional

_MOUSE_TRB = (
    "EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYQQGVLSATILYEILLGKATLYAVLVSTLVVMAMVKRKNS"
)
_MOUSE_TRA = (
    "IQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKTVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS"
)

_CONSTANTS = {
    ("beta", "mouse"): _MOUSE_TRB,
    ("alpha", "mouse"): _MOUSE_TRA,
}


def constant_aa(chain: str, species: str) -> Optional[str]:
    """Return the membrane-bound constant-region protein for a chain and species,
    or None when it is not vendored. chain is alpha or beta; species is human or
    mouse (lowercased)."""
    return _CONSTANTS.get((chain, (species or "").lower()))
