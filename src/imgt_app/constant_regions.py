"""Canonical membrane-bound TCR constant-region proteins.

The constant region is the invariant portion of a TCR chain that follows the
variable domain (after FR4). It is curated reference data per species and chain,
not derived per query, because the genomic germline does not cleanly carry the
membrane-bound protein isoform.

Provenance and validation status:
- Mouse TRBC (P01852) and TRAC (P01849) are byte-exact against UniProt and are
  validated byte-exact against the paired-chain ground-truth fixture (oracle).
- Human TRAC (P01848) is vendored verbatim from UniProt.
- Human TRBC1 (P01850) is vendored with the CH1 junction residue E restored at
  position 1: UniProt P01850 (SV=4) omits that leading E, but the mouse entry
  P01852 carries it and the mouse assembly is oracle-validated with it present,
  so restoring it keeps the human beta junction consistent (...GTRLTV|EDLNKV...).
  The human constants are NOT validated against a human ground-truth set; their
  provenance strings say so. Human beta defaults to TRBC1; TRBC2 (A0A5B9) differs
  at a few residues and cannot be selected from V, J, and CDR3 alone.
"""
from __future__ import annotations

from typing import Optional

# Mouse, UniProt P01852 (TCB1_MOUSE), oracle-validated.
_MOUSE_TRB = (
    "EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYC"
    "LSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYQQGVLSAT"
    "ILYEILLGKATLYAVLVSTLVVMAMVKRKNS"
)
# Mouse, UniProt P01849 (TRAC_MOUSE), oracle-validated.
_MOUSE_TRA = (
    "IQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKTVLDMKAMDSKSNGAIAWSNQTSFT"
    "CQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS"
)
# Human, UniProt P01848 (TRAC_HUMAN), verbatim.
_HUMAN_TRA = (
    "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFA"
    "CANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS"
)
# Human, UniProt P01850 (TRBC1_HUMAN) with the CH1 junction E restored (see docstring).
_HUMAN_TRB = "E" + (
    "DLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDS"
    "RYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVL"
    "SATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
)

# (chain, species) -> (sequence, provenance string)
_CONSTANTS: dict[tuple[str, str], tuple[str, str]] = {
    ("beta", "mouse"): (
        _MOUSE_TRB,
        "curated membrane-bound TRBC (mouse), oracle-validated; UniProt P01852",
    ),
    ("alpha", "mouse"): (
        _MOUSE_TRA,
        "curated membrane-bound TRAC (mouse), oracle-validated; UniProt P01849",
    ),
    ("beta", "human"): (
        _HUMAN_TRB,
        "curated membrane-bound TRBC1 (human), not oracle-validated; "
        "UniProt P01850 with CH1 junction E restored (present in mouse P01852)",
    ),
    ("alpha", "human"): (
        _HUMAN_TRA,
        "curated membrane-bound TRAC (human), not oracle-validated; UniProt P01848",
    ),
}


def constant_aa(chain: str, species: str) -> Optional[str]:
    """Return the membrane-bound constant-region protein for a chain and species,
    or None when it is not vendored. chain is alpha or beta; species is human or
    mouse (lowercased)."""
    entry = _CONSTANTS.get((chain, (species or "").lower()))
    return entry[0] if entry else None


def constant_source(chain: str, species: str) -> Optional[str]:
    """Return the provenance string for the vendored constant (accession and
    validation status), or None when the constant is not vendored."""
    entry = _CONSTANTS.get((chain, (species or "").lower()))
    return entry[1] if entry else None
