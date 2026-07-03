"""Germline allele catalog for the TCR sequence aligner.

Wraps the stitchr germline loader (reconstructor._load_allele_map) plus the
vendored human D-REGION table (d_regions) into one uniform, cached interface:
germline_alleles(species, chain, segment) -> list[Allele].

Honesty: a segment absent from the germline (e.g. TRA has no D, or a species
with no vendored data for that segment) yields an empty list, never a guess.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from imgt_app import d_regions
from imgt_app.reconstructor import _j_frame_and_fw, _load_allele_map, _translate

_MARKERS: dict[str, tuple[str, ...]] = {
    "V": ("VARIABLE", "V-REGION"),
    "J": ("JOINING", "J-REGION"),
    "C": ("CONSTANT", "C-REGION", "EX"),
}


@dataclass(frozen=True)
class Allele:
    name: str
    nt: str
    aa: str


def _translate_frame0_trimmed(nt: str) -> str:
    """Translate *nt* in frame 0 after trimming to a multiple of three (drop
    up to two trailing bases). Strips a single trailing stop codon if
    present, keeping any internal residues as-is."""
    trimmed = nt[: len(nt) - (len(nt) % 3)]
    aa = _translate(trimmed)
    if aa.endswith("*"):
        aa = aa[:-1]
    return aa


def _build_v_or_c_alleles(chain: str, species: str, markers: tuple[str, ...]) -> tuple[Allele, ...]:
    allele_map = _load_allele_map(chain, species, markers)
    return tuple(
        Allele(name=name, nt=nt, aa=_translate_frame0_trimmed(nt))
        for name, nt in sorted(allele_map.items())
    )


def _build_j_alleles(chain: str, species: str) -> tuple[Allele, ...]:
    allele_map = _load_allele_map(chain, species, _MARKERS["J"])
    out = []
    for name, nt in sorted(allele_map.items()):
        frame, _fw = _j_frame_and_fw(nt)
        out.append(Allele(name=name, nt=nt, aa=_translate(nt[frame:])))
    return tuple(out)


def _build_d_alleles(chain: str, species: str) -> tuple[Allele, ...]:
    # Only TRB has a D segment (TRA, TRG have none; TRD does but is not
    # vendored here). Gate on chain so an alpha-chain D query is honestly
    # empty rather than returning the vendored TRB D table by mistake.
    if chain.upper() != "TRB":
        return ()
    return tuple(
        Allele(name=name, nt=nt, aa="")
        for name, nt in sorted(d_regions.d_alleles(species).items())
    )


@lru_cache(maxsize=64)
def _cached_germline_alleles(species: str, chain: str, segment: str) -> tuple[Allele, ...]:
    if segment == "D":
        return _build_d_alleles(chain, species)
    if segment in ("V", "C"):
        return _build_v_or_c_alleles(chain, species, _MARKERS[segment])
    if segment == "J":
        return _build_j_alleles(chain, species)
    return ()


def germline_alleles(species: str, chain: str, segment: str) -> list[Allele]:
    """Return the germline allele catalog for (species, chain, segment).

    segment is one of "V", "J", "D", "C". Returns an empty list when the
    segment is not vendored for that species and chain (e.g. TRA has no D).
    Cached per (species, chain, segment); callers get a fresh list each call
    so the cached tuple cannot be mutated.
    """
    return list(_cached_germline_alleles(species, chain, segment))
