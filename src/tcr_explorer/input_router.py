"""Deterministic input router for TCR dossier queries.

Classifies a raw query string into one of: raw_nt | raw_aa | gene_name |
allele | id | unknown, honoring an explicit override and emitting
machine-readable warnings for ambiguous or unresolved input.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

_DNA = set("ACGT")
_AA = set("ACDEFGHIKLMNPQRSTVWY")
_AA_ONLY = set("EFILPQZX*")  # letters that never appear in a nucleotide alphabet
_AA_EXT = _AA | _AA_ONLY  # amino-acid alphabet widened with the signature letters
_ID_RE = re.compile(r"^(vdjdb|iedb):\w+$", re.I)
_GENE_RE = re.compile(r"^TR[ABGD][VDJC]\d+(-\d+)?(\*\d+)?$", re.I)


@dataclass
class RoutedQuery:
    detected_type: str
    normalized: str
    source: Optional[str] = None
    warnings: list[tuple[str, str]] = field(default_factory=list)


# tidytcells names species differently than the rest of the tool. Map ours to
# its convention so a mouse gene is standardized against musmusculus instead of
# the tidytcells default (homosapiens), which logs a spurious "Failed to
# standardize ... for species homosapiens" on every mouse gene.
_TT_SPECIES: dict[str, str] = {
    "human": "homosapiens",
    "mouse": "musmusculus",
    "other": "homosapiens",
}


def _normalize_gene(q: str, species: str = "human") -> str:
    try:
        import tidytcells as tt

        # on_fail="keep" returns the input unchanged on failure (no FutureWarning)
        # and is type-correct (Optional[str]), unlike the removed log_failures bool.
        tt_species = _TT_SPECIES.get(species.lower(), "homosapiens")
        fixed = tt.tr.standardize(symbol=q, species=tt_species, on_fail="keep")
        return fixed or q.upper()
    except Exception:
        return q.upper()


def route(query: str, input_type: str = "auto", species: str = "human") -> RoutedQuery:
    q = query.strip()

    if input_type != "auto":
        src = None
        m = _ID_RE.match(q)
        if input_type == "id" and m:
            src = m.group(1).lower()
        return RoutedQuery(
            input_type,
            q.upper() if input_type in ("raw_nt", "raw_aa") else q,
            src,
        )

    m = _ID_RE.match(q)
    if m:
        return RoutedQuery("id", q, source=m.group(1).lower())

    if _GENE_RE.match(q):
        norm = _normalize_gene(q, species)
        return RoutedQuery("allele" if "*" in q else "gene_name", norm)

    s = q.upper().replace(" ", "").replace("\n", "")
    if s and all(c in _AA_EXT for c in s):
        if set(s) & _AA_ONLY:
            return RoutedQuery("raw_aa", s)
        acgt_frac = sum(c in _DNA for c in s) / len(s)
        if acgt_frac >= 0.9:
            return RoutedQuery("raw_nt", s)
        return RoutedQuery(
            "raw_aa",
            s,
            warnings=[
                (
                    "ambiguous_alphabet",
                    "sequence letters are valid as both nucleotides and amino acids; defaulting to raw_aa",
                )
            ],
        )

    return RoutedQuery(
        "unknown",
        q,
        warnings=[
            (
                "unresolved_input_type",
                "could not classify the query; pass an explicit input_type",
            )
        ],
    )
