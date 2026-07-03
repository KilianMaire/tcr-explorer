"""Free-text query parsing for `/v1/records`.

`parse_query` reads a natural-language string such as "souris CASSGGTGEQYF"
or "mouse TCR with CDR3 CASSLGTEAFF" and pulls out the structured fields the
retrieval engine actually needs: species, a CDR3, a V/J gene, or a namespaced
database id. Prose words that match none of these are simply dropped.

Gene detection mirrors the segment-style check in `ask._find_gene`: a token
starting with TR[ABGD] then V or J. A CDR3 is an amino-acid token of length
8..22 starting with C that is not itself a gene token.
"""
from __future__ import annotations

import re

_ID_RE = re.compile(r"^(vdjdb|iedb|mcpas|tcr3d):")
_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
_TOKEN_SPLIT_RE = re.compile(r"[\s,;]+")

_SPECIES_WORDS = {
    "mouse": "mouse",
    "murine": "mouse",
    "souris": "mouse",
    "human": "human",
    "humain": "human",
    "homo": "human",  # first word of "homo sapiens"
}


def _is_gene_token(tok: str) -> bool:
    """Same segment-detection style as `ask._find_gene`: TR[ABGD] then V/J."""
    u = tok.upper()
    return len(u) >= 5 and u[:2] == "TR" and u[2] in "ABGD" and u[3:4] in ("V", "J")


def _gene_segment(tok: str) -> str:
    return tok.upper()[3]


def _is_cdr3_token(tok: str) -> bool:
    u = tok.upper()
    return 8 <= len(u) <= 22 and u.startswith("C") and bool(_AA_RE.match(u)) and not _is_gene_token(tok)


def _detect_species(text_lower: str) -> str | None:
    for tok in re.split(r"[^a-z]+", text_lower):
        if tok in _SPECIES_WORDS:
            return _SPECIES_WORDS[tok]
    return None


def parse_query(text: str) -> dict:
    """Parse free text into species/cdr3_aa/v_gene/j_gene/record_id.

    All keys are always present; unmatched fields are None.
    """
    result: dict = {
        "species": None,
        "cdr3_aa": None,
        "v_gene": None,
        "j_gene": None,
        "record_id": None,
    }
    if not text:
        return result

    result["species"] = _detect_species(text.lower())

    for tok in _TOKEN_SPLIT_RE.split(text.strip()):
        if not tok:
            continue
        if result["record_id"] is None and _ID_RE.match(tok.lower()):
            result["record_id"] = tok
            continue
        if _is_gene_token(tok):
            seg = _gene_segment(tok)
            if seg == "V" and result["v_gene"] is None:
                result["v_gene"] = tok
            elif seg == "J" and result["j_gene"] is None:
                result["j_gene"] = tok
            continue
        if result["cdr3_aa"] is None and _is_cdr3_token(tok):
            result["cdr3_aa"] = tok

    return result
