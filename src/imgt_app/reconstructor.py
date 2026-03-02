"""
TCR full-gene reconstruction from V gene + CDR3 + J gene.

Algorithm
---------
1. Load V-REGION nt from stitchr IMGT data (same source as cdr_enricher).
2. Load J-REGION nt from stitchr IMGT data (~JOINING segments).
3. Back-translate CDR3 amino-acid sequence using human-optimised codons.
4. Assemble junction:
   - V_nt[:-3]   : V-REGION up to but NOT including Cys104 codon
                   (Cys104 is the first residue of the VDJdb CDR3)
   - CDR3_nt     : back-translated full CDR3 (Cys104 … Phe/Trp118)
   - J_nt[fw+3:] : J-REGION AFTER the conserved Phe/Trp codon
                   (Phe/Trp118 is already the last residue of CDR3)

IMGT CDR3 boundaries used here follow VDJdb convention:
  CDR3 = Cys104 … Phe/Trp118  (boundaries inclusive)
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional

from .cdr_enricher import (
    _CODON,
    _SPECIES_STITCHR,
    _gene_to_chain,
    _stitchr_data_dir,
    _cached_v_map,
    _translate,
)

# ---------------------------------------------------------------------------
# Human-optimised codon table (most frequent codon per amino acid)
# ---------------------------------------------------------------------------
_HUMAN_CODON: dict[str, str] = {
    "A": "GCC", "R": "AGG", "N": "AAC", "D": "GAC",
    "C": "TGC", "E": "GAG", "Q": "CAG", "G": "GGC",
    "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC",
    "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


def back_translate(aa_seq: str) -> str:
    """Back-translate an amino-acid string to nucleotides using human codons."""
    return "".join(_HUMAN_CODON.get(aa, "NNN") for aa in aa_seq.upper())


# ---------------------------------------------------------------------------
# J-region loader  (mirrors _load_v_region_map from cdr_enricher)
# ---------------------------------------------------------------------------
def _load_j_region_map(chain: str, species: str) -> dict[str, str]:
    """
    Read stitchr FASTA and return {gene_name: nt_sequence} for J-REGION entries.
    Prefers *01 allele; falls back to any available allele.
    """
    data_dir = _stitchr_data_dir()
    if data_dir is None:
        return {}

    fa_path = data_dir / species.upper() / f"{chain.upper()}.fasta"
    if not fa_path.exists():
        return {}

    gene_map: dict[str, tuple[str, str]] = {}   # gene → (best_allele, seq)
    current_header = ""
    current_parts: list[str] = []

    def _commit() -> None:
        if not current_header:
            return
        seq = "".join(current_parts).upper()
        parts = current_header.split("|")
        if len(parts) < 2:
            return
        allele = parts[1]
        segment = parts[-1].strip().upper()
        # Accept both ~JOINING and J-REGION style markers
        if "JOINING" not in segment and "J-REGION" not in segment:
            return
        gene = allele.split("*")[0].upper()
        existing = gene_map.get(gene)
        if existing is None or "*01" in allele:
            gene_map[gene] = (allele, seq)

    with fa_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                _commit()
                current_header = line[1:]
                current_parts = []
            else:
                current_parts.append(line)
    _commit()

    return {gene: seq for gene, (_, seq) in gene_map.items()}


@lru_cache(maxsize=8)
def _cached_j_map(chain: str, species_str: str) -> dict[str, str]:
    return _load_j_region_map(chain, species_str)


# ---------------------------------------------------------------------------
# Junction helpers
# ---------------------------------------------------------------------------
def _find_fw_codon_offset(nt: str) -> Optional[int]:
    """
    Return the byte offset of the first Phe (F) or Trp (W) codon in *nt*.
    Returns None if no F/W codon is found.
    """
    for i in range(0, len(nt) - 2, 3):
        aa = _CODON.get(nt[i : i + 3].upper(), "?")
        if aa in ("F", "W"):
            return i
    return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def reconstruct_tcr(
    v_gene: str,
    j_gene: str,
    cdr3_aa: str,
    species: str = "human",
) -> dict:
    """
    Reconstruct a full-length TCR coding sequence from V gene, CDR3, and J gene.

    Parameters
    ----------
    v_gene   : IMGT V gene name, e.g. "TRBV19" or "TRBV19*01"
    j_gene   : IMGT J gene name, e.g. "TRBJ2-7" or "TRBJ2-7*01"
    cdr3_aa  : CDR3 amino-acid sequence (VDJdb convention: Cys…Phe/Trp)
    species  : "human" | "mouse" | "other"

    Returns
    -------
    dict with keys:
      v_gene, j_gene, cdr3_aa, species,
      full_nt, full_aa,           — assembled sequence (None if V or J not found)
      v_region_nt, cdr3_nt, j_region_nt,
      v_found, j_found,
      note
    """
    stitchr_species = _SPECIES_STITCHR.get(species.lower(), "HUMAN")
    chain = _gene_to_chain(v_gene)

    v_gene_base = v_gene.split("*")[0].strip().upper()
    j_gene_base = j_gene.split("*")[0].strip().upper()

    v_map = _cached_v_map(chain, stitchr_species)
    j_map = _cached_j_map(chain, stitchr_species)

    v_nt = v_map.get(v_gene_base, "")
    j_nt = j_map.get(j_gene_base, "")

    cdr3_nt = back_translate(cdr3_aa)

    full_nt: Optional[str] = None
    full_aa: Optional[str] = None

    if v_nt and j_nt:
        # V-REGION from stitchr ends at Cys104 (inclusive).
        # VDJdb CDR3 starts with that Cys104 → drop V's last codon to avoid duplication.
        v_prefix = v_nt[:-3] if len(v_nt) >= 3 else v_nt

        # VDJdb CDR3 ends with Phe/Trp118 = first coding codon of J segment.
        # Drop that first F/W codon from J to avoid duplication.
        fw_offset = _find_fw_codon_offset(j_nt)
        j_suffix = j_nt[fw_offset + 3 :] if fw_offset is not None else j_nt

        full_nt = v_prefix + cdr3_nt + j_suffix
        full_aa = _translate(full_nt).rstrip("*") or None

    return {
        "v_gene": v_gene_base,
        "j_gene": j_gene_base,
        "cdr3_aa": cdr3_aa,
        "species": species,
        "full_nt": full_nt,
        "full_aa": full_aa,
        "v_region_nt": v_nt or None,
        "cdr3_nt": cdr3_nt,
        "j_region_nt": j_nt or None,
        "v_found": bool(v_nt),
        "j_found": bool(j_nt),
        "note": (
            "CDR3 back-translated with human-optimised codons. "
            "V/J boundaries follow IMGT/VDJdb convention: CDR3 = Cys104…Phe/Trp118."
        ),
    }
