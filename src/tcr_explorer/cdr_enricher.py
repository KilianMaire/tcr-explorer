"""
CDR1/CDR2 enrichment for TCR records using stitchr's IMGT V-gene data.

CDR positions follow IMGT numbering (1-indexed, inclusive):
  CDR1-IMGT : positions 27–38  (up to 12 aa)
  CDR2-IMGT : positions 56–65  (up to 10 aa)

Usage
-----
    from tcr_explorer.cdr_enricher import get_cdr1_cdr2
    result = get_cdr1_cdr2("TRBV19", "HUMAN")
    # {"cdr1_aa": "LNHDAMYWYRQD", "cdr2_aa": "QKGDIAEGYS", "allele": "TRBV19*01"}
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# IMGT CDR positions (1-indexed, inclusive) in the V-REGION protein sequence
# ---------------------------------------------------------------------------
_CDR1_START = 27   # aa position (IMGT)
_CDR1_END   = 38
_CDR2_START = 56
_CDR2_END   = 65

# Minimal codon table (IUPAC standard codons only)
_CODON: dict[str, str] = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def _translate(nt: str) -> str:
    nt = nt.upper()
    return "".join(_CODON.get(nt[i:i+3], "?") for i in range(0, len(nt) - 2, 3))


def _looks_like_stitchr_data(candidate: Path) -> bool:
    """True if `candidate` actually holds stitchr germline FASTA (HUMAN/MOUSE/... subdirs),
    not just any unrelated directory that happens to be named "Data" (macOS' default
    case-insensitive filesystem can alias this app's own local "data/" cache dir)."""
    return any(
        next((candidate / species).glob("*.fasta"), None) is not None
        for species in ("HUMAN", "MOUSE")
    )


def _stitchr_data_dir() -> Optional[Path]:
    """Return stitchr Data directory, or None if stitchr is not installed."""
    # Stitchr installs to a user or system site-packages directory.
    # We probe sys.path entries for the Data folder it creates.
    for sp in sys.path:
        candidate = Path(sp) / "Data"
        if candidate.is_dir() and _looks_like_stitchr_data(candidate):
            return candidate
    # Also try the known pip --user install location
    user_lib = Path.home() / "Library" / "Python"
    for p in user_lib.glob("*/lib/python/site-packages/Data"):
        if p.is_dir() and _looks_like_stitchr_data(p):
            return p
    return None


def _load_v_region_map(chain: str, species: str) -> dict[str, str]:
    """
    Read stitchr FASTA for `chain` (TRA, TRB, …) and `species` (HUMAN, MOUSE …).
    Returns { gene_name_upper: nt_sequence } for the primary (*01) allele.
    Prefer *01 but fall back to any allele when *01 is absent.
    """
    data_dir = _stitchr_data_dir()
    if data_dir is None:
        return {}

    fa_path = data_dir / species.upper() / f"{chain.upper()}.fasta"
    if not fa_path.exists():
        return {}

    # Parse FASTA: header format is "{acc}|{allele}|{species}|…|~V-REGION" etc.
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
        allele = parts[1]          # e.g. "TRBV19*01"
        segment = parts[-1].strip()  # e.g. "~V-REGION"
        if "VARIABLE" not in segment:
            return                 # skip LEADER (~LEADER), JOINING (~JOINING), CONSTANT (~CONSTANT)
        # stitchr marks V-REGION entries as ~VARIABLE in the last pipe field
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


# ---------------------------------------------------------------------------
# Species helpers
# ---------------------------------------------------------------------------
_SPECIES_STITCHR: dict[str, str] = {
    "human": "HUMAN",
    "mouse": "MOUSE",
    "other": "HUMAN",   # fall back to human data
}

_CHAIN_FROM_GENE: dict[str, str] = {
    "TRAV": "TRA", "TRAJ": "TRA",
    "TRBV": "TRB", "TRBJ": "TRB",
    "TRGV": "TRG", "TRGJ": "TRG",
    "TRDV": "TRD", "TRDJ": "TRD",
}


def _gene_to_chain(gene: str) -> str:
    gene_up = gene.strip("*").split("*")[0].upper()
    for prefix, chain in _CHAIN_FROM_GENE.items():
        if gene_up.startswith(prefix):
            return chain
    return "TRB"   # default


# ---------------------------------------------------------------------------
# Cache: keyed by (chain, species_str)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=8)
def _cached_v_map(chain: str, species_str: str) -> dict[str, str]:
    return _load_v_region_map(chain, species_str)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def get_cdr1_cdr2(v_gene: str, species: str = "human") -> dict[str, Optional[str]]:
    """
    Return CDR1 and CDR2 amino acid sequences for a given TCR V gene.

    Parameters
    ----------
    v_gene  : IMGT V gene name, e.g. "TRBV19", "TRAV27", "TRBV19*01"
    species : "human" | "mouse" | "other"

    Returns
    -------
    dict with keys: cdr1_aa, cdr2_aa, allele, cdr1_nt, cdr2_nt
    All values are None when the gene is not found.
    """
    gene_base = v_gene.split("*")[0].strip().upper()
    chain = _gene_to_chain(gene_base)
    stitchr_species = _SPECIES_STITCHR.get(species.lower(), "HUMAN")

    v_map = _cached_v_map(chain, stitchr_species)
    nt_seq = v_map.get(gene_base)
    if not nt_seq:
        return {"cdr1_aa": None, "cdr2_aa": None, "allele": None,
                "cdr1_nt": None, "cdr2_nt": None}

    aa_seq = _translate(nt_seq)

    # Convert IMGT 1-indexed positions to 0-indexed Python slices
    cdr1_aa = aa_seq[_CDR1_START - 1 : _CDR1_END]
    cdr2_aa = aa_seq[_CDR2_START - 1 : _CDR2_END]

    # Nucleotide slices
    cdr1_nt = nt_seq[(_CDR1_START - 1) * 3 : _CDR1_END * 3]
    cdr2_nt = nt_seq[(_CDR2_START - 1) * 3 : _CDR2_END * 3]

    # Strip stop codons if they appear at the end
    cdr1_aa = cdr1_aa.rstrip("*")
    cdr2_aa = cdr2_aa.rstrip("*")

    return {
        "cdr1_aa":  cdr1_aa  or None,
        "cdr2_aa":  cdr2_aa  or None,
        "cdr1_nt":  cdr1_nt  or None,
        "cdr2_nt":  cdr2_nt  or None,
        "allele":   gene_base + "*01",
    }
