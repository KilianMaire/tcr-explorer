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


def _packaged_germline_dir() -> Path:
    """The IMGT germline bundled in the package (CC BY 4.0). Always present."""
    return Path(__file__).resolve().parent / "data" / "germline"


def _native_stitchr_data_dir() -> Optional[Path]:
    """Where stitchr/stitchrdl itself writes its Data (site-packages `Data` or a
    pip --user location). Used only by the opt-in germline refresh to find the
    freshly downloaded copy; query-time reads go through _stitchr_data_dir."""
    for sp in sys.path:
        candidate = Path(sp) / "Data"
        if candidate.is_dir() and _looks_like_stitchr_data(candidate):
            return candidate
    user_lib = Path.home() / "Library" / "Python"
    for p in user_lib.glob("*/lib/python/site-packages/Data"):
        if p.is_dir() and _looks_like_stitchr_data(p):
            return p
    return None


def _stitchr_data_dir() -> Optional[Path]:
    """Resolve the germline Data dir the enricher/reconstructor read at query time.

    Order: a user-refreshed copy (opt-in `tcr-explorer-refresh --germline`), then
    the germline bundled in the package (always present, so the tool works offline
    and without IMGT), then a native stitchr install (back-compat). Returns None
    only if none look like stitchr germline data."""
    from . import data_paths

    for c in (data_paths.germline_dir(), _packaged_germline_dir()):
        if c.is_dir() and _looks_like_stitchr_data(c):
            return c
    return _native_stitchr_data_dir()


def _iter_v_region_entries(chain: str, species: str) -> list[tuple[str, str]]:
    """Read the stitchr FASTA for `chain` (TRA, TRB, …) and `species`
    (HUMAN, MOUSE …) and return every V-REGION entry as (allele_upper, nt_seq),
    in file order. Header format is "{acc}|{allele}|{species}|…|~V-REGION".
    Empty when the germline data or file is absent."""
    data_dir = _stitchr_data_dir()
    if data_dir is None:
        return []
    fa_path = data_dir / species.upper() / f"{chain.upper()}.fasta"
    if not fa_path.exists():
        return []

    entries: list[tuple[str, str]] = []
    current_header = ""
    current_parts: list[str] = []

    def _commit() -> None:
        if not current_header:
            return
        seq = "".join(current_parts).upper()
        parts = current_header.split("|")
        if len(parts) < 2:
            return
        segment = parts[-1].strip()  # e.g. "~V-REGION"
        # stitchr marks V-REGION entries as ~VARIABLE in the last pipe field;
        # skip LEADER (~LEADER), JOINING (~JOINING), CONSTANT (~CONSTANT).
        if "VARIABLE" not in segment:
            return
        entries.append((parts[1].strip().upper(), seq))

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
    return entries


def _load_v_region_map(chain: str, species: str) -> dict[str, str]:
    """
    Return { gene_name_upper: nt_sequence } for the primary (*01) allele.
    Prefer *01 but fall back to any allele when *01 is absent. This gene-base map
    backs alignment (kmer_aligner), which scores against one representative
    germline per gene; allele-specific reads go through `_load_v_allele_map`.
    """
    gene_map: dict[str, str] = {}
    for allele, seq in _iter_v_region_entries(chain, species):
        gene = allele.split("*")[0]
        if gene not in gene_map or "*01" in allele:
            gene_map[gene] = seq
    return gene_map


def _load_v_allele_map(chain: str, species: str) -> dict[str, str]:
    """Return { full_allele_upper: nt_sequence }, keeping every allele so a caller
    can request a specific one (e.g. TRAV6-5*04). First occurrence wins on the
    rare duplicate allele id."""
    allele_map: dict[str, str] = {}
    for allele, seq in _iter_v_region_entries(chain, species):
        allele_map.setdefault(allele, seq)
    return allele_map


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


@lru_cache(maxsize=8)
def _cached_v_allele_map(chain: str, species_str: str) -> dict[str, str]:
    return _load_v_allele_map(chain, species_str)


def resolve_v_germline(v_gene: str, species: str = "human") -> tuple[Optional[str], str]:
    """Resolve a V gene name to (allele_used, v_region_nt), honoring an explicit
    allele in the name (e.g. TRAV6-5*04). A bare gene name, or an allele absent
    from germline, defaults to *01 (falling back to the lowest-numbered allele
    when *01 itself is absent). Returns (None, "") when the gene is not found.

    This mirrors reconstructor._resolve_germline so the dossier's reported allele,
    germline sequence, and CDR loops all reflect the same requested allele instead
    of silently collapsing to *01."""
    gene_base = v_gene.split("*")[0].strip().upper()
    parts = v_gene.split("*")
    requested = parts[1].strip().upper() if len(parts) > 1 and parts[1].strip() else None
    chain = _gene_to_chain(gene_base)
    stitchr_species = _SPECIES_STITCHR.get(species.lower(), "HUMAN")
    allele_map = _cached_v_allele_map(chain, stitchr_species)

    if requested and f"{gene_base}*{requested}" in allele_map:
        name = f"{gene_base}*{requested}"
        return name, allele_map[name]
    candidates = sorted(k for k in allele_map if k.split("*")[0] == gene_base)
    if not candidates:
        return None, ""
    name = f"{gene_base}*01" if f"{gene_base}*01" in allele_map else candidates[0]
    return name, allele_map[name]


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
    allele_used, nt_seq = resolve_v_germline(v_gene, species)
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
        "allele":   allele_used,
    }
