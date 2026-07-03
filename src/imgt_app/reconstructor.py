"""
TCR full-gene reconstruction from V gene + CDR3 + J gene.

Algorithm
---------
1. Load V-REGION nt from stitchr IMGT data (same source as cdr_enricher).
2. Load J-REGION nt from stitchr IMGT data (~JOINING segments).
3. Back-translate CDR3 amino-acid sequence using human-optimised codons.
4. Assemble junction, codon-aligned:
   - _v_cys_cut(V_nt)     : locates the nt offset of the conserved Cys104
                             codon in V_nt (V-REGION length is not always a
                             multiple of 3, so this is NOT simply len-3)
   - V_nt[:cut]           : V-REGION up to but NOT including Cys104 codon
                             (Cys104 is the first residue of the VDJdb CDR3)
   - CDR3_nt              : back-translated full CDR3 (Cys104 … Phe/Trp118)
   - _j_frame_and_fw(J_nt): detects J's coding frame (not always frame 0)
                             and the conserved Phe/Trp118 codon within it
   - J_nt[fw+3:]          : J-REGION AFTER the conserved Phe/Trp codon,
                             in the detected frame
                             (Phe/Trp118 is already the last residue of CDR3)

IMGT CDR3 boundaries used here follow VDJdb convention:
  CDR3 = Cys104 … Phe/Trp118  (boundaries inclusive)
"""
from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional

from .constant_regions import constant_aa
from .cdr_enricher import (
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
def _v_cys_cut(v_nt: str) -> int:
    """
    Return the nt index of the conserved Cys104 codon start in *v_nt*, so that
    v_nt[:_v_cys_cut(v_nt)] is codon-aligned and ends immediately BEFORE
    Cys104. The stitchr V-REGION sequence is not guaranteed to be a whole
    number of codons long (e.g. TRBV4-1 is 287 nt), so the naive "drop the
    last 3 nt" trick used previously is not codon-aligned in general.

    Cys104 is found by translating v_nt in frame 0 and taking the LAST "C"
    residue within the final 8 residues of the translation (the conserved
    Cys that begins the CDR3 sits right at the end of FR3/V-REGION).
    """
    aa = _translate(v_nt)
    tail_start = max(0, len(aa) - 8)
    cys_index = None
    for i in range(len(aa) - 1, tail_start - 1, -1):
        if aa[i] == "C":
            cys_index = i
            break
    if cys_index is not None:
        return cys_index * 3
    # Fallback: drop the trailing partial codon plus the last full codon
    # (codon-aligned version of the old "v_nt[:-3]" behavior).
    return len(v_nt) - (len(v_nt) % 3) - 3


def _j_frame_and_fw(j_nt: str) -> tuple[int, int]:
    """
    Determine the coding frame of a J-REGION nt sequence and locate the
    conserved Phe/Trp118 codon within it.

    stitchr J-REGION entries are not always frame 0 (e.g. TRBJ1-1 codes in
    frame 2), so each of the 3 frames is translated and scored: the frame
    whose translation has an F or W followed within 2 residues by a G-x-G
    motif, with no stop codon before that motif, is preferred.

    Returns (frame, fw_nt_index) where fw_nt_index is the nt offset (within
    j_nt) of the F/W codon's first base, or (0, -1) if no such motif is
    found in any frame.
    """
    # Phe/Trp118 sits immediately before the conserved G-x-G of FR4 (the
    # canonical F-G-x-G / W-G-x-G motif). Match that tight motif first so the
    # F/W adjacent to G-x-G is chosen, not an earlier F/W (e.g. the double F in
    # TRBJ1-1 "EAFFGQG", where Phe118 is the second F). Fall back to a looser
    # motif only when no tight one exists in any frame.
    fw_gxg_tight = re.compile(r"[FW]G.G")
    fw_gxg_loose = re.compile(r"[FW].{0,2}?G.G")

    candidates: list[tuple[int, int, int]] = []  # (frame, fw_res_index, n_stops)
    fewest_stops_frame = 0
    fewest_stops = None

    for frame in (0, 1, 2):
        translated = _translate(j_nt[frame:])
        n_stops = translated.count("*")
        if fewest_stops is None or n_stops < fewest_stops:
            fewest_stops = n_stops
            fewest_stops_frame = frame

        match = fw_gxg_tight.search(translated) or fw_gxg_loose.search(translated)
        if match is None:
            continue
        fw_res_index = match.start()
        # No stop codon before the motif.
        if "*" in translated[:fw_res_index]:
            continue
        candidates.append((frame, fw_res_index, n_stops))

    if candidates:
        # Prefer a candidate frame with no stop codons at all; otherwise
        # keep the first match found (frame order 0, 1, 2).
        clean = [c for c in candidates if c[2] == 0]
        frame, fw_res_index, _ = clean[0] if clean else candidates[0]
        return frame, fw_res_index * 3 + frame

    return fewest_stops_frame, -1


def _germline_key(gene_base: str) -> str:
    """Map an alpha/delta dual gene written with a dash (TRAV6-7-DV9) to the
    germline key convention with a slash (TRAV6-7/DV9). Idempotent for names
    already using a slash or with no DV segment."""
    return re.sub(r"-DV(\d)", r"/DV\1", gene_base)


def _j_suffix_by_overlap(j_nt: str, cdr3_aa: str) -> Optional[str]:
    """Return the J-REGION nt that contributes FR4 after the CDR3.

    The J germline 5' end overlaps the CDR3 (the CDR3 C-terminal residues are
    J encoded), so FR4 begins right after the longest CDR3 suffix that appears
    in the J N-terminal residues. This works for J genes with no canonical
    FG-x-G motif. Returns None if no overlap is found in any reading frame.
    """
    best = None  # (frame, fr4_residue_start, overlap_len)
    for frame in (0, 1, 2):
        aa = _translate(j_nt[frame:])
        head = aa[:15]  # FR4 starts early; do not match deep into the J
        for k in range(min(len(cdr3_aa), len(head)), 0, -1):
            pos = head.find(cdr3_aa[-k:])
            if pos >= 0 and "*" not in aa[: pos + k]:
                if best is None or k > best[2]:
                    best = (frame, pos + k, k)
                break
    if best is None:
        return None
    frame, fr4_start, _ = best
    return j_nt[frame + fr4_start * 3 :]


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

    v_nt = v_map.get(_germline_key(v_gene_base), "")
    j_nt = j_map.get(_germline_key(j_gene_base), "")

    cdr3_nt = back_translate(cdr3_aa)

    full_nt: Optional[str] = None
    full_aa: Optional[str] = None
    v_prefix: Optional[str] = None
    j_suffix: Optional[str] = None

    if v_nt and j_nt:
        # V-REGION from stitchr is not guaranteed codon-aligned to Cys104
        # (e.g. TRBV4-1 is 287 nt), so find the Cys104 codon boundary
        # directly rather than assuming the sequence ends 3 nt after it.
        v_prefix = v_nt[: _v_cys_cut(v_nt)]

        # FR4 begins right after the CDR3. The J germline 5' overlaps the CDR3
        # C-terminus, so align them to find where FR4 starts. This is robust for
        # J genes without a canonical FG-x-G motif (e.g. TRBJ1-6). Fall back to
        # motif-based frame detection only when no overlap is found.
        j_suffix = _j_suffix_by_overlap(j_nt, cdr3_aa)
        if j_suffix is None:
            frame, fw_nt = _j_frame_and_fw(j_nt)
            j_suffix = j_nt[fw_nt + 3 :] if fw_nt >= 0 else j_nt[frame:]

        full_nt = v_prefix + cdr3_nt + j_suffix
        full_aa = _translate(full_nt).rstrip("*") or None

    # Append the curated membrane-bound constant region (invariant per species
    # and chain) to form the full chain. Only the variable domain is derived
    # from V/J germline; the constant is vendored reference data.
    chain_name = {"TRA": "alpha", "TRB": "beta", "TRG": "gamma", "TRD": "delta"}.get(chain, "")
    constant = constant_aa(chain_name, species) if full_aa else None
    full_chain_aa = (full_aa + constant) if (full_aa and constant) else None
    constant_source = (
        f"curated membrane-bound {chain_name} constant ({species.lower()})"
        if constant else None
    )

    return {
        "v_gene": v_gene_base,
        "j_gene": j_gene_base,
        "cdr3_aa": cdr3_aa,
        "species": species,
        "full_nt": full_nt,
        "full_aa": full_aa,
        "full_chain_aa": full_chain_aa,
        "constant_source": constant_source,
        "v_region_nt": v_nt or None,
        "cdr3_nt": cdr3_nt,
        "j_region_nt": j_nt or None,
        # The actual in-frame pieces used in the assembly: the V germline up to
        # (not including) Cys104, and the J germline FR4 after Phe/Trp118. These
        # translate cleanly, unlike the raw region nt which are not frame 0.
        "v_prefix_nt": v_prefix,
        "j_suffix_nt": j_suffix,
        "v_found": bool(v_nt),
        "j_found": bool(j_nt),
        "note": (
            "CDR3 back-translated with human-optimised codons. "
            "V/J boundaries follow IMGT/VDJdb convention: CDR3 = Cys104…Phe/Trp118."
        ),
    }
