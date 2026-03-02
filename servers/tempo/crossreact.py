"""MixTCRcross: Cross-reactivity prediction for peptide-MHC variants.

Reimplements the algorithm from GfellerLab/MixTCRcross (R package)
based on Liu et al. (2025). Uses position-dependent penalties and
binding affinity to predict cross-reactivity risk.

Score scale: 0 (lowest risk) to 3 (highest risk / most cross-reactive)
Categories: Low (0-1), Medium (1-2), High (2-3)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

# MHC Class I anchor positions for 9-mers
_CLASS_I_ANCHORS_9MER = {2, 9}

# TCR-determining positions for 9-mers (non-anchor residues)
_TCR_DETERMINING_9MER = {3, 4, 5, 6, 7, 8}

# Position groups for 9-mers
_P3P5 = {3, 4, 5}
_P6P7 = {6, 7}

# Affinity thresholds
_AFFINITY_CUTOFF = 2.0  # %rank above this = non-binder
_AFFINITY_MIN = 1e-6    # log floor
_AFFINITY_MAX = 100.0   # %rank ceiling

# MHC TCR-contact residue positions (from MHC_seq.rda)
# These are the key positions in the MHC alpha chain that contact the TCR
MHC_TCR_CONTACTS: dict[str, dict[int, str] | set[int]] = {
    "HLA-A*02:01": {62: "G", 65: "Q", 66: "I", 69: "A", 72: "Q", 73: "T", 76: "V",
                     150: "A", 151: "H", 152: "V", 155: "Q", 158: "A", 159: "Y",
                     162: "T", 163: "L", 166: "E"},
    "HLA-A*02:05": {62: "G", 65: "Q", 66: "I", 69: "A", 72: "Q", 73: "T", 76: "V",
                     150: "A", 151: "H", 152: "V", 155: "Q", 156: "W", 158: "A", 159: "Y",
                     162: "T", 163: "L", 166: "E"},
    "HLA-A*01:01": {62: "G", 65: "Q", 66: "I", 69: "A", 72: "Q", 73: "T", 76: "V",
                     150: "A", 151: "H", 152: "V", 155: "Q", 158: "A", 159: "Y",
                     162: "T", 163: "L", 166: "E"},
    "HLA-B*07:02": {62: "G", 65: "Q", 66: "N", 69: "A", 72: "Q", 73: "T", 76: "E",
                     150: "A", 151: "H", 152: "V", 155: "Q", 158: "A", 159: "Y",
                     162: "T", 163: "L", 166: "E"},
    "HLA-B*08:01": {62: "G", 65: "Q", 66: "N", 69: "A", 72: "Q", 73: "T", 76: "E",
                     150: "A", 151: "H", 152: "V", 155: "Q", 158: "A", 159: "Y",
                     162: "T", 163: "L", 166: "E"},
}

# Class II MHC-TCR contact residues from published crystal structures
# Positions are 1-indexed relative to the 9-mer binding core
MHC_TCR_CONTACTS.update({
    "HLA-DRB1*01:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*03:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*04:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*04:05": {2, 3, 5, 7, 8},
    "HLA-DRB1*07:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*08:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*11:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*13:01": {2, 3, 5, 7, 8},
    "HLA-DRB1*15:01": {2, 3, 5, 7, 8},
    "HLA-DQ2.5": {2, 3, 5, 7, 8},
    "HLA-DQA1*05:01/DQB1*02:01": {2, 3, 5, 7, 8},
    "HLA-DQ8": {2, 3, 5, 7, 8},
    "HLA-DQA1*03:01/DQB1*03:02": {2, 3, 5, 7, 8},
    "HLA-DPB1*04:01": {2, 3, 5, 7, 8},
    "HLA-DP4": {2, 3, 5, 7, 8},
})

CLASS2_ANCHORS: dict[str, dict[int, str]] = {
    "DRB1*01:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DRB1*03:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DRB1*04:01": {1: "aromatic", 4: "small", 6: "negative", 9: "hydrophobic"},
    "DRB1*04:05": {1: "aromatic", 4: "small", 6: "negative", 9: "hydrophobic"},
    "DRB1*07:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DRB1*08:01": {1: "hydrophobic", 4: "negative", 6: "small", 9: "hydrophobic"},
    "DRB1*11:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DRB1*13:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DRB1*15:01": {1: "hydrophobic", 4: "positive", 6: "small", 9: "hydrophobic"},
    "DQ2.5": {1: "hydrophobic", 4: "negative", 6: "hydrophobic", 7: "small", 9: "positive"},
    "DQA1*05:01/DQB1*02:01": {1: "hydrophobic", 4: "negative", 6: "hydrophobic", 7: "small", 9: "positive"},
    "DQ8": {1: "hydrophobic", 4: "negative", 6: "hydrophobic", 9: "hydrophobic"},
    "DQA1*03:01/DQB1*03:02": {1: "hydrophobic", 4: "negative", 6: "hydrophobic", 9: "hydrophobic"},
    "DPB1*04:01": {1: "charged", 6: "hydrophobic"},
    "DP4": {1: "charged", 6: "hydrophobic"},
}
_CLASS_II_ANCHORS_GENERIC: set[int] = {1, 4, 6, 9}


def _standardize_allele(allele: str) -> str:
    """Normalize MHC allele name to standard format."""
    a = allele.strip().replace(" ", "")
    # HLA-A0201 -> HLA-A*02:01
    if a.startswith("HLA-") and "*" not in a and len(a) >= 8:
        gene = a[4]
        digits = a[5:]
        if len(digits) == 4 and digits.isdigit():
            a = f"HLA-{gene}*{digits[:2]}:{digits[2:]}"
    return a


def _tcr_contacts_match(allele1: str, allele2: str) -> bool:
    """Check if two MHC alleles have identical TCR-contact residues."""
    a1 = _standardize_allele(allele1)
    a2 = _standardize_allele(allele2)

    if a1 == a2:
        return True

    contacts1 = MHC_TCR_CONTACTS.get(a1)
    contacts2 = MHC_TCR_CONTACTS.get(a2)

    # If we don't have data for either allele, assume conserved
    if contacts1 is None or contacts2 is None:
        return True

    return contacts1 == contacts2


def _affinity_to_score(affinity_rank: float) -> float:
    """Convert binding affinity (%rank) to cross-reactivity score component.

    Formula from MixTCRcross R code:
    score = 2 - (log10(affinity) - log10(1e-6)) / (log10(100) - log10(1e-6)) + 1
    """
    if affinity_rank <= 0:
        affinity_rank = _AFFINITY_MIN
    if affinity_rank > _AFFINITY_MAX:
        affinity_rank = _AFFINITY_MAX

    log_range = math.log10(_AFFINITY_MAX) - math.log10(_AFFINITY_MIN)
    normalized = (math.log10(affinity_rank) - math.log10(_AFFINITY_MIN)) / log_range
    return 2.0 - normalized + 1.0


def classify_position(pos: int, peptide_length: int, mhc_class: str) -> str:
    """Classify a 1-indexed peptide position for a 9-mer."""
    if mhc_class == "II":
        # Class II anchors: 1, 4, 6, 9 for core 9-mer
        if pos in {1, 4, 6, 9}:
            return "anchor"
        if pos in {2, 3, 5}:
            return "nonanchor_p3p5"
        if pos in {7, 8}:
            return "nonanchor_p6p7"
        return "other"

    if peptide_length != 9:
        # Non-9-mers: basic anchor detection
        if pos in {2, peptide_length}:
            return "anchor"
        return "other"

    if pos in _CLASS_I_ANCHORS_9MER:
        return "anchor"
    if pos in _P3P5:
        return "nonanchor_p3p5"
    if pos in _P6P7:
        return "nonanchor_p6p7"
    return "other"


@dataclass
class CrossReactivityResult:
    score: float               # 0-3 scale
    category: str              # "High", "Medium", "Low"
    variant_positions: list[dict] = field(default_factory=list)
    mhc_interface_conserved: bool = True
    length_match: bool = True
    normalized_score: float = 0.0  # 0-1 scale for backward compat


def _categorize(score: float) -> str:
    """Categorize cross-reactivity score."""
    if score >= 2.0:
        return "High"
    if score >= 1.0:
        return "Medium"
    return "Low"


def predict_cross_reactivity(
    reference_peptide: str,
    variant_peptide: str,
    mhc_allele: str,
    mhc_class: str = "I",
    reference_mhc: Optional[str] = None,
    affinity_ref: Optional[float] = None,
    affinity_var: Optional[float] = None,
) -> CrossReactivityResult:
    """Predict cross-reactivity between a reference and variant peptide.

    Implements the MixTCRcross algorithm (Liu et al. 2025).

    Args:
        reference_peptide: Reference peptide sequence (9-mer recommended)
        variant_peptide: Variant peptide to test
        mhc_allele: MHC allele for the variant peptide
        mhc_class: "I" or "II"
        reference_mhc: MHC allele for reference (defaults to mhc_allele)
        affinity_ref: Binding affinity %rank for reference peptide (optional)
        affinity_var: Binding affinity %rank for variant peptide (optional)

    Returns:
        CrossReactivityResult with score (0-3), category, and details.
    """
    ref = reference_peptide.upper()
    var = variant_peptide.upper()
    ref_mhc = reference_mhc or mhc_allele
    query_mhc = mhc_allele

    length_match = len(ref) == len(var)

    # Non-9-mers: reduced confidence, basic comparison only
    if len(ref) != 9 or len(var) != 9:
        if not length_match:
            return CrossReactivityResult(
                score=0.0, category="Low",
                variant_positions=[{"position": 0, "classification": "length_change"}],
                length_match=False, normalized_score=0.0,
            )
        # For non-9-mers, count differences and score simply
        diffs = sum(1 for a, b in zip(ref, var) if a != b)
        score = max(0.0, 3.0 * (1.0 - diffs / len(ref)))
        return CrossReactivityResult(
            score=score, category=_categorize(score),
            length_match=True, normalized_score=score / 3.0,
        )

    # Identical peptides -> maximum cross-reactivity
    if ref == var:
        return CrossReactivityResult(
            score=3.0, category="High",
            mhc_interface_conserved=_tcr_contacts_match(ref_mhc, query_mhc),
            length_match=True, normalized_score=1.0,
        )

    # Find variant positions
    variant_positions = []
    tcr_det_diffs = 0
    p3p5_diff = False
    p6p7_diff = False

    for i in range(9):
        pos = i + 1
        if ref[i] != var[i]:
            pos_class = classify_position(pos, 9, mhc_class)
            variant_positions.append({
                "position": pos,
                "reference_aa": ref[i],
                "variant_aa": var[i],
                "classification": pos_class,
            })
            if pos in _TCR_DETERMINING_9MER:
                tcr_det_diffs += 1
            if pos in _P3P5:
                p3p5_diff = True
            if pos in _P6P7:
                p6p7_diff = True

    # Check MHC TCR-contact conservation
    mhc_conserved = _tcr_contacts_match(ref_mhc, query_mhc)

    # If >1 TCR-determining position differs -> Low
    if tcr_det_diffs > 1:
        return CrossReactivityResult(
            score=0.0, category="Low",
            variant_positions=variant_positions,
            mhc_interface_conserved=mhc_conserved,
            length_match=True, normalized_score=0.0,
        )

    # If MHC TCR contacts differ -> Low
    if not mhc_conserved:
        return CrossReactivityResult(
            score=0.5, category="Low",
            variant_positions=variant_positions,
            mhc_interface_conserved=False,
            length_match=True, normalized_score=0.5 / 3.0,
        )

    # Affinity-based scoring
    if affinity_var is not None:
        # Poor binder -> Low
        if affinity_var > _AFFINITY_CUTOFF:
            return CrossReactivityResult(
                score=0.0, category="Low",
                variant_positions=variant_positions,
                mhc_interface_conserved=mhc_conserved,
                length_match=True, normalized_score=0.0,
            )
        score = _affinity_to_score(affinity_var)
    else:
        # No affinity data: use position-based heuristic
        # Start at 3.0 (max) and reduce based on mutations
        score = 3.0

    # Apply position-based caps (exclusive: keep out of the next category)
    if p3p5_diff:
        score = min(score, 0.999)  # Cap below Medium boundary -> stays Low
    elif p6p7_diff:
        score = min(score, 1.999)  # Cap below High boundary -> stays at most Medium

    score = max(0.0, min(3.0, score))

    return CrossReactivityResult(
        score=score,
        category=_categorize(score),
        variant_positions=variant_positions,
        mhc_interface_conserved=mhc_conserved,
        length_match=True,
        normalized_score=score / 3.0,
    )
