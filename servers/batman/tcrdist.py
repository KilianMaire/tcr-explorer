"""TCR similarity scoring using tcrdist3.

For a query CDR3+Vgene, finds the minimum tcrdist to any known binder
of the target epitope in the reference set, and normalizes to [0, 1].
"""
from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

# tcrdist3 import wrapped to allow graceful failure if not installed
try:
    from tcrdist.repertoire import TCRrep
except ImportError:
    TCRrep = None  # type: ignore

logger = logging.getLogger(__name__)

# Class-specific normalization constants
_MAX_TCRDIST_CLASS_I: float = 96.0   # 99th percentile for human beta (Dash et al. 2017)
_MAX_TCRDIST_CLASS_II: float = 120.0  # CD4+ TCRs: longer CDR3, empirical 99th percentile
_MAX_TCRDIST = _MAX_TCRDIST_CLASS_I   # Backward compatibility
_REQUIRED_COLS = {"cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope", "count"}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def normalize_tcrdist(raw_dist: float, mhc_class: str = "I", *, max_dist: float | None = None) -> float:
    """Convert raw tcrdist to [0, 1] score.

    Args:
        raw_dist: Raw tcrdist3 distance value.
        mhc_class: "I" for Class I (CD8+), "II" for Class II (CD4+).
        max_dist: Explicit override for the normalization ceiling (backward compat).
    """
    if max_dist is not None:
        ceiling = max_dist
    else:
        ceiling = _MAX_TCRDIST_CLASS_II if mhc_class == "II" else _MAX_TCRDIST_CLASS_I
    if ceiling <= 0:
        return 1.0 if raw_dist <= 0 else 0.0
    return float(max(0.0, min(1.0, 1.0 - raw_dist / ceiling)))


def build_reference_df(vdjdb_records: list[dict], chain: str = "beta") -> pd.DataFrame:
    """Convert VDJdb GeneRecord dicts → tcrdist3-compatible DataFrame.

    Keeps only records with v_gene and j_gene present.
    Normalises V/J gene names to tcrdist3 format (TRBVXX*01 allele suffix).

    Parameters
    ----------
    vdjdb_records : list[dict]
        VDJdb-style records with ``sequence`` and ``metadata`` keys.
    chain : str
        ``"beta"`` (default) or ``"alpha"``.
    """
    rows = []
    for rec in vdjdb_records:
        meta = rec.get("metadata") or {}
        epitope = (meta.get("antigen_epitope") or rec.get("antigen_epitope") or "").strip()

        if chain == "beta":
            cdr3 = (rec.get("sequence") or "").strip()
            v_raw = (meta.get("v_segm") or "").strip()
            j_raw = (meta.get("j_segm") or "").strip()
            if not cdr3 or not v_raw or not j_raw:
                continue
            v_gene = v_raw if "*" in v_raw else f"{v_raw}*01"
            j_gene = j_raw if "*" in j_raw else f"{j_raw}*01"
            rows.append({
                "cdr3_b_aa": cdr3,
                "v_b_gene": v_gene,
                "j_b_gene": j_gene,
                "epitope": epitope,
                "count": 1,
            })
        else:
            # Alpha chain: try cdr3_a metadata first, fall back to sequence
            cdr3 = (meta.get("cdr3_a") or "").strip()
            v_raw = (meta.get("v_a_segm") or "").strip()
            j_raw = (meta.get("j_a_segm") or "").strip()
            if not cdr3 or not v_raw or not j_raw:
                cdr3 = (rec.get("sequence") or "").strip()
                v_raw = (meta.get("v_segm") or "").strip()
                j_raw = (meta.get("j_segm") or "").strip()
            if not cdr3 or not v_raw or not j_raw:
                continue
            v_gene = v_raw if "*" in v_raw else f"{v_raw}*01"
            j_gene = j_raw if "*" in j_raw else f"{j_raw}*01"
            rows.append({
                "cdr3_a_aa": cdr3,
                "v_a_gene": v_gene,
                "j_a_gene": j_gene,
                "epitope": epitope,
                "count": 1,
            })

    cols = ({"cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope", "count"} if chain == "beta"
            else {"cdr3_a_aa", "v_a_gene", "j_a_gene", "epitope", "count"})
    if not rows:
        return pd.DataFrame(columns=list(cols))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Scorer class
# ---------------------------------------------------------------------------
class TCRdistScorer:
    """Compute TCR-epitope similarity using nearest-neighbour tcrdist.

    Supports both alpha and beta chains, human and mouse organisms.

    Usage:
        scorer = TCRdistScorer()
        scorer.load_reference(vdjdb_records)                    # beta (default)
        scorer.load_reference(alpha_records, chain="alpha")     # alpha
        sim = scorer.score(cdr3_b, v_gene, j_gene, epitope)    # legacy beta
        sim = scorer.score(cdr3="X", v_gene="V", j_gene="J",
                           epitope="E", chain="alpha")          # alpha
        sim = scorer.score_paired(cdr3_a=..., v_a=..., j_a=...,
                                  cdr3_b=..., v_b=..., j_b=...,
                                  epitope=...)                  # paired
    """

    def __init__(self) -> None:
        self._ref_beta: pd.DataFrame = pd.DataFrame(columns=list(_REQUIRED_COLS))
        self._ref_alpha: pd.DataFrame = pd.DataFrame(
            columns=["cdr3_a_aa", "v_a_gene", "j_a_gene", "epitope", "count"],
        )

    # Legacy property: existing code that reads/writes scorer._ref_df
    # transparently maps to _ref_beta.
    @property
    def _ref_df(self) -> pd.DataFrame:
        return self._ref_beta

    @_ref_df.setter
    def _ref_df(self, value: pd.DataFrame) -> None:
        self._ref_beta = value

    def load_reference(self, vdjdb_records: list[dict], chain: str = "beta") -> None:
        """Load reference binders from VDJdb records."""
        if chain == "beta":
            self._ref_beta = build_reference_df(vdjdb_records, chain="beta")
            logger.info("Loaded %d beta reference binders", len(self._ref_beta))
        else:
            self._ref_alpha = build_reference_df(vdjdb_records, chain="alpha")
            logger.info("Loaded %d alpha reference binders", len(self._ref_alpha))

    def score(
        self,
        cdr3: str = "",
        v_gene: str = "",
        j_gene: str = "",
        epitope: str = "",
        chain: str = "beta",
        organism: str = "human",
        max_dist: float = _MAX_TCRDIST,
        # Legacy param for backward compatibility
        cdr3_b: str = "",
        **kwargs,
    ) -> Optional[float]:
        """Return similarity score [0, 1] or None if no reference binders exist.

        Computes tcrdist from query CDR3 to all known binders of the target
        epitope in the reference set, and returns 1 - (min_dist / max_dist).

        Parameters
        ----------
        cdr3 : str
            CDR3 amino-acid sequence (preferred parameter name).
        cdr3_b : str
            Legacy alias for *cdr3* (backward compatibility).
        v_gene, j_gene : str
            V/J gene names.
        epitope : str
            Target epitope.
        chain : str
            ``"beta"`` (default) or ``"alpha"``.
        organism : str
            ``"human"`` (default) or ``"mouse"``.
        max_dist : float
            Maximum distance for normalisation.
        """
        # Backward compat: if cdr3_b given but not cdr3
        if not cdr3 and cdr3_b:
            cdr3 = cdr3_b

        if TCRrep is None:
            logger.warning("tcrdist3 not installed; returning None")
            return None

        ref_df = self._ref_beta if chain == "beta" else self._ref_alpha
        if len(ref_df) == 0:
            return None

        cdr3_col = "cdr3_b_aa" if chain == "beta" else "cdr3_a_aa"
        v_col = "v_b_gene" if chain == "beta" else "v_a_gene"
        j_col = "j_b_gene" if chain == "beta" else "j_a_gene"
        chains_param = ["beta"] if chain == "beta" else ["alpha"]

        # Filter reference to target epitope binders
        ref_epitope = ref_df[ref_df["epitope"].str.upper() == epitope.upper()].copy()
        if len(ref_epitope) == 0:
            return None

        v_allele = v_gene if "*" in v_gene else f"{v_gene}*01"
        j_allele = j_gene if "*" in j_gene else f"{j_gene}*01"

        query_df = pd.DataFrame([{
            cdr3_col: cdr3,
            v_col: v_allele,
            j_col: j_allele,
            "count": 1,
        }])

        combined = pd.concat(
            [query_df, ref_epitope[[cdr3_col, v_col, j_col, "count"]]],
            ignore_index=True,
        )

        try:
            tr = TCRrep(
                cell_df=combined,
                organism=organism,
                chains=chains_param,
                compute_distances=True,
            )
            dist_attr = f"pw_{chain}"
            dist_matrix = getattr(tr, dist_attr)
            # Row 0 is query; columns 1..n are distances to reference
            dists_to_ref = dist_matrix[0, 1:]
            min_dist = float(np.min(dists_to_ref))
            return normalize_tcrdist(min_dist, max_dist=max_dist)
        except Exception as exc:
            logger.warning("tcrdist3 computation failed: %s", exc)
            return None

    def score_paired(
        self,
        cdr3_a: str,
        v_a: str,
        j_a: str,
        cdr3_b: str,
        v_b: str,
        j_b: str,
        epitope: str,
        organism: str = "human",
        max_dist: float = _MAX_TCRDIST,
    ) -> Optional[float]:
        """Score paired TCR alpha-beta: average of alpha and beta scores."""
        alpha = self.score(
            cdr3=cdr3_a, v_gene=v_a, j_gene=j_a,
            epitope=epitope, chain="alpha", organism=organism,
            max_dist=max_dist,
        )
        beta = self.score(
            cdr3=cdr3_b, v_gene=v_b, j_gene=j_b,
            epitope=epitope, chain="beta", organism=organism,
            max_dist=max_dist,
        )
        scores = [s for s in (alpha, beta) if s is not None]
        if not scores:
            return None
        return sum(scores) / len(scores)
