"""TEMPO scorer: log-likelihood TCR-epitope interaction predictor.

Implements the TEMPO model from Liu et al. (2025):
  S(TCR) = sum_chain [ log(P(V)/Q(V)) + log(P(J)/Q(J))
           + log(P(L)/Q_P(V,J)(L)) + sum_i log(P(CDR3_i|L)/Q_P(V,J|L)(CDR3_i|L)) ]
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .baseline import AA_INDEX, BaselineModel


@dataclass
class TempoResult:
    log_likelihood: float
    v_contribution: float = 0.0
    j_contribution: float = 0.0
    length_contribution: float = 0.0
    cdr3_contribution: float = 0.0


def _safe_log_ratio(p: float, q: float) -> float:
    """Compute log(p / q) with pseudocount regularization.

    Args:
        p: Epitope-specific frequency (numerator).
        q: Background frequency (denominator).
    Returns:
        Log ratio, regularized to avoid division by zero.
    """
    p = max(p, 1e-10)
    q = max(q, 1e-10)
    return math.log(p / q)


def _cdr3_to_indices(cdr3: str) -> list[int]:
    """Convert CDR3 amino acid string to indices in the 20-letter alphabet."""
    indices = []
    for aa in cdr3.upper():
        idx = AA_INDEX.get(aa)
        if idx is not None:
            indices.append(idx)
        else:
            indices.append(0)  # unknown AA mapped to A
    return indices


class TempoScorer:
    """TEMPO TCR-epitope interaction scorer."""

    def __init__(
        self,
        alpha_baseline: Optional[BaselineModel] = None,
        beta_baseline: Optional[BaselineModel] = None,
    ) -> None:
        self._alpha_bl = alpha_baseline
        self._beta_bl = beta_baseline

    def score_single_chain(
        self,
        v_gene: str,
        j_gene: str,
        cdr3: str,
        chain: str,
    ) -> TempoResult:
        """Score a single TCR chain against baseline."""
        bl = self._alpha_bl if chain == "alpha" else self._beta_bl
        if bl is None:
            return TempoResult(log_likelihood=0.0)

        q_v = bl.q_v(v_gene)
        v_score = _safe_log_ratio(1.0, q_v) if q_v > 0 else 0.0

        q_j = bl.q_j(j_gene)
        j_score = _safe_log_ratio(1.0, q_j) if q_j > 0 else 0.0

        cdr3_len = len(cdr3)
        q_len = bl.q_length(v_gene, j_gene, cdr3_len)
        len_score = _safe_log_ratio(1.0, q_len) if q_len > 0 else 0.0

        q_cdr3_mat = bl.q_cdr3(v_gene, j_gene, cdr3_len)
        indices = _cdr3_to_indices(cdr3)
        cdr3_score = 0.0
        for pos, aa_idx in enumerate(indices):
            if pos < q_cdr3_mat.shape[1]:
                q_aa = float(q_cdr3_mat[aa_idx, pos])
                cdr3_score += _safe_log_ratio(1.0, q_aa)

        total = v_score + j_score + len_score + cdr3_score
        return TempoResult(
            log_likelihood=total,
            v_contribution=v_score,
            j_contribution=j_score,
            length_contribution=len_score,
            cdr3_contribution=cdr3_score,
        )

    def score_paired(
        self,
        v_a: str, j_a: str, cdr3_a: str,
        v_b: str, j_b: str, cdr3_b: str,
    ) -> TempoResult:
        """Score a paired TCRab by summing alpha and beta chain log-likelihoods."""
        alpha = self.score_single_chain(v_a, j_a, cdr3_a, "alpha")
        beta = self.score_single_chain(v_b, j_b, cdr3_b, "beta")
        return TempoResult(
            log_likelihood=alpha.log_likelihood + beta.log_likelihood,
            v_contribution=alpha.v_contribution + beta.v_contribution,
            j_contribution=alpha.j_contribution + beta.j_contribution,
            length_contribution=alpha.length_contribution + beta.length_contribution,
            cdr3_contribution=alpha.cdr3_contribution + beta.cdr3_contribution,
        )
