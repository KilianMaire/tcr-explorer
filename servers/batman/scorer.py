"""BATMAN scorer: wrapper around pybatman.functions with validation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

# pybatman imports — wrapped to allow mocking in tests
try:
    from pybatman.functions import train, peptide2index  # noqa: F401
except ImportError:
    def train(*args, **kwargs): raise RuntimeError("pybatman not installed")
    def peptide2index(*args, **kwargs): raise RuntimeError("pybatman not installed")

_REQUIRED_COLS = {"tcr", "index", "peptide", "activation"}
_MAX_DIST = 20.0  # empirical upper bound for normalisation


@dataclass
class TrainingResult:
    tcr_id: str
    index_peptide: str
    weights: np.ndarray
    aa_matrix: np.ndarray


def normalize_distance(dist: float, max_dist: float = _MAX_DIST) -> float:
    """Convert raw BATMAN distance to [0, 1] activation-likelihood score."""
    return float(max(0.0, min(1.0, 1.0 - dist / max_dist)))


class BATMANScorer:
    def train(
        self,
        df: pd.DataFrame,
        mode: str = "full",
        aa_matrix: str = "blosum100",
        steps: int = 20_000,
        seed: int = 100,
    ) -> TrainingResult:
        """Validate DataFrame then call pybatman.train().

        Raises:
            ValueError: if columns missing or activation labels not consecutive.
        """
        missing = _REQUIRED_COLS - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns: {missing}")

        labels = sorted(df["activation"].unique())
        expected = list(range(len(labels)))
        if labels != expected:
            raise ValueError(
                f"activation labels must be consecutive integers from 0. "
                f"Got {labels}, expected {expected}"
            )

        tcr_id = str(df["tcr"].iloc[0])
        index_peptide = str(df["index"].iloc[0])

        weights, aa_mat = train(df, mode, aa_matrix, steps=steps, seed=seed)
        return TrainingResult(
            tcr_id=tcr_id,
            index_peptide=index_peptide,
            weights=np.asarray(weights),
            aa_matrix=np.asarray(aa_mat),
        )

    def score(
        self,
        index_peptide: str,
        candidate_peptide: str,
        weights: np.ndarray,
        aa_matrix: np.ndarray,
        max_dist: float = _MAX_DIST,
    ) -> float:
        """Score how likely candidate_peptide activates the TCR that recognises index_peptide.

        Returns: float in [0, 1], where 1 = identical to reference.
        """
        distances = peptide2index(index_peptide, [candidate_peptide], aa_matrix, weights)
        return normalize_distance(float(distances[0]), max_dist=max_dist)
