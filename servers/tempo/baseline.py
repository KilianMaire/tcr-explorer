"""TEMPO baseline TCR repertoire model.

Stores Q distributions: Q(V), Q(J), Q(L|V,J), Q(CDR3|V,J,L) derived from
iReceptor data for {human,mouse} x {alpha,beta} chains.
"""
from __future__ import annotations

import gzip
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

_AA_ALPHABET_SIZE = 20
_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "baselines"

# Standard amino acid order (IMGT convention)
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}


@dataclass
class BaselineModel:
    """Baseline TCR repertoire model for one species+chain combination."""

    species: str  # "human" or "mouse"
    chain: str    # "alpha" or "beta"
    v_freq: dict[str, float] = field(default_factory=dict)
    j_freq: dict[str, float] = field(default_factory=dict)
    v_std: dict[str, float] = field(default_factory=dict)
    j_std: dict[str, float] = field(default_factory=dict)
    length_dist: dict[tuple[str, str], dict[int, float]] = field(default_factory=dict)
    cdr3_freq: dict[tuple[str, str, int], np.ndarray] = field(default_factory=dict)

    _PSEUDOCOUNT: float = 50.0

    def q_v(self, v_gene: str) -> float:
        """Return baseline frequency Q(V) for a V gene."""
        return self.v_freq.get(v_gene, 0.0)

    def q_j(self, j_gene: str) -> float:
        """Return baseline frequency Q(J) for a J gene."""
        return self.j_freq.get(j_gene, 0.0)

    def q_length(self, v_gene: str, j_gene: str, length: int) -> float:
        """Return Q(L|V,J) -- CDR3 length probability given V,J usage."""
        vj = (v_gene, j_gene)
        if vj not in self.length_dist:
            return 0.0
        return self.length_dist[vj].get(length, 0.0)

    def q_cdr3(self, v_gene: str, j_gene: str, length: int) -> np.ndarray:
        """Return Q(CDR3|V,J,L) as a (20, L) amino acid frequency matrix.

        Returns uniform distribution if the (V,J,L) combination is unknown.
        """
        key = (v_gene, j_gene, length)
        if key in self.cdr3_freq:
            return self.cdr3_freq[key]
        return np.full((_AA_ALPHABET_SIZE, length), 1.0 / _AA_ALPHABET_SIZE)

    def v_zscore(self, v_gene: str, observed_freq: float) -> float:
        """Z-score of observed V usage vs baseline."""
        q = self.q_v(v_gene)
        std = self.v_std.get(v_gene, 0.01)
        if std <= 0:
            return 0.0
        return (observed_freq - q) / std

    def j_zscore(self, j_gene: str, observed_freq: float) -> float:
        """Z-score of observed J usage vs baseline."""
        q = self.q_j(j_gene)
        std = self.j_std.get(j_gene, 0.01)
        if std <= 0:
            return 0.0
        return (observed_freq - q) / std


def load_baseline(species: str, chain: str) -> BaselineModel:
    """Load a precomputed baseline from data/baselines/{species}_{chain}.json.gz.

    Raises FileNotFoundError if the baseline file is not present.
    """
    path = _DATA_DIR / f"{species}_{chain}.json.gz"
    if not path.exists():
        raise FileNotFoundError(f"Baseline not found: {path}")

    with gzip.open(path, "rt", encoding="utf-8") as f:
        data = json.load(f)

    cdr3_freq: dict[tuple[str, str, int], np.ndarray] = {}
    for key_str, matrix_list in data.get("cdr3_freq", {}).items():
        parts = key_str.split("|")
        v, j, length = parts[0], parts[1], int(parts[2])
        cdr3_freq[(v, j, length)] = np.array(matrix_list)

    length_dist: dict[tuple[str, str], dict[int, float]] = {}
    for key_str, dist in data.get("length_dist", {}).items():
        v, j = key_str.split("|")
        length_dist[(v, j)] = {int(k): float(v_) for k, v_ in dist.items()}

    return BaselineModel(
        species=species,
        chain=chain,
        v_freq=data.get("v_freq", {}),
        j_freq=data.get("j_freq", {}),
        v_std=data.get("v_std", {}),
        j_std=data.get("j_std", {}),
        length_dist=length_dist,
        cdr3_freq=cdr3_freq,
    )
