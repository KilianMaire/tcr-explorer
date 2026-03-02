"""TCR Specificity Profile (TSP) computation.

Compares epitope-specific TCR sets against baseline TCR repertoires
to identify enriched V/J genes, CDR3 length distributions, and CDR3 motifs.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .baseline import AA_INDEX, BaselineModel


@dataclass
class TSProfile:
    """TCR Specificity Profile for one chain."""
    v_enrichment: dict[str, dict[str, float]] = field(default_factory=dict)
    j_enrichment: dict[str, dict[str, float]] = field(default_factory=dict)
    cdr3_length_dist: dict[int, float] = field(default_factory=dict)
    cdr3_motifs: dict[int, np.ndarray] = field(default_factory=dict)
    n_tcrs: int = 0

    def to_dict(self) -> dict[str, Any]:
        """Serialize to JSON-friendly dict."""
        motifs = {}
        for length, mat in self.cdr3_motifs.items():
            motifs[str(length)] = mat.tolist()
        return {
            "v_enrichment": self.v_enrichment,
            "j_enrichment": self.j_enrichment,
            "cdr3_length_dist": {str(k): v for k, v in self.cdr3_length_dist.items()},
            "cdr3_motifs": motifs,
            "n_tcrs": self.n_tcrs,
        }


def compute_tsp(
    tcrs: list[dict[str, str]],
    baseline: BaselineModel,
) -> TSProfile:
    """Compute TSP comparing a set of epitope-specific TCRs against baseline."""
    if not tcrs:
        return TSProfile()

    n = len(tcrs)

    # V usage
    v_counts = Counter(t["v_gene"] for t in tcrs)
    v_enrichment = {}
    for v_gene, count in v_counts.items():
        p_v = count / n
        q_v = baseline.q_v(v_gene)
        fold_change = p_v / q_v if q_v > 0 else float("inf")
        z_score = baseline.v_zscore(v_gene, p_v)
        v_enrichment[v_gene] = {
            "observed_freq": p_v,
            "baseline_freq": q_v,
            "fold_change": fold_change,
            "z_score": z_score,
        }

    # J usage
    j_counts = Counter(t["j_gene"] for t in tcrs)
    j_enrichment = {}
    for j_gene, count in j_counts.items():
        p_j = count / n
        q_j = baseline.q_j(j_gene)
        fold_change = p_j / q_j if q_j > 0 else float("inf")
        z_score = baseline.j_zscore(j_gene, p_j)
        j_enrichment[j_gene] = {
            "observed_freq": p_j,
            "baseline_freq": q_j,
            "fold_change": fold_change,
            "z_score": z_score,
        }

    # CDR3 length distribution
    lengths = [len(t["cdr3"]) for t in tcrs]
    length_counts = Counter(lengths)
    cdr3_length_dist = {length: count / n for length, count in length_counts.items()}

    # CDR3 amino acid motifs per length
    cdr3_motifs: dict[int, np.ndarray] = {}
    by_length: dict[int, list[str]] = {}
    for t in tcrs:
        cdr3 = t["cdr3"]
        by_length.setdefault(len(cdr3), []).append(cdr3)

    for length, seqs in by_length.items():
        mat = np.zeros((20, length))
        for seq in seqs:
            for pos, aa in enumerate(seq):
                idx = AA_INDEX.get(aa.upper())
                if idx is not None and pos < length:
                    mat[idx, pos] += 1
        total = mat.sum(axis=0, keepdims=True)
        total[total == 0] = 1
        mat /= total
        cdr3_motifs[length] = mat

    return TSProfile(
        v_enrichment=v_enrichment,
        j_enrichment=j_enrichment,
        cdr3_length_dist=cdr3_length_dist,
        cdr3_motifs=cdr3_motifs,
        n_tcrs=n,
    )
