"""Tests for the offline tcrdist engine.

The parity tests are the honesty gate: they assert that our pwseqdist + vendored
table computation reproduces tcrdist3's own TCRrep, so we may honestly label the
engine "tcrdist". They skip when tcrdist3 (a dev only dependency) is absent.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from tcr_explorer import tcrdist_engine as te

pytest.importorskip("pwseqdist")

_HAS_TCRDIST3 = False
try:
    import tcrdist.repertoire  # noqa: F401

    _HAS_TCRDIST3 = True
except Exception:
    _HAS_TCRDIST3 = False


# Reference receptors used across parity checks (human beta).
_BETA = [
    {"v_b_gene": "TRBV19*01", "cdr3_b_aa": "CASSIRSSYEQYF"},
    {"v_b_gene": "TRBV28*01", "cdr3_b_aa": "CASSLGQAYEQYF"},
    {"v_b_gene": "TRBV7-9*01", "cdr3_b_aa": "CASSLAPGATNEKLFF"},
    {"v_b_gene": "TRBV19*01", "cdr3_b_aa": "CASSIRSAYEQYF"},
]

# Human alpha/beta paired receptors for alpha and paired parity checks.
_PAIRED = [
    {"v_a_gene": "TRAV12-1*01", "cdr3_a_aa": "CAVNFGGGKLIF", "v_b_gene": "TRBV19*01", "cdr3_b_aa": "CASSIRSSYEQYF"},
    {"v_a_gene": "TRAV1-2*01", "cdr3_a_aa": "CAVRDSNYQLIW", "v_b_gene": "TRBV28*01", "cdr3_b_aa": "CASSLGQAYEQYF"},
    {"v_a_gene": "TRAV12-1*01", "cdr3_a_aa": "CAVNFGAGKLIF", "v_b_gene": "TRBV19*01", "cdr3_b_aa": "CASSIRSAYEQYF"},
]


def _tcrrep(cells, chains):
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = pd.DataFrame(cells)
        df["count"] = 1
        return TCRrep(cell_df=df, organism="human", chains=chains,
                      db_file="alphabeta_gammadelta_db.tsv", compute_distances=True)


def test_available_true_with_pwseqdist():
    assert te.tcrdist_available() is True


def test_resolve_v_loops_returns_three_loops():
    loops = te.resolve_v_loops("TRBV19*01", "human")
    assert loops is not None
    assert len(loops) == 3
    # CDR2.5 (pMHC) loop for human TRBV19*01, gap dots preserved.
    assert loops[2] == "E.KKES"


def test_resolve_v_loops_defaults_to_star01():
    # A bare gene name resolves to the *01 allele.
    assert te.resolve_v_loops("TRBV19", "human") == te.resolve_v_loops("TRBV19*01", "human")


def test_resolve_v_loops_unknown_returns_none():
    assert te.resolve_v_loops("TRBVNONSENSE*01", "human") is None
    assert te.resolve_v_loops("", "human") is None


def test_identical_receptors_distance_zero():
    d = te.tcrdist_pair("CASSIRSSYEQYF", "TRBV19*01", "CASSIRSSYEQYF", "TRBV19*01")
    assert d == 0.0


def test_one_to_many_skips_unresolvable_v():
    dists, skipped = te.tcrdist_one_to_many(
        "CASSIRSSYEQYF", "TRBV19*01",
        ["CASSLGQAYEQYF", "CASSXXX"],
        ["TRBV28*01", "TRBVNONSENSE*01"],
    )
    assert skipped == 1
    assert dists[0] is not None
    assert dists[1] is None


def test_one_to_many_raises_when_query_v_unresolvable():
    with pytest.raises(ValueError):
        te.tcrdist_one_to_many("CASSIRSSYEQYF", "TRBVNONSENSE*01", ["CASSLGQAYEQYF"], ["TRBV28*01"])


@pytest.mark.skipif(not _HAS_TCRDIST3, reason="tcrdist3 not installed (dev only parity dependency)")
def test_parity_pair_matches_tcrdist3_beta():
    tr = _tcrrep(_BETA, ["beta"])
    clone = [(r.v_b_gene, r.cdr3_b_aa) for r in tr.clone_df.itertuples(index=False)]
    ref = np.asarray(tr.pw_beta, dtype=float)
    mine = np.array([[te.tcrdist_pair(ci, vi, cj, vj, "human") for vj, cj in clone] for vi, ci in clone])
    assert np.array_equal(mine, ref), f"\nmine=\n{mine}\nref=\n{ref}"


@pytest.mark.skipif(not _HAS_TCRDIST3, reason="tcrdist3 not installed (dev only parity dependency)")
def test_parity_one_to_many_matches_tcrdist3_beta():
    tr = _tcrrep(_BETA, ["beta"])
    clone = [(r.v_b_gene, r.cdr3_b_aa) for r in tr.clone_df.itertuples(index=False)]
    ref = np.asarray(tr.pw_beta, dtype=float)
    q_v, q_c = clone[0]
    dists, skipped = te.tcrdist_one_to_many(q_c, q_v, [c for _, c in clone], [v for v, _ in clone], "human")
    assert skipped == 0
    assert dists == list(ref[0])


@pytest.mark.skipif(not _HAS_TCRDIST3, reason="tcrdist3 not installed (dev only parity dependency)")
def test_parity_pair_matches_tcrdist3_alpha():
    # The engine is chain agnostic: the same tcrdist_pair reproduces the alpha chain.
    tr = _tcrrep(_PAIRED, ["alpha", "beta"])
    clone = [(r.v_a_gene, r.cdr3_a_aa) for r in tr.clone_df.itertuples(index=False)]
    ref = np.asarray(tr.pw_alpha, dtype=float)
    mine = np.array([[te.tcrdist_pair(ci, vi, cj, vj, "human") for vj, cj in clone] for vi, ci in clone])
    assert np.array_equal(mine, ref), f"\nmine=\n{mine}\nref=\n{ref}"


@pytest.mark.skipif(not _HAS_TCRDIST3, reason="tcrdist3 not installed (dev only parity dependency)")
def test_parity_paired_matches_tcrdist3_alpha_plus_beta():
    # Paired tcrdist is the elementwise sum of the alpha and beta single chain matrices.
    tr = _tcrrep(_PAIRED, ["alpha", "beta"])
    clone = [te.PairedTCR(r.cdr3_a_aa, r.v_a_gene, r.cdr3_b_aa, r.v_b_gene)
             for r in tr.clone_df.itertuples(index=False)]
    ref = np.asarray(tr.pw_alpha, dtype=float) + np.asarray(tr.pw_beta, dtype=float)
    mine = np.array([[te.tcrdist_paired(qi, cj, "human") for cj in clone] for qi in clone])
    assert np.array_equal(mine, ref), f"\nmine=\n{mine}\nref=\n{ref}"


@pytest.mark.skipif(not _HAS_TCRDIST3, reason="tcrdist3 not installed (dev only parity dependency)")
def test_parity_paired_one_to_many_matches_tcrdist3():
    tr = _tcrrep(_PAIRED, ["alpha", "beta"])
    clone = [te.PairedTCR(r.cdr3_a_aa, r.v_a_gene, r.cdr3_b_aa, r.v_b_gene)
             for r in tr.clone_df.itertuples(index=False)]
    ref = np.asarray(tr.pw_alpha, dtype=float) + np.asarray(tr.pw_beta, dtype=float)
    dists, skipped = te.tcrdist_paired_one_to_many(clone[0], clone, "human")
    assert skipped == 0
    assert dists == list(ref[0])


def test_paired_skips_when_a_chain_v_unresolvable():
    q = te.PairedTCR("CAVNFGGGKLIF", "TRAV12-1*01", "CASSIRSSYEQYF", "TRBV19*01")
    cands = [
        te.PairedTCR("CAVRDSNYQLIW", "TRAV1-2*01", "CASSLGQAYEQYF", "TRBV28*01"),
        te.PairedTCR("CAVRDSNYQLIW", "TRAVNONSENSE", "CASSLGQAYEQYF", "TRBV28*01"),  # bad alpha V
    ]
    dists, skipped = te.tcrdist_paired_one_to_many(q, cands, "human")
    assert skipped == 1
    assert dists[0] is not None and dists[1] is None
