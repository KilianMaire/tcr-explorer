"""Alpha single chain and alpha/beta paired similarity search."""
import pandas as pd
import pytest

from tcr_explorer import similarity as S
from tcr_explorer import tcrdist_engine

_needs_tcrdist = pytest.mark.skipif(
    not tcrdist_engine.tcrdist_available(), reason="pwseqdist not installed (tcrdist extra)")


@pytest.fixture
def alpha_idx(tmp_path):
    # New schema index (cdr3_aa / v_gene / chain), alpha rows only.
    df = pd.DataFrame([
        {"cdr3_aa": "CAVNFGGGKLIF", "v_gene": "TRAV12-1", "j_gene": "TRAJ23", "chain": "alpha",
         "species": "human", "epitope_aa": "NLVPMVATV"},
        {"cdr3_aa": "CAVNFGAGKLIF", "v_gene": "TRAV12-1", "j_gene": "TRAJ23", "chain": "alpha",
         "species": "human", "epitope_aa": "GILGFVFTL"},
        {"cdr3_aa": "CAVRDSNYQLIW", "v_gene": "TRAV1-2", "j_gene": "TRAJ33", "chain": "alpha",
         "species": "human", "epitope_aa": None},
    ])
    p = tmp_path / "alpha_idx.parquet"
    df.to_parquet(p, index=False)
    return str(p)


@pytest.fixture
def paired_idx(tmp_path):
    # One row per chain, alpha and beta linked by pairing_key.
    rows = []
    cells = [
        ("pk1", "TRAV12-1", "CAVNFGGGKLIF", "TRBV19", "CASSIRSSYEQYF", "NLVPMVATV"),
        ("pk2", "TRAV1-2", "CAVRDSNYQLIW", "TRBV28", "CASSLGQAYEQYF", "GILGFVFTL"),
        ("pk3", "TRAV12-1", "CAVNFGAGKLIF", "TRBV19", "CASSIRSAYEQYF", "GLCTLVAML"),
    ]
    for pk, va, ca, vb, cb, epi in cells:
        rows.append({"pairing_key": pk, "chain": "alpha", "cdr3_aa": ca, "v_gene": va,
                     "species": "human", "epitope_aa": epi, "antigen": None,
                     "antigen_organism": None, "mhc_class": "I", "mhc_a": None})
        rows.append({"pairing_key": pk, "chain": "beta", "cdr3_aa": cb, "v_gene": vb,
                     "species": "human", "epitope_aa": epi, "antigen": None,
                     "antigen_organism": None, "mhc_class": "I", "mhc_a": None})
    # a lone beta row (no alpha partner) must be ignored, not paired.
    rows.append({"pairing_key": "pk_lonely", "chain": "beta", "cdr3_aa": "CASSXXX", "v_gene": "TRBV19",
                 "species": "human", "epitope_aa": None, "antigen": None,
                 "antigen_organism": None, "mhc_class": None, "mhc_a": None})
    p = tmp_path / "paired_idx.parquet"
    pd.DataFrame(rows).to_parquet(p, index=False)
    return str(p)


def test_build_paired_index_pairs_and_falls_back_context():
    df = pd.DataFrame([
        # pk1: full pair, beta epitope present -> beta context wins
        {"pairing_key": "pk1", "chain": "alpha", "cdr3_aa": "CAA", "v_gene": "TRAV1-2",
         "species": "human", "epitope_aa": None, "mhc_a": "A_alpha"},
        {"pairing_key": "pk1", "chain": "beta", "cdr3_aa": "CBB", "v_gene": "TRBV19",
         "species": "human", "epitope_aa": "EPIB", "mhc_a": None},
        # pk2: beta epitope null -> falls back to the alpha row's epitope
        {"pairing_key": "pk2", "chain": "alpha", "cdr3_aa": "CAA2", "v_gene": "TRAV1-2",
         "species": "human", "epitope_aa": "EPIA", "mhc_a": None},
        {"pairing_key": "pk2", "chain": "beta", "cdr3_aa": "CBB2", "v_gene": "TRBV28",
         "species": "human", "epitope_aa": None, "mhc_a": None},
        # pk3: lone alpha, no beta -> excluded
        {"pairing_key": "pk3", "chain": "alpha", "cdr3_aa": "CAA3", "v_gene": "TRAV1-2",
         "species": "human", "epitope_aa": None, "mhc_a": None},
    ])
    out = S._build_paired_index(df, "human")
    assert len(out) == 2  # pk1, pk2; lone pk3 dropped
    row1 = out[out["cdr3_a_aa"] == "CAA"].iloc[0]
    assert row1["cdr3_b_aa"] == "CBB" and row1["epitope_aa"] == "EPIB" and row1["mhc_a"] == "A_alpha"
    row2 = out[out["cdr3_a_aa"] == "CAA2"].iloc[0]
    assert row2["epitope_aa"] == "EPIA"  # beta null -> alpha fallback


@_needs_tcrdist
def test_alpha_query_uses_tcrdist_through_find_similar(alpha_idx):
    # Alpha flows through the same single chain path; the engine is chain agnostic.
    neigh, engine, total, warns = S.find_similar_tcrs(
        "CAVNFGGGKLIF", "TRAV12-1", "TRAJ23", species="human", chain="alpha", index_path=alpha_idx)
    assert engine == "tcrdist"
    assert neigh[0].cdr3_b_aa == "CAVNFGGGKLIF"  # exact self match is nearest
    assert neigh[0].distance == 0.0
    assert not any(w.code == "tcrdist_unavailable" for w in warns)


@_needs_tcrdist
def test_paired_self_match_is_nearest(paired_idx):
    neigh, engine, total, warns = S.find_similar_paired_tcrs(
        "CAVNFGGGKLIF", "TRAV12-1", "CASSIRSSYEQYF", "TRBV19", species="human", index_path=paired_idx)
    assert engine == "tcrdist"
    assert total == 3  # three paired candidates; the lone beta row is ignored
    assert neigh[0].cdr3_a_aa == "CAVNFGGGKLIF"
    assert neigh[0].cdr3_b_aa == "CASSIRSSYEQYF"
    assert neigh[0].distance == 0.0
    assert neigh[0].similarity >= 0.99


@_needs_tcrdist
def test_paired_distance_equals_alpha_plus_beta(paired_idx):
    neigh, engine, total, warns = S.find_similar_paired_tcrs(
        "CAVNFGGGKLIF", "TRAV12-1", "CASSIRSSYEQYF", "TRBV19", species="human", index_path=paired_idx)
    hit = next(n for n in neigh if n.cdr3_a_aa == "CAVNFGAGKLIF")
    expected = tcrdist_engine.tcrdist_paired(
        tcrdist_engine.PairedTCR("CAVNFGGGKLIF", "TRAV12-1", "CASSIRSSYEQYF", "TRBV19"),
        tcrdist_engine.PairedTCR("CAVNFGAGKLIF", "TRAV12-1", "CASSIRSAYEQYF", "TRBV19"),
        "human")
    assert hit.distance == round(float(expected), 4)


def test_paired_without_tcrdist_returns_empty_and_warns(paired_idx, monkeypatch):
    monkeypatch.setattr(tcrdist_engine, "tcrdist_available", lambda: False)
    neigh, engine, total, warns = S.find_similar_paired_tcrs(
        "CAVNFGGGKLIF", "TRAV12-1", "CASSIRSSYEQYF", "TRBV19", species="human", index_path=paired_idx)
    assert neigh == []
    assert engine == "none"
    assert any(w.code == "tcrdist_unavailable" for w in warns)


def test_paired_unresolvable_query_v_returns_empty_and_warns(paired_idx):
    neigh, engine, total, warns = S.find_similar_paired_tcrs(
        "CAVNFGGGKLIF", "TRAVNONSENSE", "CASSIRSSYEQYF", "TRBV19", species="human", index_path=paired_idx)
    assert neigh == []
    assert any(w.code == "tcrdist_unavailable" for w in warns)


def test_paired_missing_index_is_graceful():
    neigh, engine, total, warns = S.find_similar_paired_tcrs(
        "CAVNFGGGKLIF", "TRAV12-1", "CASSIRSSYEQYF", "TRBV19", index_path="/no/such/idx.parquet")
    assert neigh == []
    assert any(w.code == "similarity_index_unavailable" for w in warns)
