import pandas as pd
import pytest

from tcr_explorer import similarity, tcrdist_engine
from tcr_explorer.similarity import find_similar_tcrs, cdr3_distance

FIX = "tests/fixtures/unitcr_tiny.parquet"

# The tcrdist engine needs the optional `pwseqdist` extra; guard tests that
# require it so a base install (without `tcr-explorer[tcrdist]`) skips rather
# than fails. The BLOSUM fallback tests below run either way.
_needs_tcrdist = pytest.mark.skipif(
    not tcrdist_engine.tcrdist_available(), reason="pwseqdist not installed (tcrdist extra)")


def test_self_match_is_nearest_blosum_when_tcrdist_unavailable(monkeypatch):
    # With the optional tcrdist engine absent, scoring falls back to the bundled
    # BLOSUM CDR3 distance and carries the honest downgrade warning.
    monkeypatch.setattr(tcrdist_engine, "tcrdist_available", lambda: False)
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", top_k=3, index_path=FIX)
    assert engine == "blosum_cdr3"
    assert neigh[0].cdr3_b_aa == "CASSLGTEAFF"
    assert neigh[0].similarity >= 0.99
    assert any(w.code == "tcrdist_unavailable" for w in warns)


@_needs_tcrdist
def test_self_match_is_nearest_tcrdist_when_available():
    # With pwseqdist installed (the tcrdist extra), scoring uses the authoritative
    # tcrdist metric: the engine is labelled "tcrdist", the exact self match sits at
    # distance 0, and the downgrade warning must NOT appear.
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", top_k=3, index_path=FIX)
    assert engine == "tcrdist"
    assert neigh[0].cdr3_b_aa == "CASSLGTEAFF"
    assert neigh[0].distance == 0.0
    assert not any(w.code == "tcrdist_unavailable" for w in warns)


@_needs_tcrdist
def test_query_v_unresolvable_falls_back_to_blosum():
    # A query V gene absent from the tcrdist reference cannot yield germline loops,
    # so scoring falls back to BLOSUM and says so, rather than silently mis-scoring.
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBVNONSENSE", "TRBJ1-1", top_k=3, index_path=FIX)
    assert engine == "blosum_cdr3"
    assert any(w.code == "tcrdist_unavailable" for w in warns)


def test_non_human_species_returns_empty_with_warning():
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", species="mouse", index_path=FIX)
    assert neigh == []
    assert engine == "none"
    assert ncand == 0
    assert any(w.code == "species_unsupported" for w in warns)


def test_distance_zero_for_identical():
    assert cdr3_distance("CASSLGTEAFF", "CASSLGTEAFF") == 0.0


def test_missing_index_is_graceful():
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", index_path="/no/such/index.parquet")
    assert neigh == []
    assert any(w.code == "similarity_index_unavailable" for w in warns)


def test_nan_optional_fields_become_none(tmp_path):
    # Real index rows carry NaN (numpy float) in unpopulated optional string columns;
    # pydantic rejects NaN for Optional[str]. Build a NaN-bearing index and query it.
    idx = tmp_path / "nan_index.parquet"
    pd.DataFrame(
        {
            "cdr3_b_aa": ["CASSLGTEAFF"],
            "v_b_gene": ["TRBV20-1"],
            "j_b_gene": ["TRBJ1-1"],
            "epitope_aa": [None],
            "mhc_class": [float("nan")],
            "mhc_a": [float("nan")],
            "antigen": [float("nan")],
            "antigen_organism": [None],
            "cluster_id": [float("nan")],
        }
    ).to_parquet(idx)
    neigh, engine, ncand, warns = find_similar_tcrs(
        "CASSLGTEAFF", "TRBV20-1", "TRBJ1-1", index_path=str(idx))
    assert neigh[0].antigen is None
    assert neigh[0].mhc_a is None
    assert neigh[0].antigen_organism is None
    assert neigh[0].cluster_id is None
