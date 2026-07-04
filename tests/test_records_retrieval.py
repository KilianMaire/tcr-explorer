import pandas as pd
import pytest
from tcr_explorer import records as R
from tcr_explorer.records_build import SCHEMA_COLUMNS
from tcr_explorer.dossier_models import RecordsRequest


@pytest.fixture
def idx(tmp_path):
    rows = [
        dict(source="vdjdb", source_record_id="vdjdb:1", pairing_key="vdjdb:c1", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV19", j_gene="TRBJ2-7", epitope_aa="NLVPMVATV",
             mhc_a="HLA-A*02:01", external_url="u"),
        dict(source="mcpas", source_record_id="mcpas:2", pairing_key="mcpas:2", species="human",
             chain="beta", cdr3_aa="CASSLGTEAFF", v_gene="TRBV19", j_gene="TRBJ2-7", cdr3_nt="TGT",
             epitope_aa="NLVPMVATV", external_url="u"),  # same cdr3 as vdjdb: concordance 2
        dict(source="vdjdb", source_record_id="vdjdb:3", pairing_key="vdjdb:c1", species="human",
             chain="alpha", cdr3_aa="CAVRDSNYQLIW", v_gene="TRAV12-1", j_gene="TRAJ33",
             epitope_aa="NLVPMVATV", external_url="u"),  # pairs with vdjdb:1 via vdjdb:c1
        dict(source="vdjdb", source_record_id="vdjdb:4", pairing_key="vdjdb:u4", species="mouse",
             chain="beta", cdr3_aa="CASSPGQGAETLYF", v_gene="TRBV13-1", j_gene="TRBJ2-3",
             epitope_aa="SIINFEKL", external_url="u"),
        dict(source="vdjdb", source_record_id="vdjdb:5", pairing_key="vdjdb:u5", species="human",
             chain="beta", cdr3_aa="CASSLGTEAYF", v_gene="TRBV19", j_gene="TRBJ2-7",
             epitope_aa="GILGFVFTL", external_url="u"),  # 1 substitution from CASSLGTEAFF: a neighbour
        dict(source="vdjdb", source_record_id="vdjdb:6", pairing_key="vdjdb:u6", species="mouse",
             chain="beta", cdr3_aa="CASSPWGGAETLYF", v_gene="TRBV13-2*01", j_gene="TRBJ2-3",
             epitope_aa="SIINFEKL", external_url="u"),  # mouse row, allele-suffixed raw v_gene
    ]
    df = pd.DataFrame(rows)
    for c in SCHEMA_COLUMNS:
        if c not in df.columns:
            df[c] = None
    df = df[SCHEMA_COLUMNS]
    p = tmp_path / "idx.parquet"
    df.to_parquet(p, index=False)
    return str(p)


def test_bare_cdr3_returns_exact_records_across_sources(idx):
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASSLGTEAFF"), index_path=idx)
    assert resp.total_exact == 2  # vdjdb + mcpas
    assert {r.source for r in resp.exact} == {"vdjdb", "mcpas"}
    assert all(r.match_kind == "exact" for r in resp.exact)
    assert all(r.concordance == 2 for r in resp.exact)  # 2 sources hold this cdr3


def test_neighbours_are_separate_and_never_in_exact(idx):
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASSLGTEAFF"), index_path=idx)
    exact_ids = {r.source_record_id for r in resp.exact}
    assert "vdjdb:5" not in exact_ids  # the 1-substitution variant is not exact
    assert any(r.source_record_id == "vdjdb:5" for r in resp.neighbours)
    assert all(r.match_kind == "neighbour" for r in resp.neighbours)


def test_species_filter(idx):
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASSPGQGAETLYF", species="human"), index_path=idx)
    assert resp.total_exact == 0  # that cdr3 is mouse only
    resp2 = R.retrieve_records(RecordsRequest(cdr3_aa="CASSPGQGAETLYF", species="mouse"), index_path=idx)
    assert resp2.total_exact == 1


def test_gene_query_returns_all_records_for_gene(idx):
    resp = R.retrieve_records(RecordsRequest(v_gene="TRBV19"), index_path=idx)
    assert {r.source_record_id for r in resp.exact} >= {"vdjdb:1", "mcpas:2", "vdjdb:5"}


def test_mouse_gene_query_matches_allele_suffixed_index_value(idx):
    # Task 1 binding note: mouse genes in the index keep their raw, allele-suffixed
    # source strings (e.g. "TRBV13-2*01"). A query on the bare gene base must still
    # match by comparing normalized gene BASES on both sides, not exact strings.
    resp = R.retrieve_records(RecordsRequest(v_gene="TRBV13-2", species="mouse"), index_path=idx)
    assert {r.source_record_id for r in resp.exact} == {"vdjdb:6"}


def test_pair_query_groups_alpha_beta(idx):
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASSLGTEAFF", cdr3_aa_b="CAVRDSNYQLIW"), index_path=idx)
    assert len(resp.pairs) >= 1
    pair = resp.pairs[0]
    assert pair.alpha is not None and pair.beta is not None
    assert pair.alpha.chain == "alpha" and pair.beta.chain == "beta"


def test_id_query_resolves_one_record(idx):
    resp = R.retrieve_records(RecordsRequest(query="mcpas:2"), index_path=idx)
    assert resp.total_exact == 1 and resp.exact[0].source_record_id == "mcpas:2"


def test_free_query_cdr3_is_detected(idx):
    resp = R.retrieve_records(RecordsRequest(query="CASSLGTEAFF"), index_path=idx)
    assert resp.total_exact == 2
