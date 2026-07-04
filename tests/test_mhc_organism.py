import pandas as pd
from tcr_explorer import records as R
from tcr_explorer.records_build import SCHEMA_COLUMNS
from tcr_explorer.dossier_models import RecordsRequest


def test_mhc_organism_from_allele():
    assert R.mhc_organism("HLA-A*02:01") == "human"
    assert R.mhc_organism("HLA-DQA1*03:01") == "human"
    assert R.mhc_organism("H2-Db") == "mouse"
    assert R.mhc_organism("H-2Kb") == "mouse"
    assert R.mhc_organism(None) is None
    assert R.mhc_organism("") is None


def _idx(tmp_path, rows):
    df = pd.DataFrame(rows)
    for c in SCHEMA_COLUMNS:
        if c not in df.columns:
            df[c] = None
    df = df[SCHEMA_COLUMNS]
    p = tmp_path / "idx.parquet"
    df.to_parquet(p, index=False)
    return str(p)


def test_mouse_query_hides_human_hla_by_default(tmp_path):
    idx = _idx(tmp_path, [
        dict(source="vdjdb", source_record_id="m1", pairing_key="p1", species="mouse",
             chain="beta", cdr3_aa="CASGGTGEQYF", v_gene="TRBV13-2", j_gene="TRBJ2-7",
             mhc_a="H2-Db", external_url="u"),
        dict(source="vdjdb", source_record_id="m2", pairing_key="p2", species="mouse",
             chain="beta", cdr3_aa="CASGGTGEQYF", v_gene="TRBV13-2", j_gene="TRBJ2-7",
             mhc_a="HLA-DQA1*03:01", external_url="u"),  # HLA-transgenic mouse
    ])
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASGGTGEQYF", species="mouse"), index_path=idx)
    ids = {r.source_record_id for r in resp.exact}
    assert "m1" in ids and "m2" not in ids  # human-HLA mouse record hidden by default
    resp2 = R.retrieve_records(
        RecordsRequest(cdr3_aa="CASGGTGEQYF", species="mouse", include_cross_species_mhc=True),
        index_path=idx)
    ids2 = {r.source_record_id for r in resp2.exact}
    assert "m2" in ids2
    m2 = next(r for r in resp2.exact if r.source_record_id == "m2")
    assert m2.mhc_organism == "human" and m2.mhc_is_cross_species is True


def test_no_species_does_not_apply_mhc_filter(tmp_path):
    idx = _idx(tmp_path, [
        dict(source="vdjdb", source_record_id="m2", pairing_key="p2", species="mouse",
             chain="beta", cdr3_aa="CASGGTGEQYF", mhc_a="HLA-DQA1*03:01", external_url="u"),
    ])
    resp = R.retrieve_records(RecordsRequest(cdr3_aa="CASGGTGEQYF"), index_path=idx)
    assert len(resp.exact) == 1  # no species chosen, nothing filtered
