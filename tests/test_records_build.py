from pathlib import Path
import pandas as pd
from imgt_app import records_build as rb

FIX = Path(__file__).parent / "fixtures" / "records"


def test_schema_columns_stable():
    assert rb.SCHEMA_COLUMNS[:5] == ["source", "source_record_id", "pairing_key", "species", "chain"]
    assert "cdr3_aa" in rb.SCHEMA_COLUMNS and "nt" not in rb.SCHEMA_COLUMNS[0]


def test_vdjdb_harmonizes_to_schema():
    df = rb.harmonize_vdjdb(str(FIX / "vdjdb_slim.txt"))
    assert list(df.columns) == rb.SCHEMA_COLUMNS
    assert (df["source"] == "vdjdb").all()
    assert df["chain"].isin(["alpha", "beta"]).all()
    # the fixture must exercise both chain branches, not just alpha
    assert set(df["chain"].unique()) == {"alpha", "beta"}
    # cdr3 is uppercased and stripped
    assert df["cdr3_aa"].str.strip().eq(df["cdr3_aa"]).all()
    assert df["cdr3_aa"].str.upper().eq(df["cdr3_aa"]).all()


def test_iedb_explodes_two_chains_sharing_pairing_key():
    df = rb.harmonize_iedb(str(FIX / "iedb_tcr.csv"))
    assert list(df.columns) == rb.SCHEMA_COLUMNS
    assert (df["source"] == "iedb").all()
    # a receptor with both chains present yields two rows with one pairing_key
    paired = df.groupby("pairing_key")["chain"].apply(set)
    assert any(s == {"alpha", "beta"} for s in paired), "expected at least one paired receptor"
    # IEDB carries deposited protein sequence for at least some rows
    assert df["full_aa"].notna().any()
    # the Assay-group "MHC Allele Names" column must be reachable, not just
    # the (unrelated) Receptor/Reference/Epitope columns
    assert df["mhc_a"].notna().any()
    # external_url must resolve to the receptor, never to the epitope IRI
    # that happens to share the same second-level column name
    assert not df["external_url"].str.contains("/epitope/").any()


def test_mcpas_deposited_cdr3_nt_is_kept():
    df = rb.harmonize_mcpas(str(FIX / "mcpas.csv"))
    assert (df["source"] == "mcpas").all()
    # McPAS supplies deposited CDR3 nt on at least some rows
    assert df["cdr3_nt"].notna().any()


def test_tcr3d_carries_pdb_and_cdr1_cdr2():
    df = rb.harmonize_tcr3d(str(FIX / "tcr3d_complexes.tsv"), str(FIX / "tcr3d_chain.tsv"))
    assert (df["source"] == "tcr3d").all()
    assert df["pdb_id"].notna().all()
    assert df["cdr1_aa"].notna().any() and df["cdr2_aa"].notna().any()
    assert df["external_url"].str.contains("rcsb.org").any()


def test_species_normalization():
    assert rb.normalize_species("HomoSapiens") == "human"
    assert rb.normalize_species("Homo Sapiens") == "human"
    assert rb.normalize_species("MusMusculus") == "mouse"
    assert rb.normalize_species("Mus musculus") == "mouse"
    assert rb.normalize_species("GallusGallus") == "gallusgallus"  # unmapped, lowercased, not guessed


def test_build_index_concatenates_dedups_and_writes_meta(tmp_path):
    out = tmp_path / "idx.parquet"
    meta = tmp_path / "idx.meta.json"
    result = rb.build_index(str(FIX), str(out), str(meta))
    assert out.exists() and meta.exists()
    df = pd.read_parquet(out)
    assert list(df.columns) == rb.SCHEMA_COLUMNS
    assert set(df["source"].unique()) <= {"vdjdb", "iedb", "mcpas", "tcr3d"}
    assert result["rows_total"] == len(df)
    assert set(result["per_source"]) <= {"vdjdb", "iedb", "mcpas", "tcr3d"}
    # dedup: no exact duplicate on the dedup key within a source
    key = ["source", "chain", "cdr3_aa", "v_gene", "j_gene", "epitope_aa", "mhc_a"]
    assert not df.duplicated(subset=key).any()
