import pandas as pd
from pathlib import Path
from scripts.build_unitcr_index import build_index

def _make_tiny(path):
    rows = [
        # a label==1 beta positive (kept)
        dict(cdr3_b_aa="CASSLGTEAFF", v_b_gene="TRBV20-1", j_b_gene="TRBJ1-1", epitope_aa="NLVPMVATV",
             mhc_class="I", mhc_a="HLA-A*02:01", antigen="pp65", antigen_organism="HCMV", cluster_id=1, label=1),
        # a label==0 OLGA negative (MUST be dropped)
        dict(cdr3_b_aa="CASSNEGATIVE", v_b_gene="TRBV20-1", j_b_gene="TRBJ1-1", epitope_aa="NLVPMVATV",
             mhc_class="I", mhc_a="HLA-A*02:01", antigen="pp65", antigen_organism="HCMV", cluster_id=2, label=0),
        # a label==1 but missing cdr3 (dropped)
        dict(cdr3_b_aa=None, v_b_gene="TRBV20-1", j_b_gene="TRBJ1-1", epitope_aa="NLVPMVATV",
             mhc_class="I", mhc_a="HLA-A*02:01", antigen="x", antigen_organism="y", cluster_id=3, label=1),
    ]
    pd.DataFrame(rows).to_parquet(path)

def test_build_drops_negatives_and_incomplete(tmp_path):
    raw = tmp_path / "raw.parquet"; _make_tiny(raw)
    out = tmp_path / "idx.parquet"; meta = tmp_path / "idx.meta.json"
    m = build_index(str(raw), str(out), str(meta))
    df = pd.read_parquet(out)
    assert len(df) == 1
    assert df.iloc[0]["cdr3_b_aa"] == "CASSLGTEAFF"
    assert (df["cdr3_b_aa"] == "CASSNEGATIVE").sum() == 0   # negative excluded
    assert m["rows_kept"] == 1 and m["rows_dropped_negative"] == 1
    assert set(["cdr3_b_aa","v_b_gene","j_b_gene","epitope_aa","mhc_a","antigen"]).issubset(df.columns)
