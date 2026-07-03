from __future__ import annotations
import json
import sys
import pandas as pd

_KEEP = ["cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope_aa", "mhc_class", "mhc_a",
         "antigen", "antigen_organism", "cluster_id"]
_REQUIRED = ["cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope_aa"]

def build_index(raw_parquet_path: str, out_parquet: str, out_meta: str) -> dict:
    df = pd.read_parquet(raw_parquet_path)
    total = len(df)
    neg = int((df["label"] == 0).sum()) if "label" in df.columns else 0
    pos = df[df["label"] == 1] if "label" in df.columns else df
    kept = pos.dropna(subset=_REQUIRED)
    cols = [c for c in _KEEP if c in kept.columns]
    kept = kept[cols].drop_duplicates(
        subset=[c for c in ["cdr3_b_aa", "v_b_gene", "j_b_gene", "epitope_aa", "mhc_a"] if c in cols]
    ).reset_index(drop=True)
    kept.to_parquet(out_parquet, index=False)
    meta = {
        "rows_total": total,
        "rows_kept": int(len(kept)),
        "rows_dropped_negative": neg,
        "columns": cols,
        "chain": "beta",
        "note": "label==1 positives only; OLGA negatives excluded",
    }
    with open(out_meta, "w") as fh:
        json.dump(meta, fh, indent=2)
    return meta

if __name__ == "__main__":
    raw = sys.argv[1]
    build_index(raw, "data/unitcr_beta_index.parquet", "data/unitcr_beta_index.meta.json")
