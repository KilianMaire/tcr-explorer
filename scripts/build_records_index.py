from __future__ import annotations
import sys
from tcr_explorer.records_build import build_index

if __name__ == "__main__":
    raw = sys.argv[1] if len(sys.argv) > 1 else "data/raw"
    meta = build_index(raw, "data/records_index.parquet", "data/records_index.meta.json")
    print(meta)
