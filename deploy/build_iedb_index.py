"""Build the IEDB-only records index for the public deploy.

The public instance serves IEDB only (CC BY 4.0). Datacenter build runners (such
as Hugging Face Spaces) time out downloading from iedb.org, so the index is
built once here and the resulting parquet is baked into the deploy image (copied
into the Space repo next to the Dockerfile). See deploy/DEPLOY.md.

Usage:
    python deploy/build_iedb_index.py OUT_DIR
writes OUT_DIR/records_index.parquet and OUT_DIR/records_index.meta.json .
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from tcr_explorer.data_sources import download_iedb
from tcr_explorer.records_build import build_index


def main() -> None:
    out_dir = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = Path(tempfile.mkdtemp(prefix="iedb_raw_"))

    result = download_iedb(raw)
    if not result.ok:
        raise SystemExit(f"IEDB download failed: {result.error}")

    meta = build_index(
        str(raw),
        str(out_dir / "records_index.parquet"),
        str(out_dir / "records_index.meta.json"),
    )
    print(f"IEDB index built: {meta['rows_total']} records -> {out_dir}")


if __name__ == "__main__":
    main()
