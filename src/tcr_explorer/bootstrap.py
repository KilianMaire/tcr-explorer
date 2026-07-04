"""First-run data bootstrap: download every source, build the index, fetch germline.

Downloading is only ever done here, via the explicit `tcr-explorer-refresh`
command. Queries never trigger a download; they call ensure_ready() and get a
clear message when the data is absent.
"""
from __future__ import annotations

import subprocess
import sys

from . import data_paths, data_sources
from .records_build import build_index


class DataNotReadyError(RuntimeError):
    pass


def ensure_ready() -> None:
    if not data_paths.data_present():
        raise DataNotReadyError(
            "TCR Explorer data is not downloaded yet. Run `tcr-explorer-refresh` "
            "once to download the datasets, then retry."
        )


def _run_stitchrdl() -> dict:
    out = {}
    for species in ("human", "mouse"):
        try:
            subprocess.run(["stitchrdl", "-s", species], check=True,
                           capture_output=True, text=True)
            out[species] = True
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            out[species] = f"failed: {exc}"
    return out


def refresh(force: bool = False) -> dict:
    raw = data_paths.raw_dir()
    results = {
        r.source: {"ok": r.ok, "bytes": r.bytes, "error": r.error}
        for r in (
            data_sources.download_vdjdb(raw),
            data_sources.download_iedb(raw),
            data_sources.download_mcpas(raw),
            data_sources.download_tcr3d(raw),
        )
    }
    meta = build_index(str(raw), str(data_paths.records_index_path()),
                       str(data_paths.meta_path()))
    germline = _run_stitchrdl()
    return {"sources": results, "rows_total": meta.get("rows_total", 0),
            "per_source": meta.get("per_source", {}), "germline": germline,
            "data_dir": str(data_paths.data_dir())}


def main() -> None:
    print("Downloading TCR datasets and building the index...", file=sys.stderr)
    summary = refresh()
    for src, r in summary["sources"].items():
        status = "ok" if r["ok"] else f"FAILED ({r['error']})"
        print(f"  {src}: {status}", file=sys.stderr)
    print(f"  germline: {summary['germline']}", file=sys.stderr)
    print(f"Built {summary['rows_total']} records into {summary['data_dir']}", file=sys.stderr)


if __name__ == "__main__":
    main()
