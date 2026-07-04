"""First-run data bootstrap: download every source, build the index, fetch germline.

Downloading is only ever done here, via the explicit `tcr-explorer-refresh`
command. Queries never trigger a download; they call ensure_ready() and get a
clear message when the data is absent.
"""
from __future__ import annotations

import json
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
    # Invoke stitchr's downloader through the current interpreter so it resolves
    # regardless of PATH (a venv console script is not on PATH unless activated).
    out = {}
    for species in ("human", "mouse"):
        try:
            subprocess.run([sys.executable, "-m", "Stitchr.stitchrdl", "-s", species],
                           check=True, capture_output=True, text=True)
            out[species] = True
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            out[species] = f"failed: {exc}"
    return out


def _augment_meta_with_status(meta_file, status: dict) -> None:
    """Record which sources actually downloaded, so an incomplete index is inspectable."""
    try:
        m = json.loads(meta_file.read_text())
        m["download_status"] = status
        meta_file.write_text(json.dumps(m, indent=2))
    except (OSError, ValueError):
        pass


def refresh(force: bool = False) -> dict:
    raw = data_paths.raw_dir()
    results = [
        data_sources.download_vdjdb(raw),
        data_sources.download_iedb(raw),
        data_sources.download_mcpas(raw),
        data_sources.download_tcr3d(raw),
    ]
    status = {r.source: {"ok": r.ok, "bytes": r.bytes, "error": r.error} for r in results}
    n_ok = sum(r.ok for r in results)
    all_ok = n_ok == len(results)
    had_index = data_paths.data_present()
    data_dir = str(data_paths.data_dir())

    # Never overwrite a good index with an empty or strictly-worse one.
    if n_ok == 0:
        return {"built": False, "reason": "all downloads failed",
                "sources": status, "rows_total": 0, "data_dir": data_dir}
    if not all_ok and had_index:
        return {"built": False, "reason": "partial download; kept the existing index",
                "sources": status, "rows_total": 0, "data_dir": data_dir}

    meta = build_index(str(raw), str(data_paths.records_index_path()),
                       str(data_paths.meta_path()))
    _augment_meta_with_status(data_paths.meta_path(), status)
    germline = _run_stitchrdl()
    return {"built": True, "complete": all_ok, "sources": status,
            "rows_total": meta.get("rows_total", 0), "per_source": meta.get("per_source", {}),
            "germline": germline, "data_dir": data_dir}


def main() -> None:
    print("Downloading TCR datasets and building the index...", file=sys.stderr)
    summary = refresh()
    for src, r in summary["sources"].items():
        status = "ok" if r["ok"] else f"FAILED ({r['error']})"
        print(f"  {src}: {status}", file=sys.stderr)
    if not summary.get("built"):
        print(f"Refresh aborted: {summary.get('reason')}. Existing data left unchanged.",
              file=sys.stderr)
        raise SystemExit(1)
    print(f"  germline: {summary['germline']}", file=sys.stderr)
    print(f"Built {summary['rows_total']} records into {summary['data_dir']}", file=sys.stderr)
    if not summary.get("complete", True):
        print("Warning: some sources failed; index built from the rest. Re-run to complete.",
              file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
