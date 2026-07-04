"""First-run data bootstrap: download every source, build the index, fetch germline.

Downloading is only ever done here, via the explicit `tcr-explorer-refresh`
command. Queries never trigger a download; they call ensure_ready() and get a
clear message when the data is absent.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

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


def _germline_present(species: str) -> bool:
    """True if the stitchr Data dir the enricher reads at query time actually
    holds this species' germline FASTA. This is the honest post-condition:
    stitchrdl exits 0 even when it downloaded nothing, so we verify the files
    the runtime will look for, not the child's exit code."""
    from .cdr_enricher import _stitchr_data_dir

    d = _stitchr_data_dir()
    if d is None:
        return False
    return next((d / species.upper()).glob("*.fasta"), None) is not None


def _stitchrdl_cmd() -> list[str]:
    # stitchrdl.py ships NO `if __name__ == "__main__"` guard, so
    # `python -m Stitchr.stitchrdl` imports the module and exits WITHOUT running
    # main() -- a silent no-op that never downloads. Its declared entry point is
    # the `stitchrdl` console script (Stitchr.stitchrdl:main), a sibling of the
    # current interpreter. Prefer it; fall back to invoking main() explicitly so
    # it still works if the console script is missing.
    stitchrdl = Path(sys.executable).parent / "stitchrdl"
    if stitchrdl.exists():
        return [str(stitchrdl)]
    return [sys.executable, "-c", "from Stitchr.stitchrdl import main; main()"]


# stitchrdl scrapes IMGT/GENE-DB through IMGTgeneDL, which can hang for a long
# time when IMGT is slow or unreachable. Bound it so `tcr-explorer-refresh`
# cannot freeze forever; overridable for slow links.
_STITCHRDL_TIMEOUT = int(os.environ.get("TCR_EXPLORER_STITCHRDL_TIMEOUT", "600"))


def _run_stitchrdl() -> dict:
    # stitchrdl shells out to `IMGTgeneDL` via subprocess with shell=True, which
    # resolves against PATH, so put the interpreter's bin dir (where the venv's
    # IMGTgeneDL lives) first. Without this the child's IMGTgeneDL call silently
    # fails when the venv is not activated.
    out = {}
    cmd = _stitchrdl_cmd()
    bindir = str(Path(sys.executable).parent)
    env = dict(os.environ)
    env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
    for species in ("human", "mouse"):
        try:
            subprocess.run(cmd + ["-s", species], check=True,
                           capture_output=True, text=True, env=env,
                           timeout=_STITCHRDL_TIMEOUT)
            # Do NOT trust the exit code: verify the FASTA actually landed.
            out[species] = _germline_present(species)
        except subprocess.TimeoutExpired:
            out[species] = f"failed: IMGT germline download timed out after {_STITCHRDL_TIMEOUT}s"
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
