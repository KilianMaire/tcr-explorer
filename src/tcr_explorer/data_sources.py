"""Download each TCR dataset from its official public endpoint into a raw dir.

Verified 2026-07-04. Each function fetches, validates the payload format, and
writes the raw file under the name records_build.build_index expects. The tool
never redistributes these files; it fetches them on the user's own machine.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import httpx

_TIMEOUT = 120.0
_ZIP_MAGIC = b"PK\x03\x04"

VDJDB_LATEST_API = "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest"
IEDB_URL = "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
MCPAS_URL = "https://friedmanlab.weizmann.ac.il/McPAS-TCR/McPAS-TCR.csv"
TCR3D_COMPLEXES_URL = "https://tcr3d.ibbr.umd.edu/static/download/tcr_complexes_data.tsv"
TCR3D_CHAIN_URL = "https://tcr3d.ibbr.umd.edu/static/download/tcr_chain_data.tsv"


@dataclass
class DownloadResult:
    source: str
    ok: bool
    path: Optional[str] = None
    bytes: int = 0
    error: Optional[str] = None


def _client() -> httpx.Client:
    return httpx.Client(timeout=_TIMEOUT, follow_redirects=True,
                        headers={"User-Agent": "tcr-explorer"})


def _get(url: str) -> bytes:
    with _client() as c:
        r = c.get(url)
        r.raise_for_status()
        return r.content


def download_vdjdb(raw_dir: Path) -> DownloadResult:
    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        with _client() as c:
            meta = c.get(VDJDB_LATEST_API).json()
        asset = next(a for a in meta.get("assets", [])
                     if a.get("name", "").startswith("vdjdb-") and a["name"].endswith(".zip"))
        data = _get(asset["browser_download_url"])
        if not data.startswith(_ZIP_MAGIC):
            return DownloadResult("vdjdb", False, error="payload is not a zip")
        out = raw_dir / asset["name"]
        out.write_bytes(data)
        return DownloadResult("vdjdb", True, str(out), len(data))
    except (httpx.HTTPError, StopIteration, KeyError, ValueError) as exc:
        return DownloadResult("vdjdb", False, error=str(exc))


def download_iedb(raw_dir: Path) -> DownloadResult:
    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        data = _get(IEDB_URL)
        if not data.startswith(_ZIP_MAGIC):
            return DownloadResult("iedb", False, error="payload is not a zip")
        out = raw_dir / "iedb_receptor_full_v3.zip"
        out.write_bytes(data)
        return DownloadResult("iedb", True, str(out), len(data))
    except httpx.HTTPError as exc:
        return DownloadResult("iedb", False, error=str(exc))


def download_mcpas(raw_dir: Path) -> DownloadResult:
    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        text = _get(MCPAS_URL).decode("utf-8-sig")  # utf-8-sig strips the BOM
        if "CDR3.alpha.aa" not in text.splitlines()[0]:
            return DownloadResult("mcpas", False, error="unexpected CSV header")
        out = raw_dir / "mcpas.csv"
        out.write_text(text, encoding="utf-8")
        return DownloadResult("mcpas", True, str(out), len(text))
    except (httpx.HTTPError, IndexError, UnicodeDecodeError) as exc:
        return DownloadResult("mcpas", False, error=str(exc))


def download_tcr3d(raw_dir: Path) -> DownloadResult:
    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        complexes = _get(TCR3D_COMPLEXES_URL)
        chain = _get(TCR3D_CHAIN_URL)
        if b"PDB_ID" not in complexes.splitlines()[0] or b"pdb_id" not in chain.splitlines()[0]:
            return DownloadResult("tcr3d", False, error="unexpected TSV header")
        (raw_dir / "tcr3d_complexes.tsv").write_bytes(complexes)
        (raw_dir / "tcr3d_chain.tsv").write_bytes(chain)
        return DownloadResult("tcr3d", True, str(raw_dir / "tcr3d_complexes.tsv"),
                              len(complexes) + len(chain))
    except (httpx.HTTPError, IndexError) as exc:
        return DownloadResult("tcr3d", False, error=str(exc))
