from __future__ import annotations

import csv
import os
import sys
from pathlib import Path
from typing import Any

import httpx
from fastapi import FastAPI

# Make cdr_enricher importable when running the server directly
_src = Path(__file__).resolve().parent.parent / "src"
if str(_src) not in sys.path:
    sys.path.insert(0, str(_src))

try:
    from tcr_explorer.cdr_enricher import get_cdr1_cdr2
    _CDR_AVAILABLE = True
except Exception:
    _CDR_AVAILABLE = False

app = FastAPI(title="VDJdb Tool Server", version="0.1.0")

# Local seed / override file — set VDJDB_DATA_FILE env var to point at a full TSV download.
_DEFAULT_SEED = Path(__file__).resolve().parent.parent / "examples" / "vdjdb_seed.csv"
DATA_FILE = Path(os.getenv("VDJDB_DATA_FILE", str(_DEFAULT_SEED)))

# Optional: upstream VDJdb REST API.  Set VDJDB_API_URL="" to disable.
# POST endpoint: https://vdjdb.cdr3.net/api/database/search
# Meta endpoint: https://vdjdb.cdr3.net/api/database/meta
_VDJDB_BASE = os.getenv("VDJDB_API_BASE", "https://vdjdb.cdr3.net")
VDJDB_API_URL = os.getenv("VDJDB_API_URL", f"{_VDJDB_BASE}/api/database/search")
_VDJDB_META_URL = f"{_VDJDB_BASE}/api/database/meta"

# Column order cache populated at first API call
_VDJDB_COLUMNS: list[str] = []

# ---------------------------------------------------------------------------
# Load local seed data at startup
# ---------------------------------------------------------------------------
RECORDS: list[dict[str, str]] = []
if DATA_FILE.exists():
    with DATA_FILE.open() as f:
        reader = csv.DictReader(f)
        RECORDS = [dict(r) for r in reader]


# ---------------------------------------------------------------------------
# Species name normalisation
# ---------------------------------------------------------------------------
_SPECIES_TO_IMGT: dict[str, str] = {
    "homosapiens": "human",
    "human": "human",
    "mussculmus": "mouse",
    "musmusculus": "mouse",
    "mouse": "mouse",
}

_IMGT_TO_VDJDB: dict[str, str] = {
    "human": "HomoSapiens",
    "mouse": "MusMusculus",
}


def _normalise_species(raw: str) -> str:
    return _SPECIES_TO_IMGT.get(raw.lower().replace(" ", ""), "other")


# ---------------------------------------------------------------------------
# Record builder — maps a VDJdb row to the shared GeneRecord schema
# ---------------------------------------------------------------------------
def _row_to_record(row: dict[str, Any]) -> dict[str, Any]:
    # Real API uses "v.segm" or "v.beta"/"v.alpha" depending on DB version
    v_gene = (
        row.get("v_segm") or row.get("v.segm")
        or row.get("v.beta") or row.get("v.alpha")
        or "unknown"
    )
    species = _normalise_species(row.get("species", ""))

    # Predict CDR1/CDR2 from V gene via stitchr IMGT data
    cdr_info: dict[str, Any] = {}
    if _CDR_AVAILABLE and v_gene != "unknown":
        cdr_info = get_cdr1_cdr2(v_gene, species)

    return {
        "source": "vdjdb",
        "species": species,
        "gene_name": v_gene,
        "allele_name": None,
        "region": "cdr3",
        "sequence": row.get("cdr3", ""),
        "metadata": {
            "j_segm": (
                row.get("j_segm") or row.get("j.segm")
                or row.get("j.beta") or row.get("j.alpha")
            ),
            "antigen_epitope": row.get("antigen_epitope") or row.get("antigen.epitope"),
            "antigen_gene": row.get("antigen_gene") or row.get("antigen.gene"),
            "antigen_species": row.get("antigen_species") or row.get("antigen.species"),
            "mhc_a": row.get("mhc_a") or row.get("mhc.a"),
            "mhc_b": row.get("mhc_b") or row.get("mhc.b"),
            "mhc_class": row.get("mhc_class") or row.get("mhc.class"),
            "score": row.get("score") or row.get("vdjdb.score"),
            # CDR1/CDR2 predicted from V gene germline (IMGT via stitchr)
            "cdr1_aa": cdr_info.get("cdr1_aa"),
            "cdr2_aa": cdr_info.get("cdr2_aa"),
            "backend": "vdjdb-server",
        },
    }


# ---------------------------------------------------------------------------
# Search helpers
# ---------------------------------------------------------------------------
def _search_local(
    gene: str,
    species: str,
    cdr3_contains: str,
    antigen_epitope: str,
    limit: int,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for r in RECORDS:
        v = r.get("v_segm", "")
        sp = _normalise_species(r.get("species", ""))
        cdr3 = r.get("cdr3", "").upper()
        epitope = (r.get("antigen_epitope") or "").upper()

        if gene and gene.upper() not in v.upper():
            continue
        if species and species != sp:
            continue
        if cdr3_contains and cdr3_contains.upper() not in cdr3:
            continue
        if antigen_epitope and antigen_epitope.upper() not in epitope:
            continue
        out.append(_row_to_record(r))
    return out[:limit]


async def _fetch_columns() -> list[str]:
    """Fetch column order from VDJdb metadata endpoint and cache globally."""
    global _VDJDB_COLUMNS
    if _VDJDB_COLUMNS:
        return _VDJDB_COLUMNS
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            r = await client.get(_VDJDB_META_URL)
            r.raise_for_status()
            data = r.json()
            _VDJDB_COLUMNS = [col["name"] for col in data.get("columns", [])]
    except Exception:
        pass
    return _VDJDB_COLUMNS


def _entries_to_row(entries: list[Any], columns: list[str]) -> dict[str, Any]:
    """Map a positional VDJdb entries array to a named dict."""
    return {col: val for col, val in zip(columns, entries)}


async def _search_vdjdb_api(
    gene: str,
    species: str,
    cdr3_contains: str,
    antigen_epitope: str,
    limit: int,
) -> list[dict[str, Any]]:
    if not VDJDB_API_URL:
        return []

    columns = await _fetch_columns()

    filters: list[dict[str, Any]] = []
    if gene:
        filters.append({"column": "v.segm", "filterType": "substring:set", "value": gene})
    if species and species in _IMGT_TO_VDJDB:
        filters.append({"column": "species", "filterType": "exact", "value": _IMGT_TO_VDJDB[species]})
    if cdr3_contains:
        filters.append({"column": "cdr3", "filterType": "substring:set", "value": cdr3_contains})
    if antigen_epitope:
        filters.append({"column": "antigen.epitope", "filterType": "substring:set", "value": antigen_epitope})

    payload: dict[str, Any] = {
        "filters": filters,
        "pageSize": min(limit, 500),
        "page": 0,
    }

    try:
        async with httpx.AsyncClient(timeout=12.0) as client:
            r = await client.post(VDJDB_API_URL, json=payload)
            r.raise_for_status()
            data = r.json()
    except Exception:
        return []

    api_rows = data.get("rows", [])
    if columns:
        rows = [_entries_to_row(row.get("entries", []), columns) for row in api_rows]
    else:
        # Fallback if metadata fetch failed — entries are unnamed, skip
        return []

    return [_row_to_record(row) for row in rows][:limit]


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "vdjdb"}


@app.post("/search")
async def search(req: dict) -> dict:
    limit = int(req.get("limit", 50))
    gene = (req.get("gene_name") or "").strip()
    species = (req.get("species") or "").strip().lower()
    cdr3_contains = (req.get("sequence_contains") or "").strip()
    antigen_epitope = (req.get("antigen_epitope") or "").strip()

    # Local seed first
    local = _search_local(gene, species, cdr3_contains, antigen_epitope, limit)

    # Merge with upstream API results (deduplicate by cdr3+v_segm)
    remote = await _search_vdjdb_api(gene, species, cdr3_contains, antigen_epitope, limit)

    seen: set[str] = {(r["sequence"] + "|" + r["gene_name"]) for r in local}
    for rec in remote:
        key = rec["sequence"] + "|" + rec["gene_name"]
        if key not in seen:
            local.append(rec)
            seen.add(key)

    return {"total": len(local), "records": local[:limit]}
