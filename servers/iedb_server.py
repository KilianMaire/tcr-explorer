from __future__ import annotations

import os
from typing import Any

import httpx
from fastapi import FastAPI

app = FastAPI(title="IEDB Tool Server", version="0.1.0")

# PostgREST-based IEDB Query API — T cell assay endpoint
IEDB_API_URL = os.getenv("IEDB_API_URL", "https://query-api.iedb.org/tcell_search")

# Fields we request from IEDB to keep responses lean
_SELECT = ",".join([
    "linear_sequence",
    "mhc_allele_name",
    "mhc_class",
    "source_organism_name",
    "parent_source_antigen",
    "assay_type",
    "effector_cell_type",
    "assay_qualitative_measure",
])

# ---------------------------------------------------------------------------
# Species normalisation
# ---------------------------------------------------------------------------
_SPECIES_KEYWORDS: list[tuple[str, str]] = [
    ("homo sapiens", "human"),
    ("human", "human"),
    ("mus musculus", "mouse"),
    ("mouse", "mouse"),
]

_IMGT_TO_IEDB: dict[str, str] = {
    "human": "homo sapiens",
    "mouse": "mus musculus",
}


def _normalise_species(raw: str) -> str:
    r = raw.lower()
    for keyword, norm in _SPECIES_KEYWORDS:
        if keyword in r:
            return norm
    return "other"


# ---------------------------------------------------------------------------
# Record builder — maps an IEDB tcell row to the shared GeneRecord schema
# ---------------------------------------------------------------------------
def _row_to_record(row: dict[str, Any]) -> dict[str, Any]:
    org = row.get("source_organism_name") or ""
    return {
        "source": "iedb",
        "species": _normalise_species(org),
        # MHC allele is the closest analogue to gene_name for IEDB records
        "gene_name": row.get("mhc_allele_name") or "unknown",
        "allele_name": None,
        "region": "epitope",
        "sequence": row.get("linear_sequence") or "",
        "metadata": {
            "mhc_class": row.get("mhc_class"),
            "antigen_name": row.get("parent_source_antigen"),
            "source_organism": org,
            "assay_type": row.get("assay_type"),
            "effector_cell_type": row.get("effector_cell_type"),
            "qualitative_measure": row.get("assay_qualitative_measure"),
            "backend": "iedb-server",
        },
    }


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "iedb"}


@app.post("/search")
async def search(req: dict) -> dict:
    limit = int(req.get("limit", 50))
    gene = (req.get("gene_name") or "").strip()
    species = (req.get("species") or "").strip().lower()
    antigen_epitope = (req.get("antigen_epitope") or "").strip()
    sequence_contains = (req.get("sequence_contains") or "").strip()

    # Build PostgREST query params (column=operator.value syntax)
    # httpx will percent-encode % → %25; PostgREST decodes it back to %
    params: dict[str, Any] = {
        "limit": min(limit, 500),
        "select": _SELECT,
    }

    epitope_filter = antigen_epitope or sequence_contains
    if epitope_filter:
        params["linear_sequence"] = f"ilike.%{epitope_filter}%"

    if gene:
        params["mhc_allele_name"] = f"ilike.%{gene}%"

    if species:
        organism = _IMGT_TO_IEDB.get(species, species)
        params["source_organism_name"] = f"ilike.%{organism}%"

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            r = await client.get(IEDB_API_URL, params=params)
            r.raise_for_status()
            rows: list[dict[str, Any]] = r.json() if isinstance(r.json(), list) else []
    except Exception:
        return {"total": 0, "records": []}

    records = [_row_to_record(row) for row in rows]
    return {"total": len(records), "records": records[:limit]}
