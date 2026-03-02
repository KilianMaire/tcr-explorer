"""IPD-MHC Tool Server — queries EBI IPD REST API for MHC allelic variants.

Supports both human (project=HLA) and non-human (project=MHC) species.
Non-human species include NHP (macaque, chimp, gorilla), BoLA (bovine),
DLA (canine), SLA (swine), RT1 (rat), and more.

Run:
    uvicorn servers.mhc_server:app --port 8105 --reload
"""
from __future__ import annotations

import asyncio
import os
from collections import OrderedDict
from typing import Any

import httpx
from fastapi import FastAPI

app = FastAPI(title="IPD-MHC Tool Server (EBI IPD API)", version="0.1.0")

# ---------------------------------------------------------------------------
# EBI IPD REST API
#
# Base:   https://www.ebi.ac.uk/cgi-bin/ipd/api/allele
# List:   GET /allele?project={HLA|MHC}&query=startsWith(name,{prefix})&limit=N
#         -> { "data": [{"accession":"NHP01224","name":"Mamu-A1*026:01:01:01"}, ...],
#              "meta": {"total": 685, "next": "..."} }
# Single: GET /allele/{accession}?project={HLA|MHC}
#         -> { accession, name, class, locus, organism, sequence:{coding,genomic,protein}, ... }
#
# Strategy: list -> collect accessions -> parallel individual fetches (max 10).
# ---------------------------------------------------------------------------
EBI_BASE = os.getenv("EBI_MHC_API_URL", "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele")
_MAX_PARALLEL = 10

# ---------------------------------------------------------------------------
# In-memory LRU cache
# ---------------------------------------------------------------------------
_CACHE: OrderedDict[tuple[str, str, str], list[dict[str, Any]]] = OrderedDict()
_CACHE_MAX = 64

# ---------------------------------------------------------------------------
# Species prefix → IPD project mapping
# Human HLA alleles start with locus letters (A*, B*, DRB1*, etc.)
# Non-human alleles start with species prefix (Mamu-, BoLA-, Patr-, etc.)
# ---------------------------------------------------------------------------
_KNOWN_NHP_PREFIXES = {
    "mamu", "patr", "gogo", "popy", "hymo",   # primates
    "aona", "caja", "saoe", "sasc",            # new world monkeys
    "mane", "mafa", "macy",                    # macaques
}
_KNOWN_MHC_PREFIXES = {
    "bola", "dla", "sla", "rt1", "eqca",      # bovine, canine, swine, rat, horse
    "ovar",                                     # sheep
}
_ALL_MHC_PREFIXES = _KNOWN_NHP_PREFIXES | _KNOWN_MHC_PREFIXES


def _detect_project(gene_name: str) -> str:
    """Determine EBI project (HLA or MHC) from gene name prefix."""
    lower = gene_name.lower().replace("-", "").replace("_", "")
    for prefix in _ALL_MHC_PREFIXES:
        if lower.startswith(prefix):
            return "MHC"
    # If starts with HLA or is a bare locus like A*, B*, DRB1* → HLA
    if gene_name.upper().startswith("HLA") or (len(gene_name) <= 5 and "*" in gene_name):
        return "HLA"
    # Default: try MHC first (broader database)
    return "MHC"


def _build_query_prefix(gene_name: str) -> str:
    """Build the allele name prefix for the startsWith query.

    Examples:
        'Mamu-A1' → 'Mamu-A1'
        'BoLA-2'  → 'BoLA-2'
        'HLA-A'   → 'A'  (HLA project uses bare locus)
        'DRB1'    → 'DRB1'
    """
    name = gene_name.strip()
    if name.upper().startswith("HLA-"):
        return name[4:]  # strip 'HLA-' prefix for HLA project
    return name


def _allele_to_record(item: dict[str, Any], project: str) -> dict[str, Any]:
    """Map a full single-allele EBI response to our GeneRecord schema."""
    name: str = item.get("name", "")
    locus = item.get("locus", "")
    organism = item.get("organism") or {}

    # Derive gene_name from allele name (part before '*')
    gene_name = name.split("*")[0] if "*" in name else name
    if project == "HLA" and not gene_name.upper().startswith("HLA"):
        gene_name = f"HLA-{gene_name}"

    seq_obj = item.get("sequence") or {}
    coding_seq: str = seq_obj.get("coding") or ""
    protein_seq: str = seq_obj.get("protein") or ""

    # Species mapping
    common_name = (organism.get("commonName") or "").lower()
    species = "human" if "human" in common_name or project == "HLA" else "other"

    return {
        "source": "mhc",
        "species": species,
        "gene_name": gene_name,
        "allele_name": name,
        "region": "coding",
        "sequence": coding_seq,
        "metadata": {
            "accession": item.get("accession"),
            "mhc_class": item.get("class"),
            "locus": locus,
            "organism_common": organism.get("commonName"),
            "organism_scientific": organism.get("scientificName"),
            "organism_group": organism.get("group"),
            "taxon": organism.get("taxon"),
            "genomic_seq_len": len(seq_obj.get("genomic") or ""),
            "protein_seq": protein_seq[:80],
            "date_assigned": item.get("date_assigned"),
            "date_modified": item.get("date_modified"),
            "status": item.get("status"),
            "backend": "ebi-ipd-mhc",
        },
    }


async def _fetch_single(
    client: httpx.AsyncClient, accession: str, project: str
) -> dict[str, Any] | None:
    """Fetch one allele by accession and return mapped record, or None on error."""
    try:
        r = await client.get(EBI_BASE + f"/{accession}", params={"project": project})
        r.raise_for_status()
        return _allele_to_record(r.json(), project)
    except Exception:
        return None


async def _fetch_ebi(
    gene_name: str,
    seq_contains: str,
    limit: int,
) -> list[dict[str, Any]]:
    """Query EBI IPD API with results cached by (project, prefix, seq_contains)."""
    project = _detect_project(gene_name)
    prefix = _build_query_prefix(gene_name)

    cache_key = (project, prefix.upper(), seq_contains.upper())
    if cache_key in _CACHE:
        _CACHE.move_to_end(cache_key)
        return _CACHE[cache_key][:limit]

    list_params: dict[str, Any] = {
        "project": project,
        "limit": _MAX_PARALLEL,
    }
    if prefix:
        list_params["query"] = f"startsWith(name,{prefix})"

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            r = await client.get(EBI_BASE, params=list_params)
            r.raise_for_status()
            accessions = [item["accession"] for item in r.json().get("data", [])]
            if not accessions:
                return []

            tasks = [_fetch_single(client, acc, project) for acc in accessions]
            results = await asyncio.gather(*tasks)
    except Exception:
        return []

    records: list[dict[str, Any]] = []
    for rec in results:
        if rec is None:
            continue
        if seq_contains and seq_contains.upper() not in rec["sequence"].upper():
            continue
        records.append(rec)

    # LRU cache
    _CACHE[cache_key] = records
    _CACHE.move_to_end(cache_key)
    while len(_CACHE) > _CACHE_MAX:
        _CACHE.popitem(last=False)

    return records[:limit]


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "mhc"}


@app.post("/search")
async def search(req: dict) -> dict:
    limit = int(req.get("limit", 50))
    gene = (req.get("gene_name") or "").strip()
    region = (req.get("region") or "").strip().lower()
    seq_contains = (req.get("sequence_contains") or "").strip()

    if not gene:
        return {"total": 0, "records": []}

    records = await _fetch_ebi(gene, seq_contains, limit)

    # Region post-filter
    if region:
        records = [r for r in records if region in (r.get("region") or "").lower()]

    out = records[:limit]
    return {"total": len(out), "records": out}
