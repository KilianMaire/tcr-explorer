from __future__ import annotations

import asyncio
import csv
import os
from collections import OrderedDict
from pathlib import Path
from typing import Any

import httpx
from fastapi import FastAPI

app = FastAPI(title="HLA Tool Server (IMGT/HLA via EBI)", version="0.2.0")

# ---------------------------------------------------------------------------
# Local seed fallback
# ---------------------------------------------------------------------------
_DEFAULT_SEED = Path(__file__).resolve().parent.parent / "examples" / "hla_seed.csv"
DATA_FILE = Path(os.getenv("HLA_DATA_FILE", str(_DEFAULT_SEED)))
SEED_RECORDS: list[dict[str, str]] = []
if DATA_FILE.exists():
    with DATA_FILE.open() as f:
        SEED_RECORDS = [dict(r) for r in csv.DictReader(f)]

# ---------------------------------------------------------------------------
# EBI IPD-IMGT/HLA REST API
#
# List endpoint:  GET /allele?project=HLA&query=startsWith(name,{locus})&limit=N
#   -> { "data": [{"accession":"HLA00001","name":"A*01:01:01:01"}, ...],
#        "meta": {"total": 45069, "next": "..."} }
#   Note: the list does NOT include sequence data.
#
# Single-allele:  GET /allele/{accession}
#   -> { "accession", "name", "sequence": {"coding","genomic","protein"}, ... }
#
# Strategy: list -> collect accessions -> parallel individual fetches (max 10).
# ---------------------------------------------------------------------------
EBI_BASE = os.getenv("EBI_HLA_API_URL", "https://www.ebi.ac.uk/cgi-bin/ipd/api/allele")
_MAX_PARALLEL = 10   # max individual allele fetches per request

# ---------------------------------------------------------------------------
# In-memory LRU cache — keyed by (locus, seq_contains); avoids redundant EBI
# round-trips when the same query is repeated within a server session.
# ---------------------------------------------------------------------------
_EBI_CACHE: OrderedDict[tuple[str, str], list[dict[str, Any]]] = OrderedDict()
_EBI_CACHE_MAX = 64

_LOCUS_FROM_GENE: dict[str, str] = {
    "hla-a": "A", "hla-b": "B", "hla-c": "C",
    "hla-dpa1": "DPA1", "hla-dpb1": "DPB1",
    "hla-dqa1": "DQA1", "hla-dqb1": "DQB1",
    "hla-dra": "DRA", "hla-drb1": "DRB1", "hla-drb3": "DRB3",
    "hla-drb4": "DRB4", "hla-drb5": "DRB5",
    "hla-e": "E", "hla-f": "F", "hla-g": "G",
}


def _gene_to_locus(gene: str) -> str:
    """Convert 'HLA-A' or 'A' to EBI locus prefix (e.g. 'A')."""
    return _LOCUS_FROM_GENE.get(gene.lower(), gene.upper().replace("HLA-", ""))


def _allele_item_to_record(item: dict[str, Any]) -> dict[str, Any]:
    """Map a full single-allele EBI response to our GeneRecord schema."""
    name: str = item.get("name", "")          # e.g. "A*01:01:01:01"
    locus = name.split("*")[0] if "*" in name else name
    gene_name = f"HLA-{locus}" if locus else "unknown"
    seq_obj = item.get("sequence") or {}
    coding_seq: str = seq_obj.get("coding") or ""
    return {
        "source": "hla",
        "species": "human",
        "gene_name": gene_name,
        "allele_name": name,
        "region": "coding",
        "sequence": coding_seq,
        "metadata": {
            "accession": item.get("accession"),
            "hla_class": item.get("class"),
            "genomic_seq_len": len(seq_obj.get("genomic") or ""),
            "protein_seq": (seq_obj.get("protein") or "")[:60],
            "backend": "ebi-imgt-hla",
        },
    }


async def _fetch_single(client: httpx.AsyncClient, accession: str) -> dict[str, Any] | None:
    """Fetch one allele by accession and return mapped record, or None on error."""
    try:
        r = await client.get(f"{EBI_BASE}/{accession}")
        r.raise_for_status()
        return _allele_item_to_record(r.json())
    except Exception:
        return None


async def _fetch_ebi(
    locus: str,
    seq_contains: str,
    limit: int,
) -> list[dict[str, Any]]:
    """
    Query EBI IMGT/HLA API, with results cached by (locus, seq_contains).
    Step 1: list endpoint -> get accession list for the given locus.
    Step 2: parallel individual fetches (up to _MAX_PARALLEL) -> get sequences.
    """
    cache_key = (locus.upper(), seq_contains.upper())
    if cache_key in _EBI_CACHE:
        _EBI_CACHE.move_to_end(cache_key)
        return _EBI_CACHE[cache_key][:limit]

    list_params: dict[str, Any] = {
        "project": "HLA",
        "limit": _MAX_PARALLEL,  # always fetch max so cached results cover any limit
    }
    if locus:
        # allele names start with "{locus}*" (e.g. "A*01:01"), filter by prefix
        list_params["query"] = f"startsWith(name,{locus})"

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            r = await client.get(EBI_BASE, params=list_params)
            r.raise_for_status()
            accessions = [item["accession"] for item in r.json().get("data", [])]
            if not accessions:
                return []

            # Parallel fetch of full records (sequences included)
            tasks = [_fetch_single(client, acc) for acc in accessions]
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

    # Populate LRU cache; evict least-recently-used entry when full
    _EBI_CACHE[cache_key] = records
    _EBI_CACHE.move_to_end(cache_key)
    while len(_EBI_CACHE) > _EBI_CACHE_MAX:
        _EBI_CACHE.popitem(last=False)

    return records[:limit]


# ---------------------------------------------------------------------------
# Local seed search (supplement for non-human or when EBI returns nothing)
# ---------------------------------------------------------------------------
def _search_seed(
    gene: str, species: str, region: str, seq_contains: str, limit: int
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for r in SEED_RECORDS:
        if gene and gene.lower() not in r["gene_name"].lower():
            continue
        if region and region.lower() not in (r.get("region") or "").lower():
            continue
        if species and species != (r.get("species") or "").lower():
            continue
        if seq_contains and seq_contains.upper() not in (r.get("sequence") or "").upper():
            continue
        out.append({
            "source": "hla",
            "species": r.get("species") or "other",
            "gene_name": r.get("gene_name") or "unknown",
            "allele_name": r.get("allele_name"),
            "region": r.get("region"),
            "sequence": r.get("sequence") or "",
            "metadata": {"backend": "hla-seed"},
        })
    return out[:limit]


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "hla"}


@app.post("/search")
async def search(req: dict) -> dict:
    limit = int(req.get("limit", 50))
    gene = (req.get("gene_name") or "").strip()
    region = (req.get("region") or "").strip().lower()
    species = (req.get("species") or "").strip().lower()
    seq_contains = (req.get("sequence_contains") or "").strip()

    locus = _gene_to_locus(gene) if gene else ""

    # EBI IMGT/HLA is human-only
    ebi_records: list[dict[str, Any]] = []
    if not species or species == "human":
        ebi_records = await _fetch_ebi(locus, seq_contains, limit)

    # Region post-filter
    if region:
        ebi_records = [r for r in ebi_records if region in (r.get("region") or "").lower()]

    # Supplement with seed for non-human or when EBI returned nothing
    seed_records = _search_seed(gene, species, region, seq_contains, limit)
    seen = {r["allele_name"] for r in ebi_records if r["allele_name"]}
    for sr in seed_records:
        if sr["allele_name"] not in seen:
            ebi_records.append(sr)

    out = ebi_records[:limit]
    return {"total": len(out), "records": out}
