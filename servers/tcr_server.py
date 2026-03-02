from __future__ import annotations

import asyncio
import csv
import io
import os
from collections import OrderedDict
from pathlib import Path
from typing import Any

import httpx
from fastapi import FastAPI

app = FastAPI(title="TCR Tool Server (IMGT genes via NCBI Entrez)", version="0.2.0")

# ---------------------------------------------------------------------------
# Local seed fallback
# ---------------------------------------------------------------------------
_DEFAULT_SEED = Path(__file__).resolve().parent.parent / "examples" / "tcr_seed.csv"
DATA_FILE = Path(os.getenv("TCR_DATA_FILE", str(_DEFAULT_SEED)))
SEED_RECORDS: list[dict[str, str]] = []
if DATA_FILE.exists():
    with DATA_FILE.open() as f:
        SEED_RECORDS = [dict(r) for r in csv.DictReader(f)]

# ---------------------------------------------------------------------------
# NCBI Entrez configuration
# IMGT TCR gene segments (TRAV, TRBV, TRGV, TRDV, TRAJ, TRBJ …) are deposited
# as individual GenBank/EMBL records.  We search nuccore with:
#   - [Gene Name]    → IMGT gene symbol annotation on the feature
#   - [SLEN]         → sequence length filter to exclude full chromosomes
#   - [Organism]     → species filter
# The sequences returned are germline V/J/D gene segments (~200–600 bp).
# ---------------------------------------------------------------------------
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")   # optional — raises rate limit to 10 req/s

# ---------------------------------------------------------------------------
# In-memory LRU cache — keyed by (gene_name, species, seq_contains); avoids
# redundant NCBI round-trips when the same query is repeated in a session.
# ---------------------------------------------------------------------------
_TCR_CACHE: OrderedDict[tuple[str, str, str], list[dict[str, Any]]] = OrderedDict()
_TCR_CACHE_MAX = 64

# Species mapping: internal → NCBI taxon name
_SPECIES_NCBI: dict[str, str] = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
}

# TCR V-gene segment size range (bp).  Filters out chromosomes / large contigs.
_MIN_LEN = 200
_MAX_LEN = 700


async def _entrez_esearch(term: str, retmax: int) -> list[str]:
    """Return NCBI nuccore IDs for the given search term."""
    params: dict[str, Any] = {
        "db": "nuccore",
        "term": term,
        "retmode": "json",
        "retmax": retmax,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        async with httpx.AsyncClient(timeout=12.0) as client:
            r = await client.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params)
            r.raise_for_status()
            return r.json()["esearchresult"]["idlist"]
    except Exception:
        return []


async def _entrez_efetch_fasta(ids: list[str]) -> str:
    """Fetch sequences in FASTA format for a list of nuccore IDs."""
    if not ids:
        return ""
    params: dict[str, Any] = {
        "db": "nuccore",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params)
            r.raise_for_status()
            return r.text
    except Exception:
        return ""


def _parse_fasta_text(fasta_text: str) -> list[tuple[str, str, str]]:
    """
    Parse raw FASTA text → list of (ncbi_id, description, sequence).
    Skips sequences outside the expected V-gene segment length range.
    """
    results: list[tuple[str, str, str]] = []
    current_id = ""
    current_desc = ""
    current_seq_parts: list[str] = []

    def _flush() -> None:
        if current_id:
            seq = "".join(current_seq_parts).upper()
            if _MIN_LEN <= len(seq) <= _MAX_LEN:
                results.append((current_id, current_desc, seq))

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            _flush()
            parts = line[1:].split(None, 1)
            current_id = parts[0] if parts else ""
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq_parts = []
        else:
            current_seq_parts.append(line)
    _flush()
    return results


def _infer_species(description: str) -> str:
    desc_lower = description.lower()
    if "homo sapiens" in desc_lower or "human" in desc_lower:
        return "human"
    if "mus musculus" in desc_lower or "mouse" in desc_lower:
        return "mouse"
    return "other"


def _infer_region(gene_name: str) -> str:
    g = gene_name.upper()
    if g.startswith(("TRAV", "TRBV", "TRGV", "TRDV")):
        return "v-region"
    if g.startswith(("TRAJ", "TRBJ", "TRGJ", "TRDJ")):
        return "j-region"
    if g.startswith(("TRBD", "TRDD")):
        return "d-region"
    if g.startswith(("TRAC", "TRBC", "TRGC", "TRDC")):
        return "c-region"
    return "gene-segment"


async def _fetch_ncbi_imgt(
    gene_name: str,
    species: str,
    seq_contains: str,
    limit: int,
) -> list[dict[str, Any]]:
    """
    Query NCBI nuccore for IMGT TCR gene segments, with results cached by
    (gene_name, species, seq_contains).

    Search strategy:
      "{gene_name}"[Gene Name]     → IMGT symbol on the GenBank feature
      {MIN}:{MAX}[SLEN]            → exclude chromosomes / large contigs
      "{ncbi_species}"[Organism]   → optional species filter
    """
    if not gene_name:
        return []

    cache_key = (gene_name.upper(), species.lower(), seq_contains.upper())
    if cache_key in _TCR_CACHE:
        _TCR_CACHE.move_to_end(cache_key)
        return _TCR_CACHE[cache_key][:limit]

    ncbi_species = _SPECIES_NCBI.get(species, "") if species else ""

    term_parts = [
        f'"{gene_name}"[Gene Name]',
        f"{_MIN_LEN}:{_MAX_LEN}[SLEN]",
    ]
    if ncbi_species:
        term_parts.append(f'"{ncbi_species}"[Organism]')

    term = " AND ".join(term_parts)
    ids = await _entrez_esearch(term, retmax=40)  # fetch max so cache covers any limit
    if not ids:
        return []

    fasta_text = await _entrez_efetch_fasta(ids[:20])
    parsed = _parse_fasta_text(fasta_text)

    records: list[dict[str, Any]] = []
    for ncbi_id, desc, seq in parsed:
        if seq_contains and seq_contains.upper() not in seq:
            continue
        sp = _infer_species(desc)
        if species and species != "other" and sp != species:
            continue
        records.append({
            "source": "tcr",
            "species": sp,
            "gene_name": gene_name,
            "allele_name": ncbi_id,
            "region": _infer_region(gene_name),
            "sequence": seq,
            "metadata": {
                "ncbi_id": ncbi_id,
                "description": desc[:200],
                "seq_len": len(seq),
                "backend": "ncbi-entrez-imgt",
            },
        })

    # Populate LRU cache; evict least-recently-used entry when full
    _TCR_CACHE[cache_key] = records
    _TCR_CACHE.move_to_end(cache_key)
    while len(_TCR_CACHE) > _TCR_CACHE_MAX:
        _TCR_CACHE.popitem(last=False)

    return records[:limit]


# ---------------------------------------------------------------------------
# Local seed search (always included as supplement)
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
            "source": "tcr",
            "species": r.get("species") or "other",
            "gene_name": r.get("gene_name") or "unknown",
            "allele_name": r.get("allele_name"),
            "region": r.get("region"),
            "sequence": r.get("sequence") or "",
            "metadata": {"backend": "tcr-seed"},
        })
    return out[:limit]


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "server": "tcr"}


@app.post("/search")
async def search(req: dict) -> dict:
    limit = int(req.get("limit", 50))
    gene = (req.get("gene_name") or "").strip()
    region = (req.get("region") or "").strip().lower()
    species = (req.get("species") or "").strip().lower()
    seq_contains = (req.get("sequence_contains") or "").strip()

    # Fetch from NCBI IMGT when a gene name is provided
    ncbi_records = await _fetch_ncbi_imgt(gene, species, seq_contains, limit)

    # Region post-filter
    if region:
        ncbi_records = [r for r in ncbi_records if region in (r.get("region") or "").lower()]

    # Always supplement with seed records (deduplicate by allele_name/NCBI ID)
    seed_records = _search_seed(gene, species, region, seq_contains, limit)
    seen = {r["allele_name"] for r in ncbi_records if r["allele_name"]}
    for sr in seed_records:
        if sr["allele_name"] not in seen:
            ncbi_records.append(sr)

    out = ncbi_records[:limit]
    return {"total": len(out), "records": out}
