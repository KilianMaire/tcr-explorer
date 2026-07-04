from __future__ import annotations

import json
import re

import httpx

from .config import settings
from .models import ParseQueryResult


PROMPT = """
Convert user query into JSON with keys:
source (hla|tcr|vdjdb|iedb|mhc|null), species (human|mouse|other|null), gene_name, region,
sequence_contains, antigen_epitope.
Use "mhc" when the query mentions non-human MHC species (Mamu, BoLA, Patr, DLA, SLA, RT1, macaque, bovine, canine, swine) or IPD-MHC.
Use "iedb" when the query mentions IEDB, immune epitope database, T-cell assay, MHC ligand assay, or effector cell.
Use "vdjdb" when the query mentions CDR3, antigen epitope, or T-cell receptor specificity.
Return JSON only.
""".strip()


def heuristic_parse(query: str) -> ParseQueryResult:
    q = query.lower()

    # MHC species prefixes (non-human)
    _MHC_KEYWORDS = [
        "mamu", "patr", "gogo", "popy", "bola", "dla", "sla", "rt1",
        "eqca", "ovar", "mane", "mafa", "caja", "aona",
        "macaque", "bovine", "canine", "swine", "rhesus",
        "chimpanzee", "gorilla", "ipd-mhc", "ipd mhc",
    ]

    source = None
    if "hla" in q:
        source = "hla"
    elif any(kw in q for kw in _MHC_KEYWORDS):
        source = "mhc"
    elif "iedb" in q or "immune epitope database" in q or "tcell assay" in q or "mhc ligand" in q or "effector cell" in q:
        source = "iedb"
    elif "vdjdb" in q or "cdr3" in q or "epitope" in q or "antigen" in q or "specificity" in q:
        source = "vdjdb"
    elif "tcr" in q or "t-cell receptor" in q:
        source = "tcr"

    species = None
    if "human" in q or "homo sapiens" in q:
        species = "human"
    elif "mouse" in q or "mus musculus" in q:
        species = "mouse"

    region = None
    for r in ["v-region", "d-region", "j-region", "c-region", "exon", "intron", "leader"]:
        if r in q:
            region = r
            break

    gene_name = None
    m = re.search(r"\b([A-Z]{2,6}[A-Z0-9-]*\*?[0-9:]*)\b", query)
    if m:
        gene_name = m.group(1)

    sequence_contains = None
    m2 = re.search(r"\b([ACGT]{8,})\b", query.upper())
    if m2:
        sequence_contains = m2.group(1)

    # Antigen epitope: amino-acid peptide pattern (uppercase letters, 8-15 AA)
    antigen_epitope = None
    m3 = re.search(r"\b([ACDEFGHIKLMNPQRSTVWY]{8,15})\b", query.upper())
    if m3 and m3.group(1) != sequence_contains:
        antigen_epitope = m3.group(1)

    return ParseQueryResult(
        source=source,  # type: ignore[arg-type]
        species=species,  # type: ignore[arg-type]
        gene_name=gene_name,
        region=region,
        sequence_contains=sequence_contains,
        antigen_epitope=antigen_epitope,
    )


async def lmstudio_parse(query: str) -> ParseQueryResult:
    payload = {
        "model": settings.lmstudio_model,
        "messages": [
            {"role": "system", "content": PROMPT},
            {"role": "user", "content": query},
        ],
        "temperature": 0,
    }

    try:
        async with httpx.AsyncClient(timeout=8.0) as client:
            r = await client.post(f"{settings.lmstudio_base_url}/chat/completions", json=payload)
            r.raise_for_status()
            content = r.json()["choices"][0]["message"]["content"].strip()
            parsed = json.loads(content)
            return ParseQueryResult(**parsed)
    except Exception:
        return heuristic_parse(query)
