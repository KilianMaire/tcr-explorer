from __future__ import annotations

import io
import re
from typing import Iterable, Optional, Tuple

from Bio import SeqIO

from .models import GeneRecord, Species


def _infer_species(text: str) -> Species:
    t = text.lower()
    if "homo sapiens" in t or "human" in t:
        return "human"
    if "mus musculus" in t or "mouse" in t:
        return "mouse"
    return "other"


def _parse_header(description: str) -> Tuple[str, Optional[str], Optional[str]]:
    # Common style: >ALLELE|GENE|REGION|...
    parts = [p.strip() for p in description.split("|")]
    allele = parts[0] if parts else None
    gene = None
    region = None

    if len(parts) > 1 and parts[1]:
        gene = parts[1]

    if len(parts) > 2 and parts[2]:
        region = parts[2]

    if not gene:
        m = re.search(r"\b([A-Z]{1,4}[A-Z0-9-]*\*?[0-9:]*)\b", description)
        if m:
            gene = m.group(1)

    return gene or "unknown", allele, region


def parse_cdr3_fasta(raw: bytes) -> list[tuple[str, str]]:
    """
    Parse a FASTA file where each sequence is a CDR3 amino-acid string.

    Returns a list of (record_id, cdr3_sequence) tuples.
    Sequences containing only valid amino-acid characters (A-Z) are returned
    as-is; sequences that look like nucleotides (only ACGTN) are also accepted
    — the caller decides how to use them.
    """
    handle = io.StringIO(raw.decode("utf-8", errors="ignore"))
    out: list[tuple[str, str]] = []
    for rec in SeqIO.parse(handle, "fasta"):
        seq = str(rec.seq).strip().upper()
        if seq:
            out.append((rec.id, seq))
    return out


def parse_fasta_bytes(raw: bytes, source: str, default_species: Species = "other") -> Iterable[GeneRecord]:
    handle = io.StringIO(raw.decode("utf-8", errors="ignore"))
    for rec in SeqIO.parse(handle, "fasta"):
        gene_name, allele_name, region = _parse_header(rec.description)
        inferred_species = _infer_species(rec.description)
        species = inferred_species if inferred_species != "other" else default_species
        yield GeneRecord(
            source=source,  # type: ignore[arg-type]
            species=species,
            gene_name=gene_name,
            allele_name=allele_name,
            region=region,
            sequence=str(rec.seq).upper(),
            metadata={"description": rec.description, "record_id": rec.id},
        )
