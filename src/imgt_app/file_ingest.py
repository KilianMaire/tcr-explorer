from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path
from typing import Iterable

from .fasta_parser import parse_fasta_bytes
from .models import GeneRecord, Species

DNA_RE = re.compile(r"\b[ACGT]{12,}\b", re.IGNORECASE)


def _row_to_record(row: dict, source: str, species: Species) -> GeneRecord:
    return GeneRecord(
        source=source,  # type: ignore[arg-type]
        species=(row.get("species") or species),  # type: ignore[arg-type]
        gene_name=row.get("gene_name") or row.get("gene") or row.get("v_segm") or "unknown",
        allele_name=row.get("allele_name") or row.get("allele"),
        region=row.get("region"),
        sequence=(row.get("sequence") or row.get("cdr3") or "").upper(),
        antigen_epitope=row.get("antigen_epitope") or None,
        metadata={"raw": row},
    )


def _parse_csv_like(raw: bytes, source: str, species: Species, delimiter: str) -> Iterable[GeneRecord]:
    text = raw.decode("utf-8", errors="ignore")
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    for row in reader:
        if row.get("sequence") or row.get("cdr3"):
            yield _row_to_record(row, source, species)


def _parse_json(raw: bytes, source: str, species: Species) -> Iterable[GeneRecord]:
    try:
        data = json.loads(raw.decode("utf-8", errors="ignore"))
    except Exception:
        return []

    rows = data if isinstance(data, list) else [data]
    out = []
    for row in rows:
        if isinstance(row, dict) and row.get("sequence"):
            out.append(_row_to_record(row, source, species))
    return out


def _parse_plain_text(raw: bytes, source: str, species: Species) -> Iterable[GeneRecord]:
    text = raw.decode("utf-8", errors="ignore")
    records = []
    for i, m in enumerate(DNA_RE.finditer(text), start=1):
        seq = m.group(0).upper()
        gene = f"text-seq-{i}"
        records.append(
            GeneRecord(
                source=source,  # type: ignore[arg-type]
                species=species,
                gene_name=gene,
                sequence=seq,
                metadata={"parser": "plain-text"},
            )
        )
    return records


_VDJDB_SPECIES_MAP: dict[str, Species] = {
    "homosapiens": "human",
    "homo sapiens": "human",
    "musmusculus": "mouse",
    "mus musculus": "mouse",
}


def parse_vdjdb_tsv(raw: bytes) -> list[GeneRecord]:
    """Parse a VDJdb TSV (or CSV) export and return GeneRecord list.

    Handles VDJdb-specific column names and normalises species values.
    Works with both tab-delimited (official VDJdb export) and comma-delimited files.
    """
    text = raw.decode("utf-8", errors="ignore")
    first_line = text.split("\n")[0]
    delimiter = "\t" if "\t" in first_line else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    records: list[GeneRecord] = []
    for row in reader:
        cdr3 = (row.get("cdr3") or "").strip().upper()
        if not cdr3:
            continue
        raw_species = (row.get("species") or "").strip()
        species: Species = _VDJDB_SPECIES_MAP.get(raw_species.lower().replace(" ", ""), "other")
        records.append(
            GeneRecord(
                source="vdjdb",
                species=species,
                gene_name=row.get("v_segm") or row.get("v_gene") or "unknown",
                allele_name=row.get("mhc_a") or None,
                region="CDR3",
                sequence=cdr3,
                antigen_epitope=row.get("antigen_epitope") or None,
                metadata={
                    "j_segm": row.get("j_segm"),
                    "mhc_b": row.get("mhc_b"),
                    "mhc_class": row.get("mhc_class"),
                    "antigen_gene": row.get("antigen_gene"),
                    "antigen_species": row.get("antigen_species"),
                    "score": row.get("score"),
                },
            )
        )
    return records


def parse_file(raw: bytes, filename: str, source: str, species: Species) -> list[GeneRecord]:
    ext = Path(filename).suffix.lower()
    if ext in {".fa", ".fasta", ".fna", ".faa"}:
        return list(parse_fasta_bytes(raw, source=source, default_species=species))
    if ext == ".csv":
        return list(_parse_csv_like(raw, source, species, ","))
    if ext in {".tsv", ".txtv"}:
        return list(_parse_csv_like(raw, source, species, "\t"))
    if ext == ".json":
        return list(_parse_json(raw, source, species))

    # Fallback: extract candidate sequences from plain text content.
    return list(_parse_plain_text(raw, source, species))
