from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import List, Union

from .models import GeneRecord, SearchRequest, SearchResponse


class SearchIndex:
    def __init__(self, db_path: str) -> None:
        self.db_path = db_path
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    @contextmanager
    def _conn(self):
        conn = sqlite3.connect(self.db_path)
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()

    def _init_db(self) -> None:
        with self._conn() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS gene_records (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    source TEXT NOT NULL,
                    species TEXT NOT NULL,
                    gene_name TEXT NOT NULL,
                    allele_name TEXT,
                    region TEXT,
                    sequence TEXT NOT NULL,
                    antigen_epitope TEXT,
                    metadata_json TEXT NOT NULL
                )
                """
            )
            # Migration: add antigen_epitope to databases created before this column existed
            try:
                conn.execute("ALTER TABLE gene_records ADD COLUMN antigen_epitope TEXT")
            except Exception:
                pass  # column already exists
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_records_main
                ON gene_records (source, species, gene_name, region)
                """
            )
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_records_epitope
                ON gene_records (antigen_epitope)
                """
            )

    def upsert_many(self, records: list[GeneRecord]) -> int:
        if not records:
            return 0
        rows = [
            (
                r.source,
                r.species,
                r.gene_name,
                r.allele_name,
                r.region,
                r.sequence,
                r.antigen_epitope,
                json.dumps(r.metadata, ensure_ascii=True),
            )
            for r in records
        ]
        with self._conn() as conn:
            conn.executemany(
                """
                INSERT INTO gene_records
                (source, species, gene_name, allele_name, region, sequence, antigen_epitope, metadata_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows,
            )
        return len(records)

    def search(self, req: SearchRequest) -> SearchResponse:
        where = ["1=1"]
        params: List[Union[str, int]] = []

        if req.source:
            where.append("source = ?")
            params.append(req.source)
        if req.species:
            where.append("species = ?")
            params.append(req.species)
        if req.gene_name:
            where.append("gene_name LIKE ?")
            params.append(f"%{req.gene_name}%")
        if req.region:
            where.append("region LIKE ?")
            params.append(f"%{req.region}%")
        if req.sequence_contains:
            where.append("sequence LIKE ?")
            params.append(f"%{req.sequence_contains.upper()}%")
        if req.antigen_epitope:
            where.append("antigen_epitope LIKE ?")
            params.append(f"%{req.antigen_epitope}%")

        # Fetch limit + offset rows so the API layer has enough records to apply
        # the offset after merging local and remote results.
        fetch_limit = req.limit + req.offset
        sql = f"""
            SELECT source, species, gene_name, allele_name, region, sequence, antigen_epitope, metadata_json
            FROM gene_records
            WHERE {' AND '.join(where)}
            LIMIT ?
        """
        params.append(fetch_limit)

        count_sql = f"SELECT COUNT(*) FROM gene_records WHERE {' AND '.join(where)}"

        with self._conn() as conn:
            total = conn.execute(count_sql, params[:-1]).fetchone()[0]
            rows = conn.execute(sql, params).fetchall()

        records = [
            GeneRecord(
                source=row[0],
                species=row[1],
                gene_name=row[2],
                allele_name=row[3],
                region=row[4],
                sequence=row[5],
                antigen_epitope=row[6],
                metadata=json.loads(row[7]),
            )
            for row in rows
        ]
        return SearchResponse(total=total, records=records, limit=req.limit, offset=req.offset)
