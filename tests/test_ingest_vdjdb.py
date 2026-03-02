"""Tests for /ingest/vdjdb endpoint and parse_vdjdb_tsv parser."""
from __future__ import annotations

import io

import pytest
from fastapi.testclient import TestClient

# Use PYTHONPATH=src convention — tests are run from project root.
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from imgt_app.file_ingest import parse_vdjdb_tsv

SAMPLE_TSV = (
    "cdr3\tv_segm\tj_segm\tspecies\tmhc_a\tmhc_b\tmhc_class\tantigen_epitope\tantigen_gene\tantigen_species\tscore\n"
    "CASSIRSSYEQYF\tTRBV19\tTRBJ2-7\tHomoSapiens\tHLA-A*02:01\tB2M\tI\tGILGFVFTL\tM1\tInfluenzaA\t3\n"
    "CASGDSSYEQYF\tTRBV13-2\tTRBJ2-7\tMusMusculus\tH-2Kb\tB2M\tI\tSIINFEKL\tOVA\tHomoSapiens\t3\n"
    "CAVLDSNYQLIW\tTRAV41\tTRAJ33\tUnknownSpecies\tHLA-A*02:01\tB2M\tI\tNLVPMVATV\tpp65\tCMV\t2\n"
)

SAMPLE_CSV = SAMPLE_TSV.replace("\t", ",")


class TestParseVdjdbTsv:
    def test_parses_tsv_rows(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert len(records) == 3

    def test_parses_csv_rows(self):
        records = parse_vdjdb_tsv(SAMPLE_CSV.encode())
        assert len(records) == 3

    def test_source_is_vdjdb(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert all(r.source == "vdjdb" for r in records)

    def test_species_normalisation(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].species == "human"
        assert records[1].species == "mouse"
        assert records[2].species == "other"

    def test_sequence_is_cdr3_uppercased(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].sequence == "CASSIRSSYEQYF"

    def test_gene_name_from_v_segm(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].gene_name == "TRBV19"

    def test_antigen_epitope_captured(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].antigen_epitope == "GILGFVFTL"

    def test_allele_name_from_mhc_a(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].allele_name == "HLA-A*02:01"

    def test_region_is_cdr3(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert all(r.region == "CDR3" for r in records)

    def test_empty_file_returns_no_records(self):
        records = parse_vdjdb_tsv(b"cdr3\tv_segm\tspecies\n")
        assert records == []

    def test_metadata_includes_score(self):
        records = parse_vdjdb_tsv(SAMPLE_TSV.encode())
        assert records[0].metadata["score"] == "3"


class TestIngestVdjdbEndpoint:
    @pytest.fixture
    def client(self, tmp_path):
        import os
        os.environ.setdefault("DATABASE_PATH", str(tmp_path / "test.db"))
        # Re-import to pick up temp db path
        from imgt_app.config import settings
        settings.database_path = str(tmp_path / "test.db")
        from imgt_app.api import app, index
        index.db_path = settings.database_path
        index._init_db()
        return TestClient(app)

    def test_ingest_returns_inserted_count(self, client):
        response = client.post(
            "/ingest/vdjdb",
            files={"file": ("vdjdb.tsv", SAMPLE_TSV.encode(), "text/tab-separated-values")},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["inserted"] == 3
        assert data["source"] == "vdjdb"

    def test_ingest_makes_records_searchable(self, client):
        client.post(
            "/ingest/vdjdb",
            files={"file": ("vdjdb.tsv", SAMPLE_TSV.encode(), "text/tab-separated-values")},
        )
        search_resp = client.post(
            "/search",
            json={"source": "vdjdb", "antigen_epitope": "GILGFVFTL", "limit": 10},
        )
        assert search_resp.status_code == 200
        results = search_resp.json()
        assert results["total"] >= 1
        assert any(r["antigen_epitope"] == "GILGFVFTL" for r in results["records"])
