"""Tests for the _enrich_with_iedb helper in api.py."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from imgt_app.models import GeneRecord, IEDBHit, SearchResponse


# ---------------------------------------------------------------------------
# Helpers — build fake VDJdb records and IEDB SearchResponse
# ---------------------------------------------------------------------------

def _vdjdb_record(cdr3: str, epitope: str | None) -> GeneRecord:
    return GeneRecord(
        source="vdjdb",
        gene_name="TRBV19",
        sequence=cdr3,
        antigen_epitope=epitope,
    )


def _iedb_response(epitopes: list[str]) -> SearchResponse:
    records = [
        GeneRecord(
            source="iedb",
            gene_name="HLA-A*02:01",
            sequence=ep,
            region="epitope",
            metadata={
                "mhc_class": "I",
                "antigen_name": "Matrix protein",
                "source_organism": "Influenza A virus",
                "assay_type": "cytokine release",
                "effector_cell_type": "CD8+",
                "qualitative_measure": "Positive",
            },
        )
        for ep in epitopes
    ]
    return SearchResponse(total=len(records), records=records)


# ---------------------------------------------------------------------------
# Import target
# ---------------------------------------------------------------------------

from imgt_app.api import _enrich_with_iedb  # noqa: E402


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestEnrichWithIEDB:
    def test_attaches_hits_to_matching_vdjdb_record(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", "GILGFVFTL")]
        iedb = _iedb_response(["GILGFVFTL"])

        result = _enrich_with_iedb(vdjdb, iedb)

        assert len(result) == 1
        assert result[0].iedb_hits is not None
        assert len(result[0].iedb_hits) == 1
        assert result[0].iedb_hits[0].epitope_sequence == "GILGFVFTL"

    def test_hit_fields_mapped_correctly(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", "GILGFVFTL")]
        iedb = _iedb_response(["GILGFVFTL"])

        result = _enrich_with_iedb(vdjdb, iedb)
        hit = result[0].iedb_hits[0]

        assert hit.mhc_allele == "HLA-A*02:01"
        assert hit.mhc_class == "I"
        assert hit.source_organism == "Influenza A virus"
        assert hit.assay_type == "cytokine release"
        assert hit.effector_cell_type == "CD8+"
        assert hit.qualitative_measure == "Positive"

    def test_caps_hits_at_five_per_record(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", "GILGFVFTL")]
        # build 8 IEDB records all matching the same epitope
        iedb_records = [
            GeneRecord(source="iedb", gene_name=f"HLA-A*02:0{i}", sequence="GILGFVFTL",
                       region="epitope", metadata={})
            for i in range(8)
        ]
        iedb = SearchResponse(total=8, records=iedb_records)

        result = _enrich_with_iedb(vdjdb, iedb)
        assert len(result[0].iedb_hits) == 5

    def test_orphan_iedb_hit_creates_phantom_record(self):
        vdjdb: list[GeneRecord] = []  # no VDJdb results
        iedb = _iedb_response(["GILGFVFTL"])

        result = _enrich_with_iedb(vdjdb, iedb)

        assert len(result) == 1
        phantom = result[0]
        assert phantom.source == "vdjdb"
        assert phantom.sequence == ""
        assert phantom.gene_name == "GILGFVFTL"
        assert phantom.iedb_hits is not None
        assert len(phantom.iedb_hits) == 1

    def test_no_epitope_on_vdjdb_record_leaves_iedb_hits_none(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", None)]
        iedb = _iedb_response(["GILGFVFTL"])

        result = _enrich_with_iedb(vdjdb, iedb)

        assert result[0].iedb_hits is None
        # phantom created because GILGFVFTL is orphan
        assert any(r.sequence == "" for r in result)

    def test_empty_vdjdb_and_iedb_returns_empty(self):
        result = _enrich_with_iedb([], SearchResponse(total=0, records=[]))
        assert result == []

    def test_vdjdb_record_with_no_matching_iedb_hit(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", "NLVPMVATV")]
        iedb = _iedb_response(["GILGFVFTL"])  # different epitope

        result = _enrich_with_iedb(vdjdb, iedb)

        # VDJdb record has no matching hit
        vdjdb_result = next(r for r in result if r.sequence == "CASSIRSSYEQYF")
        assert vdjdb_result.iedb_hits is None
        # GILGFVFTL is orphan → phantom
        assert any(r.sequence == "" and r.gene_name == "GILGFVFTL" for r in result)

    def test_case_insensitive_epitope_matching(self):
        vdjdb = [_vdjdb_record("CASSIRSSYEQYF", "gilgfvftl")]  # lowercase
        iedb_records = [
            GeneRecord(source="iedb", gene_name="HLA-A*02:01",
                       sequence="GILGFVFTL",  # uppercase in IEDB
                       region="epitope", metadata={})
        ]
        iedb = SearchResponse(total=1, records=iedb_records)

        result = _enrich_with_iedb(vdjdb, iedb)
        vdjdb_result = next(r for r in result if r.sequence == "CASSIRSSYEQYF")
        assert vdjdb_result.iedb_hits is not None
        assert len(vdjdb_result.iedb_hits) == 1
