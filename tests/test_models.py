"""Tests for Pydantic model validation in models.py."""
from __future__ import annotations

import sys
import pathlib

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from pydantic import ValidationError

from imgt_app.models import (
    GeneRecord,
    IEDBHit,        # ← add this
    SearchRequest,
    SearchResponse,
    NLQueryRequest,
    ParseQueryResult,
    IngestResponse,
    CDRPredictResponse,
)


class TestGeneRecord:
    def test_minimal_valid_record(self):
        r = GeneRecord(source="hla", gene_name="HLA-A", sequence="ATCG")
        assert r.source == "hla"
        assert r.species == "other"  # default
        assert r.gene_name == "HLA-A"
        assert r.sequence == "ATCG"

    def test_all_sources_accepted(self):
        for src in ("hla", "tcr", "vdjdb", "iedb", "mhc"):
            r = GeneRecord(source=src, gene_name="X", sequence="ATCG")
            assert r.source == src

    def test_invalid_source_raises(self):
        with pytest.raises(ValidationError):
            GeneRecord(source="unknown", gene_name="X", sequence="ATCG")

    def test_all_species_accepted(self):
        for sp in ("human", "mouse", "other"):
            r = GeneRecord(source="hla", gene_name="X", sequence="ATCG", species=sp)
            assert r.species == sp

    def test_invalid_species_raises(self):
        with pytest.raises(ValidationError):
            GeneRecord(source="hla", gene_name="X", sequence="ATCG", species="rat")

    def test_optional_fields_default_to_none(self):
        r = GeneRecord(source="hla", gene_name="X", sequence="ATCG")
        assert r.allele_name is None
        assert r.region is None
        assert r.antigen_epitope is None

    def test_metadata_defaults_to_empty_dict(self):
        r = GeneRecord(source="hla", gene_name="X", sequence="ATCG")
        assert r.metadata == {}

    def test_antigen_epitope_stored(self):
        r = GeneRecord(source="vdjdb", gene_name="TRBV19", sequence="CASSIRSSYEQYF",
                       antigen_epitope="GILGFVFTL")
        assert r.antigen_epitope == "GILGFVFTL"


class TestSearchRequest:
    def test_defaults(self):
        req = SearchRequest()
        assert req.source is None
        assert req.species is None
        assert req.gene_name is None
        assert req.region is None
        assert req.sequence_contains is None
        assert req.antigen_epitope is None
        assert req.limit == 50

    def test_limit_boundary_min(self):
        req = SearchRequest(limit=1)
        assert req.limit == 1

    def test_limit_boundary_max(self):
        req = SearchRequest(limit=500)
        assert req.limit == 500

    def test_limit_below_min_raises(self):
        with pytest.raises(ValidationError):
            SearchRequest(limit=0)

    def test_limit_above_max_raises(self):
        with pytest.raises(ValidationError):
            SearchRequest(limit=501)

    def test_valid_source_accepted(self):
        req = SearchRequest(source="vdjdb")
        assert req.source == "vdjdb"

    def test_invalid_source_raises(self):
        with pytest.raises(ValidationError):
            SearchRequest(source="bad")

    def test_valid_species_accepted(self):
        req = SearchRequest(species="human")
        assert req.species == "human"

    def test_invalid_species_raises(self):
        with pytest.raises(ValidationError):
            SearchRequest(species="cat")

    def test_antigen_epitope_field(self):
        req = SearchRequest(antigen_epitope="GILGFVFTL")
        assert req.antigen_epitope == "GILGFVFTL"


class TestSearchResponse:
    def test_valid_response(self):
        rec = GeneRecord(source="hla", gene_name="HLA-A", sequence="ATCG")
        resp = SearchResponse(total=1, records=[rec])
        assert resp.total == 1
        assert len(resp.records) == 1

    def test_empty_records(self):
        resp = SearchResponse(total=0, records=[])
        assert resp.total == 0
        assert resp.records == []


class TestNLQueryRequest:
    def test_minimal(self):
        req = NLQueryRequest(query="find HLA-A sequences")
        assert req.query == "find HLA-A sequences"
        assert req.limit == 50

    def test_custom_limit(self):
        req = NLQueryRequest(query="test", limit=10)
        assert req.limit == 10

    def test_empty_query_accepted(self):
        req = NLQueryRequest(query="")
        assert req.query == ""

    def test_limit_below_min_raises(self):
        with pytest.raises(ValidationError):
            NLQueryRequest(query="test", limit=0)


class TestParseQueryResult:
    def test_all_none_defaults(self):
        r = ParseQueryResult()
        assert r.source is None
        assert r.species is None
        assert r.gene_name is None
        assert r.region is None
        assert r.sequence_contains is None
        assert r.antigen_epitope is None

    def test_valid_source(self):
        r = ParseQueryResult(source="iedb")
        assert r.source == "iedb"

    def test_invalid_source_raises(self):
        with pytest.raises(ValidationError):
            ParseQueryResult(source="bad")

    def test_all_fields_set(self):
        r = ParseQueryResult(
            source="vdjdb",
            species="human",
            gene_name="TRBV19",
            region="v-region",
            sequence_contains="ATCGATCG",
            antigen_epitope="GILGFVFTL",
        )
        assert r.source == "vdjdb"
        assert r.species == "human"
        assert r.gene_name == "TRBV19"


class TestIngestResponse:
    def test_valid(self):
        r = IngestResponse(inserted=5, source="vdjdb", species="human")
        assert r.inserted == 5
        assert r.source == "vdjdb"
        assert r.species == "human"

    def test_invalid_source_raises(self):
        with pytest.raises(ValidationError):
            IngestResponse(inserted=0, source="bad", species="human")

    def test_invalid_species_raises(self):
        with pytest.raises(ValidationError):
            IngestResponse(inserted=0, source="vdjdb", species="rat")


class TestCDRPredictResponse:
    def test_full_response(self):
        r = CDRPredictResponse(
            v_gene="TRBV19",
            species="human",
            allele="TRBV19*01",
            cdr1_aa="MNHEYMSW",
            cdr2_aa="SVGAGITD",
            cdr1_nt="ATGATGCAT",
            cdr2_nt="TCTGTTGGT",
        )
        assert r.v_gene == "TRBV19"
        assert r.allele == "TRBV19*01"

    def test_optional_fields_none(self):
        r = CDRPredictResponse(
            v_gene="TRBV9999",
            species="other",
            allele=None,
            cdr1_aa=None,
            cdr2_aa=None,
            cdr1_nt=None,
            cdr2_nt=None,
        )
        assert r.cdr1_aa is None
        assert r.allele is None

    def test_invalid_species_raises(self):
        with pytest.raises(ValidationError):
            CDRPredictResponse(
                v_gene="TRBV19",
                species="dog",
                allele=None,
                cdr1_aa=None,
                cdr2_aa=None,
                cdr1_nt=None,
                cdr2_nt=None,
            )


class TestIEDBHit:
    def test_minimal(self):
        h = IEDBHit(epitope_sequence="GILGFVFTL")
        assert h.epitope_sequence == "GILGFVFTL"
        assert h.mhc_allele is None
        assert h.mhc_class is None
        assert h.source_organism is None
        assert h.antigen_name is None
        assert h.assay_type is None
        assert h.effector_cell_type is None
        assert h.qualitative_measure is None

    def test_full(self):
        h = IEDBHit(
            epitope_sequence="GILGFVFTL",
            mhc_allele="HLA-A*02:01",
            mhc_class="I",
            source_organism="Influenza A virus",
            antigen_name="Matrix protein 1",
            assay_type="cytokine release",
            effector_cell_type="CD8+",
            qualitative_measure="Positive",
        )
        assert h.mhc_allele == "HLA-A*02:01"
        assert h.qualitative_measure == "Positive"

    def test_missing_epitope_raises(self):
        with pytest.raises(ValidationError):
            IEDBHit()


class TestGeneRecordIedbHits:
    def test_iedb_hits_defaults_to_none(self):
        r = GeneRecord(source="hla", gene_name="HLA-A", sequence="ATCG")
        assert r.iedb_hits is None

    def test_iedb_hits_accepts_list(self):
        hit = IEDBHit(epitope_sequence="GILGFVFTL", mhc_allele="HLA-A*02:01")
        r = GeneRecord(
            source="vdjdb", gene_name="TRBV19", sequence="CASSIRSSYEQYF",
            iedb_hits=[hit],
        )
        assert len(r.iedb_hits) == 1
        assert r.iedb_hits[0].mhc_allele == "HLA-A*02:01"

    def test_iedb_hits_accepts_empty_list(self):
        r = GeneRecord(source="vdjdb", gene_name="TRBV19", sequence="CASSIRSSYEQYF",
                       iedb_hits=[])
        assert r.iedb_hits == []


class TestGeneRecordScoringFields:
    def test_scoring_fields_default_to_none(self):
        r = GeneRecord(source="vdjdb", gene_name="TRBV19", sequence="CASSIRSSYEQYF")
        assert r.batman_score is None
        assert r.pmhc_score is None
        assert r.tcrdist_score is None
        assert r.composite_score is None

    def test_scoring_fields_accept_float(self):
        r = GeneRecord(
            source="vdjdb", gene_name="TRBV19", sequence="CASSIRSSYEQYF",
            batman_score=0.72,
            pmhc_score=0.88,
            tcrdist_score=0.65,
            composite_score=0.74,
        )
        assert r.batman_score == pytest.approx(0.72)
        assert r.composite_score == pytest.approx(0.74)

    def test_scoring_fields_reject_out_of_range(self):
        with pytest.raises(ValidationError):
            GeneRecord(source="vdjdb", gene_name="X", sequence="Y", batman_score=1.5)
        with pytest.raises(ValidationError):
            GeneRecord(source="vdjdb", gene_name="X", sequence="Y", pmhc_score=-0.1)
