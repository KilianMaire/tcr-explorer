"""Tests for heuristic_parse in nl_query.py."""
from __future__ import annotations

import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from imgt_app.nl_query import heuristic_parse


class TestSourceDetection:
    def test_hla_keyword(self):
        result = heuristic_parse("Show me HLA-A alleles")
        assert result.source == "hla"

    def test_iedb_keyword(self):
        result = heuristic_parse("Search IEDB for epitopes")
        assert result.source == "iedb"

    def test_immune_epitope_database(self):
        result = heuristic_parse("query the immune epitope database")
        assert result.source == "iedb"

    def test_tcell_assay(self):
        result = heuristic_parse("tcell assay results for influenza")
        assert result.source == "iedb"

    def test_mhc_ligand(self):
        result = heuristic_parse("mhc ligand binding data")
        assert result.source == "iedb"

    def test_effector_cell(self):
        result = heuristic_parse("data from effector cell experiments")
        assert result.source == "iedb"

    def test_vdjdb_keyword(self):
        result = heuristic_parse("vdjdb CDR3 sequences")
        assert result.source == "vdjdb"

    def test_cdr3_keyword(self):
        result = heuristic_parse("find CDR3 sequences for influenza")
        assert result.source == "vdjdb"

    def test_epitope_keyword(self):
        result = heuristic_parse("what epitope does this TCR bind")
        assert result.source == "vdjdb"

    def test_antigen_keyword(self):
        result = heuristic_parse("antigen specificity of TRBV19")
        assert result.source == "vdjdb"

    def test_tcr_keyword(self):
        result = heuristic_parse("show TCR V gene sequences")
        assert result.source == "tcr"

    def test_t_cell_receptor(self):
        result = heuristic_parse("t-cell receptor beta chain")
        assert result.source == "tcr"

    def test_no_keywords_returns_none(self):
        result = heuristic_parse("list all sequences")
        assert result.source is None

    def test_hla_takes_priority_over_vdjdb(self):
        # "hla" appears in query along with "epitope" — hla is checked first
        result = heuristic_parse("HLA epitope binding")
        assert result.source == "hla"

    # MHC species detection
    def test_mamu_keyword(self):
        result = heuristic_parse("Mamu-A1 alleles in rhesus macaque")
        assert result.source == "mhc"

    def test_bola_keyword(self):
        result = heuristic_parse("BoLA-2 bovine MHC class I")
        assert result.source == "mhc"

    def test_patr_keyword(self):
        result = heuristic_parse("Patr-A chimpanzee alleles")
        assert result.source == "mhc"

    def test_dla_keyword(self):
        result = heuristic_parse("DLA-88 canine alleles")
        assert result.source == "mhc"

    def test_sla_keyword(self):
        result = heuristic_parse("SLA-1 swine MHC")
        assert result.source == "mhc"

    def test_rt1_keyword(self):
        result = heuristic_parse("RT1 rat MHC alleles")
        assert result.source == "mhc"

    def test_macaque_keyword(self):
        result = heuristic_parse("macaque MHC class I sequences")
        assert result.source == "mhc"

    def test_bovine_keyword(self):
        result = heuristic_parse("bovine MHC allele variants")
        assert result.source == "mhc"

    def test_rhesus_keyword(self):
        result = heuristic_parse("rhesus monkey MHC alleles")
        assert result.source == "mhc"

    def test_ipd_mhc_keyword(self):
        result = heuristic_parse("search ipd-mhc database")
        assert result.source == "mhc"

    def test_hla_takes_priority_over_mhc(self):
        # "hla" is checked before MHC keywords
        result = heuristic_parse("HLA macaque alleles")
        assert result.source == "hla"


class TestSpeciesDetection:
    def test_human_keyword(self):
        result = heuristic_parse("human TRBV19 sequences")
        assert result.species == "human"

    def test_homo_sapiens(self):
        result = heuristic_parse("Homo sapiens HLA alleles")
        assert result.species == "human"

    def test_mouse_keyword(self):
        result = heuristic_parse("mouse TCR sequences")
        assert result.species == "mouse"

    def test_mus_musculus(self):
        result = heuristic_parse("Mus musculus immune repertoire")
        assert result.species == "mouse"

    def test_no_species_returns_none(self):
        result = heuristic_parse("TRBV19 sequences")
        assert result.species is None


class TestRegionDetection:
    def test_v_region(self):
        result = heuristic_parse("find the v-region of TRBV19")
        assert result.region == "v-region"

    def test_j_region(self):
        result = heuristic_parse("j-region of TRAJ33")
        assert result.region == "j-region"

    def test_no_region_returns_none(self):
        result = heuristic_parse("TRBV19 human")
        assert result.region is None


class TestGeneNameExtraction:
    def test_extracts_trbv_gene(self):
        result = heuristic_parse("show TRBV19 sequences")
        assert result.gene_name == "TRBV19"

    def test_extracts_hla_allele(self):
        result = heuristic_parse("search HLA-A*02:01 alleles")
        # gene_name regex finds first all-caps token; HLA will match
        assert result.gene_name is not None

    def test_no_gene_name_returns_none(self):
        result = heuristic_parse("show all human sequences")
        assert result.gene_name is None


class TestSequenceExtraction:
    def test_extracts_dna_sequence(self):
        result = heuristic_parse("sequence ATCGATCGATCG in human")
        assert result.sequence_contains == "ATCGATCGATCG"

    def test_short_sequence_not_extracted(self):
        # fewer than 8 nt — should not be captured
        result = heuristic_parse("ATCG fragment")
        assert result.sequence_contains is None


class TestAntigenEpitopeExtraction:
    def test_extracts_peptide(self):
        # GILGFVFTL is 9 aa — valid epitope pattern
        result = heuristic_parse("CDR3 against GILGFVFTL epitope")
        assert result.antigen_epitope == "GILGFVFTL"

    def test_short_peptide_not_extracted(self):
        result = heuristic_parse("peptide MGIV")
        assert result.antigen_epitope is None


class TestReturnType:
    def test_returns_parse_query_result(self):
        from imgt_app.models import ParseQueryResult
        result = heuristic_parse("HLA-A sequences")
        assert isinstance(result, ParseQueryResult)

    def test_all_fields_present(self):
        result = heuristic_parse("")
        assert hasattr(result, "source")
        assert hasattr(result, "species")
        assert hasattr(result, "gene_name")
        assert hasattr(result, "region")
        assert hasattr(result, "sequence_contains")
        assert hasattr(result, "antigen_epitope")
