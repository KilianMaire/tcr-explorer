"""Tests for cdr_enricher.py."""
from __future__ import annotations

import sys
import pathlib
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from tcr_explorer.cdr_enricher import (
    _translate,
    _gene_to_chain,
    get_cdr1_cdr2,
    _CDR1_START,
    _CDR1_END,
    _CDR2_START,
    _CDR2_END,
)


class TestTranslate:
    def test_start_codon(self):
        assert _translate("ATG") == "M"

    def test_stop_codon(self):
        assert _translate("TAA") == "*"

    def test_multiple_codons(self):
        # ATG (M) + GCT (A) + TAT (Y)
        assert _translate("ATGGCTTAT") == "MAY"

    def test_unknown_codon_returns_question_mark(self):
        # NNN is not in the codon table
        assert _translate("NNN") == "?"

    def test_empty_string(self):
        assert _translate("") == ""

    def test_incomplete_codon_ignored(self):
        # 4 nt: first codon ATG → M; last 1 nt ignored
        assert _translate("ATGG") == "M"

    def test_case_insensitive(self):
        assert _translate("atg") == "M"

    def test_phenylalanine_codons(self):
        assert _translate("TTTTTC") == "FF"

    def test_leucine_codons(self):
        assert _translate("TTACTA") == "LL"


class TestGeneToChain:
    def test_trbv_to_trb(self):
        assert _gene_to_chain("TRBV19") == "TRB"

    def test_trav_to_tra(self):
        assert _gene_to_chain("TRAV27") == "TRA"

    def test_trgv_to_trg(self):
        assert _gene_to_chain("TRGV2") == "TRG"

    def test_trdv_to_trd(self):
        assert _gene_to_chain("TRDV1") == "TRD"

    def test_unknown_gene_defaults_to_trb(self):
        assert _gene_to_chain("UNKNOWNGENE") == "TRB"

    def test_allele_stripped(self):
        # Gene with allele suffix — the *01 part should be stripped before lookup
        assert _gene_to_chain("TRBV19*01") == "TRB"

    def test_lowercase_input(self):
        # _gene_to_chain uppercases internally
        assert _gene_to_chain("trbv19") == "TRB"


class TestCdrPositions:
    def test_cdr1_positions(self):
        assert _CDR1_START == 27
        assert _CDR1_END == 38

    def test_cdr2_positions(self):
        assert _CDR2_START == 56
        assert _CDR2_END == 65


class TestGetCdr1Cdr2MissingData:
    """Tests for behavior when stitchr data is unavailable."""

    def test_returns_none_when_gene_not_found(self):
        # A gene that definitely won't be in any data
        result = get_cdr1_cdr2("TRBV9999")
        assert result["cdr1_aa"] is None
        assert result["cdr2_aa"] is None
        assert result["allele"] is None
        assert result["cdr1_nt"] is None
        assert result["cdr2_nt"] is None

    def test_returns_dict_with_expected_keys(self):
        result = get_cdr1_cdr2("TRBV9999")
        assert set(result.keys()) == {"cdr1_aa", "cdr2_aa", "allele", "cdr1_nt", "cdr2_nt"}

    def test_missing_stitchr_dir_returns_nones(self):
        # Patch _stitchr_data_dir to return None (simulates stitchr not installed)
        with patch("tcr_explorer.cdr_enricher._stitchr_data_dir", return_value=None):
            # Clear cache so mock takes effect
            from tcr_explorer.cdr_enricher import _cached_v_map
            _cached_v_map.cache_clear()
            result = get_cdr1_cdr2("TRBV19", "human")
            _cached_v_map.cache_clear()
        assert result["cdr1_aa"] is None
        assert result["cdr2_aa"] is None


class TestGetCdr1Cdr2WithMockData:
    """Tests using injected mock V-gene sequences."""

    # Synthetic 300-nt V-region covering CDR1 (pos 27-38) and CDR2 (pos 56-65).
    # The sequence must be long enough: CDR2_END=65 aa → 65*3=195 nt minimum.
    # We build exactly 200 codons (600 nt) using ATG (M) for every codon.
    _MOCK_NT = "ATG" * 200  # 600 nt → 200 Met residues

    def _make_mock_v_map(self):
        return {"TRBV1": self._MOCK_NT}

    def test_cdr1_extracted_as_met_string(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        expected_len = _CDR1_END - _CDR1_START + 1  # 12 aa
        assert result["cdr1_aa"] == "M" * expected_len

    def test_cdr2_extracted_as_met_string(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        expected_len = _CDR2_END - _CDR2_START + 1  # 10 aa
        assert result["cdr2_aa"] == "M" * expected_len

    def test_cdr1_nt_length(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        expected_nt_len = (_CDR1_END - _CDR1_START + 1) * 3  # 36 nt
        assert len(result["cdr1_nt"]) == expected_nt_len

    def test_cdr2_nt_length(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        expected_nt_len = (_CDR2_END - _CDR2_START + 1) * 3  # 30 nt
        assert len(result["cdr2_nt"]) == expected_nt_len

    def test_allele_set_to_star01(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        assert result["allele"] == "TRBV1*01"

    def test_allele_suffix_stripped_before_lookup(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1*02", "human")
        # Should find TRBV1 by stripping *02
        assert result["cdr1_aa"] is not None

    def test_stop_codon_stripped_from_cdr(self):
        # Build a sequence where CDR1 ends in a stop codon (TAA = *)
        # CDR1 aa positions 27-38 → nt 78-114 (0-indexed: [78:114])
        # Make all ATG then overwrite CDR1_END codon with TAA
        nt = list("ATG" * 200)
        stop_pos = (_CDR1_END - 1) * 3  # 0-indexed start of codon at aa pos 38
        nt[stop_pos : stop_pos + 3] = list("TAA")
        mock_map = {"TRBV1": "".join(nt)}
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        assert result["cdr1_aa"] is not None
        assert not result["cdr1_aa"].endswith("*")

    def test_species_human_mapped(self):
        """human species maps to HUMAN stitchr dir; no crash."""
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "human")
        assert isinstance(result, dict)

    def test_species_mouse_mapped(self):
        mock_map = self._make_mock_v_map()
        with patch("tcr_explorer.cdr_enricher._cached_v_map", return_value=mock_map):
            result = get_cdr1_cdr2("TRBV1", "mouse")
        assert isinstance(result, dict)
