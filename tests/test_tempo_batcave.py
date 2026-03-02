"""Tests for BATCAVE database client."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
from unittest.mock import patch, AsyncMock
from tempo.batcave import BatcaveClient, BatcaveVariant


class TestBatcaveClient:
    def test_parse_variants_from_csv_row(self):
        client = BatcaveClient()
        row = {
            "peptide": "GILGFVFTL",
            "variant": "AILGFVFTL",
            "activation_score": "0.85",
            "mhc_allele": "HLA-A*02:01",
            "mhc_class": "I",
            "position": "1",
            "original_aa": "G",
            "mutant_aa": "A",
        }
        variant = client._parse_variant(row)
        assert isinstance(variant, BatcaveVariant)
        assert variant.variant_peptide == "AILGFVFTL"
        assert variant.activation_score == pytest.approx(0.85)

    def test_filter_by_mhc_class(self):
        client = BatcaveClient()
        variants = [
            BatcaveVariant("GILGFVFTL", "AILGFVFTL", 0.85, "HLA-A*02:01", "I", 1, "G", "A"),
            BatcaveVariant("XXXYYYZZZ", "AAAYYYZZZ", 0.60, "HLA-DRB1*04:01", "II", 1, "X", "A"),
        ]
        filtered = client._filter_variants(variants, mhc_class="I")
        assert len(filtered) == 1
        assert filtered[0].mhc_class == "I"

    def test_filter_by_reference_peptide(self):
        client = BatcaveClient()
        variants = [
            BatcaveVariant("GILGFVFTL", "AILGFVFTL", 0.85, "HLA-A*02:01", "I", 1, "G", "A"),
            BatcaveVariant("LLWNGPMAV", "ALWNGPMAV", 0.70, "HLA-A*02:01", "I", 1, "L", "A"),
        ]
        filtered = client._filter_variants(variants, reference_peptide="GILGFVFTL")
        assert len(filtered) == 1
        assert filtered[0].reference_peptide == "GILGFVFTL"
