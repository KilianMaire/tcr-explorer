"""Tests for IEDB → BATMAN training data extraction."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import pandas as pd
from batman.training_data import iedb_records_to_batman_df, merge_user_data


class TestIedbRecordsToBatmanDf:
    def _make_iedb_records(self):
        return [
            {"sequence": "GILGFVFTL", "metadata": {"qualitative_measure": "Positive"}},
            {"sequence": "AILGFVFTL", "metadata": {"qualitative_measure": "Negative"}},
            {"sequence": "NLVPMVATV", "metadata": {"qualitative_measure": "Positive"}},
        ]

    def test_returns_dataframe_with_required_columns(self):
        df = iedb_records_to_batman_df(
            self._make_iedb_records(), tcr_id="TCR1", index_peptide="GILGFVFTL"
        )
        assert set(df.columns) >= {"tcr", "index", "peptide", "activation"}

    def test_positive_maps_to_2(self):
        records = [{"sequence": "GILGFVFTL", "metadata": {"qualitative_measure": "Positive"}}]
        df = iedb_records_to_batman_df(records, tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert df["activation"].iloc[0] == 2

    def test_negative_maps_to_0(self):
        records = [{"sequence": "AILGFVFTL", "metadata": {"qualitative_measure": "Negative"}}]
        df = iedb_records_to_batman_df(records, tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert df["activation"].iloc[0] == 0

    def test_intermediate_maps_to_1(self):
        records = [{"sequence": "BILGFVFTL", "metadata": {"qualitative_measure": "Intermediate"}}]
        df = iedb_records_to_batman_df(records, tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert df["activation"].iloc[0] == 1

    def test_empty_records_returns_empty_df(self):
        df = iedb_records_to_batman_df([], tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert len(df) == 0

    def test_deduplication_keeps_highest_activation(self):
        records = [
            {"sequence": "GILGFVFTL", "metadata": {"qualitative_measure": "Positive"}},
            {"sequence": "GILGFVFTL", "metadata": {"qualitative_measure": "Negative"}},
        ]
        df = iedb_records_to_batman_df(records, tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert len(df) == 1
        assert df["activation"].iloc[0] == 2


class TestMergeUserData:
    def test_user_data_overrides_iedb(self):
        iedb_df = pd.DataFrame({
            "tcr": ["TCR1"], "index": ["GILGFVFTL"],
            "peptide": ["GILGFVFTL"], "activation": [2],
        })
        user_data = [{"peptide": "GILGFVFTL", "activation": 0}]
        merged = merge_user_data(iedb_df, user_data, tcr_id="TCR1", index_peptide="GILGFVFTL")
        row = merged[merged["peptide"] == "GILGFVFTL"]
        assert row["activation"].iloc[0] == 0

    def test_user_data_adds_new_peptides(self):
        iedb_df = pd.DataFrame({
            "tcr": ["TCR1"], "index": ["GILGFVFTL"],
            "peptide": ["GILGFVFTL"], "activation": [2],
        })
        user_data = [{"peptide": "NLVPMVATV", "activation": 1}]
        merged = merge_user_data(iedb_df, user_data, tcr_id="TCR1", index_peptide="GILGFVFTL")
        assert len(merged) == 2
        assert "NLVPMVATV" in merged["peptide"].values
