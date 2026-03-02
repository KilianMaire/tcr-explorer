"""Tests for src/imgt_app/frontend.py — data-shaping logic only (no UI rendering)."""
from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Ensure src/ is on the path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pd = pytest.importorskip("pandas")

# ---------------------------------------------------------------------------
# Stub out streamlit so frontend.py can be imported without it installed
# ---------------------------------------------------------------------------

def _stub_streamlit() -> None:
    """Insert a minimal streamlit stub into sys.modules if not already present."""
    if "streamlit" in sys.modules:
        return
    st_stub = types.ModuleType("streamlit")
    # Add every attribute used at module level in frontend.py as a no-op mock.
    for attr in (
        "set_page_config", "title", "caption", "header", "subheader", "info",
        "error", "success", "spinner", "sidebar", "button", "text_input",
        "selectbox", "slider", "radio", "dataframe", "download_button",
        "columns", "metric", "text_area",
    ):
        setattr(st_stub, attr, MagicMock())
    # sidebar needs the same attributes
    st_stub.sidebar = MagicMock()
    sys.modules["streamlit"] = st_stub


_stub_streamlit()

from imgt_app.frontend import _records_to_df  # noqa: E402  (must come after stub)


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

SAMPLE_RECORDS = [
    {
        "source": "vdjdb",
        "species": "human",
        "gene_name": "TRBV12-3",
        "allele_name": None,
        "region": "cdr3",
        "sequence": "CASSQDRGNTGELFF",
        "antigen_epitope": None,
        "metadata": {
            "cdr1_aa": "SGHDN",
            "cdr2_aa": "FNNNVP",
            "antigen_epitope": "GILGFVFTL",
            "mhc_a": "HLA-A*02:01",
            "mhc_class": "MHCI",
            "j_segm": "TRBJ2-2",
        },
    },
    {
        "source": "hla",
        "species": "human",
        "gene_name": "HLA-A",
        "allele_name": "HLA-A*02:01",
        "region": "exon2",
        "sequence": "ATGGCC",
        "antigen_epitope": None,
        "metadata": {},
    },
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_records_to_df_shape():
    df = _records_to_df(SAMPLE_RECORDS)
    assert len(df) == 2
    assert "CDR3" in df.columns
    assert "CDR1" in df.columns
    assert "CDR2" in df.columns
    assert "antigen_epitope" in df.columns


def test_records_to_df_vdjdb_row():
    df = _records_to_df(SAMPLE_RECORDS)
    row = df[df["source"] == "vdjdb"].iloc[0]
    assert row["CDR3"] == "CASSQDRGNTGELFF"
    assert row["CDR1"] == "SGHDN"
    assert row["CDR2"] == "FNNNVP"
    assert row["antigen_epitope"] == "GILGFVFTL"
    assert row["gene_name"] == "TRBV12-3"
    assert row["j_segm"] == "TRBJ2-2"


def test_records_to_df_hla_row():
    df = _records_to_df(SAMPLE_RECORDS)
    row = df[df["source"] == "hla"].iloc[0]
    assert row["allele_name"] == "HLA-A*02:01"
    assert row["CDR1"] == ""
    assert row["CDR2"] == ""
    assert row["antigen_epitope"] == ""


def test_records_to_df_empty():
    df = _records_to_df([])
    assert len(df) == 0


def test_records_to_df_antigen_from_top_level():
    """antigen_epitope on the record itself should take priority over metadata."""
    rec = {
        "source": "vdjdb",
        "species": "human",
        "gene_name": "TRBV7",
        "allele_name": None,
        "region": "cdr3",
        "sequence": "CASSX",
        "antigen_epitope": "DIRECT",
        "metadata": {"antigen_epitope": "META_EPITOPE"},
    }
    df = _records_to_df([rec])
    assert df.iloc[0]["antigen_epitope"] == "DIRECT"


def test_records_to_df_missing_metadata_key():
    """Records without CDR keys in metadata should produce empty strings."""
    rec = {
        "source": "tcr",
        "species": "mouse",
        "gene_name": "TRAV14D-3",
        "allele_name": None,
        "region": "v",
        "sequence": "ATGCTA",
        "antigen_epitope": None,
        "metadata": {},
    }
    df = _records_to_df([rec])
    assert df.iloc[0]["CDR1"] == ""
    assert df.iloc[0]["CDR2"] == ""
    assert df.iloc[0]["mhc_a"] == ""
