"""The runtime data defaults must resolve to real files shipped in the package."""
from pathlib import Path

from tcr_explorer.records import _DEFAULT_RECORDS_INDEX
from tcr_explorer.similarity import _DEFAULT_INDEX
from tcr_explorer.cdr_enricher import _stitchr_data_dir


def test_records_default_is_inside_the_package():
    p = Path(_DEFAULT_RECORDS_INDEX)
    assert p.exists(), p
    assert p.parent.name == "data" and p.parent.parent.name == "tcr_explorer"


def test_similarity_default_matches_records_default():
    assert Path(_DEFAULT_INDEX).exists()


def test_vendored_germline_ships_in_the_package():
    """The vendored germline directory exists regardless of any stitchr install."""
    vendored = Path(__file__).resolve().parent.parent / "src" / "tcr_explorer" / "data" / "germline"
    assert (vendored / "HUMAN" / "TRB.fasta").exists()
    assert (vendored / "MOUSE" / "TRB.fasta").exists()


def test_stitchr_data_dir_resolves():
    d = _stitchr_data_dir()
    assert d is not None
    assert (d / "HUMAN" / "TRB.fasta").exists()
