"""The package ships the static blosum matrix and the IMGT germline (CC BY 4.0),
but NO record datasets (those are download-on-first-run)."""
from pathlib import Path

import tcr_explorer


def test_package_data_dir_has_blosum_germline_and_tcrdist():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    names = sorted(p.name for p in data.iterdir()) if data.exists() else []
    assert names == ["blosum62.json", "germline", "tcrdist"], names


def test_tcrdist_reference_table_is_vendored():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    # tcrdist V/J CDR reference table IS vendored (tcrdist3, MIT; CDRs from IMGT CC BY 4.0).
    assert (data / "tcrdist" / "alphabeta_gammadelta_db.tsv").exists()
    assert (data / "tcrdist" / "ATTRIBUTION.md").exists()


def test_germline_is_vendored_records_are_not():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    # germline IS vendored (IMGT CC BY 4.0), with attribution.
    assert (data / "germline" / "HUMAN" / "TRB.fasta").exists()
    assert (data / "germline" / "MOUSE" / "TRB.fasta").exists()
    assert (data / "germline" / "ATTRIBUTION.md").exists()
    # record datasets are NOT vendored (VDJdb / McPAS licensing): downloaded only.
    assert not (data / "records_index.parquet").exists()
