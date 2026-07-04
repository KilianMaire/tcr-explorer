"""The package ships the static blosum matrix and the IMGT germline (CC BY 4.0),
but NO record datasets (those are download-on-first-run)."""
from pathlib import Path

import tcr_explorer


def test_package_data_dir_has_blosum_and_germline():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    names = sorted(p.name for p in data.iterdir()) if data.exists() else []
    assert names == ["blosum62.json", "germline"], names


def test_germline_is_vendored_records_are_not():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    # germline IS vendored (IMGT CC BY 4.0), with attribution.
    assert (data / "germline" / "HUMAN" / "TRB.fasta").exists()
    assert (data / "germline" / "MOUSE" / "TRB.fasta").exists()
    assert (data / "germline" / "ATTRIBUTION.md").exists()
    # record datasets are NOT vendored (VDJdb / McPAS licensing): downloaded only.
    assert not (data / "records_index.parquet").exists()
