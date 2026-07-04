"""The package must ship only the static blosum matrix, no vendored datasets."""
from pathlib import Path

import tcr_explorer


def test_package_data_dir_has_only_blosum():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    names = sorted(p.name for p in data.iterdir()) if data.exists() else []
    assert names == ["blosum62.json"], names


def test_no_vendored_germline_or_parquet():
    data = Path(tcr_explorer.__file__).resolve().parent / "data"
    assert not (data / "records_index.parquet").exists()
    assert not (data / "germline").exists()
