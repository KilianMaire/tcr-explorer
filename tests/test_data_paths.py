import json

from tcr_explorer import data_paths


def test_data_dir_honors_env(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    assert data_paths.data_dir() == tmp_path
    assert data_paths.raw_dir() == tmp_path / "raw"
    assert data_paths.records_index_path() == tmp_path / "records_index.parquet"


def test_data_present_false_when_absent(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    assert data_paths.data_present() is False


def test_is_stale_from_meta(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    (tmp_path / "records_index.parquet").write_bytes(b"x")
    (tmp_path / "records_index.meta.json").write_text(json.dumps({"built_at": "2000-01-01"}))
    assert data_paths.data_present() is True
    assert data_paths.is_stale(30) is True
    assert (data_paths.index_age_days() or 0) > 1000


def test_default_dir_uses_platformdirs(monkeypatch):
    monkeypatch.delenv("TCR_EXPLORER_DATA", raising=False)
    d = data_paths.data_dir()
    assert "tcr-explorer" in str(d).lower()
