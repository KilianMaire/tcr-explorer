from tcr_explorer import data_paths, records


def test_records_default_points_at_user_data_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    monkeypatch.delenv("RECORDS_INDEX_PATH", raising=False)
    assert records._default_records_index_path() == str(tmp_path / "records_index.parquet")
    assert str(data_paths.records_index_path()) == str(tmp_path / "records_index.parquet")


def test_retrieve_records_absent_data_is_clear(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    monkeypatch.delenv("RECORDS_INDEX_PATH", raising=False)
    records.load_records_index.cache_clear()
    # No parquet present: load returns None (absent), not a raw error.
    assert records.load_records_index() is None
