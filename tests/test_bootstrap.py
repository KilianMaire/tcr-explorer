import pytest

from tcr_explorer import bootstrap, data_paths, data_sources


def test_ensure_ready_raises_when_absent(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    with pytest.raises(bootstrap.DataNotReadyError) as e:
        bootstrap.ensure_ready()
    assert "tcr-explorer-refresh" in str(e.value)


def test_refresh_builds_index_from_downloads(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))

    def fake_mcpas(raw):
        raw.mkdir(parents=True, exist_ok=True)
        (raw / "mcpas.csv").write_text(
            "CDR3.alpha.aa,CDR3.beta.aa,Species,TRBV,TRBJ\n,CASSF,Human,TRBV19,TRBJ2-1\n")
        return data_sources.DownloadResult("mcpas", True, str(raw / "mcpas.csv"), 1)

    monkeypatch.setattr(data_sources, "download_vdjdb",
                        lambda raw: data_sources.DownloadResult("vdjdb", False, error="skip"))
    monkeypatch.setattr(data_sources, "download_iedb",
                        lambda raw: data_sources.DownloadResult("iedb", False, error="skip"))
    monkeypatch.setattr(data_sources, "download_tcr3d",
                        lambda raw: data_sources.DownloadResult("tcr3d", False, error="skip"))
    monkeypatch.setattr(data_sources, "download_mcpas", fake_mcpas)
    monkeypatch.setattr(bootstrap, "_run_stitchrdl", lambda: {"ok": True})

    summary = bootstrap.refresh()
    assert data_paths.data_present()
    assert summary["rows_total"] >= 1
    assert summary["sources"]["mcpas"]["ok"] is True
    bootstrap.ensure_ready()  # now does not raise
