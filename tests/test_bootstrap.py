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


def test_refresh_all_downloads_fail_does_not_write_index(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    for name in ("download_vdjdb", "download_iedb", "download_mcpas", "download_tcr3d"):
        monkeypatch.setattr(data_sources, name,
                            lambda raw, s=name: data_sources.DownloadResult(s, False, error="net"))
    summary = bootstrap.refresh()
    assert summary["built"] is False
    assert not data_paths.data_present()  # good/absent index not clobbered by an empty build


def test_refresh_partial_keeps_existing_index(monkeypatch, tmp_path):
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    # simulate a pre-existing good index
    data_paths.records_index_path().parent.mkdir(parents=True, exist_ok=True)
    data_paths.records_index_path().write_bytes(b"PRIOR")
    monkeypatch.setattr(data_sources, "download_vdjdb",
                        lambda raw: data_sources.DownloadResult("vdjdb", True, bytes=1))
    for name in ("download_iedb", "download_mcpas", "download_tcr3d"):
        monkeypatch.setattr(data_sources, name,
                            lambda raw, s=name: data_sources.DownloadResult(s, False, error="net"))
    summary = bootstrap.refresh()
    assert summary["built"] is False and "partial" in summary["reason"]
    assert data_paths.records_index_path().read_bytes() == b"PRIOR"  # untouched


def test_stitchrdl_cmd_is_not_the_no_op_module_form():
    # `python -m Stitchr.stitchrdl` is a silent no-op: stitchrdl.py has no
    # `if __name__ == "__main__"` guard, so importing it never runs main() and
    # nothing downloads. Guard against regressing to that form.
    cmd = bootstrap._stitchrdl_cmd()
    assert "-m" not in cmd
    assert not (len(cmd) >= 2 and cmd[1] == "-m"), "must not use the no-op -m form"


def test_run_stitchrdl_does_not_trust_exit_code(monkeypatch, tmp_path):
    # stitchrdl exits 0 even when it downloaded nothing. _run_stitchrdl must
    # report success only when the species FASTA actually landed in stitchr's
    # native Data dir, not on exit 0.
    from tcr_explorer import cdr_enricher
    monkeypatch.setattr(bootstrap.subprocess, "run",
                        lambda *a, **k: None)  # "succeeds", writes nothing
    empty = tmp_path / "Data"
    empty.mkdir()
    monkeypatch.setattr(cdr_enricher, "_native_stitchr_data_dir", lambda: empty)
    out = bootstrap._run_stitchrdl()
    assert out["human"] is False and out["mouse"] is False


def test_run_stitchrdl_true_when_fasta_present(monkeypatch, tmp_path):
    from tcr_explorer import cdr_enricher
    data = tmp_path / "Data"
    for sp in ("HUMAN", "MOUSE"):
        (data / sp).mkdir(parents=True)
        (data / sp / "TRB.fasta").write_text(">x\nACGT\n")
    monkeypatch.setattr(bootstrap.subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(cdr_enricher, "_native_stitchr_data_dir", lambda: data)
    out = bootstrap._run_stitchrdl()
    assert out["human"] is True and out["mouse"] is True


def test_packaged_germline_is_bundled():
    # The germline ships with the package (CC BY 4.0), so germline features work
    # offline and without IMGT.
    from tcr_explorer.cdr_enricher import _packaged_germline_dir
    d = _packaged_germline_dir()
    for sp in ("HUMAN", "MOUSE"):
        trb = d / sp / "TRB.fasta"
        assert trb.exists() and trb.stat().st_size > 0, f"missing bundled {sp}/TRB.fasta"


def test_stitchr_data_dir_resolves_to_bundled_offline(monkeypatch, tmp_path):
    # With no user germline and no native stitchr install, the resolver still
    # finds germline: the bundled copy.
    from tcr_explorer import cdr_enricher
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))  # no germline/ subdir
    monkeypatch.setattr(cdr_enricher, "_native_stitchr_data_dir", lambda: None)
    resolved = cdr_enricher._stitchr_data_dir()
    assert resolved == cdr_enricher._packaged_germline_dir()


def test_user_germline_overrides_bundled(monkeypatch, tmp_path):
    from tcr_explorer import cdr_enricher, data_paths
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    g = data_paths.germline_dir()
    (g / "HUMAN").mkdir(parents=True)
    (g / "HUMAN" / "TRB.fasta").write_text(">x\nACGT\n")
    (g / "MOUSE").mkdir(parents=True)
    (g / "MOUSE" / "TRB.fasta").write_text(">x\nACGT\n")
    assert cdr_enricher._stitchr_data_dir() == g


def test_default_refresh_never_touches_imgt(monkeypatch, tmp_path):
    # Germline is bundled, so the default refresh must not invoke stitchrdl.
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))

    def boom():
        raise AssertionError("default refresh must not run stitchrdl")

    monkeypatch.setattr(bootstrap, "_run_stitchrdl", boom)

    def fake_mcpas(raw):
        raw.mkdir(parents=True, exist_ok=True)
        (raw / "mcpas.csv").write_text(
            "CDR3.alpha.aa,CDR3.beta.aa,Species,TRBV,TRBJ\n,CASSF,Human,TRBV19,TRBJ2-1\n")
        return data_sources.DownloadResult("mcpas", True, str(raw / "mcpas.csv"), 1)

    for name in ("download_vdjdb", "download_iedb", "download_tcr3d"):
        monkeypatch.setattr(data_sources, name,
                            lambda raw, s=name: data_sources.DownloadResult(s, False, error="skip"))
    monkeypatch.setattr(data_sources, "download_mcpas", fake_mcpas)

    summary = bootstrap.refresh()
    assert summary["built"] is True
    assert summary["germline"]["source"] == "bundled"
    assert summary["germline"]["release"] != "unknown"


def test_load_records_index_does_not_cache_absent(monkeypatch, tmp_path):
    from tcr_explorer import records as R
    monkeypatch.setenv("TCR_EXPLORER_DATA", str(tmp_path))
    monkeypatch.delenv("RECORDS_INDEX_PATH", raising=False)
    R.load_records_index.cache_clear()
    p = str(tmp_path / "records_index.parquet")
    assert R.load_records_index(p) is None
    import pandas as pd
    from tcr_explorer.records_build import SCHEMA_COLUMNS
    pd.DataFrame({c: ["x"] for c in SCHEMA_COLUMNS}).to_parquet(p, index=False)
    # after the file appears, the same call must see it (absent result was not cached)
    assert R.load_records_index(p) is not None
