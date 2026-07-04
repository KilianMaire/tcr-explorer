import io
import zipfile

import httpx

from tcr_explorer import data_sources


def _zip_bytes(name: str, content: bytes) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr(name, content)
    return buf.getvalue()


def test_download_mcpas_strips_bom_and_writes(monkeypatch, tmp_path):
    body = "﻿CDR3.alpha.aa,CDR3.beta.aa\nCAA,CASS\n".encode("utf-8")

    def handler(request):
        return httpx.Response(200, content=body)

    monkeypatch.setattr(data_sources, "_client",
                        lambda: httpx.Client(transport=httpx.MockTransport(handler)))
    res = data_sources.download_mcpas(tmp_path)
    assert res.ok
    out = (tmp_path / "mcpas.csv").read_text(encoding="utf-8")
    assert out.startswith("CDR3.alpha.aa")  # BOM stripped


def test_download_iedb_rejects_non_zip(monkeypatch, tmp_path):
    def handler(request):
        return httpx.Response(200, content=b"<html>not a zip</html>")

    monkeypatch.setattr(data_sources, "_client",
                        lambda: httpx.Client(transport=httpx.MockTransport(handler)))
    res = data_sources.download_iedb(tmp_path)
    assert not res.ok
    assert not (tmp_path / "iedb_receptor_full_v3.zip").exists()


def test_download_iedb_accepts_zip(monkeypatch, tmp_path):
    payload = _zip_bytes("tcr_full_v3.csv", b"a,b\n1,2\n")

    def handler(request):
        return httpx.Response(200, content=payload)

    monkeypatch.setattr(data_sources, "_client",
                        lambda: httpx.Client(transport=httpx.MockTransport(handler)))
    res = data_sources.download_iedb(tmp_path)
    assert res.ok
    assert (tmp_path / "iedb_receptor_full_v3.zip").exists()


def test_download_tcr3d_writes_two_tsvs(monkeypatch, tmp_path):
    def handler(request):
        if "complexes" in str(request.url):
            return httpx.Response(200, content=b"PDB_ID\tTCR_name\n1abc\tx\n")
        return httpx.Response(200, content=b"pdb_id\ttcr_name\n1abc\tx\n")

    monkeypatch.setattr(data_sources, "_client",
                        lambda: httpx.Client(transport=httpx.MockTransport(handler)))
    res = data_sources.download_tcr3d(tmp_path)
    assert res.ok
    assert (tmp_path / "tcr3d_complexes.tsv").exists()
    assert (tmp_path / "tcr3d_chain.tsv").exists()
