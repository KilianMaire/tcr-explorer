from imgt_app import mcp_server

def test_get_tcr_dossier_tool_returns_dossier_dict():
    out = mcp_server.get_tcr_dossier("TRBV20-1", species="human")
    assert out["schema_version"] == "1.0"
    assert "summary" in out and "genes" in out
    assert set(out["genes"].keys()) == {"v", "d", "j", "c"}

def test_ask_tcr_tool_returns_routed_dict():
    out = mcp_server.ask_tcr("TRBV20-1")
    assert out["intent"] in ("dossier", "search", "similar")
    assert "llm_used" in out

def test_find_similar_tool_returns_engine(monkeypatch):
    # avoid depending on the full vendored index in unit tests
    import imgt_app.mcp_server as m
    monkeypatch.setattr(m, "find_similar_tcrs_fn",
        lambda *a, **k: ([], "blosum_cdr3", 0, []))
    out = mcp_server.find_similar_tcrs("CASSLGTEAFF", "TRBV20-1", "TRBJ1-1")
    assert out["engine"] == "blosum_cdr3"
    assert "neighbours" in out

def test_mcp_object_exists():
    from mcp.server.fastmcp import FastMCP
    assert isinstance(mcp_server.mcp, FastMCP)
