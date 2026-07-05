"""The MCP-over-HTTP mount is gated by TCR_EXPLORER_MCP_HTTP.

Local runs and the rest of the suite must see no /mcp mount; a deployed instance
sets the env var to expose the same MCP tools to remote chat clients.
"""
from __future__ import annotations

import importlib


def _has_mcp_mount(app) -> bool:
    return any(getattr(r, "path", None) == "/mcp" for r in app.routes)


def test_no_remote_mcp_mount_by_default():
    from tcr_explorer.api import app
    assert not _has_mcp_mount(app)


def test_remote_mcp_mount_when_enabled(monkeypatch):
    monkeypatch.setenv("TCR_EXPLORER_MCP_HTTP", "1")
    from tcr_explorer import api as api_mod
    importlib.reload(api_mod)
    try:
        assert _has_mcp_mount(api_mod.app)
    finally:
        # Restore the default (unmounted) app so later tests are unaffected.
        monkeypatch.delenv("TCR_EXPLORER_MCP_HTTP", raising=False)
        importlib.reload(api_mod)
