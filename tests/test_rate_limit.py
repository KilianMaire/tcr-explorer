"""The public-instance rate limit is gated by TCR_EXPLORER_RATE_LIMIT_PER_MIN.

Off by default (local runs and the suite are unaffected); when set, a per-IP
sliding window returns 429 past the limit.
"""
from __future__ import annotations

import importlib

from fastapi.testclient import TestClient


def test_rate_limit_off_by_default():
    from tcr_explorer.api import app
    client = TestClient(app)
    codes = {client.get("/health").status_code for _ in range(6)}
    assert codes == {200}


def test_rate_limit_returns_429_past_the_limit(monkeypatch):
    monkeypatch.setenv("TCR_EXPLORER_RATE_LIMIT_PER_MIN", "3")
    from tcr_explorer import api as api_mod
    importlib.reload(api_mod)
    try:
        client = TestClient(api_mod.app)
        codes = [client.get("/health").status_code for _ in range(5)]
        assert codes[:3] == [200, 200, 200]
        assert codes[3] == 429 and codes[4] == 429
    finally:
        monkeypatch.delenv("TCR_EXPLORER_RATE_LIMIT_PER_MIN", raising=False)
        importlib.reload(api_mod)
