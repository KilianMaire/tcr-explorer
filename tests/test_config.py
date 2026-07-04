"""Tests for Settings configuration."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

import pytest
from tcr_explorer.config import Settings


class TestBatmanSettings:
    def test_batman_server_url_default(self):
        s = Settings()
        assert s.batman_server_url == "http://127.0.0.1:8105"

    def test_batman_server_url_from_env(self, monkeypatch):
        monkeypatch.setenv("BATMAN_SERVER_URL", "http://batman:8105")
        import importlib, tcr_explorer.config as cfg
        importlib.reload(cfg)
        s = cfg.Settings()
        assert s.batman_server_url == "http://batman:8105"

    def test_batman_timeout_default(self):
        s = Settings()
        assert s.batman_timeout == pytest.approx(30.0)

    def test_batman_enable_default_true(self):
        s = Settings()
        assert s.batman_enable is True

    def test_batman_enable_from_env(self, monkeypatch):
        monkeypatch.setenv("BATMAN_ENABLE", "false")
        import importlib, tcr_explorer.config as cfg
        importlib.reload(cfg)
        s = cfg.Settings()
        assert s.batman_enable is False
