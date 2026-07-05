"""Tests for Settings configuration."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from tcr_explorer.config import Settings


class TestServerSettings:
    def test_hla_server_url_default(self):
        s = Settings()
        assert s.hla_server_url == "http://127.0.0.1:8101"

    def test_mhc_server_url_default(self):
        s = Settings()
        assert s.mhc_server_url == "http://127.0.0.1:8105"

    def test_hla_server_url_from_env(self, monkeypatch):
        monkeypatch.setenv("HLA_SERVER_URL", "http://hla:8101")
        import importlib, tcr_explorer.config as cfg
        importlib.reload(cfg)
        s = cfg.Settings()
        assert s.hla_server_url == "http://hla:8101"
