"""Application settings loaded from environment variables."""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class Settings:
    database_path: str = os.getenv("DATABASE_PATH", "data/imgt.db")
    hla_server_url: str = os.getenv("HLA_SERVER_URL", "http://127.0.0.1:8101")
    tcr_server_url: str = os.getenv("TCR_SERVER_URL", "http://127.0.0.1:8102")
    vdjdb_server_url: str = os.getenv("VDJDB_SERVER_URL", "http://127.0.0.1:8103")
    iedb_server_url: str = os.getenv("IEDB_SERVER_URL", "http://127.0.0.1:8104")
    mhc_server_url: str = os.getenv("MHC_SERVER_URL", "http://127.0.0.1:8105")
    batman_server_url: str = os.getenv("BATMAN_SERVER_URL", "http://127.0.0.1:8105")
    batman_timeout: float = float(os.getenv("BATMAN_TIMEOUT", "30.0"))
    batman_enable: bool = os.getenv("BATMAN_ENABLE", "true").lower() != "false"
    tempo_server_url: str = os.getenv("TEMPO_SERVER_URL", "http://127.0.0.1:8106")
    tempo_timeout: float = float(os.getenv("TEMPO_TIMEOUT", "30.0"))
    tempo_enable: bool = os.getenv("TEMPO_ENABLE", "true").lower() != "false"
    llm_base_url: str = os.getenv("LLM_BASE_URL", "")
    llm_model: str = os.getenv("LLM_MODEL", "local-model")
    llm_api_key: str = os.getenv("LLM_API_KEY", "")
    llm_timeout: float = float(os.getenv("LLM_TIMEOUT", "8.0"))
    llm_enable: bool = os.getenv("LLM_ENABLE", "true").lower() != "false"


settings = Settings()
