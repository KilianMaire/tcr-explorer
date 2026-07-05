"""Application settings loaded from environment variables."""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class Settings:
    hla_server_url: str = os.getenv("HLA_SERVER_URL", "http://127.0.0.1:8101")
    mhc_server_url: str = os.getenv("MHC_SERVER_URL", "http://127.0.0.1:8105")
    llm_base_url: str = os.getenv("LLM_BASE_URL", "")
    llm_model: str = os.getenv("LLM_MODEL", "local-model")
    llm_api_key: str = os.getenv("LLM_API_KEY", "")
    llm_timeout: float = float(os.getenv("LLM_TIMEOUT", "8.0"))
    llm_enable: bool = os.getenv("LLM_ENABLE", "true").lower() != "false"


settings = Settings()
