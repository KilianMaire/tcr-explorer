"""Agnostic OpenAI-compatible LLM client (optional, graceful).

Talks to any OpenAI-compatible `/chat/completions` endpoint (e.g. LM Studio,
Ollama, vLLM, or a hosted OpenAI-compatible gateway). Disabled by default
unless `LLM_BASE_URL` is set. Never raises: any failure (network error,
malformed response, disabled) results in `None` so callers can fall back to
non-LLM behavior.
"""
from __future__ import annotations

import json
import os
from typing import Optional

import httpx


def _base_url() -> str:
    return os.environ.get("LLM_BASE_URL", "").rstrip("/")


def _model() -> str:
    return os.environ.get("LLM_MODEL", "local-model")


def _timeout() -> float:
    return float(os.environ.get("LLM_TIMEOUT", "8.0"))


def llm_available() -> bool:
    return bool(_base_url()) and os.environ.get("LLM_ENABLE", "true").lower() != "false"


def llm_json(system: str, user: str) -> Optional[dict]:
    if not llm_available():
        return None
    payload = {
        "model": _model(),
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "temperature": 0,
        "response_format": {"type": "json_object"},
    }
    headers = {}
    key = os.environ.get("LLM_API_KEY")
    if key:
        headers["Authorization"] = f"Bearer {key}"
    try:
        with httpx.Client(timeout=_timeout()) as client:
            r = client.post(f"{_base_url()}/chat/completions", json=payload, headers=headers)
            r.raise_for_status()
            content = r.json()["choices"][0]["message"]["content"]
            parsed = json.loads(content)
            return parsed if isinstance(parsed, dict) else None
    except Exception:
        return None
