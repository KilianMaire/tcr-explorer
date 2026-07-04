"""Resolve the local user data directory and index staleness.

Data lives in a user directory (platformdirs), overridable with TCR_EXPLORER_DATA.
Nothing here downloads; it only computes paths and freshness.
"""
from __future__ import annotations

import json
import os
from datetime import date
from pathlib import Path
from typing import Optional

from platformdirs import user_data_dir

_DEFAULT_MAX_AGE_DAYS = 30.0


def data_dir() -> Path:
    override = os.environ.get("TCR_EXPLORER_DATA")
    return Path(override) if override else Path(user_data_dir("tcr-explorer"))


def raw_dir() -> Path:
    return data_dir() / "raw"


def records_index_path() -> Path:
    return data_dir() / "records_index.parquet"


def meta_path() -> Path:
    return data_dir() / "records_index.meta.json"


def data_present() -> bool:
    return records_index_path().exists()


def index_age_days() -> Optional[float]:
    p = meta_path()
    if not p.exists():
        return None
    try:
        built = json.loads(p.read_text()).get("built_at")
        if not built:
            return None
        return float((date.today() - date.fromisoformat(built)).days)
    except (ValueError, OSError):
        return None


def _max_age() -> float:
    raw = os.environ.get("TCR_EXPLORER_MAX_AGE_DAYS")
    if raw:
        try:
            return float(raw)
        except ValueError:
            pass
    return _DEFAULT_MAX_AGE_DAYS


def is_stale(max_age_days: Optional[float] = None) -> bool:
    age = index_age_days()
    if age is None:
        return False
    return age > (max_age_days if max_age_days is not None else _max_age())
