"""BATMAN model cache: in-memory LRU + disk (pickle) persistence."""
from __future__ import annotations

import hashlib
import logging
import pickle
import threading
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class CachedModel:
    tcr_id: str
    index_peptide: str
    weights: np.ndarray
    aa_matrix: np.ndarray
    trained_at: float = field(default_factory=time.time)


class ModelCache:
    """Two-tier cache: RAM (LRU OrderedDict) + disk (pickle files)."""

    def __init__(self, cache_dir: str | Path = "./batman_cache", max_ram: int = 50) -> None:
        self._dir = Path(cache_dir)
        self._dir.mkdir(parents=True, exist_ok=True)
        self._max_ram = max_ram
        self._ram: OrderedDict[str, CachedModel] = OrderedDict()
        self._lock = threading.Lock()

    def get(self, tcr_id: str) -> Optional[CachedModel]:
        """Return cached model or None. RAM hit promoted to MRU."""
        with self._lock:
            if tcr_id in self._ram:
                self._ram.move_to_end(tcr_id)
                return self._ram[tcr_id]
            path = self._disk_path(tcr_id)
            if path.exists():
                try:
                    with path.open("rb") as f:
                        model: CachedModel = pickle.load(f)
                    self._store_ram(tcr_id, model)
                    return model
                except (pickle.UnpicklingError, EOFError, OSError, ValueError):
                    path.unlink(missing_ok=True)
            return None

    def put(self, model: CachedModel) -> None:
        """Store model in RAM + disk."""
        with self._lock:
            self._store_ram(model.tcr_id, model)
            path = self._disk_path(model.tcr_id)
            try:
                with path.open("wb") as f:
                    pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)
            except OSError as exc:
                logger.warning("Failed to persist model %s to disk: %s", model.tcr_id, exc)

    @property
    def ram_size(self) -> int:
        return len(self._ram)

    def _store_ram(self, tcr_id: str, model: CachedModel) -> None:
        self._ram[tcr_id] = model
        self._ram.move_to_end(tcr_id)
        while len(self._ram) > self._max_ram:
            self._ram.popitem(last=False)

    def clear(self) -> None:
        """Remove all models from RAM and disk. Intended for testing."""
        with self._lock:
            # Remove disk files
            for path in self._dir.glob("batman_*.pkl"):
                path.unlink(missing_ok=True)
            self._ram.clear()

    def _disk_path(self, tcr_id: str) -> Path:
        safe = hashlib.md5(tcr_id.encode()).hexdigest()
        return self._dir / f"batman_{safe}.pkl"
