"""Tests for BATMAN model cache (RAM LRU + disk persistence)."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import numpy as np
import tempfile
from batman.cache import ModelCache, CachedModel


class TestCachedModel:
    def test_cached_model_stores_weights_and_matrix(self):
        weights = np.array([[0.1, 0.2, 0.3]])
        matrix = np.eye(20)
        m = CachedModel(tcr_id="TCR1", index_peptide="GILGFVFTL",
                        weights=weights, aa_matrix=matrix)
        assert m.tcr_id == "TCR1"
        assert m.index_peptide == "GILGFVFTL"
        np.testing.assert_array_equal(m.weights, weights)
        np.testing.assert_array_equal(m.aa_matrix, matrix)


class TestModelCache:
    def test_cache_miss_returns_none(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = ModelCache(cache_dir=tmpdir, max_ram=10)
            assert cache.get("nonexistent_tcr") is None

    def test_cache_store_and_retrieve_ram(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = ModelCache(cache_dir=tmpdir, max_ram=10)
            model = CachedModel("TCR1", "GILGFVFTL", np.ones((1, 9)), np.eye(20))
            cache.put(model)
            result = cache.get("TCR1")
            assert result is not None
            assert result.tcr_id == "TCR1"

    def test_cache_persists_to_disk(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache1 = ModelCache(cache_dir=tmpdir, max_ram=10)
            model = CachedModel("TCR1", "GILGFVFTL", np.ones((1, 9)), np.eye(20))
            cache1.put(model)
            cache2 = ModelCache(cache_dir=tmpdir, max_ram=10)
            result = cache2.get("TCR1")
            assert result is not None
            assert result.tcr_id == "TCR1"

    def test_ram_lru_eviction(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = ModelCache(cache_dir=tmpdir, max_ram=2)
            for i in range(3):
                m = CachedModel(f"TCR{i}", "GILGFVFTL", np.ones((1, 9)), np.eye(20))
                cache.put(m)
            assert cache.get("TCR0") is not None  # loaded from disk
            assert cache.ram_size <= 2
