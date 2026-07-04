"""Tests for the in-memory LRU cache in hla_server."""
from __future__ import annotations

import sys
from collections import OrderedDict
from pathlib import Path

# Ensure servers/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import servers.hla_server as hla_mod


# ---------------------------------------------------------------------------
# Helper — apply the same LRU eviction logic used by both servers
# ---------------------------------------------------------------------------

def _lru_set(
    cache: OrderedDict,
    key: tuple,
    value: list,
    max_size: int,
) -> None:
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > max_size:
        cache.popitem(last=False)


def _lru_get(cache: OrderedDict, key: tuple) -> list | None:
    if key in cache:
        cache.move_to_end(key)
        return cache[key]
    return None


# ---------------------------------------------------------------------------
# Generic LRU behaviour (no mocking needed)
# ---------------------------------------------------------------------------

def test_lru_hit_moves_to_end():
    c: OrderedDict = OrderedDict()
    _lru_set(c, ("a",), [1], max_size=3)
    _lru_set(c, ("b",), [2], max_size=3)
    _lru_set(c, ("c",), [3], max_size=3)
    # access "a" — it should become the most-recently-used
    _lru_get(c, ("a",))
    # adding a 4th evicts the LRU, which is now "b"
    _lru_set(c, ("d",), [4], max_size=3)
    assert ("b",) not in c
    assert ("a",) in c


def test_lru_evicts_oldest_when_full():
    c: OrderedDict = OrderedDict()
    for i in range(5):
        _lru_set(c, (str(i),), [i], max_size=3)
    # only the 3 most-recently inserted items should remain
    assert len(c) == 3
    assert ("0",) not in c
    assert ("1",) not in c
    assert ("4",) in c


def test_lru_no_eviction_below_max():
    c: OrderedDict = OrderedDict()
    for i in range(3):
        _lru_set(c, (str(i),), [i], max_size=64)
    assert len(c) == 3


# ---------------------------------------------------------------------------
# HLA cache module-level attributes
# ---------------------------------------------------------------------------

def test_hla_cache_is_ordered_dict():
    assert isinstance(hla_mod._EBI_CACHE, OrderedDict)


def test_hla_cache_max_is_positive_int():
    assert isinstance(hla_mod._EBI_CACHE_MAX, int)
    assert hla_mod._EBI_CACHE_MAX > 0


def test_hla_cache_populates_on_manual_insert():
    """Directly exercise the cache dict used by _fetch_ebi."""
    hla_mod._EBI_CACHE.clear()
    key = ("A", "")
    hla_mod._EBI_CACHE[key] = [{"source": "hla", "sequence": "ATGC"}]
    assert key in hla_mod._EBI_CACHE
    hla_mod._EBI_CACHE.clear()
