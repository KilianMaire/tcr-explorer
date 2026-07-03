"""Known-epitope lookup and id resolution for the dossier, via the existing search path.

Thin and fully defensive: the VDJdb/IEDB tool servers are not guaranteed to be
running, so any failure (network error, missing server, unexpected shape)
degrades to an empty result rather than raising. Callers never need to guard
against exceptions from this module.
"""
from __future__ import annotations

import asyncio
from typing import Optional

from .api import _IEDB_HITS_CAP, search
from .models import IEDBHit, SearchRequest


def _run_search(req: SearchRequest):
    """Run the async `search()` coroutine from sync code.

    Deliberately avoids `asyncio.run()`: it tears down the loop it creates
    (via `set_event_loop(None)`) on exit, which permanently breaks any other
    test or caller in the same thread that later relies on
    `asyncio.get_event_loop()`. Instead this reuses (or lazily creates and
    keeps) a loop on the current thread, mirroring what `asyncio.run` does
    minus the destructive teardown.

    Returns `None` (treated as "no result") if a loop is already running in
    this thread (e.g. we were called from inside an async context) or if
    anything else goes wrong obtaining/using a loop.
    """
    try:
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            loop = None

        if loop is not None and loop.is_running():
            # Can't block synchronously on a loop that's already running.
            return None

        if loop is None or loop.is_closed():
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)

        return loop.run_until_complete(search(req))
    except Exception:
        return None


def lookup_known_epitopes(
    gene: Optional[str], cdr3_aa: Optional[str], species: str
) -> tuple[list[IEDBHit], int]:
    """Look up known epitopes for a V/J gene via the vdjdb/iedb search path.

    Returns `([], 0)` whenever `gene` is falsy, the search layer is
    unavailable, or anything unexpected happens. Never raises.
    """
    if not gene:
        return [], 0

    try:
        req = SearchRequest(source="vdjdb", gene_name=gene, species=species, limit=50)
        resp = _run_search(req)
        if resp is None:
            return [], 0

        hits: list[IEDBHit] = []
        for rec in resp.records:
            if rec.iedb_hits:
                hits.extend(rec.iedb_hits)
            elif rec.antigen_epitope:
                hits.append(IEDBHit(epitope_sequence=rec.antigen_epitope))

        return hits[:_IEDB_HITS_CAP], len(hits)
    except Exception:
        return [], 0


def resolve_id(source: str, ident: str) -> dict:
    """Best-effort resolve a `vdjdb:`/`iedb:` id via the same search layer.

    Returns `{}` on any failure; the caller is responsible for surfacing a
    `source_unavailable` warning in that case.
    """
    try:
        gene_source = source if source in ("vdjdb", "iedb", "hla", "tcr", "mhc") else None
        if gene_source is None:
            return {}
        req = SearchRequest(source=gene_source, limit=1)
        resp = _run_search(req)
        if resp is None or not resp.records:
            return {}
        for rec in resp.records:
            if rec.metadata.get("id") == ident or ident in (rec.gene_name, rec.allele_name):
                return rec.model_dump()
        return {}
    except Exception:
        return {}
