"""Known-epitope lookup and id resolution for the dossier, via the existing search path.

Thin and fully defensive: the VDJdb/IEDB tool servers are not guaranteed to be
running, so any failure (network error, missing server, unexpected shape)
degrades to an empty result rather than raising. Callers never need to guard
against exceptions from this module.
"""
from __future__ import annotations

import asyncio
from typing import Optional

from .api import search
from .models import IEDBHit, SearchRequest

# Cap on the number of known-epitope hits returned per dossier lookup.
_IEDB_HITS_CAP = 5


def _run_search(req: SearchRequest):
    """Run the async `search()` coroutine from sync code.

    Never touches `asyncio.set_event_loop()` (which is what poisoned the
    thread's global loop and broke unrelated tests earlier); always closes any
    loop it creates. Returns `None` (treated as "no result") when we're already
    inside a running loop, since we cannot block on it synchronously.
    """
    # If we're already inside a running loop, we can't block; degrade gracefully.
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        pass  # no running loop -> safe to create one below
    else:
        return None  # caller treats None as "no result"

    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(search(req))
    finally:
        loop.close()


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

    NOTE: `SearchRequest`/`SearchIndex` have no id filter field, so this cannot
    do true server-side id resolution within sub-project A (we must not modify
    `search()`/models here). It pulls a page of records and strictly matches
    `ident` client-side; the common no-live-server case returns `{}` and the
    dossier degrades to `source_unavailable` + partial. Reliable id resolution
    needs a dedicated id filter on SearchRequest/SearchIndex, deferred beyond
    sub-project A. Never fabricates a record. Returns `{}` on any failure.
    """
    try:
        # The router only ever produces vdjdb/iedb sources for id inputs.
        if source not in ("vdjdb", "iedb"):
            return {}
        req = SearchRequest(source=source, limit=50)
        resp = _run_search(req)
        if resp is None or not resp.records:
            return {}
        for rec in resp.records:
            if rec.metadata.get("id") == ident or ident in (rec.gene_name, rec.allele_name):
                return rec.model_dump()
        return {}
    except Exception:
        return {}
