"""conftest.py — install lightweight stubs for optional tcr_explorer dependencies.

The worktree branch only ships the core source files needed for the D-series
batman tasks.  Modules that exist in the main branch (cdr_enricher, fasta_parser,
mcp_clients, nl_query, reconstructor) are stubbed here so that ``api.py`` can be
imported without error.

On the main branch, pure-logic modules (nl_query, cdr_enricher, fasta_parser,
reconstructor) use their real implementations.  I/O-dependent modules
(mcp_clients, search_index, file_ingest) are always stubbed so tests don't
hit external services or local databases.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

# Ensure src/ is importable so real modules are preferred over stubs.
_src = str(Path(__file__).resolve().parent.parent / "src")
if _src not in sys.path:
    sys.path.insert(0, _src)


def _stub(module_name: str, **attrs):
    """Create and register a minimal stub module with given attributes."""
    mod = types.ModuleType(module_name)
    for k, v in attrs.items():
        setattr(mod, k, v)
    sys.modules[module_name] = mod
    return mod


def _stub_if_missing(module_name: str, **attrs):
    """Only stub the module if it cannot be imported from the real source."""
    try:
        __import__(module_name)
    except ImportError:
        _stub(module_name, **attrs)


# ---- Pure-logic modules: use real if available, stub otherwise ---------------

_stub_if_missing("tcr_explorer.cdr_enricher", get_cdr1_cdr2=lambda v_gene, species="human": {
    "allele": None, "cdr1_aa": None, "cdr2_aa": None,
    "cdr1_nt": None, "cdr2_nt": None,
})

_stub_if_missing("tcr_explorer.fasta_parser",
      parse_cdr3_fasta=lambda raw: [],
      parse_fasta_bytes=lambda raw, source=None, default_species="other": [])

async def _fake_lmstudio_parse(query: str):
    from tcr_explorer.models import ParseQueryResult
    return ParseQueryResult()

def _fake_heuristic_parse(query: str):
    from tcr_explorer.models import ParseQueryResult
    return ParseQueryResult()

_stub_if_missing("tcr_explorer.nl_query",
      lmstudio_parse=_fake_lmstudio_parse,
      heuristic_parse=_fake_heuristic_parse)

_stub_if_missing("tcr_explorer.reconstructor", reconstruct_tcr=lambda *a, **kw: {
    "v_gene": "", "j_gene": "", "cdr3_aa": "", "species": "human",
    "full_nt": None, "full_aa": None, "v_region_nt": None,
    "cdr3_nt": "", "j_region_nt": None, "v_found": False, "j_found": False,
    "note": "stub",
})

# ---- I/O-dependent modules: always stub to avoid external deps ---------------

class _FakeToolServerClient:
    def __init__(self, url: str) -> None:
        self.url = url

    async def search(self, req):
        from tcr_explorer.models import SearchResponse
        return SearchResponse(total=0, records=[], limit=50, offset=0)


_stub("tcr_explorer.mcp_clients", ToolServerClient=_FakeToolServerClient)

_stub_if_missing("tcr_explorer.file_ingest",
      parse_file=lambda raw, filename, source=None, species="other": [],
      parse_vdjdb_tsv=lambda raw: [])
