from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

from .dossier_models import DossierWarning, Neighbour

_DEFAULT_INDEX = str(Path(__file__).resolve().parent.parent.parent / "data" / "records_index.parquet")
_BLOSUM_PATH = Path(__file__).resolve().parent / "data" / "blosum62.json"
_GAP = 8
_NTRIM, _CTRIM = 3, 2  # trim conserved CDR3 ends (tcrdist convention)


def tcrdist3_available() -> bool:
    try:
        import tcrdist  # noqa: F401

        return True
    except Exception:
        return False


@lru_cache(maxsize=1)
def _blosum() -> dict:
    return json.loads(_BLOSUM_PATH.read_text())


def _sub(bl: dict, a: str, b: str) -> int:
    try:
        return bl[a][b]
    except KeyError:
        return -4


def cdr3_distance(a: str, b: str) -> float:
    """Bundled BLOSUM CDR3 distance. 0 for identical; grows with substitutions and length gaps."""
    bl = _blosum()
    a2, b2 = a[_NTRIM : len(a) - _CTRIM] or a, b[_NTRIM : len(b) - _CTRIM] or b
    if len(a2) == len(b2):
        return float(sum(_sub(bl, x, x) - _sub(bl, x, y) for x, y in zip(a2, b2)))
    short, long = (a2, b2) if len(a2) <= len(b2) else (b2, a2)
    best = None
    for off in range(len(long) - len(short) + 1):
        window = long[off : off + len(short)]
        d = sum(_sub(bl, x, x) - _sub(bl, x, y) for x, y in zip(short, window))
        best = d if best is None else min(best, d)
    return float((best or 0) + _GAP * (len(long) - len(short)))


@lru_cache(maxsize=4)
def _load_index(path: str) -> Optional[pd.DataFrame]:
    p = Path(path)
    if not p.exists():
        return None
    return pd.read_parquet(p)


def _v_family(v: str) -> str:
    return (v or "").split("*")[0].split("-")[0]


def _clean(v):
    """Coerce NaN (numpy float) in optional string columns to None. Real index
    rows carry NaN in unpopulated optional fields; pydantic rejects NaN for Optional[str]."""
    return None if v is None or (isinstance(v, float) and pd.isna(v)) else v


def distance_to_similarity(d: float, max_d: float = 1.0) -> float:
    """Within-query relative similarity, normalised to the candidate maximum
    distance for this query. Only similarity==1.0 (identical) is absolute;
    cross-query comparisons must threshold on the raw distance instead."""
    max_d = max(float(max_d), 1.0)
    return 1.0 - float(d) / max_d


def find_similar_tcrs(
    cdr3: str,
    v_gene: str,
    j_gene: str,
    species: str = "human",
    top_k: int = 10,
    min_similarity: float = 0.0,
    index_path: Optional[str] = None,
    chain: Optional[str] = None,
) -> tuple[list[Neighbour], str, int, list[DossierWarning]]:
    if not isinstance(cdr3, str) or not isinstance(v_gene, str) or not isinstance(j_gene, str):
        raise TypeError("cdr3, v_gene, and j_gene must be strings")
    warnings: list[DossierWarning] = []
    path = index_path or os.environ.get("UNITCR_INDEX_PATH") or _DEFAULT_INDEX
    df = _load_index(path)
    if df is None:
        warnings.append(
            DossierWarning(
                code="similarity_index_unavailable",
                block="neighbours",
                message=f"vendored index not found at {path}",
            )
        )
        return [], "none", 0, warnings

    is_legacy = "cdr3_b_aa" in df.columns

    if is_legacy:
        # The legacy vendored index is human beta only. We never score a
        # non-human query against it; returning human neighbours for a mouse
        # query would be dishonest.
        if species != "human":
            warnings.append(
                DossierWarning(
                    code="species_unsupported",
                    block="neighbours",
                    message=f"similarity index is human beta only; species '{species}' not supported",
                )
            )
            return [], "none", 0, warnings
        cdr3_col, v_col, j_col = "cdr3_b_aa", "v_b_gene", "j_b_gene"
        cand = df
    else:
        # The new multi-source index is multi-species and multi-chain.
        # Species is a filter here, never a gate.
        cdr3_col, v_col, j_col = "cdr3_aa", "v_gene", "j_gene"
        cand = df
        if species:
            cand = cand[cand["species"].fillna("").str.strip().str.lower() == species.strip().lower()]
        if chain and "chain" in cand.columns:
            cand = cand[cand["chain"] == chain]

    fam = _v_family(v_gene)
    fam_cand = cand[cand[v_col].map(_v_family) == fam]
    if len(fam_cand) >= 5:
        cand = fam_cand
    total = len(cand)
    if total == 0:
        warnings.append(
            DossierWarning(
                code="no_reference_candidates",
                block="neighbours",
                message="no reference candidates for this query",
            )
        )
        return [], "none", 0, warnings

    # NOTE: the tcrdist3-authoritative scoring path is NOT implemented yet (deferred
    # follow-up). Only the bundled BLOSUM CDR3 distance is ever used to score
    # candidates below, so `engine` always reports "blosum_cdr3" regardless of
    # whether tcrdist3 happens to be importable -- the label must match the distance
    # actually used, not the distance that would ideally be used.
    engine = "blosum_cdr3"
    # Emit the downgrade warning whenever the authoritative tcrdist engine was NOT
    # the one used to score -- key on the engine actually used, not on whether
    # tcrdist3 merely happens to be importable. While `engine` is hardcoded to
    # "blosum_cdr3" the tcrdist path is never authoritative, so this always fires.
    if engine != "tcrdist":
        warnings.append(
            DossierWarning(
                code="tcrdist_unavailable",
                block="neighbours",
                message="tcrdist3 authoritative scoring not wired; used the bundled BLOSUM CDR3 distance",
            )
        )

    dists = cand[cdr3_col].map(lambda s: cdr3_distance(cdr3, s))
    # NOTE: similarity below is within-query relative -- normalised to this
    # query's candidate maximum distance, so only similarity==1.0 (identical) is
    # absolute. Cross-query comparisons must threshold on the `distance` field.
    max_d = max(float(dists.max()), 1.0)
    # NOTE: itertuples() below needs plain identifier column names; a leading
    # underscore (e.g. "_dist") is not a valid namedtuple field, so pandas
    # silently falls back to positional access. Use "dist_score"/"sim_score" instead.
    scored = cand.assign(dist_score=dists, sim_score=dists.map(lambda d: distance_to_similarity(d, max_d)))
    scored = scored[scored["sim_score"] >= min_similarity].nsmallest(top_k, "dist_score")

    neigh = [
        Neighbour(
            cdr3_b_aa=getattr(r, cdr3_col),
            v_b_gene=getattr(r, v_col),
            j_b_gene=getattr(r, j_col),
            similarity=round(float(r.sim_score), 4),
            distance=round(float(r.dist_score), 4),
            epitope_aa=_clean(getattr(r, "epitope_aa", None)),
            mhc_class=_clean(getattr(r, "mhc_class", None)),
            mhc_a=_clean(getattr(r, "mhc_a", None)),
            antigen=_clean(getattr(r, "antigen", None)),
            antigen_organism=_clean(getattr(r, "antigen_organism", None)),
            cluster_id=(int(r.cluster_id) if pd.notna(getattr(r, "cluster_id", None)) else None),
        )
        for r in scored.itertuples(index=False)
    ]
    return neigh, engine, total, warnings
