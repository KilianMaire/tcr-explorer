from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

from .dossier_models import DossierWarning, Neighbour, PairedNeighbour

_BLOSUM_PATH = Path(__file__).resolve().parent / "data" / "blosum62.json"
_GAP = 8
_NTRIM, _CTRIM = 3, 2  # trim conserved CDR3 ends (tcrdist convention)


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
def _read_index_cached(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def _load_index(path: str) -> Optional[pd.DataFrame]:
    # Do not cache the absent case, so a query after tcr-explorer-refresh sees
    # the newly built index in a long-running process.
    p = Path(path)
    if not p.exists():
        return None
    return _read_index_cached(path)


_load_index.cache_clear = _read_index_cached.cache_clear


def _v_family(v: object) -> str:
    # The multi-source index carries NaN (float) gene values; a NaN is truthy so
    # "v or ''" would not guard it. Coerce anything non-string to empty.
    if not isinstance(v, str):
        return ""
    return v.split("*")[0].split("-")[0]


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
    from .data_paths import records_index_path
    path = index_path or os.environ.get("UNITCR_INDEX_PATH") or str(records_index_path())
    df = _load_index(path)
    if df is None:
        warnings.append(
            DossierWarning(
                code="similarity_index_unavailable",
                block="neighbours",
                message="records data is not downloaded yet. Run `tcr-explorer-refresh` once.",
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

    # Score with the authoritative tcrdist metric when the optional `pwseqdist`
    # engine is installed, the species is supported, and the query V gene resolves
    # to germline loops. Otherwise fall back to the bundled BLOSUM CDR3 distance.
    # The label reported below must always match the engine actually used.
    from . import tcrdist_engine

    engine = "blosum_cdr3"
    dists = None
    organism = (species or "human").strip().lower()
    use_tcrdist = (
        tcrdist_engine.tcrdist_available()
        and organism in ("human", "mouse")
        and tcrdist_engine.resolve_v_loops(v_gene, organism) is not None
    )
    if use_tcrdist:
        raw, n_skipped = tcrdist_engine.tcrdist_one_to_many(
            cdr3, v_gene, cand[cdr3_col].tolist(), cand[v_col].tolist(), organism=organism
        )
        dser = pd.Series(raw, index=cand.index, dtype="float64").dropna()
        if not dser.empty:
            engine = "tcrdist"
            cand = cand.loc[dser.index]
            dists = dser
            if n_skipped:
                warnings.append(
                    DossierWarning(
                        code="tcrdist_candidates_skipped",
                        block="neighbours",
                        message=f"{n_skipped} candidate(s) skipped: V gene not in the tcrdist reference table",
                    )
                )
        # dser empty -> no candidate scorable under tcrdist; fall through to BLOSUM honestly.

    # Emit the downgrade warning whenever the authoritative tcrdist engine was NOT
    # the one used to score -- key on the engine actually used, never on whether
    # pwseqdist merely happens to be importable.
    if engine != "tcrdist":
        warnings.append(
            DossierWarning(
                code="tcrdist_unavailable",
                block="neighbours",
                message="tcrdist scoring unavailable; used the bundled BLOSUM CDR3 distance. "
                "Install the tcrdist extra (`pip install tcr-explorer[tcrdist]`) for authoritative scoring.",
            )
        )

    if dists is None:
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
            v_b_gene=_clean(getattr(r, v_col)),
            j_b_gene=_clean(getattr(r, j_col)),
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


_PAIRED_CTX_COLS = ("epitope_aa", "antigen", "antigen_organism", "mhc_class", "mhc_a")


def _build_paired_index(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """Reconstruct alpha/beta paired references from the per chain index, vectorised.

    The index carries one row per chain, linked by `pairing_key`. A pairing key
    that has both an alpha row and a beta row yields one paired candidate. The
    alpha chain contributes `cdr3_a`/`v_a`, the beta chain `cdr3_b`/`v_b`, and the
    epitope and MHC context is taken from the beta row (falling back to the alpha
    row when the beta row leaves it null).

    This is a single merge of the alpha and beta sub frames on `pairing_key`, not a
    Python loop over groups, so it scales to the full multi hundred thousand row
    index. The first row per (pairing_key, chain) is kept, matching a first wins
    reduction.
    """
    need = {"pairing_key", "chain", "cdr3_aa", "v_gene"}
    if not need.issubset(df.columns):
        return pd.DataFrame(columns=["cdr3_a_aa", "v_a_gene", "cdr3_b_aa", "v_b_gene", *_PAIRED_CTX_COLS])
    sub = df[df["chain"].isin(["alpha", "beta"])]
    if species:
        sub = sub[sub["species"].fillna("").str.strip().str.lower() == species.strip().lower()]
    sub = sub[sub["pairing_key"].notna() & sub["cdr3_aa"].notna()]

    ctx_present = [c for c in _PAIRED_CTX_COLS if c in sub.columns]
    keep = ["pairing_key", "cdr3_aa", "v_gene", *ctx_present]
    alpha = sub[sub["chain"] == "alpha"][keep].drop_duplicates("pairing_key", keep="first")
    beta = sub[sub["chain"] == "beta"][keep].drop_duplicates("pairing_key", keep="first")
    merged = alpha.merge(beta, on="pairing_key", suffixes=("_a", "_b"), how="inner")

    out = pd.DataFrame(
        {
            "cdr3_a_aa": merged["cdr3_aa_a"].to_numpy(),
            "v_a_gene": merged["v_gene_a"].to_numpy(),
            "cdr3_b_aa": merged["cdr3_aa_b"].to_numpy(),
            "v_b_gene": merged["v_gene_b"].to_numpy(),
        }
    )
    for c in _PAIRED_CTX_COLS:
        if c in ctx_present:
            b, a = merged[f"{c}_b"], merged[f"{c}_a"]
            out[c] = b.where(b.notna(), a).to_numpy()
        else:
            out[c] = None
    return out.reset_index(drop=True)


@lru_cache(maxsize=4)
def _paired_index_cached(path: str, species: str) -> pd.DataFrame:
    """Reconstructed paired index for a records index path, built once per process.

    Keyed on (path, species) so repeated paired queries reuse the reconstruction.
    Invalidated together with the record index cache via `_load_index.cache_clear`.
    """
    return _build_paired_index(_read_index_cached(path), species)


def _clear_index_caches() -> None:
    _read_index_cached.cache_clear()
    _paired_index_cached.cache_clear()


# Override the earlier record-only clear so a refresh also drops the paired index.
_load_index.cache_clear = _clear_index_caches


def find_similar_paired_tcrs(
    cdr3_a: str,
    v_a: str,
    cdr3_b: str,
    v_b: str,
    species: str = "human",
    top_k: int = 10,
    min_similarity: float = 0.0,
    index_path: Optional[str] = None,
) -> tuple[list[PairedNeighbour], str, int, list[DossierWarning]]:
    """Find paired (alpha/beta) neighbours by authoritative paired tcrdist.

    Paired scoring is tcrdist only: the paired distance is the sum of the alpha and
    beta single chain tcrdist, which has no BLOSUM analogue here. When the tcrdist
    extra is not installed, or a query V gene is unresolvable, this returns no
    neighbours with an explanatory warning rather than a misleading fallback.
    """
    if not all(isinstance(x, str) for x in (cdr3_a, v_a, cdr3_b, v_b)):
        raise TypeError("cdr3_a, v_a, cdr3_b, and v_b must be strings")
    warnings: list[DossierWarning] = []
    from . import tcrdist_engine
    from .data_paths import records_index_path

    path = index_path or os.environ.get("UNITCR_INDEX_PATH") or str(records_index_path())
    df = _load_index(path)
    if df is None:
        warnings.append(DossierWarning(
            code="similarity_index_unavailable", block="neighbours",
            message="records data is not downloaded yet. Run `tcr-explorer-refresh` once."))
        return [], "none", 0, warnings

    organism = (species or "human").strip().lower()
    if not tcrdist_engine.tcrdist_available() or organism not in ("human", "mouse"):
        warnings.append(DossierWarning(
            code="tcrdist_unavailable", block="neighbours",
            message="paired similarity requires the tcrdist extra (`pip install tcr-explorer[tcrdist]`) "
            "and a supported species (human or mouse)."))
        return [], "none", 0, warnings
    if tcrdist_engine.resolve_v_loops(v_a, organism) is None or tcrdist_engine.resolve_v_loops(v_b, organism) is None:
        warnings.append(DossierWarning(
            code="tcrdist_unavailable", block="neighbours",
            message="a query V gene is not in the tcrdist reference table, so a paired tcrdist "
            "cannot be computed. Provide alpha and beta V genes present in the reference."))
        return [], "none", 0, warnings

    cand = _paired_index_cached(path, organism)
    total = len(cand)
    if total == 0:
        warnings.append(DossierWarning(
            code="no_reference_candidates", block="neighbours",
            message="no paired (alpha and beta) reference candidates for this species"))
        return [], "none", 0, warnings

    query = tcrdist_engine.PairedTCR(cdr3_a, v_a, cdr3_b, v_b)
    cand_tcrs = [
        tcrdist_engine.PairedTCR(r.cdr3_a_aa, r.v_a_gene or "", r.cdr3_b_aa, r.v_b_gene or "")
        for r in cand.itertuples(index=False)
    ]
    raw, n_skipped = tcrdist_engine.tcrdist_paired_one_to_many(query, cand_tcrs, organism=organism)
    if n_skipped:
        warnings.append(DossierWarning(
            code="tcrdist_candidates_skipped", block="neighbours",
            message=f"{n_skipped} paired candidate(s) skipped: a V gene not in the tcrdist reference table"))

    dser = pd.Series(raw, index=cand.index, dtype="float64").dropna()
    if dser.empty:
        warnings.append(DossierWarning(
            code="no_reference_candidates", block="neighbours",
            message="no paired candidate could be scored under tcrdist for this query"))
        return [], "tcrdist", total, warnings

    cand = cand.loc[dser.index]
    max_d = max(float(dser.max()), 1.0)
    scored = cand.assign(
        dist_score=dser,
        sim_score=dser.map(lambda d: distance_to_similarity(d, max_d)),
    )
    scored = scored[scored["sim_score"] >= min_similarity].nsmallest(top_k, "dist_score")

    neigh = [
        PairedNeighbour(
            cdr3_a_aa=r.cdr3_a_aa, v_a_gene=_clean(r.v_a_gene),
            cdr3_b_aa=r.cdr3_b_aa, v_b_gene=_clean(r.v_b_gene),
            distance=round(float(r.dist_score), 4),
            similarity=round(float(r.sim_score), 4),
            epitope_aa=_clean(r.epitope_aa), antigen=_clean(r.antigen),
            antigen_organism=_clean(r.antigen_organism),
            mhc_class=_clean(r.mhc_class), mhc_a=_clean(r.mhc_a),
        )
        for r in scored.itertuples(index=False)
    ]
    return neigh, "tcrdist", total, warnings
