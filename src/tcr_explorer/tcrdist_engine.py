"""Authoritative tcrdist metric, computed offline.

This reproduces the tcrdist distance used by tcrdist3, byte for byte on the
integer distances, without depending on tcrdist3 at runtime. It combines two
vendored or optional pieces:

- the germline CDR1, CDR2, and CDR2.5 (pMHC) loops for each V gene, read from
  the vendored reference table (`data/tcrdist/alphabeta_gammadelta_db.tsv`, MIT,
  from tcrdist3; CDR sequences derive from IMGT, CC BY 4.0), and
- the `pwseqdist` distance engine (the same engine tcrdist3 delegates to),
  installed only with the optional extra `pip install tcr-explorer[tcrdist]`.

When `pwseqdist` is not installed, `tcrdist_available()` returns False and the
caller falls back to the bundled BLOSUM CDR3 distance. The label reported to the
user must always match the engine actually used, never the one merely available.

The per region parameters below are the tcrdist3 defaults, confirmed against
tcrdist3's own `TCRrep` (parity test in `tests/test_tcrdist_engine.py`):
CDR3 carries weight 3 and trims its conserved ends (ntrim 3, ctrim 2,
fixed_gappos False); CDR1, CDR2, and CDR2.5 carry weight 1 and are compared
untrimmed with a fixed gap position. The gap penalty is 4 throughout.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional

import pandas as pd

_TABLE_PATH = Path(__file__).resolve().parent / "data" / "tcrdist" / "alphabeta_gammadelta_db.tsv"

# tcrdist3 default weights and per region metric parameters.
_CDR3_WEIGHT = 3
_LOOP_WEIGHT = 1
_GAP_PENALTY = 4


def tcrdist_available() -> bool:
    """True when the optional `pwseqdist` engine is importable."""
    try:
        import pwseqdist  # noqa: F401

        return True
    except Exception:
        return False


@lru_cache(maxsize=1)
def _vgene_cdr_table() -> dict[tuple[str, str], tuple[str, str, str]]:
    """Map (organism, V gene id) -> (cdr1, cdr2, cdr2.5), IMGT gap dots preserved."""
    df = pd.read_csv(_TABLE_PATH, sep="\t")
    df = df[df["region"] == "V"]
    out: dict[tuple[str, str], tuple[str, str, str]] = {}
    for row in df.itertuples(index=False):
        parts = str(row.cdrs).split(";")
        if len(parts) >= 3:
            out[(str(row.organism).lower(), str(row.id))] = (parts[0], parts[1], parts[2])
    return out


def resolve_v_loops(v_gene: str, organism: str = "human") -> Optional[tuple[str, str, str]]:
    """Germline CDR1, CDR2, CDR2.5 loops for a V gene, or None if not in the table.

    Tries the exact allele id first, then the `*01` default allele of the gene,
    matching tcrdist3's behaviour of defaulting an unspecified allele to `*01`.
    """
    if not isinstance(v_gene, str) or not v_gene:
        return None
    table = _vgene_cdr_table()
    org = organism.strip().lower()
    if (org, v_gene) in table:
        return table[(org, v_gene)]
    gene = v_gene.split("*")[0]
    return table.get((org, f"{gene}*01"))


def _metric_kwargs(is_cdr3: bool) -> dict:
    import pwseqdist as pw

    return dict(
        distance_matrix=pw.matrices.tcr_nb_distance_matrix,
        dist_weight=1,
        gap_penalty=_GAP_PENALTY,
        ntrim=3 if is_cdr3 else 0,
        ctrim=2 if is_cdr3 else 0,
        fixed_gappos=not is_cdr3,
    )


def _rect(seqs1: list[str], seqs2: list[str], is_cdr3: bool):
    import pwseqdist as pw

    return pw.apply_pairwise_rect(
        metric=pw.metrics.nb_vector_tcrdist,
        seqs1=seqs1,
        seqs2=seqs2,
        ncpus=1,
        uniqify=True,
        use_numba=True,
        **_metric_kwargs(is_cdr3),
    )


def tcrdist_pair(
    cdr3_a: str,
    v_a: str,
    cdr3_b: str,
    v_b: str,
    organism: str = "human",
) -> Optional[float]:
    """Full single chain tcrdist between two receptors, or None if either V is unresolvable."""
    la = resolve_v_loops(v_a, organism)
    lb = resolve_v_loops(v_b, organism)
    if la is None or lb is None:
        return None
    total = 0.0
    for xa, xb in zip(la, lb):
        total += _LOOP_WEIGHT * float(_rect([xa], [xb], is_cdr3=False)[0, 0])
    total += _CDR3_WEIGHT * float(_rect([cdr3_a], [cdr3_b], is_cdr3=True)[0, 0])
    return total


def tcrdist_one_to_many(
    query_cdr3: str,
    query_v: str,
    cand_cdr3s: list[str],
    cand_vs: list[str],
    organism: str = "human",
) -> tuple[list[Optional[float]], int]:
    """Full single chain tcrdist from one query to many candidates.

    Returns (distances, n_skipped). A candidate whose V gene is not in the
    reference table gets distance None and is counted in n_skipped: its germline
    loops are unknown, so a full tcrdist cannot be computed honestly.
    """
    ql = resolve_v_loops(query_v, organism)
    if ql is None:
        # The query V is required for the germline half; without it there is no
        # full tcrdist. The caller decides whether to fall back to BLOSUM.
        raise ValueError(f"query V gene {query_v!r} not in tcrdist reference for organism {organism!r}")

    n = len(cand_cdr3s)
    cand_loops = [resolve_v_loops(v, organism) for v in cand_vs]
    resolvable = [i for i, cl in enumerate(cand_loops) if cl is not None]
    n_skipped = n - len(resolvable)

    distances: list[Optional[float]] = [None] * n
    if not resolvable:
        return distances, n_skipped

    # CDR3 half: one vectorised 1xN rectangle over the resolvable candidates.
    cdr3_seqs = [cand_cdr3s[i] for i in resolvable]
    cdr3_d = _rect([query_cdr3], cdr3_seqs, is_cdr3=True)[0]

    # Germline half: three 1xN rectangles. pwseqdist uniqifies repeated loops,
    # so the many candidates sharing a V gene cost little.
    loop_d = [None, None, None]
    for k in range(3):
        cand_k = [cand_loops[i][k] for i in resolvable]
        loop_d[k] = _rect([ql[k]], cand_k, is_cdr3=False)[0]

    for pos, i in enumerate(resolvable):
        total = _CDR3_WEIGHT * float(cdr3_d[pos])
        for k in range(3):
            total += _LOOP_WEIGHT * float(loop_d[k][pos])
        distances[i] = total
    return distances, n_skipped


class PairedTCR(NamedTuple):
    """One alpha/beta paired receptor: alpha CDR3 and V gene, beta CDR3 and V gene."""

    cdr3_a: str
    v_a: str
    cdr3_b: str
    v_b: str


def tcrdist_paired(query: PairedTCR, cand: PairedTCR, organism: str = "human") -> Optional[float]:
    """Paired tcrdist: the sum of the per chain single chain tcrdist for alpha and beta.

    This is tcrdist3's paired convention (`pw_alpha + pw_beta`, each chain full
    weight). Returns None if any of the four V genes is unresolvable, since a full
    paired tcrdist cannot then be computed honestly.
    """
    da = tcrdist_pair(query.cdr3_a, query.v_a, cand.cdr3_a, cand.v_a, organism)
    db = tcrdist_pair(query.cdr3_b, query.v_b, cand.cdr3_b, cand.v_b, organism)
    if da is None or db is None:
        return None
    return da + db


def tcrdist_paired_one_to_many(
    query: PairedTCR,
    cands: list[PairedTCR],
    organism: str = "human",
) -> tuple[list[Optional[float]], int]:
    """Paired tcrdist from one paired query to many paired candidates.

    Returns (distances, n_skipped). A candidate is skipped (distance None) when
    either its alpha or beta V gene is unresolvable. Raises ValueError if either
    of the query's own V genes is unresolvable (both chains are required).
    """
    a_dists, _ = tcrdist_one_to_many(
        query.cdr3_a, query.v_a, [c.cdr3_a for c in cands], [c.v_a for c in cands], organism
    )
    b_dists, _ = tcrdist_one_to_many(
        query.cdr3_b, query.v_b, [c.cdr3_b for c in cands], [c.v_b for c in cands], organism
    )
    distances: list[Optional[float]] = []
    n_skipped = 0
    for da, db in zip(a_dists, b_dists):
        if da is None or db is None:
            distances.append(None)
            n_skipped += 1
        else:
            distances.append(da + db)
    return distances, n_skipped
