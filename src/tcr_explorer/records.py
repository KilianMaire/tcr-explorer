"""Per-row record builder and the federated retrieval engine.

Turns one harmonized index row (see `records_build.SCHEMA_COLUMNS`) into a
`TCRRecord`. Deposited sequences already present on the row are kept
verbatim and tagged "deposited". A value is only reconstructed when it is
missing on the row and both V and J genes are present; reconstructed
sequences are tagged "reconstructed" and flip `nt_is_synthetic`. Nothing is
ever fabricated when the genes are absent.

`retrieve_records` routes a `RecordsRequest` against the vendored harmonized
index to exact hits, BLOSUM neighbours, gene matches, alpha/beta pairs, or a
namespaced database id.
"""
from __future__ import annotations

import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd

from . import input_router
from .cdr_enricher import _gene_to_chain, _translate
from .dossier_models import (
    DossierWarning,
    PairedRecord,
    RecordsRequest,
    RecordsResponse,
    TCRRecord,
    Composition,
)
from .query_nl import parse_query
from .reconstructor import reconstruct_tcr
from .similarity import cdr3_distance, distance_to_similarity


def mhc_organism(mhc_a: Optional[str]) -> Optional[str]:
    """Organism implied by an MHC allele name: 'human' for HLA*, 'mouse' for
    H2*/H-2*, else None. Independent of the record's own species field, so it
    can flag HLA-transgenic mice and similar cross-species constructs."""
    if not mhc_a or (isinstance(mhc_a, float) and pd.isna(mhc_a)):
        return None
    s = str(mhc_a).strip().upper()
    if not s:
        return None
    if s.startswith("HLA"):
        return "human"
    if s.startswith("H2") or s.startswith("H-2"):
        return "mouse"
    return None


def build_record(
    row: dict,
    match_kind: str = "exact",
    similarity: Optional[float] = None,
    concordance: int = 1,
) -> TCRRecord:
    rec = TCRRecord(
        source=row["source"],
        source_record_id=row["source_record_id"],
        external_url=row["external_url"],
        pairing_key=row["pairing_key"],
        chain=row["chain"],
        species=row.get("species") or "",
        cdr3_aa=row["cdr3_aa"],
        v_gene=row.get("v_gene"),
        d_gene=row.get("d_gene"),
        j_gene=row.get("j_gene"),
        cdr1_aa=row.get("cdr1_aa"),
        cdr2_aa=row.get("cdr2_aa"),
        epitope_aa=row.get("epitope_aa"),
        antigen=row.get("antigen"),
        antigen_organism=row.get("antigen_organism"),
        mhc_class=row.get("mhc_class"),
        mhc_a=row.get("mhc_a"),
        mhc_b=row.get("mhc_b"),
        pdb_id=row.get("pdb_id"),
        reference_pmid=row.get("reference_pmid"),
        score=row.get("score"),
        match_kind=match_kind,
        similarity=similarity,
        concordance=concordance,
    )
    rec.mhc_organism = mhc_organism(row.get("mhc_a"))
    rec.mhc_is_cross_species = bool(rec.mhc_organism and rec.mhc_organism != (rec.species or "").strip().lower())

    # Deposited sequences first, kept verbatim.
    if row.get("cdr3_nt"):
        rec.cdr3_nt = row["cdr3_nt"]
        rec.cdr3_nt_kind = "deposited"
    if row.get("full_aa"):
        rec.full_aa = row["full_aa"]
        rec.full_aa_kind = "deposited"
    if row.get("full_nt"):
        rec.full_nt = row["full_nt"]
        rec.full_nt_kind = "deposited"

    # Reconstruct only when no nucleotide value was deposited at all, and
    # only when both V and J genes are present. Mixing a deposited nt
    # fragment with a back translated reconstruction would be dishonest
    # (the two could disagree), so any deposited cdr3_nt or full_nt blocks
    # reconstruction outright. No genes means no reconstruction and no
    # fabrication.
    has_deposited_nt = bool(row.get("cdr3_nt")) or bool(row.get("full_nt"))
    v, j = row.get("v_gene"), row.get("j_gene")
    if v and j and not has_deposited_nt:
        rc = reconstruct_tcr(v, j, row["cdr3_aa"], row.get("species") or "human")
        if rc.get("full_nt"):
            if rec.full_nt is None:
                rec.full_nt = rc["full_nt"]
                rec.full_nt_kind = "reconstructed"
                rec.nt_is_synthetic = True
            if rec.full_aa is None and rc.get("full_aa"):
                rec.full_aa = rc["full_aa"]
                rec.full_aa_kind = "reconstructed"
            if rec.cdr3_nt is None and rc.get("cdr3_nt"):
                rec.cdr3_nt = rc["cdr3_nt"]
                rec.cdr3_nt_kind = "reconstructed"
            # Use the in-frame pieces actually assembled (V up to Cys104, J FR4
            # after Phe/Trp118), not the raw region nt, which are not frame 0 and
            # would translate to garbage with stop codons.
            v_piece = rc.get("v_prefix_nt")
            j_piece = rc.get("j_suffix_nt")
            rec.composition = Composition(
                cdr3_aa=row["cdr3_aa"],
                v_region_nt=v_piece,
                cdr3_nt=rc.get("cdr3_nt"),
                j_region_nt=j_piece,
                v_germline_aa=_translate(v_piece) if v_piece else None,
                j_germline_aa=_translate(j_piece) if j_piece else None,
                note="reconstructed from V and J germline plus back translated CDR3",
            )
    return rec


# ---------------------------------------------------------------------------
# Retrieval engine
# ---------------------------------------------------------------------------
_ID_RE = re.compile(r"^(vdjdb|iedb|mcpas|tcr3d):")
_MAX_NEIGHBOURS = 25
_CHAIN_LABELS = {"TRA": "alpha", "TRB": "beta", "TRG": "gamma", "TRD": "delta"}


def _default_records_index_path() -> str:
    from .data_paths import records_index_path
    return os.environ.get("RECORDS_INDEX_PATH") or str(records_index_path())


@lru_cache(maxsize=4)
def _read_index_cached(path: str) -> pd.DataFrame:
    return pd.read_parquet(path)


def load_records_index(path: Optional[str] = None) -> Optional[pd.DataFrame]:
    """Load the harmonized index parquet, or None if it is absent.

    `path` defaults to the `RECORDS_INDEX_PATH` env var, else the user data dir
    index. The absent case is NOT cached, so a query issued after
    `tcr-explorer-refresh` picks up the newly built index in a long-running
    process (only successful loads are memoized).
    """
    p = Path(path or _default_records_index_path())
    if not p.exists():
        return None
    return _read_index_cached(str(p))


load_records_index.cache_clear = _read_index_cached.cache_clear


def infer_vj_from_cdr3(
    cdr3_aa: str,
    species: Optional[str] = None,
    *,
    index_path: Optional[str] = None,
    top: int = 5,
) -> list[dict]:
    """Infer likely V and J genes for a CDR3 by tallying database records that
    carry the exact same CDR3. Returns pairings ranked by supporting record
    count: ``[{"chain", "v_gene", "j_gene", "count"}]`` (gene names normalized
    to their base, allele stripped). Empty when the CDR3 matches no record, the
    index is unavailable, or matched rows lack a V or J gene. This is a
    frequency inference over deposited records, not a germline assignment: the
    caller must label a reconstruction built on it as inferred.
    """
    if not cdr3_aa or not cdr3_aa.strip():
        return []
    frame = load_records_index(index_path)
    if frame is None or frame.empty:
        return []
    work = frame[_cdr3_mask(frame, cdr3_aa)]
    if species:
        work = _species_filter(work, species)
    if work.empty:
        return []
    vb = _gene_base_series(work["v_gene"])
    jb = _gene_base_series(work["j_gene"])
    work = work.assign(_vb=vb, _jb=jb)
    work = work[(work["_vb"] != "") & (work["_jb"] != "")]
    if work.empty:
        return []
    counts = work.groupby(["_vb", "_jb"]).size().sort_values(ascending=False)
    out: list[dict] = []
    for (v, j), n in counts.items():
        chain = _segment_and_chain(v)[1] or _segment_and_chain(j)[1]
        out.append({"chain": chain, "v_gene": v, "j_gene": j, "count": int(n)})
        if len(out) >= top:
            break
    return out


def _segment_and_chain(gene: str) -> tuple[Optional[str], Optional[str]]:
    """Segment letter (V/D/J/C) and chain label for a gene, guarded against
    `_gene_to_chain`'s default-to-TRB behaviour on unrecognized input: only
    call it once the raw string is confirmed to start with a real chain
    prefix (TRA/TRB/TRG/TRD)."""
    s = (gene or "").strip().upper()
    if len(s) < 4 or not (s.startswith("TRA") or s.startswith("TRB") or s.startswith("TRG") or s.startswith("TRD")):
        return None, None
    chain_code = _gene_to_chain(s)
    return s[3], _CHAIN_LABELS.get(chain_code)


def _safe_normalize_gene(g: str) -> str:
    try:
        return input_router._normalize_gene(g)
    except Exception:
        return g


def _gene_base(g: Optional[str]) -> str:
    """Normalized gene BASE: strip any '*allele' suffix, uppercase. Mouse
    genes in the index keep their raw allele-suffixed source strings, so
    gene matching must never rely on exact string equality (see Task 1
    binding note)."""
    if not g:
        return ""
    return str(g).strip().upper().split("*")[0]


def _gene_base_series(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip().str.upper().str.split("*").str[0]


def _score_num(v) -> float:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return float("-inf")
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("-inf")


def _parse_request(request: RecordsRequest) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Resolve (cdr3_aa, cdr3_aa_b, v_gene, j_gene, id_lookup) from a request.

    When only `query` is set, classify it with `input_router.route`: a
    namespaced id (vdjdb:/iedb:/mcpas:/tcr3d:) is an id lookup; a raw_aa of
    length 8..22 becomes a cdr3_aa search; a gene_name/allele becomes v_gene
    or j_gene depending on its segment letter.
    """
    cdr3_aa = request.cdr3_aa
    cdr3_aa_b = request.cdr3_aa_b
    v_gene = request.v_gene
    j_gene = request.j_gene
    id_lookup: Optional[str] = None

    if request.query and not (cdr3_aa or cdr3_aa_b or v_gene or j_gene):
        q = request.query.strip()
        if _ID_RE.match(q):
            id_lookup = q
        else:
            routed = input_router.route(q, "auto")
            if routed.detected_type == "id":
                id_lookup = routed.normalized
            elif routed.detected_type == "raw_aa" and 8 <= len(routed.normalized) <= 22:
                cdr3_aa = routed.normalized
            elif routed.detected_type in ("gene_name", "allele"):
                gene = routed.normalized
                seg, _chain = _segment_and_chain(gene)
                if seg == "V":
                    v_gene = gene
                elif seg == "J":
                    j_gene = gene

    return cdr3_aa, cdr3_aa_b, v_gene, j_gene, id_lookup


def _apply_nl_query(request: RecordsRequest) -> RecordsRequest:
    """Layer natural-language parsing on top of `request.query` before
    `_parse_request` runs its single-token classification.

    `parse_query` is the richer front end for free text (prose, mixed
    species/gene/CDR3 phrases, French or English species words). Its parsed
    species always overrides `request.species` when present, per the task
    contract. Its parsed cdr3/gene/id only fill fields the caller did not
    already set explicitly, so an explicit `cdr3_aa=`/`v_gene=`/`j_gene=`
    request field is never clobbered. When a record id is parsed, `query` is
    rewritten to that id so the existing `_ID_RE` id-lookup path in
    `_parse_request` fires unchanged; otherwise `query` is left as-is so the
    existing single-token id/gene/cdr3 classification keeps working for
    queries `parse_query` does not touch (e.g. a bare gene name or sequence
    it cannot confidently classify)."""
    if not request.query:
        return request

    parsed = parse_query(request.query)
    updates: dict = {}

    if parsed["species"] is not None:
        updates["species"] = parsed["species"]

    if parsed["record_id"] is not None:
        updates["query"] = parsed["record_id"]
    else:
        if request.cdr3_aa is None and parsed["cdr3_aa"] is not None:
            updates["cdr3_aa"] = parsed["cdr3_aa"]
        if request.v_gene is None and parsed["v_gene"] is not None:
            updates["v_gene"] = parsed["v_gene"]
        if request.j_gene is None and parsed["j_gene"] is not None:
            updates["j_gene"] = parsed["j_gene"]

    if not updates:
        return request
    return request.model_copy(update=updates)


def _species_filter(frame: pd.DataFrame, species: Optional[str]) -> pd.DataFrame:
    if not species:
        return frame
    norm = species.strip().lower()
    return frame[frame["species"].fillna("").astype(str).str.strip().str.lower() == norm]


def _mhc_species_filter(frame: pd.DataFrame, species: str) -> pd.DataFrame:
    """Drop rows whose MHC allele implies a different organism than the
    query species (the HLA-transgenic mouse case). Rows with no inferrable
    MHC organism are kept, since absence of a known allele is not evidence
    of a cross-species record."""
    norm = species.strip().lower()
    organisms = frame["mhc_a"].map(mhc_organism)
    keep = organisms.isna() | (organisms == norm)
    return frame[keep]


def _cdr3_mask(frame: pd.DataFrame, cdr3: str) -> pd.Series:
    return frame["cdr3_aa"].fillna("").astype(str).str.upper() == cdr3.strip().upper()


def _concordance_for(frame: pd.DataFrame, cdr3: Optional[str]) -> int:
    if not cdr3:
        return 1
    matches = frame[_cdr3_mask(frame, cdr3)]
    n = matches["source"].nunique()
    return max(int(n), 1)


def _sanitize_row(row: dict) -> dict:
    """Coerce NaN (numpy float, from unpopulated parquet columns) to None.
    Real index rows carry NaN in unpopulated optional fields; pydantic
    rejects NaN for Optional[str]."""
    return {k: (None if isinstance(v, float) and pd.isna(v) else v) for k, v in row.items()}


def _build_exact_records(frame: pd.DataFrame, species_frame: pd.DataFrame) -> list[TCRRecord]:
    return [
        build_record(
            _sanitize_row(row),
            match_kind="exact",
            concordance=_concordance_for(species_frame, row.get("cdr3_aa")),
        )
        for row in frame.to_dict("records")
    ]


def _find_exact(work: pd.DataFrame, cdr3_aa: Optional[str], v_gene: Optional[str], j_gene: Optional[str]) -> pd.DataFrame:
    mask = pd.Series(True, index=work.index)
    if cdr3_aa:
        mask &= _cdr3_mask(work, cdr3_aa)
    if v_gene:
        qbase = _gene_base(_safe_normalize_gene(v_gene))
        mask &= _gene_base_series(work["v_gene"]) == qbase
    if j_gene:
        qbase = _gene_base(_safe_normalize_gene(j_gene))
        mask &= _gene_base_series(work["j_gene"]) == qbase
    return work[mask]


def _find_neighbours(work: pd.DataFrame, cdr3_aa: str, chains: Optional[set]) -> list[TCRRecord]:
    cand = work
    if chains:
        cand = cand[cand["chain"].isin(chains)]
    cand = cand[~_cdr3_mask(cand, cdr3_aa)]
    if cand.empty:
        return []
    query = cdr3_aa.strip().upper()
    dists = cand["cdr3_aa"].map(lambda s: cdr3_distance(query, str(s or "")))
    max_d = max(float(dists.max()), 1.0)
    scored = cand.assign(_dist=dists).nsmallest(_MAX_NEIGHBOURS, "_dist")
    out = []
    for row, d in zip(scored.to_dict("records"), scored["_dist"]):
        if str(row.get("cdr3_aa") or "").strip().upper() == query:
            continue  # never emit a neighbour whose cdr3 equals the query
        sim = round(distance_to_similarity(float(d), max_d), 4)
        out.append(build_record(_sanitize_row(row), match_kind="neighbour", similarity=sim))
    return out


_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _nonstandard_residues(seq: Optional[str]) -> list[str]:
    """Sorted unique non-standard residues in an amino acid query, or [] when the
    query is empty or all-standard. Lets a CDR3 with e.g. 'Z' still search but be
    flagged, rather than silently returning fuzzy neighbours for gibberish."""
    if not seq:
        return []
    return sorted({c for c in seq.upper() if c not in _STANDARD_AA})


def retrieve_records(request: RecordsRequest, index_path: Optional[str] = None) -> RecordsResponse:
    warnings: list[DossierWarning] = []
    request = _apply_nl_query(request)
    query_echo = request.model_dump()
    resolved_path = index_path or _default_records_index_path()
    df = load_records_index(resolved_path)

    if df is None or df.empty:
        present = Path(resolved_path).exists()
        message = (
            "records index is empty (a refresh may have failed). Run `tcr-explorer-refresh`."
            if present else
            "records data is not downloaded yet. Run `tcr-explorer-refresh` once."
        )
        warnings.append(
            DossierWarning(code="records_index_unavailable", block="records", message=message)
        )
        return RecordsResponse(query_echo=query_echo, warnings=warnings)

    from .data_paths import index_age_days, is_stale
    if is_stale():
        warnings.append(
            DossierWarning(
                code="records_index_stale",
                block="records",
                message=(f"records index is {int(index_age_days() or 0)} days old; "
                         "run `tcr-explorer-refresh` to update."),
            )
        )

    sources_searched = sorted(str(s) for s in df["source"].dropna().unique().tolist())

    cdr3_aa, cdr3_aa_b, v_gene, j_gene, id_lookup = _parse_request(request)

    for label, seq in (("cdr3", cdr3_aa), ("cdr3_b", cdr3_aa_b)):
        bad = _nonstandard_residues(seq)
        if bad:
            warnings.append(
                DossierWarning(
                    code="nonstandard_residues",
                    block="records",
                    message=(f"{label} query contains non-standard amino acid residue(s): "
                             f"{', '.join(bad)}. Any results are best-effort."),
                )
            )

    # Id lookup: filter by source_record_id verbatim, no neighbours, no species gate.
    if id_lookup is not None:
        hit_rows = df[df["source_record_id"].astype(str) == id_lookup]
        exact = _build_exact_records(hit_rows, df)
        return RecordsResponse(
            query_echo=query_echo,
            exact=exact,
            total_exact=len(exact),
            sources_searched=sources_searched,
            warnings=warnings,
        )

    work = _species_filter(df, request.species)
    if request.species and not request.include_cross_species_mhc:
        work = _mhc_species_filter(work, request.species)

    if not (cdr3_aa or v_gene or j_gene):
        # Nothing to search on (e.g. a bare pair query with no primary field).
        return RecordsResponse(query_echo=query_echo, sources_searched=sources_searched, warnings=warnings)

    exact_frame = _find_exact(work, cdr3_aa, v_gene, j_gene)
    exact_frame = exact_frame.assign(_score_sort=exact_frame["score"].map(_score_num))
    exact_frame = exact_frame.sort_values(["_score_sort", "source"], ascending=[False, True])
    capped = exact_frame.head(request.top_k)
    exact_records = _build_exact_records(capped, work)

    neighbours: list[TCRRecord] = []
    if cdr3_aa and request.include_neighbours:
        chains_present = set(exact_frame["chain"].dropna().unique().tolist()) or None
        neighbours = _find_neighbours(work, cdr3_aa, chains_present)

    pairs: list[PairedRecord] = []
    if cdr3_aa_b:
        b_frame = _find_exact(work, cdr3_aa_b, None, None)
        b_records = _build_exact_records(b_frame, work)

        a_by_key: dict[str, list[TCRRecord]] = {}
        for r in exact_records:
            a_by_key.setdefault(r.pairing_key, []).append(r)
        b_by_key: dict[str, list[TCRRecord]] = {}
        for r in b_records:
            b_by_key.setdefault(r.pairing_key, []).append(r)

        for key in sorted(set(a_by_key) & set(b_by_key)):
            for ra in a_by_key[key]:
                for rb in b_by_key[key]:
                    if ra.chain == rb.chain:
                        continue  # two same-chain rows are never a pair
                    alpha = ra if ra.chain == "alpha" else (rb if rb.chain == "alpha" else None)
                    beta = ra if ra.chain == "beta" else (rb if rb.chain == "beta" else None)
                    if alpha is None or beta is None:
                        continue
                    pairs.append(PairedRecord(pairing_key=key, source=alpha.source, alpha=alpha, beta=beta))

        if not pairs:
            exact_records = exact_records + b_records
            warnings.append(
                DossierWarning(
                    code="no_pairing_found",
                    block="pairs",
                    message="no shared pairing_key between the two chains",
                )
            )

    return RecordsResponse(
        query_echo=query_echo,
        exact=exact_records,
        neighbours=neighbours,
        pairs=pairs,
        total_exact=len(exact_records),
        sources_searched=sources_searched,
        warnings=warnings,
    )
