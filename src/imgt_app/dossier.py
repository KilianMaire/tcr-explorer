"""Dossier assembler: `build_dossier`.

Routes a query, runs annotation or germline lookup, and assembles one
`TCRDossier` with honest provenance, warnings, projection gating, and status.

Honesty invariants enforced here:
  * every populated block carries a matching `Provenance` entry;
  * any back-translated `cdr3_nt` is flagged `cdr3_nt_is_synthetic=True` and
    accompanied by a `back_translated` provenance record;
  * projection (`include`) gates long `nt` fields and `germline_nt`;
  * unresolved / unknown routes yield a `partial` dossier, never an exception;
  * `status="complete"` only when a V call exists and there are no warnings.
"""
from __future__ import annotations

from hashlib import sha256
from typing import Callable, Optional

from .dossier_models import (
    DossierRequest,
    DossierWarning,
    GeneCall,
    Junction,
    Provenance,
    RegionSeq,
    TCRDossier,
)
from .input_router import route
from .annotator import annotate
from .cdr_enricher import get_cdr1_cdr2, _gene_to_chain, _cached_v_map, _translate, _SPECIES_STITCHR
from .reconstructor import reconstruct_tcr
from .models import IEDBHit
from .dossier_epitopes import lookup_known_epitopes, resolve_id
from .similarity import find_similar_tcrs


# Map cdr_enricher's chain codes (TRA/TRB/TRG/TRD) to dossier chain names.
_CHAIN_NAME: dict[str, str] = {
    "TRA": "alpha",
    "TRB": "beta",
    "TRG": "gamma",
    "TRD": "delta",
}


def _echo(value: str, detected: str) -> dict:
    """Echo the query compactly; hash long inputs instead of embedding them."""
    if len(value) > 40:
        return {
            "length": len(value),
            "sha256": sha256(value.encode()).hexdigest(),
            "detected_type": detected,
        }
    return {"value": value, "detected_type": detected}


def _chain_name_for_gene(gene: str) -> str:
    return _CHAIN_NAME.get(_gene_to_chain(gene), "unknown")


def _segment_letter(gene: str) -> Optional[str]:
    """IMGT segment letter (V/D/J/C) — the 4th char of a TR gene, e.g. TRBV20-1 -> 'V'."""
    g = gene.strip().upper()
    return g[3] if len(g) >= 4 and g[:2] == "TR" and g[3] in "VDJC" else None


def _gene_path(
    gene: str,
    request: DossierRequest,
    genes: dict[str, Optional[GeneCall]],
    regions: dict[str, Optional[RegionSeq]],
    provenance: list[Provenance],
    warnings: list[DossierWarning],
    want_seq: bool,
    want_germ: bool,
) -> str:
    """Germline lookup path for a gene-name / allele query.

    Only V genes are enriched in this build. A V gene that resolves in germline
    populates the V `GeneCall`, CDR1/CDR2 regions, and (when requested) the
    germline nucleotide sequence, each with matching provenance. An unresolved
    V gene or a non-V (D/J/C) gene is NOT slotted or fabricated — it degrades to
    a warning. The chain, derived from the gene name, is honest to set either
    way. Returns the dossier chain name.
    """
    chain = _chain_name_for_gene(gene)
    call = gene.split("*")[0]
    segment = _segment_letter(gene)

    # Non-V genes: recognized but not enriched. Never slot as V, never fabricate.
    if segment in ("D", "J", "C"):
        warnings.append(
            DossierWarning(code="partial_annotation", block="genes",
                    message=(f"germline enrichment is only available for V genes in "
                             f"this build; {gene} recognized but not enriched"))
        )
        return chain

    res = get_cdr1_cdr2(gene, request.species)

    # Unresolved V gene (valid pattern, absent from germline): no fake call.
    if not res["allele"]:
        warnings.append(
            DossierWarning(code="ambiguous_gene", block="genes",
                    message=f"V gene {gene} not found in germline")
        )
        return chain

    v = GeneCall(call=call, allele=res["allele"], score_method="kmer")

    if res["cdr1_aa"] or res["cdr2_aa"]:
        if res["cdr1_aa"]:
            regions["cdr1"] = RegionSeq(
                aa=res["cdr1_aa"],
                nt=res["cdr1_nt"] if want_seq else None,
            )
        if res["cdr2_aa"]:
            regions["cdr2"] = RegionSeq(
                aa=res["cdr2_aa"],
                nt=res["cdr2_nt"] if want_seq else None,
            )
        provenance.append(
            Provenance(block="regions", source="cdr_enricher",
                       confidence="high", kind="germline_lookup")
        )

    # Germline sequence is honest only when the stitchr V-map actually holds it.
    stitchr_species = _SPECIES_STITCHR.get(request.species.lower(), "HUMAN")
    v_nt = _cached_v_map(_gene_to_chain(gene), stitchr_species).get(call.upper())
    if v_nt:
        v.germline_aa = _translate(v_nt).rstrip("*") or None
        if want_germ:
            v.germline_nt = v_nt
        provenance.append(
            Provenance(block="germline", source="cdr_enricher",
                       confidence="high", kind="germline_lookup")
        )

    genes["v"] = v
    return chain


def _seq_path(
    seq: str,
    request: DossierRequest,
    dt: str,
    genes: dict[str, Optional[GeneCall]],
    provenance: list[Provenance],
    warnings: list[DossierWarning],
) -> str:
    """Annotation path for a raw nucleotide / amino-acid query.

    Runs `annotate`, maps the result into V/J/D `GeneCall`s, forwards annotator
    warnings, and records annotation provenance. Returns the annotated chain.
    """
    ann = annotate(seq, request.species, is_protein=(dt == "raw_aa"), mode=request.mode)
    method = "igblast_pident" if ann.source == "igblast" else "kmer"

    if ann.v_call:
        genes["v"] = GeneCall(call=ann.v_call, score=ann.v_score, score_method=method)
    if ann.j_call:
        genes["j"] = GeneCall(call=ann.j_call, score=ann.j_score, score_method=method)
    if ann.d_call:
        genes["d"] = GeneCall(call=ann.d_call, score_method=method)

    for code, msg in ann.warnings:
        warnings.append(DossierWarning(code=code, block="annotation", message=msg))

    if ann.v_call or ann.j_call or ann.d_call:
        provenance.append(
            Provenance(block="annotation", source=ann.source,
                       confidence=ann.confidence, kind="observed")
        )
    return ann.chain


def _summarize(chain: str, genes: dict[str, Optional[GeneCall]], hits: list[IEDBHit]) -> str:
    v = genes["v"].call if genes["v"] else "NA"
    j = genes["j"].call if genes["j"] else "NA"
    return f"{chain} TCR, V={v}, J={j}, {len(hits)} known epitope(s)"


def build_dossier(
    request: DossierRequest,
    epitope_lookup: Callable[..., tuple[list[IEDBHit], int]] = lookup_known_epitopes,
) -> TCRDossier:
    warnings: list[DossierWarning] = []
    provenance: list[Provenance] = []
    want_seq = "sequences" in request.include
    want_germ = "germline" in request.include

    routed = route(request.query, request.input_type)
    for code, msg in routed.warnings:
        warnings.append(DossierWarning(code=code, block=None, message=msg))

    genes: dict[str, Optional[GeneCall]] = {"v": None, "d": None, "j": None, "c": None}
    # Pre-seed every region key so the top-level sub-object shape is stable for
    # strict JSON-schema / function-calling consumers (keys present, values None).
    regions: dict[str, Optional[RegionSeq]] = {
        "fr1": None, "fr2": None, "fr3": None, "fr4": None,
        "cdr1": None, "cdr2": None, "cdr3": None,
    }
    junction: Optional[Junction] = None
    full_sequence: Optional[RegionSeq] = None
    chain = "unknown"

    dt = routed.detected_type
    # V + J + CDR3 reconstruction path. Reachable when the caller supplies v_gene
    # and j_gene alongside a CDR3 (explicit cdr3_aa, or a raw_aa query used as the
    # CDR3). Everything it emits is synthetic/derived, never marked observed.
    cdr3_for_recon = request.cdr3_aa or (routed.normalized if dt == "raw_aa" else None)
    if request.v_gene and request.j_gene and cdr3_for_recon:
        chain, junction, full_sequence = _reconstruction_path(
            request, cdr3_for_recon, genes, regions, provenance,
            warnings, want_seq, want_germ,
        )
    elif dt in ("gene_name", "allele"):
        chain = _gene_path(routed.normalized, request, genes, regions,
                           provenance, warnings, want_seq, want_germ)
    elif dt in ("raw_nt", "raw_aa"):
        chain = _seq_path(routed.normalized, request, dt, genes, provenance, warnings)
        # A raw amino-acid TCR query is treated as a CDR3 candidate. We can only
        # reconstruct (and thus emit a SYNTHETIC cdr3_nt) when both a V and J
        # call were annotated; otherwise the junction carries the aa only.
        if dt == "raw_aa":
            junction = _build_junction_from_cdr3(
                routed.normalized, genes, request, provenance
            )

    elif dt == "id":
        resolved = {}
        if routed.source:
            _, _, ident = routed.normalized.partition(":")
            resolved = resolve_id(routed.source, ident)
        if not resolved:
            warnings.append(
                DossierWarning(code="source_unavailable", block="annotation",
                        message=f"could not resolve id {routed.normalized!r}")
            )
    else:
        warnings.append(
            DossierWarning(code="unresolved_input_type", block=None,
                    message="query could not be routed")
        )

    # Known epitopes via the real vdjdb/iedb search path (fully mockable).
    v_for_ep = genes["v"].call if genes["v"] else None
    cdr3_for_ep = junction.cdr3_aa if junction else None
    hits, total = epitope_lookup(v_for_ep, cdr3_for_ep, request.species)
    hits = hits[:5]
    if hits:
        provenance.append(
            Provenance(block="known_epitopes", source="iedb",
                       confidence="high", kind="observed")
        )

    # Neighbours: an inferred, weaker signal, kept strictly separate from
    # known_epitopes. Never copy a neighbour's epitope into known_epitopes.
    dossier_neighbours: Optional[list] = None
    want_neighbours = ("neighbours" in request.include) or (request.mode == "full")
    # Fall back to annotated values: a raw query that annotates to V/J/CDR3 should
    # also drive neighbours, not just an explicit v_gene/j_gene/cdr3_aa request.
    v_q = request.v_gene or (genes["v"].call if genes.get("v") else None)
    j_q = request.j_gene or (genes["j"].call if genes.get("j") else None)
    cdr3_q = (
        request.cdr3_aa
        or (junction.cdr3_aa if junction else None)
        or (routed.normalized if routed.detected_type == "raw_aa" else None)
    )
    if want_neighbours and cdr3_q and v_q and j_q:
        neigh, engine, total_n, nwarn = find_similar_tcrs(
            cdr3_q, v_q, j_q, species=request.species)
        for w in nwarn:
            warnings.append(w)
        if neigh:
            dossier_neighbours = neigh
            provenance.append(Provenance(block="neighbours", source="unitcr",
                confidence="low", kind="neighbor_inferred"))

    status = "complete" if (genes["v"] and not warnings) else "partial"
    summary = _summarize(chain, genes, hits)
    return TCRDossier(
        status=status,
        summary=summary,
        query_echo=_echo(request.query, dt),
        chain=chain,
        species=request.species,
        genes=genes,
        regions=regions,
        junction=junction,
        full_sequence=full_sequence,
        known_epitopes=hits,
        known_epitopes_total=max(total, len(hits)),
        provenance=provenance,
        warnings=warnings,
        neighbours=dossier_neighbours,
    )


def _build_junction_from_cdr3(
    cdr3_aa: str,
    genes: dict[str, Optional[GeneCall]],
    request: DossierRequest,
    provenance: list[Provenance],
) -> Junction:
    """Build a `Junction` from a bare CDR3 amino-acid query.

    Always records the CDR3 aa (verbatim). When both a V and J call are known,
    calls `reconstruct_tcr` to obtain a back-translated cdr3_nt, which is ALWAYS
    flagged synthetic and given a `back_translated` provenance record. Without
    V and J, cdr3_nt stays None (nothing is fabricated).
    """
    v = genes["v"].call if genes["v"] else None
    j = genes["j"].call if genes["j"] else None
    if v and j:
        res = reconstruct_tcr(v, j, cdr3_aa, request.species)
        cdr3_nt = res.get("cdr3_nt")
        if cdr3_nt:
            provenance.append(
                Provenance(block="junction", source="reconstructor",
                           confidence="medium", kind="back_translated")
            )
            return Junction(cdr3_aa=cdr3_aa, cdr3_nt=cdr3_nt, cdr3_nt_is_synthetic=True)
    return Junction(cdr3_aa=cdr3_aa)


def _reconstruction_path(
    request: DossierRequest,
    cdr3_aa: str,
    genes: dict[str, Optional[GeneCall]],
    regions: dict[str, Optional[RegionSeq]],
    provenance: list[Provenance],
    warnings: list[DossierWarning],
    want_seq: bool,
    want_germ: bool,
) -> tuple[str, Junction, Optional[RegionSeq]]:
    """Assemble a full dossier from a supplied V gene, J gene, and CDR3.

    Everything produced here is synthetic/derived, never observed:
      * the V call and its CDR1/CDR2/germline come from the germline lookup
        (`_gene_path`), which also fixes the chain;
      * the J gene is recorded verbatim as the caller supplied it;
      * the junction cdr3_nt is back-translated (flagged synthetic) via the
        shared `_build_junction_from_cdr3` machinery;
      * `full_sequence` is the reconstructor's assembled aa (always) plus its
        nt (only when `include` requests sequences), tagged `back_translated`.
    Never marks any nt as observed.
    """
    # V germline enrichment (also sets chain, CDR1/CDR2, germline_aa/nt).
    chain = _gene_path(request.v_gene, request, genes, regions,
                       provenance, warnings, want_seq, want_germ)
    # Record the caller's J gene call verbatim (not scored, not fabricated).
    j_call = request.j_gene.split("*")[0]
    genes["j"] = GeneCall(call=j_call)
    # Ensure the junction reconstruction can proceed even if the V gene was not
    # found in germline: back-translation of the CDR3 does not need germline.
    if genes["v"] is None:
        genes["v"] = GeneCall(call=request.v_gene.split("*")[0])

    junction = _build_junction_from_cdr3(cdr3_aa, genes, request, provenance)

    full_sequence: Optional[RegionSeq] = None
    res = reconstruct_tcr(request.v_gene, request.j_gene, cdr3_aa, request.species)
    if res.get("full_aa"):
        full_sequence = RegionSeq(
            aa=res["full_aa"],
            nt=res["full_nt"] if want_seq else None,
        )
        provenance.append(
            Provenance(block="full_sequence", source="reconstructor",
                       confidence="medium", kind="back_translated")
        )
    return chain, junction, full_sequence
