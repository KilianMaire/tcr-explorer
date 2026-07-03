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
    GeneCall,
    Junction,
    Provenance,
    RegionSeq,
    TCRDossier,
    Warning,
)
from .input_router import route
from .annotator import annotate
from .cdr_enricher import get_cdr1_cdr2, _gene_to_chain, _cached_v_map, _translate, _SPECIES_STITCHR
from .reconstructor import reconstruct_tcr
from .models import IEDBHit


# Map cdr_enricher's chain codes (TRA/TRB/TRG/TRD) to dossier chain names.
_CHAIN_NAME: dict[str, str] = {
    "TRA": "alpha",
    "TRB": "beta",
    "TRG": "gamma",
    "TRD": "delta",
}


def _noop_lookup(gene, cdr3_aa, species):
    """Default epitope lookup: a no-op. The real one is injected in Task 6."""
    return [], 0


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


def _gene_path(
    gene: str,
    request: DossierRequest,
    genes: dict[str, Optional[GeneCall]],
    regions: dict[str, Optional[RegionSeq]],
    provenance: list[Provenance],
    want_seq: bool,
    want_germ: bool,
) -> str:
    """Germline lookup path for a gene-name / allele query.

    Populates the V `GeneCall`, CDR1/CDR2 regions from the germline, and (when
    requested) the germline nucleotide sequence, each with matching provenance.
    Returns the dossier chain name.
    """
    chain = _chain_name_for_gene(gene)
    call = gene.split("*")[0]

    res = get_cdr1_cdr2(gene, request.species)
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
    if res["allele"]:
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
    warnings: list[Warning],
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
        warnings.append(Warning(code=code, block="annotation", message=msg))

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
    epitope_lookup: Callable[..., tuple[list[IEDBHit], int]] = _noop_lookup,
) -> TCRDossier:
    warnings: list[Warning] = []
    provenance: list[Provenance] = []
    want_seq = "sequences" in request.include
    want_germ = "germline" in request.include

    routed = route(request.query, request.input_type)
    for code, msg in routed.warnings:
        warnings.append(Warning(code=code, block=None, message=msg))

    genes: dict[str, Optional[GeneCall]] = {"v": None, "d": None, "j": None, "c": None}
    regions: dict[str, Optional[RegionSeq]] = {}
    junction: Optional[Junction] = None
    chain = "unknown"

    dt = routed.detected_type
    if dt in ("gene_name", "allele"):
        chain = _gene_path(routed.normalized, request, genes, regions,
                           provenance, want_seq, want_germ)
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
        # Real id resolution (vdjdb:/iedb:) is wired in Task 6.
        warnings.append(
            Warning(code="source_unavailable", block="annotation",
                    message="id resolution is wired in a later task")
        )
    else:
        warnings.append(
            Warning(code="unresolved_input_type", block=None,
                    message="query could not be routed")
        )

    # Known epitopes (no-op by default; the real lookup is injected in Task 6).
    v_for_ep = genes["v"].call if genes["v"] else None
    cdr3_for_ep = junction.cdr3_aa if junction else None
    hits, total = epitope_lookup(v_for_ep, cdr3_for_ep, request.species)
    hits = hits[:5]
    if hits:
        provenance.append(
            Provenance(block="known_epitopes", source="iedb",
                       confidence="high", kind="observed")
        )

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
        full_sequence=None,
        known_epitopes=hits,
        known_epitopes_total=max(total, len(hits)),
        provenance=provenance,
        warnings=warnings,
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
