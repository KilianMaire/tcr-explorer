"""Per-row record builder.

Turns one harmonized index row (see `records_build.SCHEMA_COLUMNS`) into a
`TCRRecord`. Deposited sequences already present on the row are kept
verbatim and tagged "deposited". A value is only reconstructed when it is
missing on the row and both V and J genes are present; reconstructed
sequences are tagged "reconstructed" and flip `nt_is_synthetic`. Nothing is
ever fabricated when the genes are absent.
"""
from __future__ import annotations

from typing import Optional

from .cdr_enricher import _translate
from .dossier_models import Composition, TCRRecord
from .reconstructor import reconstruct_tcr


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
        species=row["species"],
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
            rec.composition = Composition(
                cdr3_aa=row["cdr3_aa"],
                v_region_nt=rc.get("v_region_nt"),
                cdr3_nt=rc.get("cdr3_nt"),
                j_region_nt=rc.get("j_region_nt"),
                v_germline_aa=_translate(rc["v_region_nt"]) if rc.get("v_region_nt") else None,
                j_germline_aa=_translate(rc["j_region_nt"]) if rc.get("j_region_nt") else None,
                note="reconstructed from V and J germline plus back translated CDR3",
            )
    return rec
