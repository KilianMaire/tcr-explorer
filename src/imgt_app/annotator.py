"""
Annotator interface: selects between IgBLAST (authoritative, deferred) and the
always-available k-mer aligner fallback (Task 3).

Selection rule: IgBLAST is only attempted when mode=="full", the input is
nucleotide (not protein), and the `igblastn` binary is present on PATH. The
actual subprocess invocation is deferred (`_run_igblast` is a stub returning
None), so this module currently always falls through to the k-mer backend;
the point of this module is the selection/fallback logic and the
`igblast_unavailable` warning, not the IgBLAST call itself.
"""
from __future__ import annotations
import shutil
from dataclasses import dataclass, field
from typing import Optional
from .kmer_aligner import annotate_sequence, KmerAnnotation


@dataclass
class Annotation:
    v_call: Optional[str] = None
    j_call: Optional[str] = None
    d_call: Optional[str] = None
    v_score: Optional[float] = None
    j_score: Optional[float] = None
    chain: str = "unknown"
    source: str = "kmer_align"
    confidence: str = "low"
    warnings: list[tuple[str, str]] = field(default_factory=list)


def igblast_available() -> bool:
    return shutil.which("igblastn") is not None


def _run_igblast(seq: str, species: str) -> Optional[Annotation]:
    # Deferred: real igblastn -outfmt 19 invocation + AIRR TSV parse.
    # Returning None signals "not implemented / failed" so callers fall back.
    return None


def _from_kmer(k: KmerAnnotation, source: str, confidence: str) -> Annotation:
    return Annotation(
        v_call=k.v_call, j_call=k.j_call, d_call=k.d_call,
        v_score=k.v_score, j_score=k.j_score, chain=k.chain,
        source=source, confidence=confidence, warnings=list(k.warnings),
    )


def annotate(seq: str, species: str, is_protein: bool, mode: str) -> Annotation:
    use_igblast = (mode == "full") and (not is_protein) and igblast_available()
    if use_igblast:
        res = _run_igblast(seq, species)
        if res is not None:
            res.source, res.confidence = "igblast", "high"
            return res
        # fell through: igblast present but failed

    k = annotate_sequence(seq, species, is_protein)
    ann = _from_kmer(k, "kmer_align", "medium" if not is_protein else "low")
    if mode == "full" and not is_protein and not igblast_available():
        ann.warnings.append((
            "igblast_unavailable",
            "igblastn not found on PATH; used the k-mer backend",
        ))
    return ann
