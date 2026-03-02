"""BATCAVE database client for cross-reactivity reference data.

BATCAVE (Bhatt et al., 2024) provides activation scores for epitope
variants across multiple TCRs. Used as ground truth for cross-reactivity
prediction validation.
"""
from __future__ import annotations

import csv
import io
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "batcave"


@dataclass
class BatcaveVariant:
    reference_peptide: str
    variant_peptide: str
    activation_score: float
    mhc_allele: str
    mhc_class: str
    mutation_position: int
    original_aa: str
    mutant_aa: str


class BatcaveClient:
    """Client for querying local BATCAVE CSV data."""

    def __init__(self, data_dir: Optional[Path] = None) -> None:
        self._data_dir = data_dir or _DATA_DIR
        self._variants: list[BatcaveVariant] = []
        self._loaded = False

    def _ensure_loaded(self) -> None:
        if self._loaded:
            return
        self._variants = []
        if not self._data_dir.exists():
            logger.warning("BATCAVE data dir not found: %s", self._data_dir)
            self._loaded = True
            return

        for csv_path in sorted(self._data_dir.glob("*.csv")):
            try:
                with csv_path.open("r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        variant = self._parse_variant(row)
                        if variant is not None:
                            self._variants.append(variant)
            except Exception as exc:
                logger.warning("Failed to read %s: %s", csv_path, exc)
        self._loaded = True
        logger.info("Loaded %d BATCAVE variants", len(self._variants))

    def _parse_variant(self, row: dict) -> Optional[BatcaveVariant]:
        try:
            return BatcaveVariant(
                reference_peptide=row.get("peptide", "").upper().strip(),
                variant_peptide=row.get("variant", "").upper().strip(),
                activation_score=float(row.get("activation_score", 0)),
                mhc_allele=row.get("mhc_allele", ""),
                mhc_class=row.get("mhc_class", "I"),
                mutation_position=int(row.get("position", 0)),
                original_aa=row.get("original_aa", ""),
                mutant_aa=row.get("mutant_aa", ""),
            )
        except (ValueError, TypeError):
            return None

    def _filter_variants(
        self,
        variants: list[BatcaveVariant],
        reference_peptide: Optional[str] = None,
        mhc_class: Optional[str] = None,
        mhc_allele: Optional[str] = None,
    ) -> list[BatcaveVariant]:
        result = variants
        if reference_peptide:
            ref = reference_peptide.upper()
            result = [v for v in result if v.reference_peptide == ref]
        if mhc_class:
            result = [v for v in result if v.mhc_class == mhc_class]
        if mhc_allele:
            allele = mhc_allele.upper()
            result = [v for v in result if v.mhc_allele.upper() == allele]
        return result

    def lookup(
        self,
        reference_peptide: Optional[str] = None,
        mhc_class: Optional[str] = None,
        mhc_allele: Optional[str] = None,
    ) -> list[BatcaveVariant]:
        self._ensure_loaded()
        return self._filter_variants(
            self._variants,
            reference_peptide=reference_peptide,
            mhc_class=mhc_class,
            mhc_allele=mhc_allele,
        )
