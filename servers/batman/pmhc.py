"""pMHC binding prediction: Class I (MHCflurry) and Class II (IEDB API)."""
from __future__ import annotations

import csv
import io
import logging
import math
from dataclasses import dataclass
from typing import Optional

import httpx

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Binding core result for Class II MHC
# ---------------------------------------------------------------------------
@dataclass
class BindingCoreResult:
    """Result of Class II binding core extraction.

    The TCR only 'sees' the 9-mer binding core presented by Class II MHC,
    so downstream cross-reactivity and scoring need this extracted core
    rather than the full 15-mer peptide.
    """

    peptide: str
    core: str
    core_start: int
    core_end: int
    score: float
    rank: float

# ---------------------------------------------------------------------------
# HLA class detection
# ---------------------------------------------------------------------------
_CLASS_II_PREFIXES = ("HLA-D", "DRA", "DRB", "DQA", "DQB", "DPA", "DPB")


def detect_mhc_class(allele: str) -> str:
    """Return 'I' or 'II' based on HLA allele string."""
    a = allele.strip().upper()
    for prefix in _CLASS_II_PREFIXES:
        if a.startswith(prefix):
            return "II"
    return "I"


# ---------------------------------------------------------------------------
# Score normalization
# ---------------------------------------------------------------------------
def rank_to_score(percentile_rank: float) -> float:
    """Convert IEDB/NetMHCpan percentile_rank [0, 100] to [0, 1].

    Logistic sigmoid: %Rank 0.5→~0.98, 2.0→~0.88, 10→0.50, 50→~0.02.
    """
    return 1.0 / (1.0 + math.exp(2.0 * (percentile_rank - 10.0)))


def ic50_to_score(ic50_nm: float) -> float:
    """Convert IC50 (nM) to presentation probability [0, 1].

    Uses a square-root transform so that IC50=50 nM (strong binder) maps to
    ~0.80 and IC50=50000 nM (threshold) maps exactly to 0.0.
    """
    if ic50_nm <= 0:
        return 1.0
    raw = max(0.0, 1.0 - math.log10(ic50_nm) / math.log10(50_000))
    return min(1.0, math.sqrt(raw))


# Lazy import — mhcflurry optional dependency
try:
    from mhcflurry import Class1PresentationPredictor
    _MHCFLURRY_AVAILABLE = True
except ImportError:
    _MHCFLURRY_AVAILABLE = False
    Class1PresentationPredictor = None  # type: ignore


class MHCflurryPredictor:
    """Wraps mhcflurry Class1PresentationPredictor as a lazy singleton."""

    def __init__(self) -> None:
        self._predictor = None

    def _load(self) -> bool:
        if self._predictor is not None:
            return True
        if Class1PresentationPredictor is None:
            return False
        try:
            self._predictor = Class1PresentationPredictor.load()
            return True
        except Exception as exc:
            logger.warning("MHCflurry model load failed: %s", exc)
            return False

    def predict_class_i(self, peptide: str, allele: str) -> Optional[float]:
        """Return MHCflurry presentation_score in [0,1] or None on failure."""
        if not self._load():
            return None
        try:
            result = self._predictor.predict(
                peptides=[peptide],
                alleles=[allele],
            )
            return float(result["presentation_score"].iloc[0])
        except Exception as exc:
            logger.warning("MHCflurry predict failed for %s/%s: %s", peptide, allele, exc)
            return None


# NetMHCIIpan — optional dependency for Class II local prediction
try:
    from mhctools import NetMHCIIpan
    _NETMHCIIPAN_AVAILABLE = True
except ImportError:
    _NETMHCIIPAN_AVAILABLE = False
    NetMHCIIpan = None  # type: ignore


def _netmhciipan_predict(peptide: str, allele: str) -> Optional[float]:
    """Run NetMHCIIpan prediction and return score [0,1]."""
    if NetMHCIIpan is None:
        return None
    try:
        predictor = NetMHCIIpan(alleles=[allele])
        binding_predictions = predictor.predict_subsequences({0: peptide})
        if not binding_predictions:
            return None
        best = min(binding_predictions, key=lambda bp: bp.percentile_rank)
        return rank_to_score(best.percentile_rank)
    except Exception as exc:
        logger.warning("NetMHCIIpan failed for %s/%s: %s", peptide, allele, exc)
        return None


class NetMHCIIpanPredictor:
    """Local NetMHCIIpan predictor for MHC Class II binding."""

    def predict(self, peptide: str, allele: str) -> Optional[float]:
        """Return binding score [0,1] via NetMHCIIpan, or None if unavailable."""
        if not _NETMHCIIPAN_AVAILABLE:
            return None
        return _netmhciipan_predict(peptide, allele)


_IEDB_MHCI_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
_IEDB_MHCII_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"


class IEDBpMHCPredictor:
    """Async IEDB pMHC API client for both Class I and Class II."""

    async def predict(self, peptide: str, allele: str) -> Optional[float]:
        """Return presentation score [0,1] via IEDB API or None on failure."""
        mhc_class = detect_mhc_class(allele)
        url = _IEDB_MHCI_URL if mhc_class == "I" else _IEDB_MHCII_URL
        peptide_length = str(len(peptide))

        payload = {
            "method": "recommended",
            "sequence_text": f">query\n{peptide}",
            "allele": allele,
            "length": peptide_length,
        }

        client = httpx.AsyncClient(timeout=15.0)
        try:
            # Enter async context manager
            entered = client.__aenter__()
            if hasattr(entered, "__await__"):
                entered = await entered
            r = entered.post(url, data=payload)
            if hasattr(r, "__await__"):
                r = await r
            r.raise_for_status()
            result = self._parse_tsv(r.text)
        except Exception as exc:
            logger.warning("IEDB pMHC API failed for %s/%s: %s", peptide, allele, exc)
            result = None
        finally:
            try:
                ex = client.__aexit__(None, None, None)
                if hasattr(ex, "__await__"):
                    await ex
            except Exception:
                pass
        return result

    def _parse_tsv(self, text: str) -> Optional[float]:
        """Parse IEDB TSV response to presentation score."""
        lines = [l for l in text.strip().splitlines() if not l.startswith("#")]
        if len(lines) < 2:
            return None
        reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
        for row in reader:
            try:
                rank = float(row.get("percentile_rank", 100))
                return rank_to_score(rank)
            except (ValueError, TypeError):
                try:
                    ic50 = float(row.get("ic50", 50_000))
                    return ic50_to_score(ic50)
                except (ValueError, TypeError):
                    pass
        return None


# Module-level singletons (created once at import time)
_mhcflurry = MHCflurryPredictor()
_netmhciipan = NetMHCIIpanPredictor()
_iedb = IEDBpMHCPredictor()


async def predict_pmhc(peptide: str, allele: str) -> Optional[float]:
    """Return P(peptide presented | HLA allele) in [0, 1] or None.

    Routing:
        Class I  -> MHCflurry (primary) -> IEDB /mhci/ (fallback)
        Class II -> NetMHCIIpan (primary) -> IEDB /mhcii/ (fallback)
    """
    mhc_class = detect_mhc_class(allele)

    if mhc_class == "I":
        score = _mhcflurry.predict_class_i(peptide, allele)
        if score is not None:
            return score
        return await _iedb.predict(peptide, allele)
    else:
        # Class II: NetMHCIIpan primary, IEDB fallback
        score = _netmhciipan.predict(peptide, allele)
        if score is not None:
            return score
        return await _iedb.predict(peptide, allele)
