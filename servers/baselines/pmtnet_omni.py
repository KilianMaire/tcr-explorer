"""pMTnet Omni API client for pan-MHC TCR-epitope binding prediction.
AUROC: 0.888. Handles both Class I and II.
"""
from __future__ import annotations
import logging
from dataclasses import dataclass
import httpx

logger = logging.getLogger(__name__)
_DEFAULT_URL = "https://dbai.biohpc.swmed.edu/pmtnet"

@dataclass
class PmtnetOmniResult:
    score: float
    rank: float = 0.0
    tcr_representation: str = "beta_only"

class PmtnetOmniClient:
    def __init__(self, base_url: str = _DEFAULT_URL, timeout: float = 60.0):
        self.base_url = base_url
        self.timeout = timeout

    async def predict(self, cdr3_beta: str, v_beta: str, peptide: str,
                      mhc_allele: str, cdr3_alpha: str = "", v_alpha: str = "") -> PmtnetOmniResult:
        payload = {"cdr3b": cdr3_beta, "vb": v_beta, "peptide": peptide, "mhc": mhc_allele}
        if cdr3_alpha:
            payload["cdr3a"] = cdr3_alpha
        if v_alpha:
            payload["va"] = v_alpha
        tcr_repr = "paired" if cdr3_alpha else "beta_only"
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                resp = await client.post(f"{self.base_url}/api/predict", json=payload)
                resp.raise_for_status()
                data = resp.json()
                return PmtnetOmniResult(score=float(data.get("score", 0.0)),
                                        rank=float(data.get("rank", 0.0)),
                                        tcr_representation=tcr_repr)
        except (httpx.HTTPError, KeyError, ValueError) as e:
            logger.error("pMTnet Omni prediction failed: %s", e)
            return PmtnetOmniResult(score=0.0, tcr_representation=tcr_repr)
