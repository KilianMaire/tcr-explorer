"""MixTCRpred client for epitope-specific TCR prediction.
Nature Communications 2024. Uses transformer encoder on curated 17,715 paired alpha-beta TCRs.
"""
from __future__ import annotations
import logging
from dataclasses import dataclass
import httpx

logger = logging.getLogger(__name__)
_DEFAULT_URL = "https://mixtcrpred.gfellerlab.org"

@dataclass
class MixTCRpredResult:
    score: float
    epitope: str = ""
    rank: float = 0.0

class MixTCRpredClient:
    def __init__(self, base_url: str = _DEFAULT_URL, timeout: float = 60.0):
        self.base_url = base_url
        self.timeout = timeout

    async def predict(self, cdr3_beta: str, v_beta: str, epitope: str,
                      cdr3_alpha: str = "", v_alpha: str = "") -> MixTCRpredResult:
        payload = {"cdr3b": cdr3_beta, "trbv": v_beta, "epitope": epitope}
        if cdr3_alpha:
            payload["cdr3a"] = cdr3_alpha
        if v_alpha:
            payload["trav"] = v_alpha
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                resp = await client.post(f"{self.base_url}/api/predict", json=payload)
                resp.raise_for_status()
                data = resp.json()
                return MixTCRpredResult(score=float(data.get("score", 0.0)),
                                        epitope=epitope,
                                        rank=float(data.get("rank", 0.0)))
        except (httpx.HTTPError, KeyError, ValueError) as e:
            logger.error("MixTCRpred prediction failed: %s", e)
            return MixTCRpredResult(score=0.0, epitope=epitope)
