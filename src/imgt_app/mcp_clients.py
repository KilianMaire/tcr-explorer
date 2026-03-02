from __future__ import annotations

import httpx

from .models import SearchRequest, SearchResponse


class ToolServerClient:
    def __init__(self, base_url: str) -> None:
        self.base_url = base_url.rstrip("/")

    async def search(self, req: SearchRequest) -> SearchResponse:
        async with httpx.AsyncClient(timeout=8.0) as client:
            r = await client.post(f"{self.base_url}/search", json=req.model_dump())
            r.raise_for_status()
            return SearchResponse(**r.json())
