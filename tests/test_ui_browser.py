"""Real-browser end-to-end test for the /ui page.

This guards the failure mode that shipped once and slipped every other test: a
JavaScript syntax error in the inline UI script that broke all page interactivity
while the HTML-only unit test stayed green. It launches uvicorn and drives the page
with a headless Chromium, asserting there are no console errors and that the forms
actually fetch and render.

Skipped automatically when Playwright or the Chromium browser is not installed, so
it never breaks environments that lack them.
"""
from __future__ import annotations

import os
import socket
import subprocess
import time
from pathlib import Path

import httpx
import pytest

sync_playwright = None
try:  # pragma: no cover - import guard
    from playwright.sync_api import sync_playwright as _spw
    sync_playwright = _spw
except Exception:  # pragma: no cover
    sync_playwright = None


def _chromium_available() -> bool:
    if sync_playwright is None:
        return False
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch()
            browser.close()
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _chromium_available(), reason="playwright chromium not installed"
)

_ROOT = Path(__file__).resolve().parent.parent


def _free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


@pytest.fixture(scope="module")
def server():
    port = _free_port()
    env = dict(os.environ)
    env["PYTHONPATH"] = str(_ROOT / "src")
    idx = _ROOT / "data" / "unitcr_beta_index.parquet"
    if idx.exists():
        env["UNITCR_INDEX_PATH"] = str(idx)
    proc = subprocess.Popen(
        [str(_ROOT / ".venv" / "bin" / "uvicorn"), "imgt_app.api:app",
         "--port", str(port), "--host", "127.0.0.1"],
        cwd=str(_ROOT), env=env,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    base = f"http://127.0.0.1:{port}"
    try:
        for _ in range(50):
            try:
                if httpx.get(f"{base}/health", timeout=1).status_code == 200:
                    break
            except Exception:
                time.sleep(0.2)
        else:
            raise RuntimeError("uvicorn did not start")
        yield base
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            proc.kill()


def _run(base):
    errors = []
    result = {}
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{base}/ui")
        # a gene query renders a real dossier
        page.fill("#q", "TRBV20-1")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('TRBV20-1')", timeout=8000)
        result["gene"] = page.locator("#out").inner_text()
        # a bare CDR3 renders the guidance hint (not a blank result)
        page.fill("#q", "CASSLGTEAFF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('intent')", timeout=8000)
        result["cdr3"] = page.locator("#out").inner_text()
        browser.close()
    return errors, result


def test_ui_has_no_console_errors_and_forms_work(server):
    errors, result = _run(server)
    assert errors == [], f"UI console errors: {errors}"
    assert "chain: beta" in result["gene"] and "TRBV20-1" in result["gene"]
    assert "cannot identify V/D/J" in result["cdr3"], \
        "bare CDR3 should render the guidance hint, not a blank result"


def test_ui_align_form_renders_colored_aa_nt(server):
    if not (_ROOT / "data" / "unitcr_beta_index.parquet").exists():
        pass  # index not needed for align; germline is
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.select_option("#a_sp", "mouse")
        page.fill("#a_chain", "TRB")
        page.select_option("#a_seg", "J")
        page.check("#a_translate")
        page.click("#af button")
        page.wait_for_function(
            "document.querySelector('#a_out').innerText.includes('engine')", timeout=8000)
        txt = page.locator("#a_out").inner_text()
        html = page.locator("#a_out").inner_html()
        browser.close()
    assert errors == [], f"align UI console errors: {errors}"
    assert "view aa_nt" in txt, "translated germline set should render the codon-aware aa+nt view"
    assert " aa  " in txt and " nt  " in txt, "both an aa row and an nt row should render"
    assert "background:#" in html, "conservation coloring spans should be present"
