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
    records_idx = _ROOT / "data" / "records_index.parquet"
    if records_idx.exists():
        env["RECORDS_INDEX_PATH"] = str(records_idx)
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
        # a bare CDR3 is looked up against the known-TCR reference (similarity),
        # not a dead-end germline annotation
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
    assert "intent: similar" in result["cdr3"], \
        "a bare CDR3 should route to a similarity lookup, not germline annotation"
    assert "Similar TCRs" in result["cdr3"], \
        "the CDR3 lookup should render matching known TCRs"


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
        page.click("details.advanced summary")  # align panel is collapsed by default
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


def test_ui_renders_record_cards(server):
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.fill("#q", "CASSLGTEAFF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('Exact records')", timeout=8000)
        text = page.locator("#out").inner_text()
        html = page.locator("#out").inner_html()
        browser.close()
    assert errors == [], f"record card console errors: {errors}"
    assert "Exact records" in text
    assert "Near neighbours" in text
    assert any(badge in text for badge in ("VDJDB", "IEDB", "MCPAS", "TCR3D")), \
        "a source badge should be rendered"
    assert 'href="http' in html, "a record card should link back to its source"
    if "nt:" in text:
        assert ("deposited" in text) or ("reconstructed" in text), \
            "a shown nt line should be labeled with its provenance"


def test_ui_bare_cdr3_query_shows_exact_records_heading(server):
    """A bare CDR3 typed into the single query box renders cards under an
    'Exact records' heading (this is the natural-language single-box path,
    not a dedicated CDR3 form)."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.fill("#q", "CASSLGTEAFF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('Exact records')", timeout=8000)
        text = page.locator("#out").inner_text()
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert "Exact records" in text


def test_ui_reconstructed_card_contains_queried_cdr3_no_stop_codon(server):
    """A reconstructed full_aa card must contain the queried CDR3 verbatim and
    never a stop codon marker, proving the in-frame reconstruction fix (Task 1)
    reaches the UI."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.fill("#q", "CASSLGTEAFF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('Exact records')", timeout=8000)
        reconstructed_aa = page.evaluate("""() => {
            const cards = [...document.querySelectorAll('#out .rec')];
            for (const c of cards) {
                const codes = [...c.querySelectorAll('code')];
                for (const code of codes) {
                    const next = code.nextElementSibling;
                    if (next && next.classList.contains('kind')
                        && next.innerText.includes('reconstructed')) {
                        return code.innerText;
                    }
                }
            }
            return null;
        }""")
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert reconstructed_aa is not None, "expected at least one reconstructed record card"
    assert "CASSLGTEAFF" in reconstructed_aa
    assert "*" not in reconstructed_aa, \
        f"the reconstructed full_aa must contain no internal stop codon, got: {reconstructed_aa}"


def test_ui_natural_language_mouse_query_hides_hla_until_toggled(server):
    """Typing natural language ("mouse <CDR3>") into the single query box
    drives a mouse-scoped search: the page echoes the detected species, and
    with the cross-species-MHC toggle off no rendered record shows an HLA
    allele (a human-organism MHC on a mouse record is filtered server-side
    unless the toggle is checked)."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        assert page.locator("#xmhc").is_checked() is False, \
            "the HLA-transgenic toggle must default to off"
        page.fill("#q", "mouse CASGGTGEQYF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#echo').innerText.length > 0", timeout=8000)
        echo_text = page.locator("#echo").inner_text()
        out_text = page.locator("#out").inner_text()
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert "mouse" in echo_text.lower(), \
        f"expected the echoed understanding to name the detected species mouse, got: {echo_text}"
    assert "HLA" not in out_text, \
        "no record card should show an HLA allele while the transgenic toggle is off"


def test_ui_reconstruction_panel_builds_full_chain(server):
    """When both V and J overrides are filled, the panel still calls the
    direct /reconstruct path unchanged: full membrane-bound chain, the
    germline allele used, and the constant-region provenance (an
    oracle-validated mouse example)."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.select_option("#rc_sp", "mouse")
        page.fill("#rc_v", "TRBV19")
        page.fill("#rc_j", "TRBJ1-4")
        page.fill("#rc_seq", "CASSMADRKFF")
        page.click("#rcf button")
        page.wait_for_function(
            "document.querySelector('#rc_out').innerText.includes('full chain')", timeout=8000)
        text = page.locator("#rc_out").inner_text()
        browser.close()
    assert errors == [], f"reconstruction UI console errors: {errors}"
    assert "CASSMADRKFF" in text, "the CDR3 should appear in the reconstructed chain"
    assert "EDLRNVTPP" in text, "the mouse constant region should be appended"
    assert "*01" in text, "the germline allele used should be reported"
    assert "oracle-validated" in text, "the constant-region provenance should be shown"


def test_ui_assign_panel_identifies_full_chain_with_blank_overrides(server):
    """Pasting a full chain with V and J left blank routes to
    POST /v1/tcr/assign (not /reconstruct): the panel reports the assigned V
    allele, the extracted CDR3, and a per-region identity row."""
    from imgt_app.reconstructor import reconstruct_tcr

    full_chain_aa = reconstruct_tcr(
        "TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse"
    )["full_chain_aa"]

    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.select_option("#rc_sp", "mouse")
        page.fill("#rc_v", "")
        page.fill("#rc_j", "")
        page.fill("#rc_seq", full_chain_aa)
        page.click("#rcf button")
        page.wait_for_function(
            "document.querySelector('#rc_out').innerText.includes('CASSMADRKFF')", timeout=8000)
        text = page.locator("#rc_out").inner_text()
        browser.close()
    assert errors == [], f"assign UI console errors: {errors}"
    assert "TRBV19" in text, "the assigned V allele should be reported"
    assert "CASSMADRKFF" in text, "the extracted CDR3 should be reported"
    assert "FR1" in text or "CDR1" in text, "a per region identity row should render"


def test_ui_assign_panel_refuses_bare_cdr3_with_db_inference(server):
    """A bare CDR3 (no framework) cannot be assigned a V allele: the panel
    shows the refusal reason plus the database frequency inference, under a
    heading clearly separate from a real germline call."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.select_option("#rc_sp", "human")
        page.fill("#rc_v", "")
        page.fill("#rc_j", "")
        page.fill("#rc_seq", "CASSLGTEAFF")
        page.click("#rcf button")
        page.wait_for_function(
            "document.querySelector('#rc_out').innerText.includes('database frequency inference')",
            timeout=8000)
        text = page.locator("#rc_out").inner_text()
        browser.close()
    assert errors == [], f"assign refusal UI console errors: {errors}"
    assert "not determinable" in text.lower(), "the V refusal reason should be shown"
    assert "database frequency inference" in text, \
        "the weaker database inference should be under its own heading"
    assert "TRBV4-1" in text, "the inferred V/J candidates should be listed"
