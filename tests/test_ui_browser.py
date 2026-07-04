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
        [str(_ROOT / ".venv" / "bin" / "uvicorn"), "tcr_explorer.api:app",
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
        # a gene query routes to the records tool through the central query box
        page.fill("#q", "TRBV20-1")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('understood as')", timeout=8000)
        result["gene"] = page.locator("#out").inner_text()
        # a bare CDR3 routes to both known records and a germline assignment
        page.fill("#q", "CASSLGTEAFF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('Exact records')", timeout=8000)
        result["cdr3"] = page.locator("#out").inner_text()
        browser.close()
    return errors, result


def test_ui_has_no_console_errors_and_forms_work(server):
    errors, result = _run(server)
    assert errors == [], f"UI console errors: {errors}"
    assert "understood as: gene_name -> records" in result["gene"] \
        or "understood as: allele -> records" in result["gene"]
    assert "Exact records" in result["gene"] and "Near neighbours" in result["gene"]
    assert "understood as: raw_aa -> records, assign" in result["cdr3"], \
        "a bare CDR3 should route to both the records and the assign tools"
    assert "Exact records" in result["cdr3"], \
        "the CDR3 lookup should render the known-records card"


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
        page.click('[data-tool="align"]')  # the align form is hidden until its chip is clicked
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


def test_ui_natural_language_mouse_query_echoes_detected_species(server):
    """Typing natural language ("mouse <CDR3>") into the single query box
    drives a mouse-scoped search: the transparency line names the species
    the router detected from free text, alongside the detected type and the
    routed tool(s)."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.fill("#q", "mouse CASGGTGEQYF")
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('understood as')", timeout=8000)
        out_text = page.locator("#out").inner_text()
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert "understood as:" in out_text
    assert "species: mouse" in out_text.lower(), \
        f"expected the transparency line to name the detected species mouse, got: {out_text}"


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
        page.click('[data-tool="reconstruct"]')  # the reconstruct form is hidden until its chip is clicked
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
    from tcr_explorer.reconstructor import reconstruct_tcr

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
        page.click('[data-tool="reconstruct"]')  # the reconstruct form is hidden until its chip is clicked
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
        page.click('[data-tool="reconstruct"]')  # the reconstruct form is hidden until its chip is clicked
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


def test_ui_sequence_query_routes_to_assignment_card(server):
    """A full chain pasted into the single query box (Task 3's central box)
    routes to the assign tool, and its germline-assignment card renders in
    the adaptive result area."""
    from tcr_explorer.reconstructor import reconstruct_tcr

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
        page.select_option("#sp", "mouse")
        page.fill("#q", full_chain_aa)
        page.click("#f button")
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('CASSMADRKFF')", timeout=8000)
        text = page.locator("#out").inner_text()
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert "understood as: raw_aa -> assign" in text, \
        "a full chain should route to the assign tool alone"
    assert "TRBV19" in text, "the assigned V allele should render"
    assert "CASSMADRKFF" in text, "the extracted CDR3 should render"


def test_ui_cdr3_query_renders_two_cards(server):
    """A bare CDR3 typed into the single query box yields two blocks
    (records then assign) and both must render in the adaptive result area."""
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
    assert "understood as: raw_aa -> records, assign" in text
    assert "Exact records" in text, "the records card should render"
    assert "not determinable" in text.lower() or "CASSLGTEAFF" in text, \
        "the assign card should also render (a bare CDR3 refuses a V call)"


def test_ui_override_chip_reruns_query_forcing_tool(server):
    """Clicking the records override chip re runs the current query forcing
    the records tool, and the forced block renders."""
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.on("console", lambda m: errors.append(m.text)
                if m.type == "error" and "favicon" not in m.text else None)
        page.goto(f"{server}/ui")
        page.fill("#q", "CASSLGTEAFF")
        page.click('[data-tool="records"]')
        page.wait_for_function(
            "document.querySelector('#out').innerText.includes('Exact records')", timeout=8000)
        text = page.locator("#out").inner_text()
        browser.close()
    assert errors == [], f"UI console errors: {errors}"
    assert "understood as: raw_aa -> records" in text, \
        "the forced tool should be reported in the transparency line"
    assert "Exact records" in text


def test_ui_reconstruct_chip_reveals_hidden_form(server):
    """The reconstruct chip reveals the (previously hidden) #rcf form
    instead of routing through /v1/tcr/query."""
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(f"{server}/ui")
        assert not page.locator("#rcf").is_visible(), \
            "the reconstruct form should be hidden before its chip is clicked"
        page.click('[data-tool="reconstruct"]')
        assert page.locator("#rcf").is_visible(), \
            "the reconstruct form should be revealed after its chip is clicked"
        browser.close()


def test_ui_onboarding_block_renders_artifacts_and_copy_button(server):
    """The assistant onboarding block renders both copy-paste artifacts (the
    MCP config JSON and the install prompt) and at least one copy button."""
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(f"{server}/ui")
        text = page.locator("body").inner_text()
        assert "Ask in plain English, use your own AI assistant" in text
        assert '"mcpServers"' in text and "tcr-explorer-mcp" in text, \
            "the MCP config JSON artifact should render"
        assert "retrieve_tcr_records" in text and "align_tcr_genes" in text, \
            "the install prompt artifact should render"
        assert page.locator(".copy-btn").count() >= 1, \
            "at least one copy button should be present"
        browser.close()
