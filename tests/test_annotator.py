from imgt_app import annotator
from imgt_app.annotator import annotate


def test_absent_igblast_falls_back_to_kmer(monkeypatch):
    monkeypatch.setattr(annotator, "igblast_available", lambda: False)
    a = annotate("ACGT" * 40, "human", is_protein=False, mode="full")
    assert a.source == "kmer_align"
    assert any(c == "igblast_unavailable" for c, _ in a.warnings)


def test_fast_mode_never_calls_igblast(monkeypatch):
    called = {"v": False}
    monkeypatch.setattr(annotator, "igblast_available", lambda: True)
    monkeypatch.setattr(annotator, "_run_igblast", lambda *a, **k: called.__setitem__("v", True) or None)
    annotate("ACGT" * 40, "human", is_protein=False, mode="fast")
    assert called["v"] is False  # fast mode uses k-mer, not igblast
