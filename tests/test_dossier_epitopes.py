from tcr_explorer.dossier_models import DossierRequest
from tcr_explorer.dossier import build_dossier
from tcr_explorer.models import IEDBHit
import tcr_explorer.dossier_epitopes as dossier_epitopes


def test_epitopes_flow_into_dossier():
    def fake_lookup(gene, cdr3_aa, species):
        return [IEDBHit(epitope_sequence="NLVPMVATV", mhc_allele="HLA-A*02:01", antigen_name="pp65")], 1

    d = build_dossier(DossierRequest(query="TRBV20-1"), epitope_lookup=fake_lookup)
    assert d.known_epitopes and d.known_epitopes[0].epitope_sequence == "NLVPMVATV"
    assert d.known_epitopes_total == 1
    assert any(p.block == "known_epitopes" for p in d.provenance)


def test_lookup_known_epitopes_no_gene_returns_empty():
    hits, total = dossier_epitopes.lookup_known_epitopes(None, None, "human")
    assert hits == []
    assert total == 0


def test_lookup_known_epitopes_defensive_when_search_raises(monkeypatch):
    def boom(req):
        raise ValueError("tool server unreachable")

    monkeypatch.setattr(dossier_epitopes, "search", boom)
    hits, total = dossier_epitopes.lookup_known_epitopes("TRBV20-1", None, "human")
    assert hits == []
    assert total == 0


def test_lookup_known_epitopes_defensive_when_no_event_loop_path_fails(monkeypatch):
    def fake_run_search(req):
        return None

    monkeypatch.setattr(dossier_epitopes, "_run_search", fake_run_search)
    hits, total = dossier_epitopes.lookup_known_epitopes("TRBV20-1", None, "human")
    assert hits == []
    assert total == 0


def test_resolve_id_defensive_on_failure(monkeypatch):
    def boom(req):
        raise ValueError("tool server unreachable")

    monkeypatch.setattr(dossier_epitopes, "search", boom)
    assert dossier_epitopes.resolve_id("vdjdb", "abc123") == {}


def test_resolve_id_unknown_source_returns_empty():
    assert dossier_epitopes.resolve_id("bogus", "abc123") == {}
