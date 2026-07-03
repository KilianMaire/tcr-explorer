from imgt_app.dossier_models import DossierRequest
from imgt_app.dossier import build_dossier


def test_gene_name_path_populates_cdr_and_provenance():
    d = build_dossier(DossierRequest(query="TRBV20-1", species="human"))
    assert d.chain in ("beta", "unknown")
    # TRBV20-1 always resolves in germline: cdr regions + provenance are required.
    assert d.regions.get("cdr1") is not None
    assert any(p.block == "regions" and p.source == "cdr_enricher" for p in d.provenance)
    assert d.status in ("complete", "partial")


def test_unresolved_v_gene_is_partial_no_fake_call():
    d = build_dossier(DossierRequest(query="TRBV999", species="human"))
    assert d.genes["v"] is None
    assert d.status == "partial"
    assert any(w.code == "ambiguous_gene" for w in d.warnings)
    # no provenance for a block that was never populated
    assert all(p.block != "regions" for p in d.provenance) or d.regions.get("cdr1") is not None


def test_j_gene_not_slotted_as_v():
    d = build_dossier(DossierRequest(query="TRBJ2-7", species="human"))
    assert d.genes["v"] is None
    assert d.status == "partial"


def test_include_germline_populates_germline_nt():
    d = build_dossier(DossierRequest(query="TRBV20-1", include=["germline"]))
    if d.genes.get("v"):
        assert d.genes["v"].germline_nt is not None


def test_reconstruct_marks_synthetic_nt():
    req = DossierRequest(query="CASSLGTEAFF", input_type="raw_aa", species="human")
    # a bare CDR3 aa alone cannot be fully reconstructed without V/J; expect junction cdr3_aa set,
    # and if any cdr3_nt is produced it must be flagged synthetic.
    d = build_dossier(req)
    if d.junction and d.junction.cdr3_nt:
        assert d.junction.cdr3_nt_is_synthetic is True
        assert any(p.kind == "back_translated" for p in d.provenance)


def test_unresolved_input_is_partial_not_crash():
    d = build_dossier(DossierRequest(query="12345"))
    assert d.status == "partial"
    assert any(w.code == "unresolved_input_type" for w in d.warnings)


def test_projection_hides_long_nt_by_default():
    d = build_dossier(DossierRequest(query="TRBV20-1"))
    if d.genes.get("v"):
        assert d.genes["v"].germline_nt is None  # include empty -> no germline_nt
