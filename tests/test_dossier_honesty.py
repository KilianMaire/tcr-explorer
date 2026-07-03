from imgt_app.dossier_models import DossierRequest
from imgt_app.dossier import build_dossier


def test_every_populated_block_has_provenance():
    d = build_dossier(DossierRequest(query="TRBV20-1"))
    if d.regions.get("cdr1"):
        assert any(p.block == "regions" for p in d.provenance)
    if d.genes.get("v") and d.genes["v"].call:
        assert any(p.block in ("annotation", "regions") for p in d.provenance)


def test_include_sequences_toggles_projection():
    d0 = build_dossier(DossierRequest(query="TRBV20-1"))
    d1 = build_dossier(DossierRequest(query="TRBV20-1", include=["sequences", "germline"]))
    if d1.genes.get("v"):
        assert d1.genes["v"].germline_nt is not None or d0.genes["v"] is None


def test_no_unlisted_warning_codes():
    d = build_dossier(DossierRequest(query="12345"))
    valid = {"igblast_unavailable","source_unavailable","ambiguous_gene","ambiguous_alphabet",
             "unresolved_input_type","d_segment_unresolved","aa_annotation_limited",
             "back_translated_nt","partial_annotation","timeout"}
    assert all(w.code in valid for w in d.warnings)
