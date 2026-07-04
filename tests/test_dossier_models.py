from tcr_explorer.dossier_models import (
    DossierRequest, TCRDossier, GeneCall, DossierWarning, WarningCode,
)

def test_request_defaults():
    r = DossierRequest(query="TRBV20-1")
    assert r.input_type == "auto"
    assert r.species == "human"
    assert r.mode == "fast"
    assert r.include == []

def test_dossier_minimal_valid():
    d = TCRDossier(
        schema_version="1.0", status="partial", summary="x",
        query_echo={"detected_type": "gene_name", "value": "TRBV20-1"},
        chain="beta", species="human",
        genes={"v": None, "d": None, "j": None, "c": None},
        regions={}, junction=None, full_sequence=None,
        known_epitopes=[], known_epitopes_total=0,
        provenance=[], warnings=[DossierWarning(code="partial_annotation", block=None, message="m")],
    )
    assert d.warnings[0].code == "partial_annotation"

def test_warning_code_is_constrained():
    import pydantic, pytest
    with pytest.raises(pydantic.ValidationError):
        DossierWarning(code="not_a_real_code", message="m")
