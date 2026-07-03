# tests/test_input_router.py
from imgt_app.input_router import route

def test_override_wins():
    r = route("ACGTACGT", input_type="raw_aa")
    assert r.detected_type == "raw_aa"

def test_gene_and_allele():
    assert route("TRBV20-1", "auto").detected_type == "gene_name"
    assert route("TRBV20-1*01", "auto").detected_type == "allele"

def test_namespaced_id():
    r = route("vdjdb:12345", "auto")
    assert r.detected_type == "id" and r.source == "vdjdb"

def test_bare_integer_is_unresolved():
    r = route("12345", "auto")
    assert r.detected_type == "unknown"
    assert any(c == "unresolved_input_type" for c, _ in r.warnings)

def test_protein_with_signature_letter():
    # contains E and F -> unambiguously amino acids
    assert route("CASSEFGYT", "auto").detected_type == "raw_aa"

def test_clean_dna_is_nt():
    assert route("ACGTACGTACGTACGT", "auto").detected_type == "raw_nt"

def test_iupac_ambiguous_flags():
    # all letters in the aa alphabet, no aa-only signature letter, low ACGT fraction
    r = route("RYSWKMRYSWKM", "auto")
    assert any(c == "ambiguous_alphabet" for c, _ in r.warnings)

def test_unresolved_residue_x_is_protein():
    assert route("CASSXYEQYF", "auto").detected_type == "raw_aa"

def test_stop_marker_is_protein():
    assert route("CASS*", "auto").detected_type == "raw_aa"
