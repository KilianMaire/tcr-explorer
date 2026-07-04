from tcr_explorer.query_nl import parse_query


def test_parses_species_words_en_and_fr():
    assert parse_query("mouse CASSLGTEAFF")["species"] == "mouse"
    assert parse_query("souris CASSLGTEAFF")["species"] == "mouse"
    assert parse_query("human TRBV20-1")["species"] == "human"
    assert parse_query("humain CASSLGTEAFF")["species"] == "human"
    assert parse_query("CASSLGTEAFF")["species"] is None


def test_extracts_cdr3_amid_prose():
    p = parse_query("mouse TCR with CDR3 CASSGGTGEQYF please")
    assert p["species"] == "mouse" and p["cdr3_aa"] == "CASSGGTGEQYF"


def test_extracts_genes_and_id():
    p = parse_query("show TRBV20-1 and TRBJ2-7 CASSLGTEAFF")
    assert p["v_gene"] == "TRBV20-1" and p["j_gene"] == "TRBJ2-7" and p["cdr3_aa"] == "CASSLGTEAFF"
    assert parse_query("vdjdb:12345")["record_id"] == "vdjdb:12345"


def test_species_word_not_confused_with_sequence():
    # a bare CDR3 is not misread as a species; prose words are dropped
    p = parse_query("find similar to CASSLGTEAFF")
    assert p["cdr3_aa"] == "CASSLGTEAFF" and p["species"] is None
