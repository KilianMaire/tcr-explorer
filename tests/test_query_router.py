from tcr_explorer.query_router import route_query


def test_gene_routes_to_records():
    r = route_query("TRBV20-1")
    assert r.tools == ["records"] and r.blocks[0].tool == "records"


def test_id_routes_to_records():
    r = route_query("vdjdb:c1")
    assert r.tools == ["records"]


def test_cdr3_routes_to_records_and_assign():
    r = route_query("CASSLGTEAFF", species="human")
    assert r.tools == ["records", "assign"]
    assert [b.tool for b in r.blocks] == ["records", "assign"]
    # the assign block is the real assignment output (bare CDR3 refuses V)
    assert r.blocks[1].data["v_determinable"] is False


def test_full_chain_aa_routes_to_assign_only():
    from tcr_explorer.reconstructor import reconstruct_tcr
    chain = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")["full_chain_aa"]
    r = route_query(chain, species="mouse")
    assert r.tools == ["assign"]
    assert r.blocks[0].data["cdr3_aa"] == "CASSMADRKFF"


def test_nt_routes_to_assign():
    r = route_query("ACGTACGTACGTACGTACGTACGT")
    assert r.tools == ["assign"] and r.detected_type == "raw_nt"


def test_phrase_routes_to_records_with_species():
    r = route_query("mouse CASSLGTEAFF TRBJ2-7")
    assert r.tools == ["records"] and r.species == "mouse"


def test_phrase_explicit_species_overrides_contradictory_word_in_records():
    # dropdown species="human" but the phrase itself says "mouse": the
    # explicit dropdown must win both in the understood banner AND in what
    # records actually filters by (not just the displayed species).
    r = route_query("mouse CASSGGTGEQYF", species="human")
    assert r.tools == ["records"] and r.species == "human"
    assert r.blocks[0].data["query_echo"]["species"] == "human"


def test_free_text_falls_to_ask():
    r = route_query("which TCRs recognize influenza")
    assert r.tools == ["ask"]


def test_force_restricts_to_one_tool():
    r = route_query("CASSLGTEAFF", species="human", force="records")
    assert r.tools == ["records"] and len(r.blocks) == 1


def test_force_similar_without_genes_warns():
    r = route_query("CASSLGTEAFF", force="similar")
    assert r.blocks == [] and any("similarity needs" in w for w in r.warnings)
