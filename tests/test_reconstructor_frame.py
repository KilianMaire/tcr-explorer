from tcr_explorer.reconstructor import reconstruct_tcr


def test_reconstruction_is_in_frame_and_contains_cdr3():
    r = reconstruct_tcr("TRBV4-1", "TRBJ1-1", "CASSLGTEAFF", "human")
    aa = r["full_aa"]
    assert aa is not None
    # the CDR3 reads cleanly inside the full-length protein
    assert "CASSLGTEAFF" in aa, aa
    # no internal stop codons (a clean coding sequence)
    assert "*" not in aa, aa
    # the J framework (FR4) follows the CDR3: TRBJ1-1 FR4 is GQGTRLTVV-ish
    assert "GTRLTVV" in aa, aa


def test_reconstruction_frame_holds_for_a_second_vj():
    r = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    aa = r["full_aa"]
    assert aa and "CASSIRSSYEQYF" in aa and "*" not in aa


def test_mouse_reconstruction_in_frame():
    r = reconstruct_tcr("TRBV13-2", "TRBJ2-7", "CASGGTGEQYF", "mouse")
    aa = r["full_aa"]
    # mouse germline may or may not be present; if reconstructed, it must be clean
    if aa is not None:
        assert "CASGGTGEQYF" in aa and "*" not in aa


def test_no_spurious_residue_at_cdr3_j_junction():
    # TRBJ1-1 FR4 is F-G-x-G where that F is Phe118, the CDR3's own last F.
    # The junction must not duplicate it into a triple F (CASSLGTEAF-F-F-GQG).
    r = reconstruct_tcr("TRBV4-1", "TRBJ1-1", "CASSLGTEAFF", "human")
    aa = r["full_aa"]
    assert "CASSLGTEAFFGQGTRLTVV" in aa, aa
    assert "AFFF" not in aa, aa


def test_composition_pieces_are_in_frame():
    # The composition germline pieces must translate cleanly (no stop codons),
    # not from the raw region nt which are not frame 0.
    from tcr_explorer.records import build_record
    row = dict(source="vdjdb", source_record_id="x", pairing_key="p", external_url="u",
               chain="beta", species="human", cdr3_aa="CASSLGTEAFF",
               v_gene="TRBV4-1", j_gene="TRBJ1-1")
    c = build_record(row).composition
    assert c is not None
    assert "*" not in (c.v_germline_aa or "") and "*" not in (c.j_germline_aa or "")
    assert c.j_germline_aa and c.j_germline_aa.startswith("GQGTRLTVV")
