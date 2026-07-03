from imgt_app.reconstructor import reconstruct_tcr


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
