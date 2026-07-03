import pytest
from imgt_app.germline_db import germline_alleles
from imgt_app.tcr_align import detect_alphabet, best_alleles, detect_chain


def test_detect_alphabet():
    assert detect_alphabet("ACGTACGTNN") == "nucleotide"
    assert detect_alphabet("CASSLGTEAFF") == "amino_acid"


def test_self_identity_is_100_percent_and_names_the_allele():
    vs = germline_alleles("human", "TRA", "V")
    target = next(v for v in vs if v.name == "TRAV1-1*01")
    call = best_alleles(target.nt, vs, "nt")
    assert "TRAV1-1*01" in call.alleles
    assert call.identity == 1.0


def test_point_mutation_keeps_call_below_100():
    vs = germline_alleles("human", "TRB", "V")
    target = next(v for v in vs if v.name.startswith("TRBV19"))
    mutated = target.nt[:30] + ("A" if target.nt[30] != "A" else "C") + target.nt[31:]
    call = best_alleles(mutated, vs, "nt")
    assert target.name.split("*")[0] in call.alleles[0]
    assert 0.9 < call.identity < 1.0


def test_ties_are_all_reported():
    # two synthetic alleles identical over the query span both survive
    from imgt_app.germline_db import Allele
    q = "ACGTACGTACGTACGTACGT"
    alleles = [
        Allele("FAKEV1*01", q + "AAA", ""),
        Allele("FAKEV1*02", q + "GGG", ""),  # differs only outside the query span
    ]
    call = best_alleles(q, alleles, "nt")
    assert set(call.alleles) == {"FAKEV1*01", "FAKEV1*02"}


def test_chain_detection_does_not_misclassify_alpha_as_beta():
    vs = germline_alleles("human", "TRA", "V")
    target = next(v for v in vs if v.name == "TRAV1-1*01")
    assert detect_chain(target.nt, "human", "nt") == "TRA"


def test_full_chain_roundtrip_recovers_v_j_and_cdr3():
    # a full reconstructed mouse beta chain from the oracle fixture (row 2)
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    chain_aa = built["full_chain_aa"]
    a = assign(chain_aa, species="mouse")
    assert a.chain == "TRB"
    assert a.v_determinable and any(n.startswith("TRBV19") for n in a.v_call["alleles"])
    assert any(n.startswith("TRBJ1-4") for n in a.j_call["alleles"])
    assert a.cdr3_aa == "CASSMADRKFF"
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])


def test_bare_cdr3_refuses_v_and_attaches_db_inference():
    from imgt_app.tcr_align import assign
    a = assign("CASSLGTEAFF", species="human")
    assert a.v_determinable is False
    assert a.v_call is None and a.v_reason
    assert a.j_call and any(n.startswith("TRBJ") for n in a.j_call["alleles"])
    # database inference is attached separately and labeled
    assert a.v_db_inference is not None


def test_nucleotide_input_roundtrip():
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    a = assign(built["full_nt"], species="mouse")
    assert a.input_kind == "nucleotide"
    assert a.v_determinable and a.cdr3_aa == "CASSMADRKFF"


def test_point_mutation_shows_in_region_breakdown():
    from imgt_app.germline_db import germline_alleles
    from imgt_app.tcr_align import assign
    v = next(x for x in germline_alleles("human", "TRB", "V") if x.name.startswith("TRBV19"))
    # mutate a CDR1 position (IMGT aa 27..38 -> nt 78..114) and assign
    nt = v.nt
    pos = 27 * 3  # into CDR1
    mutated = nt[:pos] + ("A" if nt[pos] != "A" else "C") + nt[pos + 1:]
    a = assign(mutated, species="human", chain="TRB")
    assert a.regions.get("CDR1", 1.0) < 1.0
    assert a.regions.get("FR1", 0.0) == 1.0


def test_partial_input_reading_frame_is_derived_from_target_offset():
    # A nucleotide fragment that starts a few codons into the V region aligns
    # at target start tstart > 0, with tstart % 3 != 0. The correct query
    # reading frame is (qstart - tstart) % 3, not qstart % 3; the old code used
    # qstart % 3 and translated the fragment out of frame, so the CDR3 could not
    # be mapped. Dropping the first 4 nt shifts the natural frame by 2.
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    partial = built["full_nt"][4:]  # tstart == 4 -> tstart % 3 == 1
    a = assign(partial, species="human", chain="TRB")
    # Amino acid content is recovered only with the corrected frame; under
    # qstart % 3 the translation is out of frame and cdr3_aa is None.
    assert a.cdr3_aa == "CASSIRSSYEQYF"
    assert a.v_determinable
    # Covered regions read near identity (nt-level, so the mutation-free
    # fragment stays high for every region it actually covers).
    for name, ident in a.regions.items():
        assert ident > 0.9, (name, ident)


def test_uncovered_regions_are_omitted_not_zeroed():
    # A fragment starting inside FR2 leaves FR1 and CDR1 with no covered
    # positions. They must be omitted from the regions dict, not reported as
    # 0.0 (which would read as fully divergent rather than not covered).
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    frag = built["full_nt"][120:]  # begins past FR1/CDR1
    a = assign(frag, species="human", chain="TRB")
    assert "FR1" not in a.regions
    assert "CDR1" not in a.regions
    # The covered regions are still present and high identity.
    assert a.regions.get("FR3", 0.0) > 0.9


def test_v_only_input_emits_no_confident_constant_call():
    # A V-region-only amino acid input contains no constant region. The 3'
    # remainder past the (spurious) J hit is the V-region body, so a naive
    # constant call aligns a tiny island at ~66.7% identity: a spurious
    # confident constant call with zero warnings under the old code. After the
    # fix the call must never be emitted silently as confident: either it is
    # gated out (short remainder) or, when emitted, it carries the low constant
    # identity warning. This asserts the honesty invariant (no silent confident
    # spurious call), which fails under the old code and passes after the fix.
    from imgt_app.tcr_align import assign
    v = germline_alleles("human", "TRB", "V")[0]
    a = assign(v.aa, chain="TRB", species="human")
    assert a.constant_call is None or "low constant identity" in a.warnings


def test_full_chain_still_yields_a_trbc_constant_call():
    # The real full-chain path is unaffected in that it still identifies the
    # constant region as TRBC (the length gate keeps a large remainder).
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    a = assign(built["full_chain_aa"], species="mouse")
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])


def test_clean_full_chain_constant_is_not_spuriously_warned():
    # A full-chain constant call that aligns at high identity carries no low
    # constant identity warning (the warning fires only on genuinely weak
    # evidence). The human TRBC reference resolves at identity 1.0 here.
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    a = assign(built["full_chain_aa"], species="human")
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])
    assert "low constant identity" not in a.warnings


def test_absent_d_for_alpha_is_explained_when_requested():
    # Alpha has no D segment. When the caller asks for D and none is available,
    # the result must say why rather than silently returning d_call None.
    from imgt_app.tcr_align import assign
    v = germline_alleles("human", "TRA", "V")[0]
    a = assign(v.aa, chain="TRA", species="human", want_d=True)
    assert a.d_call is None
    assert any("no D segment" in w for w in a.warnings)


def test_d_call_is_always_low_confidence_on_human_beta():
    # A human beta nucleotide chain with want_d must either return a D call
    # carrying the always-on low_confidence flag, or, when no D germline is
    # vendored, warn honestly about the absence. Never a confident D guess.
    from imgt_app.reconstructor import reconstruct_tcr
    from imgt_app.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    a = assign(built["full_nt"], species="human", chain="TRB", want_d=True)
    if a.d_call is not None:
        assert a.d_call["low_confidence"] is True
    else:
        assert any("D" in w for w in a.warnings)
