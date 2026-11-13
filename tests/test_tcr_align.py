from tcr_explorer.germline_db import germline_alleles
from tcr_explorer.tcr_align import best_alleles, detect_alphabet, detect_chain


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
    from tcr_explorer.germline_db import Allele
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
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    chain_aa = built["full_chain_aa"]
    a = assign(chain_aa, species="mouse")
    assert a.chain == "TRB"
    assert a.v_determinable and any(n.startswith("TRBV19") for n in a.v_call["alleles"])
    assert any(n.startswith("TRBJ1-4") for n in a.j_call["alleles"])
    assert a.cdr3_aa == "CASSMADRKFF"
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])


def test_bare_cdr3_refuses_v_and_attaches_db_inference():
    from tcr_explorer.tcr_align import assign
    a = assign("CASSLGTEAFF", species="human")
    assert a.v_determinable is False
    assert a.v_call is None and a.v_reason
    assert a.j_call and any(n.startswith("TRBJ") for n in a.j_call["alleles"])
    # database inference is attached separately and labeled
    assert a.v_db_inference is not None


def test_nucleotide_input_roundtrip():
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    a = assign(built["full_nt"], species="mouse")
    assert a.input_kind == "nucleotide"
    assert a.v_determinable and a.cdr3_aa == "CASSMADRKFF"


def test_point_mutation_shows_in_region_breakdown():
    from tcr_explorer.germline_db import germline_alleles
    from tcr_explorer.tcr_align import assign
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
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
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
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
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
    from tcr_explorer.tcr_align import assign
    v = germline_alleles("human", "TRB", "V")[0]
    a = assign(v.aa, chain="TRB", species="human")
    assert a.constant_call is None or "low constant identity" in a.warnings


def test_full_chain_still_yields_a_trbc_constant_call():
    # The real full-chain path is unaffected in that it still identifies the
    # constant region as TRBC (the length gate keeps a large remainder). With
    # the mouse TRBC germline translated in its true reading frame, the call
    # is now a clean high-identity match, not a low-identity one propped up by
    # a tiny stop-free island.
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ1-4", "CASSMADRKFF", "mouse")
    a = assign(built["full_chain_aa"], species="mouse")
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])
    assert a.constant_call["identity"] >= 0.90, a.constant_call["identity"]
    assert "low constant identity" not in a.warnings


def test_clean_full_chain_constant_is_not_spuriously_warned():
    # A full-chain constant call that aligns at high identity carries no low
    # constant identity warning (the warning fires only on genuinely weak
    # evidence). The human TRBC reference resolves at identity 1.0 here.
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    a = assign(built["full_chain_aa"], species="human")
    assert a.constant_call and any(n.startswith("TRBC") for n in a.constant_call["alleles"])
    assert "low constant identity" not in a.warnings


def test_absent_d_for_alpha_is_explained_when_requested():
    # Alpha has no D segment. When the caller asks for D and none is available,
    # the result must say why rather than silently returning d_call None.
    from tcr_explorer.tcr_align import assign
    v = germline_alleles("human", "TRA", "V")[0]
    a = assign(v.aa, chain="TRA", species="human", want_d=True)
    assert a.d_call is None
    assert any("no D segment" in w for w in a.warnings)


def test_d_call_is_always_low_confidence_on_human_beta():
    # A human beta nucleotide chain with want_d must either return a D call
    # carrying the always-on low_confidence flag, or, when no D germline is
    # vendored, warn honestly about the absence. Never a confident D guess.
    from tcr_explorer.reconstructor import reconstruct_tcr
    from tcr_explorer.tcr_align import assign
    built = reconstruct_tcr("TRBV19", "TRBJ2-7", "CASSIRSSYEQYF", "human")
    a = assign(built["full_nt"], species="human", chain="TRB", want_d=True)
    if a.d_call is not None:
        assert a.d_call["low_confidence"] is True
    else:
        assert any("D" in w for w in a.warnings)


def test_chain_synonym_beta_is_normalized_not_silent_nulls():
    from tcr_explorer.tcr_align import assign
    # "beta" (a natural synonym) must resolve to TRB and produce a real call,
    # not the old all-null silent result.
    a = assign("DTAVYYCASSLGTEAFFGQGTRLTVV", species="human", chain="beta")
    assert a.chain == "TRB"
    assert (a.j_call is not None) or (a.v_call is not None)


def test_unrecognized_chain_warns_and_autodetects():
    from tcr_explorer.tcr_align import assign
    a = assign("DTAVYYCASSLGTEAFFGQGTRLTVV", species="human", chain="zeta")
    # falls back to auto-detection with an explicit warning, never silent nulls
    assert any("unrecognized chain" in w for w in a.warnings)


def test_normalize_chain_maps_codes_and_names():
    from tcr_explorer.tcr_align import _normalize_chain
    assert _normalize_chain("beta") == "TRB"
    assert _normalize_chain("TRA") == "TRA"
    assert _normalize_chain("delta") == "TRD"
    assert _normalize_chain("nonsense") is None


# --- J joined to the V, and CDR3 anchors, on synthetic rearrangements --------
# Every sequence below is assembled from the bundled germline plus a made-up
# junction, so the expected call and CDR3 are known by construction.

def _human(chain, segment, name):
    return next(a for a in germline_alleles("human", chain, segment) if a.name == name)


def _rearranged(v, j, junction_nt, j_trim_codons=1):
    """V through position 104 plus two codons, a made-up junction, then the J
    with its first codons trimmed. Returns (nt, offset of the codon at 104)."""
    from tcr_explorer.reconstructor import _v_cys_cut
    cut = _v_cys_cut(v.nt)
    j_frame = (len(j.nt) - 1) % 3  # a J-REGION ends one base into a codon
    return v.nt[: cut + 9] + junction_nt + j.nt[j_frame + 3 * j_trim_codons:], cut


_SPACER = "CTAGGTCACTCTGATTCCAGAGTCTAACATGCTAGGTCACT"  # made-up, not a J


def test_j_call_is_the_j_joined_to_the_v_not_a_downstream_j():
    # A transcript that reads through into a second, untrimmed J: the second J
    # scores higher than the recombined (5' trimmed) one, but the J that closes
    # the CDR3 is the one joined to the V.
    from tcr_explorer.tcr_align import assign
    v = _human("TRA", "V", "TRAV12-2*01")
    joined = _human("TRA", "J", "TRAJ31*01")
    downstream = _human("TRA", "J", "TRAJ30*01")
    vj, _cut = _rearranged(v, joined, "GGAGCC", j_trim_codons=2)
    seq = vj + _SPACER + downstream.nt + _human("TRA", "C", "TRAC*01").nt[:90]
    a = assign(seq, species="human", chain="TRA")
    assert a.j_call["alleles"] == ["TRAJ31*01"]
    assert a.j_call["identity"] == 1.0
    assert a.j_call["aligned_span"][1] <= len(vj)
    assert a.cdr3_aa == "CAVGANARLMF"
    assert not any("CDR3" in w for w in a.warnings)


def test_j_call_ignores_a_j_upstream_of_the_v():
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV9*01")
    joined = _human("TRB", "J", "TRBJ1-5*01")
    upstream = _human("TRB", "J", "TRBJ2-1*01")
    vj, _cut = _rearranged(v, joined, "GGGACAGGC", j_trim_codons=2)
    seq = upstream.nt + _SPACER + vj
    a = assign(seq, species="human", chain="TRB")
    assert a.j_call["alleles"] == ["TRBJ1-5*01"]
    assert a.cdr3_aa == "CASGTGQPQHF"


def test_j_call_on_an_input_that_stops_inside_the_j():
    # Only 17 nt of the J are present. The J joined to the V is still called,
    # and the CDR3 is refused because its closing anchor is not in the input.
    from tcr_explorer.tcr_align import assign
    v = _human("TRA", "V", "TRAV12-2*01")
    j = _human("TRA", "J", "TRAJ10*01")
    vj, cut = _rearranged(v, j, "GGAGCC", j_trim_codons=0)
    a = assign(vj[: cut + 15 + 17], species="human", chain="TRA")
    assert a.j_call["alleles"] == ["TRAJ10*01"]
    assert a.j_call["identity"] == 1.0
    assert a.cdr3_aa is None
    assert "CDR3 anchor could not be mapped" in a.warnings


def test_cdr3_closes_on_non_canonical_j_anchors():
    # TRAJ16 (F-A-R-G), TRAJ35 (Cys instead of Phe/Trp) and TRAJ58 (germline 5'
    # stop codon ahead of a canonical F-G-E-G) have no F/W-G-x-G motif behind a
    # stop-free prefix; the J coding frame comes from the J-REGION length.
    from tcr_explorer.reconstructor import _j_frame_and_fw, _translate
    from tcr_explorer.tcr_align import assign
    v = _human("TRA", "V", "TRAV12-2*01")
    for name, anchor, cdr3 in (
        ("TRAJ16*01", "F", "CAVGASDGQKLLF"),
        ("TRAJ35*01", "C", "CAVGAGFGNVLHC"),
        ("TRAJ58*01", "F", "CAVGAETSGSRLTF"),
    ):
        j = _human("TRA", "J", name)
        frame, fw = _j_frame_and_fw(j.nt)
        assert frame == (len(j.nt) - 1) % 3, name
        assert j.aa[(fw - frame) // 3] == anchor, name
        seq, _cut = _rearranged(v, j, "GGAGCC")
        a = assign(seq, species="human", chain="TRA")
        assert a.j_call["alleles"][0] == name
        assert a.cdr3_aa == cdr3, name
        # the amino acid level takes the same path
        assert assign(_translate(seq), species="human", chain="TRA").cdr3_aa == cdr3, name


def test_cdr3_survives_a_frameshift_upstream_in_the_v():
    # One inserted base in FR1 shifts the reading frame of everything upstream
    # of it. The junction itself is intact, so the CDR3 is still read.
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV9*01")
    seq, _cut = _rearranged(v, _human("TRB", "J", "TRBJ1-5*01"), "GGGACAGGC", j_trim_codons=2)
    shifted = seq[:40] + "A" + seq[40:]
    a = assign(shifted, species="human", chain="TRB")
    assert a.cdr3_aa == "CASGTGQPQHF"


def test_out_of_frame_junction_yields_no_cdr3():
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV9*01")
    seq, _cut = _rearranged(v, _human("TRB", "J", "TRBJ1-5*01"), "GGGACAGGCA", j_trim_codons=2)
    a = assign(seq, species="human", chain="TRB")
    assert a.cdr3_aa is None
    assert "V and J are out of frame; no CDR3" in a.warnings


def test_cdr3_start_for_a_v_allele_without_the_conserved_cys():
    # TRBV7-3*03 carries an Arg at position 104. The CDR3 starts at that same
    # position (taken from the alleles of the gene that kept the Cys), not at
    # the last codon of the V-REGION.
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV7-3*03")
    j = _human("TRB", "J", "TRBJ1-5*01")
    j_frame = (len(j.nt) - 1) % 3
    seq = v.nt[:285] + "GGGACAGGC" + j.nt[j_frame + 6:]
    a = assign(seq, species="human", chain="TRB")
    assert a.v_call["alleles"] == ["TRBV7-3*03"]
    assert a.cdr3_aa == "RASSGTGQPQHF"


def test_implausibly_long_cdr3_is_flagged():
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV9*01")
    seq, _cut = _rearranged(v, _human("TRB", "J", "TRBJ1-5*01"), "GGTACC" * 20, j_trim_codons=2)
    a = assign(seq, species="human", chain="TRB")
    assert len(a.cdr3_aa) == 48
    assert "implausible CDR3 length (48 residues)" in a.warnings


def test_indel_in_the_j_at_the_anchor_keeps_the_cdr3_and_warns():
    # One base inserted in the J just ahead of the Phe codon (T run): the
    # junction still reads Cys to Phe in the V frame, FR4 is shifted.
    from tcr_explorer.reconstructor import _j_frame_and_fw
    from tcr_explorer.tcr_align import assign
    v = _human("TRB", "V", "TRBV9*01")
    j = _human("TRB", "J", "TRBJ2-3*01")
    _frame, fw = _j_frame_and_fw(j.nt)
    seq, _cut = _rearranged(v, j, "GGGACAGGC", j_trim_codons=2)
    at = seq.index(j.nt[fw - 9 : fw + 6]) + 9
    a = assign(seq[:at] + "T" + seq[at:], species="human", chain="TRB")
    assert a.cdr3_aa == "CASGTGDTQYF"
    assert "frameshift in the J at the CDR3 anchor" in a.warnings
