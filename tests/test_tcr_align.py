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
