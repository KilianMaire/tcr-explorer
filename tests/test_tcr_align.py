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
