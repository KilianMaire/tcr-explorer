from tcr_explorer.germline_db import germline_alleles, Allele
from tcr_explorer import d_regions


def test_v_catalog_has_alleles_with_nt_and_aa():
    vs = germline_alleles("human", "TRA", "V")
    assert len(vs) >= 100
    a = next(x for x in vs if x.name == "TRAV1-1*01")
    assert set(a.nt) <= set("ACGTN") and len(a.nt) > 200
    assert a.aa and set(a.aa) <= set("ACDEFGHIKLMNPQRSTVWY*X")
    assert isinstance(a, Allele)


def test_j_aa_is_in_coding_frame():
    js = germline_alleles("human", "TRB", "J")
    j = next(x for x in js if x.name.startswith("TRBJ2-7"))
    # the conserved FGxG motif must appear in the coding-frame translation
    assert "FG" in j.aa, j.aa


def test_d_alleles_human_only():
    ds = germline_alleles("human", "TRB", "D")
    assert {d.name for d in ds} >= {"TRBD1*01", "TRBD2*01"}
    assert germline_alleles("mouse", "TRB", "D") == []
    assert d_regions.d_alleles("mouse") == {}


def test_constant_catalog_distinguishes_trbc1_trbc2():
    cs = germline_alleles("human", "TRB", "C")
    names = {c.name for c in cs}
    assert any(n.startswith("TRBC1") for n in names)
    assert any(n.startswith("TRBC2") for n in names)


def test_mouse_trbc1_aa_is_translated_in_its_true_reading_frame():
    # The stitchr mouse TRBC1 C-REGION nucleotide sequence is not in frame 0;
    # translating naively in frame 0 yields a stop-riddled garbage protein.
    # The true open reading frame (frame 2) yields the real constant protein.
    cs = germline_alleles("mouse", "TRB", "C")
    c = next(x for x in cs if x.name.startswith("TRBC1"))
    assert c.aa.count("*") == 0, c.aa
    assert c.aa.startswith("DLRNVTPP"), c.aa


def test_human_trbc1_aa_is_translated_in_its_true_reading_frame():
    cs = germline_alleles("human", "TRB", "C")
    c = next(x for x in cs if x.name.startswith("TRBC1"))
    assert c.aa.count("*") == 0, c.aa
    assert "DLNKVFPP" in c.aa, c.aa


def test_absent_segment_returns_empty():
    assert germline_alleles("human", "TRA", "D") == []  # alpha has no D
