import csv
from pathlib import Path
from tcr_explorer.reconstructor import reconstruct_tcr

FIX = Path(__file__).parent / "fixtures" / "reconstruction" / "mouse_paired_tcrs.tsv"
TRBC = "EDLRNVTPP"


def _rows():
    return list(csv.reader(open(FIX), delimiter="\t"))


def _var(full, motif):
    i = full.find(motif)
    return full[:i] if i >= 0 else full


def test_alpha_delta_dv_gene_names_resolve():
    r = reconstruct_tcr("TRAV6-7-DV9", "TRAJ40", "CALFNTGNYKYVF", "mouse")
    assert r["v_found"] is True and r["full_aa"]
    r2 = reconstruct_tcr("TRAV14D-3-DV8", "TRAJ21", "CAASAGANYNVLYF", "mouse")
    assert r2["v_found"] is True and r2["full_aa"]


def test_trbj1_6_no_overlap_duplication():
    r = reconstruct_tcr("TRBV13-3", "TRBJ1-6", "CASSDRYNSPLYF", "mouse")
    aa = r["full_aa"]
    assert "CASSDRYNSPLYFAAGTRLTVT" in aa, aa
    assert "PLYFSYNSPLYF" not in aa, aa


def test_all_beta_variable_domains_reproduce_ground_truth():
    hits = 0
    for row in _rows():
        trbv, trbj, trbcdr3, trb_full = row[0], row[1], row[2], row[6]
        got = reconstruct_tcr(trbv, trbj, trbcdr3, "mouse")["full_aa"] or ""
        if got == _var(trb_full, TRBC):
            hits += 1
    assert hits == 5, f"only {hits}/5 beta variable domains matched"


def test_full_beta_chains_reproduce_ground_truth():
    # variable domain (V+CDR3+J) plus the vendored constant reproduces the
    # complete membrane-bound beta chains byte-exact for all 5 rows.
    hits = 0
    for row in _rows():
        r = reconstruct_tcr(row[0], row[1], row[2], "mouse")
        if r.get("full_chain_aa") == row[6]:
            hits += 1
    assert hits == 5, f"only {hits}/5 full beta chains matched"


def test_alpha_full_chain_constant_exact_diffs_confined_to_variable():
    # alpha full chains are length-exact; the vendored constant is reproduced
    # exactly, so any residual diffs are confined to the variable domain, where
    # they are germline allele polymorphism (not inferable from V+J+CDR3).
    from tcr_explorer.constant_regions import constant_aa
    tra_const = constant_aa("alpha", "mouse")
    for row in _rows():
        r = reconstruct_tcr(row[3], row[4], row[5], "mouse")
        fc = r.get("full_chain_aa")
        assert fc is not None and len(fc) == len(row[7])
        assert fc.endswith(tra_const) and row[7].endswith(tra_const)


def test_explicit_allele_is_honored_and_reported():
    r = reconstruct_tcr("TRAV7-4*02", "TRAJ16", "CAASLTSSGQKLVF", "mouse")
    assert r["v_allele_used"] and r["v_allele_used"].endswith("*02")
    d = reconstruct_tcr("TRAV7-4", "TRAJ16", "CAASLTSSGQKLVF", "mouse")
    assert d["v_allele_used"] and d["v_allele_used"].endswith("*01")


def test_default_allele_keeps_beta_full_chains_exact():
    # switching to allele-aware lookup must not change the default (*01) result
    hits = sum(1 for row in _rows()
               if reconstruct_tcr(row[0], row[1], row[2], "mouse").get("full_chain_aa") == row[6])
    assert hits == 5
