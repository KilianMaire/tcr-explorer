import csv
from pathlib import Path
from imgt_app.reconstructor import reconstruct_tcr

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
