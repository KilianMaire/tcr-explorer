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


def test_full_beta_chains_reproduce_ground_truth():
    # variable domain (V+CDR3+J) plus the vendored constant reproduces the
    # complete membrane-bound beta chains byte-exact for all 5 rows.
    hits = 0
    for row in _rows():
        r = reconstruct_tcr(row[0], row[1], row[2], "mouse")
        if r.get("full_chain_aa") == row[6]:
            hits += 1
    assert hits == 5, f"only {hits}/5 full beta chains matched"


def test_alpha_full_chain_is_length_exact_allele_limited():
    # alpha full chains reproduce at exact length; residual diffs are germline
    # allele polymorphism (not inferable from V+J+CDR3), so at most a few residues.
    for row in _rows():
        r = reconstruct_tcr(row[3], row[4], row[5], "mouse")
        fc = r.get("full_chain_aa")
        assert fc is not None and len(fc) == len(row[7])
        diffs = sum(1 for a, b in zip(fc, row[7]) if a != b)
        assert diffs <= 2, f"{diffs} diffs, expected allele-limited"
