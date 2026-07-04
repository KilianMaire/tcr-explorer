import pytest
from tcr_explorer.dossier_models import AlignRequest
from tcr_explorer.msa import align
from tcr_explorer.cdr_enricher import _stitchr_data_dir


def test_provided_nt_codon_aware_registers_aa_over_codons():
    # two in-frame nt sequences (same length, one substitution) -> codon-aware view
    req = AlignRequest(
        sequences=[{"name": "a", "seq": "TGTGCCAGC"}, {"name": "b", "seq": "TGTGGCAGC"}],
        seq_type="nt", translate=True)
    r = align(req)
    assert r.view == "aa_nt"
    for rec in r.records:
        assert rec.aligned_aa is not None and rec.aligned_nt is not None
        # nt row is exactly three times the aa row (a residue over its codon)
        assert len(rec.aligned_nt) == 3 * len(rec.aligned_aa)
    # conservation list has one entry per aa column, in 0..1
    assert len(r.conservation) == len(r.records[0].aligned_aa)
    assert all(0.0 <= c <= 1.0 for c in r.conservation)


def test_gap_codon_is_triple_dash():
    req = AlignRequest(
        sequences=[{"name": "a", "seq": "TGTGCCAGCTTT"}, {"name": "b", "seq": "TGTGCCAGC"}],
        seq_type="nt", translate=True)
    r = align(req)
    # the shorter sequence must carry a '---' gap codon where the aa is '-'
    short = next(rec for rec in r.records if rec.name == "b")
    for i, ch in enumerate(short.aligned_aa):
        codon = short.aligned_nt[3 * i:3 * i + 3]
        if ch == "-":
            assert codon == "---"
        else:
            assert "-" not in codon


def test_backward_compat_nt_view_unchanged():
    req = AlignRequest(
        sequences=[{"name": "a", "seq": "TGTGCCAGC"}, {"name": "b", "seq": "TGTGGCAGC"}],
        seq_type="nt", translate=False)
    r = align(req)
    assert r.view == "nt"
    assert r.records[0].aligned  # primary row still present
    assert r.records[0].aligned_aa is None


def test_conservation_identical_column_is_one():
    req = AlignRequest(
        sequences=[{"name": "a", "seq": "TGTGCCAGC"}, {"name": "b", "seq": "TGTGCCAGC"}],
        seq_type="nt", translate=True)
    r = align(req)
    # identical sequences -> every column fully conserved
    assert all(c == 1.0 for c in r.conservation)


@pytest.mark.skipif(_stitchr_data_dir() is None, reason="stitchr germline not installed")
def test_mouse_trbj_codon_aware_nt_is_triple_aa():
    r = align(AlignRequest(species="mouse", chain="TRB", segment="J",
                           seq_type="nt", translate=True))
    assert r.view == "aa_nt" and r.n_sequences == 14
    for rec in r.records:
        assert len(rec.aligned_nt) == 3 * len(rec.aligned_aa)
    # the conserved FGxG motif should appear in the aa consensus/rows
    joined = "".join(r.consensus)
    assert "FG" in joined
