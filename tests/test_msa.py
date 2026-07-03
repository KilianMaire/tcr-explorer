import pytest

from imgt_app import msa
from imgt_app.msa import align, center_star_align, to_fasta
from imgt_app.dossier_models import AlignRequest
from imgt_app.cdr_enricher import _stitchr_data_dir


def test_identical_sequences_no_gaps_full_identity():
    r = align(AlignRequest(sequences=[{"name": "a", "seq": "CASSLGTEAFF"}, {"name": "b", "seq": "CASSLGTEAFF"}], seq_type="aa"))
    assert r.alignment_length == len("CASSLGTEAFF")
    assert all("-" not in rec.aligned for rec in r.records)
    assert r.mean_pct_identity == 100.0


def test_equal_length_after_alignment(monkeypatch):
    monkeypatch.setattr(msa, "clustalo_available", lambda: False)
    r = align(AlignRequest(sequences=[
        {"name": "a", "seq": "CASSLGTEAFF"}, {"name": "b", "seq": "CASSLGTEAF"}, {"name": "c", "seq": "CASSLGTEAFFF"}], seq_type="aa"))
    lens = {len(rec.aligned) for rec in r.records}
    assert len(lens) == 1  # all padded to one alignment length
    assert r.engine == "center_star"


def test_too_few_sequences_warns():
    r = align(AlignRequest(sequences=[{"name": "a", "seq": "CASS"}], seq_type="aa"))
    assert r.records == [] and any(w.code == "too_few_sequences" for w in r.warnings)


def test_to_fasta_roundtrips_names():
    r = align(AlignRequest(sequences=[{"name": "x", "seq": "CASS"}, {"name": "y", "seq": "CASF"}], seq_type="aa"))
    fa = to_fasta(r)
    assert ">x" in fa and ">y" in fa


def test_center_star_merge_preserves_column_correspondence(monkeypatch):
    """A naive left-justified pad puts a different number of trailing gaps on each
    row depending on its own pairwise alignment to the center, so a conserved
    residue that should line up in the same column across all rows instead lands
    in different columns. This test uses a center sequence with an internal
    deletion relative to one of the other two sequences, so the correct merge
    must insert a gap into the OTHER already-aligned rows (not just pad at the
    end) to keep the trailing conserved motif in a single shared column.

    center: CASS-LGTEAFF   (this is the longest/most-central by pairwise score)
    seq b : CASSLLGTEAFF   (has an extra L that the center lacks -> center gets
                             a gap inserted relative to b, at an INTERNAL position)
    seq a : CASSLGTEAFF    (identical to center minus the gap)

    Correct center-star merge: once the a-vs-center and b-vs-center pairwise
    alignments are merged, ALL rows must be re-threaded onto the union of gap
    columns, so the trailing "GTEAFF" motif shared by all three sequences ends
    up in the same set of columns for every row. A naive right-pad would instead
    leave the motif in different columns for the row that didn't need a gap.
    """
    monkeypatch.setattr(msa, "clustalo_available", lambda: False)
    seqs = [
        ("center", "CASSLGTEAFF"),
        ("a", "CASSLGTEAFF"),
        ("b", "CASSLLGTEAFF"),
    ]
    aligned = center_star_align(seqs, "aa")
    by_name = dict(aligned)
    L = len(by_name["b"])
    assert all(len(s) == L for s in by_name.values())

    # The trailing conserved motif "GTEAFF" must occupy identical columns in
    # every row (allowing for the row's own internal gaps elsewhere).
    def motif_start(gapped: str, motif: str = "GTEAFF") -> int:
        ungapped_positions = [i for i, c in enumerate(gapped) if c != "-"]
        residues = "".join(gapped[i] for i in ungapped_positions)
        idx = residues.index(motif)
        return ungapped_positions[idx]

    starts = {name: motif_start(g) for name, g in aligned}
    assert len(set(starts.values())) == 1, f"motif columns diverge: {starts} / alignment={aligned}"


@pytest.mark.skipif(_stitchr_data_dir() is None, reason="stitchr germline data not present")
def test_translated_germline_set_aligns(monkeypatch):
    """Translating a germline J-REGION from frame 0 yields proteins containing a
    stop '*' and other letters near the matrix-alphabet edges. The aa aligner
    must tolerate the whole germline set (frame-0 translation is acceptable for
    aligning the set to itself). Regression: previously returned engine='none'
    with an alignment_failed warning because the growing profile's gapped first
    row was fed back into PairwiseAligner (which rejects '-')."""
    monkeypatch.setattr(msa, "clustalo_available", lambda: False)
    r = align(AlignRequest(species="mouse", chain="TRB", segment="J", translate=True, seq_type="aa"))
    assert r.engine == "center_star"
    assert r.n_sequences == 14
    assert r.alignment_length > 0
    assert not any(w.code == "alignment_failed" for w in r.warnings)
