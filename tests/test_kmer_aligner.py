import pytest
from imgt_app.cdr_enricher import _stitchr_data_dir
from imgt_app.kmer_aligner import annotate_sequence

pytestmark = pytest.mark.skipif(_stitchr_data_dir() is None, reason="stitchr germline not installed")

def test_v_call_from_v_region_nt():
    # Grab a real TRBV V-REGION nt from the germline map and confirm the aligner recovers its gene.
    from imgt_app.cdr_enricher import _cached_v_map
    vmap = _cached_v_map("TRB", "HUMAN")
    gene, nt = next(iter(vmap.items()))
    ann = annotate_sequence(nt, "human", is_protein=False)
    assert ann.v_call is not None
    assert ann.chain == "beta"
    # best V hit should be the source gene (allow allele/base equivalence)
    assert ann.v_call.split("*")[0] == gene.split("*")[0]


def test_v_call_from_translated_v_region_aa():
    # Protein mode must score the aa query against TRANSLATED germline, not nt germline.
    # Translate a real TRBV V-REGION nt and confirm the aligner recovers the source V gene.
    from imgt_app.cdr_enricher import _cached_v_map, _translate
    vmap = _cached_v_map("TRB", "HUMAN")
    gene, nt = next(iter(vmap.items()))
    aa = _translate(nt).rstrip("*")
    ann = annotate_sequence(aa, "human", is_protein=True)
    assert ann.v_call is not None
    assert ann.chain == "beta"
    assert ann.v_call.split("*")[0] == gene.split("*")[0]
