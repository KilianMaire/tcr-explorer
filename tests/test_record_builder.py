from imgt_app import records
from imgt_app.dossier_models import TCRRecord


def _row(**kw):
    base = dict(
        source="vdjdb", source_record_id="vdjdb:1", pairing_key="vdjdb:u1",
        species="human", chain="beta", cdr3_aa="CASSLGTEAFF", cdr3_nt=None,
        full_aa=None, full_nt=None, v_gene="TRBV20-1", d_gene=None, j_gene="TRBJ2-7",
        cdr1_aa=None, cdr2_aa=None, epitope_aa="NLVPMVATV", antigen=None,
        antigen_organism=None, mhc_class="I", mhc_a="HLA-A*02:01", mhc_b=None,
        pdb_id=None, reference_pmid="12345", external_url="https://vdjdb.cdr3.net/search",
        score=1.0,
    )
    base.update(kw)
    return base


def test_deposited_nt_is_not_reconstructed():
    r = records.build_record(_row(source="mcpas", cdr3_nt="TGTGCC", full_aa=None))
    assert r.cdr3_nt == "TGTGCC"
    assert r.cdr3_nt_kind == "deposited"
    assert r.nt_is_synthetic is False


def test_missing_nt_is_reconstructed_and_flagged():
    r = records.build_record(_row())  # vdjdb: no deposited nt, has V+J+CDR3
    assert r.full_nt is not None and r.full_nt_kind == "reconstructed"
    assert r.nt_is_synthetic is True
    # composition is populated from germline pieces
    assert r.composition is not None
    assert r.composition.cdr3_aa == "CASSLGTEAFF"
    assert r.composition.v_region_nt and r.composition.j_region_nt


def test_no_genes_yields_no_reconstruction_no_fabrication():
    r = records.build_record(_row(v_gene=None, j_gene=None))
    assert r.full_nt is None and r.nt_is_synthetic is False
    assert r.composition is None or r.composition.v_germline_aa is None


def test_iedb_deposited_full_aa_kept_verbatim():
    r = records.build_record(_row(source="iedb", full_aa="MGAAAVILQCS...", cdr3_nt=None))
    assert r.full_aa == "MGAAAVILQCS..."
    assert r.full_aa_kind == "deposited"


def test_match_kind_and_similarity_threaded():
    r = records.build_record(_row(), match_kind="neighbour", similarity=0.83, concordance=3)
    assert r.match_kind == "neighbour" and r.similarity == 0.83 and r.concordance == 3
