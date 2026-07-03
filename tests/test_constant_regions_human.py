"""Human membrane-bound constant regions and full human-chain assembly.

The human constants are byte-exact from UniProt (TRAC P01848, TRBC1 P01850)
with the beta CH1 junction residue E restored (UniProt P01850 omits it, but
the oracle-validated mouse entry P01852 carries it). They are not oracle
validated against a human ground-truth set, so their provenance string says so.
"""
from imgt_app.constant_regions import constant_aa, constant_source
from imgt_app.reconstructor import reconstruct_tcr

# UniProt P01848 (TRAC_HUMAN), verbatim.
HUMAN_TRA = (
    "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFA"
    "CANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS"
)
# UniProt P01850 (TRBC1_HUMAN) with the CH1 junction E restored at position 1.
HUMAN_TRB = (
    "EDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALND"
    "SRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGV"
    "LSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
)


def test_human_constants_are_vendored_verbatim():
    assert constant_aa("alpha", "human") == HUMAN_TRA
    assert constant_aa("beta", "human") == HUMAN_TRB
    # beta begins with the E junction residue, mirroring the mouse convention
    assert constant_aa("beta", "human").startswith("EDLNKVFPP")


def test_human_provenance_flags_not_oracle_validated():
    src = constant_source("beta", "human")
    assert src and "human" in src.lower()
    assert "not oracle" in src.lower()
    assert "P01850" in src


def test_mouse_constants_unchanged():
    # regression: mouse stays oracle-validated and byte-identical
    assert constant_aa("beta", "mouse").startswith("EDLRNVTPP")
    assert constant_aa("alpha", "mouse").startswith("IQNPEPAVY")
    assert "oracle" in (constant_source("beta", "mouse") or "").lower()


def test_human_full_chain_appends_human_constant():
    r = reconstruct_tcr("TRBV4-1", "TRBJ1-1", "CASSLGTEAFF", "human")
    assert r["full_aa"] and r["full_aa"].endswith("GQGTRLTVV")
    assert r["full_chain_aa"] == r["full_aa"] + HUMAN_TRB
    assert "not oracle" in r["constant_source"].lower()


def test_unknown_species_constant_is_none():
    assert constant_aa("beta", "rat") is None
    assert constant_source("beta", "rat") is None
