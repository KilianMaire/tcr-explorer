"""Tests for external baseline prediction clients."""
import pytest
from servers.baselines.pmtnet_omni import PmtnetOmniClient, PmtnetOmniResult
from servers.baselines.mixtcrpred import MixTCRpredClient, MixTCRpredResult

def test_pmtnet_result():
    r = PmtnetOmniResult(score=0.888, rank=1.5, tcr_representation="paired")
    assert 0.0 <= r.score <= 1.0

def test_pmtnet_client_creates():
    client = PmtnetOmniClient()
    assert client.base_url is not None

def test_mixtcrpred_result():
    r = MixTCRpredResult(score=0.75, epitope="GILGFVFTL", rank=2.0)
    assert r.epitope == "GILGFVFTL"

def test_mixtcrpred_client_creates():
    client = MixTCRpredClient()
    assert client.base_url is not None
