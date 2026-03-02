"""Tests for pMHC binding predictor (Class I and II)."""
from __future__ import annotations
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "servers"))

import pytest
import math
from batman.pmhc import detect_mhc_class, rank_to_score, ic50_to_score


class TestDetectMhcClass:
    @pytest.mark.parametrize("allele,expected", [
        ("HLA-A*02:01", "I"),
        ("HLA-B*07:02", "I"),
        ("HLA-C*03:04", "I"),
        ("HLA-E*01:01", "I"),
        ("HLA-DRA1*01:01", "II"),
        ("HLA-DRB1*04:01", "II"),
        ("HLA-DQA1*01:02", "II"),
        ("HLA-DQB1*06:02", "II"),
        ("HLA-DPA1*01:03", "II"),
        ("HLA-DPB1*04:01", "II"),
        ("DRB1*04:01", "II"),
        ("DQB1*06:02", "II"),
        ("", "I"),
        ("unknown", "I"),
    ])
    def test_detect_class(self, allele, expected):
        assert detect_mhc_class(allele) == expected


class TestRankToScore:
    def test_strong_binder_rank_0_5_near_1(self):
        score = rank_to_score(0.5)
        assert score > 0.95

    def test_weak_binder_rank_2_above_half(self):
        score = rank_to_score(2.0)
        assert score > 0.80

    def test_rank_10_near_half(self):
        score = rank_to_score(10.0)
        assert 0.45 <= score <= 0.55

    def test_rank_50_near_zero(self):
        score = rank_to_score(50.0)
        assert score < 0.10

    def test_output_in_0_1(self):
        for r in [0, 0.5, 1, 2, 5, 10, 25, 50, 100]:
            assert 0.0 <= rank_to_score(r) <= 1.0


class TestIc50ToScore:
    def test_very_low_ic50_near_1(self):
        assert ic50_to_score(1.0) > 0.95

    def test_strong_binder_50nm(self):
        score = ic50_to_score(50.0)
        assert score > 0.75

    def test_threshold_50000nm_is_zero(self):
        assert ic50_to_score(50_000.0) == pytest.approx(0.0)

    def test_above_threshold_clipped_to_zero(self):
        assert ic50_to_score(100_000.0) == pytest.approx(0.0)

    def test_output_always_in_0_1(self):
        for ic50 in [1, 50, 500, 5000, 50000, 100000]:
            assert 0.0 <= ic50_to_score(ic50) <= 1.0


from unittest.mock import patch, MagicMock
import pandas as pd
from batman.pmhc import MHCflurryPredictor


class TestMHCflurryPredictor:
    @patch("batman.pmhc.Class1PresentationPredictor")
    def test_predict_class_i_returns_presentation_score(self, mock_cls):
        mock_result = pd.DataFrame({
            "peptide": ["GILGFVFTL"],
            "allele": ["HLA-A*02:01"],
            "presentation_score": [0.87],
            "affinity": [45.0],
            "affinity_percentile": [0.12],
        })
        mock_predictor = MagicMock()
        mock_predictor.predict.return_value = mock_result
        mock_cls.load.return_value = mock_predictor

        p = MHCflurryPredictor()
        score = p.predict_class_i("GILGFVFTL", "HLA-A*02:01")

        assert score == pytest.approx(0.87)
        mock_predictor.predict.assert_called_once_with(
            peptides=["GILGFVFTL"],
            alleles=["HLA-A*02:01"],
        )

    @patch("batman.pmhc.Class1PresentationPredictor")
    def test_predict_returns_none_on_mhcflurry_error(self, mock_cls):
        mock_cls.load.side_effect = RuntimeError("model not found")
        p = MHCflurryPredictor()
        score = p.predict_class_i("GILGFVFTL", "HLA-A*02:01")
        assert score is None

    @patch("batman.pmhc.Class1PresentationPredictor")
    def test_predictor_loaded_once(self, mock_cls):
        mock_cls.load.return_value = MagicMock(
            predict=MagicMock(return_value=pd.DataFrame({
                "peptide": ["GILGFVFTL"], "allele": ["HLA-A*02:01"],
                "presentation_score": [0.7], "affinity": [100.0],
                "affinity_percentile": [0.5],
            }))
        )
        p = MHCflurryPredictor()
        p.predict_class_i("GILGFVFTL", "HLA-A*02:01")
        p.predict_class_i("NLVPMVATV", "HLA-A*02:01")
        assert mock_cls.load.call_count == 1


import asyncio
import httpx
from batman.pmhc import IEDBpMHCPredictor


class TestIEDBpMHCPredictor:
    def _run(self, coro):
        return asyncio.get_event_loop().run_until_complete(coro)

    def _tsv_class_i(self) -> str:
        return (
            "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tic50\tpercentile_rank\n"
            "HLA-A*02:01\t1\t1\t9\t9\tGILGFVFTL\tnetmhcpan_el\t31.3\t0.12\n"
        )

    def _tsv_class_ii(self) -> str:
        return (
            "allele\tseq_num\tstart\tend\tlength\tpeptide\tmethod\tic50\tpercentile_rank\n"
            "HLA-DRB1*04:01\t1\t1\t15\t15\tGILGFVFTLMTCGRT\tnetmhciipan\t245.0\t1.8\n"
        )

    @patch("batman.pmhc.httpx.AsyncClient")
    def test_predict_class_i_returns_score(self, mock_client_cls):
        mock_resp = MagicMock()
        mock_resp.text = self._tsv_class_i()
        mock_resp.raise_for_status = MagicMock()

        mock_client = MagicMock()
        mock_client.__aenter__ = MagicMock(return_value=mock_client)
        mock_client.__aexit__ = MagicMock(return_value=False)
        mock_client.post = MagicMock(return_value=mock_resp)
        mock_client_cls.return_value = mock_client

        predictor = IEDBpMHCPredictor()
        score = self._run(predictor.predict("GILGFVFTL", "HLA-A*02:01"))
        assert score is not None
        assert 0.0 < score <= 1.0
        assert score > 0.95  # rank=0.12 is strong binder

    @patch("batman.pmhc.httpx.AsyncClient")
    def test_predict_class_ii_uses_mhcii_endpoint(self, mock_client_cls):
        mock_resp = MagicMock()
        mock_resp.text = self._tsv_class_ii()
        mock_resp.raise_for_status = MagicMock()

        mock_client = MagicMock()
        mock_client.__aenter__ = MagicMock(return_value=mock_client)
        mock_client.__aexit__ = MagicMock(return_value=False)
        mock_client.post = MagicMock(return_value=mock_resp)
        mock_client_cls.return_value = mock_client

        predictor = IEDBpMHCPredictor()
        self._run(predictor.predict("GILGFVFTLMTCGRT", "HLA-DRB1*04:01"))

        call_args = mock_client.post.call_args
        assert "mhcii" in call_args[0][0]

    @patch("batman.pmhc.httpx.AsyncClient")
    def test_predict_returns_none_on_http_error(self, mock_client_cls):
        mock_client = MagicMock()
        mock_client.__aenter__ = MagicMock(return_value=mock_client)
        mock_client.__aexit__ = MagicMock(return_value=False)
        mock_client.post = MagicMock(side_effect=httpx.HTTPError("timeout"))
        mock_client_cls.return_value = mock_client

        predictor = IEDBpMHCPredictor()
        score = self._run(predictor.predict("GILGFVFTL", "HLA-A*02:01"))
        assert score is None

    @patch("batman.pmhc.httpx.AsyncClient")
    def test_class_ii_peptide_length_uses_actual_length(self, mock_client_cls):
        mock_resp = MagicMock()
        mock_resp.text = self._tsv_class_ii()
        mock_resp.raise_for_status = MagicMock()

        mock_client = MagicMock()
        mock_client.__aenter__ = MagicMock(return_value=mock_client)
        mock_client.__aexit__ = MagicMock(return_value=False)
        mock_client.post = MagicMock(return_value=mock_resp)
        mock_client_cls.return_value = mock_client

        predictor = IEDBpMHCPredictor()
        self._run(predictor.predict("GILGFVFTLMTCGRT", "HLA-DRB1*04:01"))

        posted_data = mock_client.post.call_args[1]["data"]
        assert posted_data["length"] == "15"


from batman.pmhc import predict_pmhc


class TestPredictPmhcUnified:
    def _run(self, coro):
        return asyncio.get_event_loop().run_until_complete(coro)

    @patch("batman.pmhc._mhcflurry")
    def test_class_i_uses_mhcflurry_primary(self, mock_mhcflurry):
        mock_mhcflurry.predict_class_i.return_value = 0.87
        score = self._run(predict_pmhc("GILGFVFTL", "HLA-A*02:01"))
        assert score == pytest.approx(0.87)
        mock_mhcflurry.predict_class_i.assert_called_once()

    @patch("batman.pmhc._iedb")
    @patch("batman.pmhc._mhcflurry")
    def test_class_i_falls_back_to_iedb_when_mhcflurry_fails(self, mock_mhcflurry, mock_iedb):
        mock_mhcflurry.predict_class_i.return_value = None
        async def fake_iedb(peptide, allele):
            return 0.72
        mock_iedb.predict = fake_iedb

        score = self._run(predict_pmhc("GILGFVFTL", "HLA-A*02:01"))
        assert score == pytest.approx(0.72)

    @patch("batman.pmhc._iedb")
    @patch("batman.pmhc._mhcflurry")
    def test_class_ii_skips_mhcflurry_uses_iedb(self, mock_mhcflurry, mock_iedb):
        async def fake_iedb(peptide, allele):
            return 0.65
        mock_iedb.predict = fake_iedb

        score = self._run(predict_pmhc("GILGFVFTLMTCGRT", "HLA-DRB1*04:01"))
        assert score == pytest.approx(0.65)
        mock_mhcflurry.predict_class_i.assert_not_called()

    @patch("batman.pmhc._iedb")
    @patch("batman.pmhc._mhcflurry")
    def test_both_fail_returns_none(self, mock_mhcflurry, mock_iedb):
        mock_mhcflurry.predict_class_i.return_value = None
        async def fake_iedb(peptide, allele):
            return None
        mock_iedb.predict = fake_iedb

        score = self._run(predict_pmhc("GILGFVFTL", "HLA-A*02:01"))
        assert score is None


from batman.pmhc import NetMHCIIpanPredictor


class TestNetMHCIIpanPredictor:
    @patch("batman.pmhc._NETMHCIIPAN_AVAILABLE", True)
    @patch("batman.pmhc._netmhciipan_predict")
    def test_predict_class_ii_returns_score(self, mock_predict):
        mock_predict.return_value = 0.72
        p = NetMHCIIpanPredictor()
        score = p.predict("AGFKGEQGPKGEPG", "HLA-DRB1*04:01")
        assert score == pytest.approx(0.72)

    @patch("batman.pmhc._NETMHCIIPAN_AVAILABLE", False)
    def test_predict_returns_none_when_not_installed(self):
        p = NetMHCIIpanPredictor()
        score = p.predict("AGFKGEQGPKGEPG", "HLA-DRB1*04:01")
        assert score is None

    @patch("batman.pmhc._NETMHCIIPAN_AVAILABLE", True)
    @patch("batman.pmhc._netmhciipan_predict")
    def test_predict_mouse_class_ii(self, mock_predict):
        mock_predict.return_value = 0.65
        p = NetMHCIIpanPredictor()
        score = p.predict("AGFKGEQGPKGEPG", "H-2-IAb")
        assert score == pytest.approx(0.65)
        mock_predict.assert_called_once()


class TestPredictPmhcClassIIRouting:
    def _run(self, coro):
        return asyncio.get_event_loop().run_until_complete(coro)

    @patch("batman.pmhc._netmhciipan")
    @patch("batman.pmhc._iedb")
    def test_class_ii_uses_netmhciipan_primary(self, mock_iedb, mock_netmhciipan):
        mock_netmhciipan.predict.return_value = 0.72
        score = self._run(predict_pmhc("AGFKGEQGPKGEPG", "HLA-DRB1*04:01"))
        assert score == pytest.approx(0.72)
        mock_netmhciipan.predict.assert_called_once()

    @patch("batman.pmhc._netmhciipan")
    @patch("batman.pmhc._iedb")
    def test_class_ii_falls_back_to_iedb(self, mock_iedb, mock_netmhciipan):
        mock_netmhciipan.predict.return_value = None
        async def fake_iedb(peptide, allele):
            return 0.58
        mock_iedb.predict = fake_iedb
        score = self._run(predict_pmhc("AGFKGEQGPKGEPG", "HLA-DRB1*04:01"))
        assert score == pytest.approx(0.58)
