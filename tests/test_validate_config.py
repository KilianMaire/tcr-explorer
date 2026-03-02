"""Tests for scripts/validate_config.py — pre-flight configuration validation."""
from __future__ import annotations

import json
import os
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Ensure project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.validate_config import (
    EMBEDDING_DIRS,
    ENV_VARS,
    MODEL_CHECKPOINT_PATH,
    PIPELINE_SCRIPTS,
    REQUIRED_DIRS,
    SERVICES,
    TRAINING_DATA_PATTERNS,
    CheckResult,
    SectionResult,
    ValidationReport,
    build_parser,
    check_compute,
    check_environment,
    check_file_system,
    check_model_readiness,
    check_services,
    format_json,
    format_text,
    main,
    run_all_checks,
)


# ---- Helpers ---------------------------------------------------------------

def _make_root(tmp_path: Path, *, scripts: bool = True, dirs: bool = True,
               dvc: bool = True) -> Path:
    """Build a fake project root with controllable contents."""
    root = tmp_path / "project"
    root.mkdir()
    if scripts:
        for s in PIPELINE_SCRIPTS:
            p = root / s
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("# stub\n")
    if dirs:
        for d in REQUIRED_DIRS:
            (root / d).mkdir(parents=True, exist_ok=True)
    if dvc:
        (root / ".dvc").mkdir()
    return root


# ---- CheckResult / SectionResult / ValidationReport -----------------------

class TestDataClasses:
    def test_check_result_symbols(self):
        assert CheckResult("ok", passed=True).symbol == "\u2713"
        assert CheckResult("warn", passed=False, warning=True).symbol == "~"
        assert CheckResult("fail", passed=False).symbol == "\u2717"

    def test_check_result_to_dict(self):
        cr = CheckResult("x", passed=True, message="hi", detail="d")
        d = cr.to_dict()
        assert d["name"] == "x"
        assert d["passed"] is True
        assert d["message"] == "hi"
        assert d["detail"] == "d"
        assert d["warning"] is False

    def test_section_result_counters(self):
        section = SectionResult(name="s", checks=[
            CheckResult("a", passed=True),
            CheckResult("b", passed=False, warning=True),
            CheckResult("c", passed=False),
        ])
        assert section.passed == 1
        assert section.warnings == 1
        assert section.failures == 1

    def test_section_overall_symbol_pass(self):
        s = SectionResult(name="s", checks=[
            CheckResult("a", passed=True),
            CheckResult("b", passed=True),
        ])
        assert s.overall_symbol == "\u2713"

    def test_section_overall_symbol_warn(self):
        s = SectionResult(name="s", checks=[
            CheckResult("a", passed=True),
            CheckResult("b", passed=False, warning=True),
        ])
        assert s.overall_symbol == "~"

    def test_section_overall_symbol_fail(self):
        s = SectionResult(name="s", checks=[
            CheckResult("a", passed=True),
            CheckResult("b", passed=False),
        ])
        assert s.overall_symbol == "\u2717"

    def test_section_to_dict(self):
        s = SectionResult(name="test", checks=[
            CheckResult("a", passed=True, message="ok"),
        ])
        d = s.to_dict()
        assert d["name"] == "test"
        assert d["passed"] == 1
        assert d["warnings"] == 0
        assert d["failures"] == 0
        assert len(d["checks"]) == 1

    def test_report_totals(self):
        r = ValidationReport(sections=[
            SectionResult(name="s1", checks=[
                CheckResult("a", passed=True),
                CheckResult("b", passed=False, warning=True),
            ]),
            SectionResult(name="s2", checks=[
                CheckResult("c", passed=False),
            ]),
        ])
        assert r.total_checks == 3
        assert r.total_passed == 1
        assert r.total_warnings == 1
        assert r.total_failures == 1

    def test_report_to_dict(self):
        r = ValidationReport(sections=[
            SectionResult(name="s", checks=[
                CheckResult("a", passed=True, message="ok"),
            ]),
        ])
        d = r.to_dict()
        assert d["total_checks"] == 1
        assert d["passed"] == 1
        assert len(d["sections"]) == 1


# ---- File system checks ----------------------------------------------------

class TestFileSystem:
    def test_all_present(self, tmp_path):
        root = _make_root(tmp_path)
        result = check_file_system(root)
        assert result.name == "File System"
        scripts_check = next(c for c in result.checks if c.name == "pipeline_scripts")
        assert scripts_check.passed is True
        dirs_check = next(c for c in result.checks if c.name == "directories")
        assert dirs_check.passed is True

    def test_missing_scripts(self, tmp_path):
        root = _make_root(tmp_path, scripts=False)
        result = check_file_system(root)
        scripts_check = next(c for c in result.checks if c.name == "pipeline_scripts")
        assert scripts_check.passed is False
        assert "missing" in scripts_check.message

    def test_missing_directories(self, tmp_path):
        root = _make_root(tmp_path, dirs=False)
        # Create script parents
        for s in PIPELINE_SCRIPTS:
            p = root / s
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("")
        result = check_file_system(root)
        dirs_check = next(c for c in result.checks if c.name == "directories")
        assert dirs_check.passed is False

    def test_fix_creates_directories(self, tmp_path):
        root = _make_root(tmp_path, dirs=False)
        result = check_file_system(root, fix=True)
        dirs_check = next(c for c in result.checks if c.name == "directories")
        assert dirs_check.passed is True
        assert "Created" in dirs_check.message
        # Verify directories actually exist now
        for d in REQUIRED_DIRS:
            assert (root / d).is_dir()

    def test_dvc_missing(self, tmp_path):
        root = _make_root(tmp_path, dvc=False)
        result = check_file_system(root)
        dvc_check = next(c for c in result.checks if c.name == "dvc_initialized")
        assert dvc_check.passed is False
        assert dvc_check.warning is True

    def test_dvc_present(self, tmp_path):
        root = _make_root(tmp_path)
        result = check_file_system(root)
        dvc_check = next(c for c in result.checks if c.name == "dvc_initialized")
        assert dvc_check.passed is True

    def test_import_failure_is_warning(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module",
                    side_effect=ImportError("no")):
            result = check_file_system(root)
        imp_check = next(c for c in result.checks if c.name == "python_imports")
        assert imp_check.passed is False
        assert imp_check.warning is True

    def test_import_success(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            result = check_file_system(root)
        imp_check = next(c for c in result.checks if c.name == "python_imports")
        assert imp_check.passed is True


# ---- Environment checks ----------------------------------------------------

class TestEnvironment:
    def test_all_set(self, tmp_path):
        root = _make_root(tmp_path)
        env_patch = {name: "http://example.com" for name, _, _ in ENV_VARS}
        with patch.dict(os.environ, env_patch, clear=False):
            result = check_environment(root)
        assert all(c.passed for c in result.checks)

    def test_required_missing(self, tmp_path):
        root = _make_root(tmp_path)
        # Clear all the env vars we care about
        clean_env = {k: v for k, v in os.environ.items()
                     if k not in {name for name, _, _ in ENV_VARS}}
        with patch.dict(os.environ, clean_env, clear=True):
            result = check_environment(root)
        required_names = {name for name, req, _ in ENV_VARS if req}
        for c in result.checks:
            if c.name in required_names:
                assert c.passed is False
                assert c.warning is False, f"{c.name} should be critical, not warning"

    def test_optional_missing_is_warning(self, tmp_path):
        root = _make_root(tmp_path)
        clean_env = {k: v for k, v in os.environ.items()
                     if k not in {name for name, _, _ in ENV_VARS}}
        with patch.dict(os.environ, clean_env, clear=True):
            result = check_environment(root)
        optional_names = {name for name, req, _ in ENV_VARS if not req}
        for c in result.checks:
            if c.name in optional_names:
                assert c.warning is True, f"{c.name} should be a warning"

    def test_dotenv_loading(self, tmp_path):
        root = _make_root(tmp_path)
        dotenv = root / ".env"
        dotenv.write_text('NCBI_API_KEY="test_key_123"\n')
        # No env vars set
        clean_env = {k: v for k, v in os.environ.items()
                     if k not in {name for name, _, _ in ENV_VARS}}
        with patch.dict(os.environ, clean_env, clear=True):
            result = check_environment(root)
        ncbi_check = next(c for c in result.checks if c.name == "NCBI_API_KEY")
        assert ncbi_check.passed is True
        assert "test_key_123" in ncbi_check.message

    def test_env_overrides_dotenv(self, tmp_path):
        root = _make_root(tmp_path)
        dotenv = root / ".env"
        dotenv.write_text("NCBI_API_KEY=from_file\n")
        clean_env = {k: v for k, v in os.environ.items()
                     if k not in {name for name, _, _ in ENV_VARS}}
        clean_env["NCBI_API_KEY"] = "from_env"
        with patch.dict(os.environ, clean_env, clear=True):
            result = check_environment(root)
        ncbi_check = next(c for c in result.checks if c.name == "NCBI_API_KEY")
        assert "from_env" in ncbi_check.message


# ---- Service health checks -------------------------------------------------

class TestServices:
    def test_service_up(self):
        mock_resp = MagicMock()
        mock_resp.status_code = 200

        fake_httpx = types.ModuleType("httpx")
        fake_httpx.get = MagicMock(return_value=mock_resp)

        with patch.dict(sys.modules, {"httpx": fake_httpx}):
            with patch.dict(os.environ, {}, clear=False):
                result = check_services()

        # All services should show as up
        for check in result.checks:
            assert check.passed is True
            assert "up" in check.message

    def test_service_down(self):
        fake_httpx = types.ModuleType("httpx")
        fake_httpx.get = MagicMock(side_effect=ConnectionError("refused"))

        with patch.dict(sys.modules, {"httpx": fake_httpx}):
            result = check_services()

        for check in result.checks:
            assert check.passed is False

    def test_service_unhealthy_status(self):
        mock_resp = MagicMock()
        mock_resp.status_code = 503

        fake_httpx = types.ModuleType("httpx")
        fake_httpx.get = MagicMock(return_value=mock_resp)

        with patch.dict(sys.modules, {"httpx": fake_httpx}):
            result = check_services()

        for check in result.checks:
            assert check.passed is False
            assert "unhealthy" in check.message

    def test_httpx_not_installed(self):
        # Remove httpx from modules temporarily
        with patch.dict(sys.modules, {"httpx": None}):
            # Force reimport check_services to see httpx as unavailable
            # Instead we patch the import inside the function
            def _raise_import(*a, **kw):
                raise ImportError("no httpx")

            with patch("scripts.validate_config.httpx",
                       new_callable=lambda: property(lambda s: (_ for _ in ()).throw(ImportError)),
                       create=True):
                # The function tries `import httpx` so mock that
                import builtins
                original_import = builtins.__import__

                def mock_import(name, *args, **kwargs):
                    if name == "httpx":
                        raise ImportError("no httpx")
                    return original_import(name, *args, **kwargs)

                with patch("builtins.__import__", side_effect=mock_import):
                    result = check_services()

        for check in result.checks:
            assert check.passed is False
            assert check.warning is True
            assert "httpx" in check.message

    def test_configured_service_url_from_env(self):
        mock_resp = MagicMock()
        mock_resp.status_code = 200

        fake_httpx = types.ModuleType("httpx")
        fake_httpx.get = MagicMock(return_value=mock_resp)

        with patch.dict(sys.modules, {"httpx": fake_httpx}):
            with patch.dict(os.environ,
                            {"TIER1_SERVER_URL": "http://custom:9999"},
                            clear=False):
                result = check_services()

        tier1 = next(c for c in result.checks if c.name == "Tier 1 ML")
        assert "http://custom:9999" in tier1.detail


# ---- Model readiness -------------------------------------------------------

class TestModelReadiness:
    def test_checkpoint_missing(self, tmp_path):
        root = _make_root(tmp_path)
        result = check_model_readiness(root)
        ckpt_check = next(c for c in result.checks if c.name == "model_checkpoint")
        assert ckpt_check.passed is False
        assert "not yet trained" in ckpt_check.message

    def test_checkpoint_present(self, tmp_path):
        root = _make_root(tmp_path)
        ckpt = root / MODEL_CHECKPOINT_PATH
        ckpt.parent.mkdir(parents=True, exist_ok=True)
        ckpt.write_bytes(b"\x00")
        result = check_model_readiness(root)
        ckpt_check = next(c for c in result.checks if c.name == "model_checkpoint")
        assert ckpt_check.passed is True

    def test_training_data_present(self, tmp_path):
        root = _make_root(tmp_path)
        # Create first training data pattern
        td = root / TRAINING_DATA_PATTERNS[0]
        td.parent.mkdir(parents=True, exist_ok=True)
        td.write_bytes(b"\x00")
        result = check_model_readiness(root)
        td_check = next(c for c in result.checks if c.name == "training_data")
        assert td_check.passed is True

    def test_training_data_missing(self, tmp_path):
        root = _make_root(tmp_path)
        result = check_model_readiness(root)
        td_check = next(c for c in result.checks if c.name == "training_data")
        assert td_check.passed is False
        assert td_check.warning is True

    def test_embeddings_present(self, tmp_path):
        root = _make_root(tmp_path)
        emb_dir = root / EMBEDDING_DIRS[0]
        emb_dir.mkdir(parents=True, exist_ok=True)
        (emb_dir / "sample.pt").write_bytes(b"\x00")
        result = check_model_readiness(root)
        emb_check = next(c for c in result.checks if c.name == "embeddings")
        assert emb_check.passed is True

    def test_embeddings_missing(self, tmp_path):
        root = _make_root(tmp_path)
        result = check_model_readiness(root)
        emb_check = next(c for c in result.checks if c.name == "embeddings")
        assert emb_check.passed is False
        assert emb_check.warning is True

    def test_embeddings_dir_exists_but_empty(self, tmp_path):
        root = _make_root(tmp_path)
        emb_dir = root / EMBEDDING_DIRS[0]
        emb_dir.mkdir(parents=True, exist_ok=True)
        # Directory exists but is empty
        result = check_model_readiness(root)
        emb_check = next(c for c in result.checks if c.name == "embeddings")
        assert emb_check.passed is False


# ---- Compute checks --------------------------------------------------------

class TestCompute:
    def test_cuda_device(self):
        mock_torch = MagicMock()
        mock_torch.__version__ = "2.2.0"
        mock_torch.cuda.is_available.return_value = True
        mock_torch.cuda.get_device_name.return_value = "NVIDIA A100"
        props = MagicMock()
        props.total_mem = 40 * 1024 * 1024 * 1024  # 40GB
        mock_torch.cuda.get_device_properties.return_value = props
        mock_torch.backends.mps.is_available.return_value = False

        with patch.dict(sys.modules, {"torch": mock_torch}):
            result = check_compute()

        pt_check = next(c for c in result.checks if c.name == "pytorch")
        assert pt_check.passed is True
        assert "2.2.0" in pt_check.message

        dev_check = next(c for c in result.checks if c.name == "device")
        assert dev_check.passed is True
        assert "CUDA" in dev_check.message
        assert "NVIDIA A100" in dev_check.message

    def test_mps_device(self):
        mock_torch = MagicMock()
        mock_torch.__version__ = "2.3.0"
        mock_torch.cuda.is_available.return_value = False
        mock_torch.backends.mps.is_available.return_value = True

        with patch.dict(sys.modules, {"torch": mock_torch}):
            result = check_compute()

        dev_check = next(c for c in result.checks if c.name == "device")
        assert dev_check.passed is True
        assert "MPS" in dev_check.message

    def test_cpu_only(self):
        mock_torch = MagicMock()
        mock_torch.__version__ = "2.2.0"
        mock_torch.cuda.is_available.return_value = False
        mock_torch.backends.mps.is_available.return_value = False

        with patch.dict(sys.modules, {"torch": mock_torch}):
            result = check_compute()

        dev_check = next(c for c in result.checks if c.name == "device")
        assert dev_check.passed is True
        assert "CPU" in dev_check.message

    def test_pytorch_not_installed(self):
        import builtins
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "torch":
                raise ImportError("no torch")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=mock_import):
            result = check_compute()

        pt_check = next(c for c in result.checks if c.name == "pytorch")
        assert pt_check.passed is False
        assert "not installed" in pt_check.message
        # Should only have pytorch check (device skipped)
        assert len(result.checks) == 1


# ---- JSON output -----------------------------------------------------------

class TestJSONOutput:
    def test_json_format_valid(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        output = format_json(report)
        data = json.loads(output)
        assert "total_checks" in data
        assert "passed" in data
        assert "warnings" in data
        assert "failures" in data
        assert "sections" in data
        assert isinstance(data["sections"], list)

    def test_json_sections_structure(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        data = json.loads(format_json(report))
        for section in data["sections"]:
            assert "name" in section
            assert "overall" in section
            assert "checks" in section
            for check in section["checks"]:
                assert "name" in check
                assert "passed" in check
                assert "warning" in check
                assert "message" in check

    def test_json_round_trip_totals(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        data = json.loads(format_json(report))
        assert data["total_checks"] == report.total_checks
        assert data["passed"] == report.total_passed
        assert data["warnings"] == report.total_warnings
        assert data["failures"] == report.total_failures


# ---- Text output -----------------------------------------------------------

class TestTextOutput:
    def test_text_format_has_header(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        output = format_text(report)
        assert "TCRpredictor Configuration Check" in output
        assert "=" * 35 in output

    def test_text_format_has_summary(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        output = format_text(report)
        assert "Summary:" in output
        assert "checks passed" in output

    def test_verbose_shows_details(self, tmp_path):
        root = _make_root(tmp_path)
        ckpt = root / MODEL_CHECKPOINT_PATH
        ckpt.parent.mkdir(parents=True, exist_ok=True)
        ckpt.write_bytes(b"\x00")
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        verbose_output = format_text(report, verbose=True)
        plain_output = format_text(report, verbose=False)
        # Verbose output should be longer due to detail lines
        assert len(verbose_output) > len(plain_output)


# ---- --fix mode ------------------------------------------------------------

class TestFixMode:
    def test_fix_creates_missing_directories(self, tmp_path):
        root = _make_root(tmp_path, dirs=False)
        result = check_file_system(root, fix=True)
        dirs_check = next(c for c in result.checks if c.name == "directories")
        assert dirs_check.passed is True
        for d in REQUIRED_DIRS:
            assert (root / d).is_dir()

    def test_fix_idempotent(self, tmp_path):
        root = _make_root(tmp_path)
        # All dirs already exist
        result = check_file_system(root, fix=True)
        dirs_check = next(c for c in result.checks if c.name == "directories")
        assert dirs_check.passed is True
        assert "Created" not in dirs_check.message


# ---- CLI argument parsing --------------------------------------------------

class TestCLI:
    def test_default_args(self):
        parser = build_parser()
        args = parser.parse_args([])
        assert args.verbose is False
        assert args.json_output is False
        assert args.fix is False

    def test_verbose_flag(self):
        parser = build_parser()
        args = parser.parse_args(["--verbose"])
        assert args.verbose is True

    def test_verbose_short(self):
        parser = build_parser()
        args = parser.parse_args(["-v"])
        assert args.verbose is True

    def test_json_flag(self):
        parser = build_parser()
        args = parser.parse_args(["--json"])
        assert args.json_output is True

    def test_fix_flag(self):
        parser = build_parser()
        args = parser.parse_args(["--fix"])
        assert args.fix is True

    def test_all_flags(self):
        parser = build_parser()
        args = parser.parse_args(["--verbose", "--json", "--fix"])
        assert args.verbose is True
        assert args.json_output is True
        assert args.fix is True

    def test_main_returns_zero_on_success(self, tmp_path, capsys):
        root = _make_root(tmp_path)
        # Mock everything to pass
        mock_report = ValidationReport(sections=[
            SectionResult(name="test", checks=[
                CheckResult("a", passed=True, message="ok"),
            ]),
        ])
        with patch("scripts.validate_config.run_all_checks", return_value=mock_report):
            exit_code = main([])
        assert exit_code == 0

    def test_main_returns_one_on_failures(self, capsys):
        mock_report = ValidationReport(sections=[
            SectionResult(name="test", checks=[
                CheckResult("a", passed=False, message="fail"),
            ]),
        ])
        with patch("scripts.validate_config.run_all_checks", return_value=mock_report):
            exit_code = main([])
        assert exit_code == 1

    def test_main_json_output(self, capsys):
        mock_report = ValidationReport(sections=[
            SectionResult(name="test", checks=[
                CheckResult("a", passed=True, message="ok"),
            ]),
        ])
        with patch("scripts.validate_config.run_all_checks", return_value=mock_report):
            main(["--json"])
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "total_checks" in data

    def test_main_text_output(self, capsys):
        mock_report = ValidationReport(sections=[
            SectionResult(name="test", checks=[
                CheckResult("a", passed=True, message="ok"),
            ]),
        ])
        with patch("scripts.validate_config.run_all_checks", return_value=mock_report):
            main([])
        captured = capsys.readouterr()
        assert "TCRpredictor Configuration Check" in captured.out


# ---- Integration: run_all_checks ------------------------------------------

class TestRunAllChecks:
    def test_all_sections_present(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        names = {s.name for s in report.sections}
        assert "File System" in names
        assert "Environment" in names
        assert "Model Readiness" in names
        assert "Compute" in names

    def test_services_included_by_default(self, tmp_path):
        root = _make_root(tmp_path)

        mock_resp = MagicMock()
        mock_resp.status_code = 200
        fake_httpx = types.ModuleType("httpx")
        fake_httpx.get = MagicMock(return_value=mock_resp)

        with patch.dict(sys.modules, {"httpx": fake_httpx}):
            with patch("scripts.validate_config.importlib.import_module",
                       return_value=True):
                report = run_all_checks(root, skip_services=False)

        names = {s.name for s in report.sections}
        assert "Services" in names

    def test_skip_services(self, tmp_path):
        root = _make_root(tmp_path)
        with patch("scripts.validate_config.importlib.import_module", return_value=True):
            report = run_all_checks(root, skip_services=True)
        names = {s.name for s in report.sections}
        assert "Services" not in names
