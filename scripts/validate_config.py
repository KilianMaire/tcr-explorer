#!/usr/bin/env python3
"""Pre-flight configuration validation for TCRpredictor.

Checks file system layout, environment variables, service health,
model readiness, and compute devices before running the pipeline.

Usage:
    python scripts/validate_config.py [--verbose] [--json] [--fix]

Options:
    --verbose   Show detailed output for each check
    --json      Output results as JSON (for CI integration)
    --fix       Attempt to fix issues (create missing dirs, etc.)
"""
from __future__ import annotations

import argparse
import importlib
import json as json_mod
import logging
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Pipeline scripts that must exist
PIPELINE_SCRIPTS = [
    "scripts/download_class2_data.py",
    "scripts/generate_esm2_embeddings.py",
    "scripts/prepare_training_data.py",
    "scripts/train_model.py",
    "scripts/run_pipeline.py",
    "scripts/build_mhc_pseudoseq.py",
]

# Required directories
REQUIRED_DIRS = [
    "data",
    "data/models",
    "data/models/tier1",
    "data/processed",
    "data/features",
    "data/raw",
    "models",
    "pipeline",
    "servers",
    "mlops",
    "structural",
]

# Python modules that must be importable
IMPORTABLE_MODULES = [
    "pipeline.orchestrator",
    "models.esm2_finetune",
    "models.transformer",
    "mlops.trainer",
]

# Environment variable registry: (name, required, description)
ENV_VARS = [
    ("TIER1_SERVER_URL", True, "Tier 1 ML screening server URL"),
    ("STRUCTURAL_SERVER_URL", False, "Tier 3 structural analysis server URL"),
    ("HLA_SERVER_URL", False, "HLA microservice URL"),
    ("TCR_SERVER_URL", False, "TCR microservice URL"),
    ("VDJDB_SERVER_URL", False, "VDJdb microservice URL"),
    ("IEDB_SERVER_URL", False, "IEDB microservice URL"),
    ("BATMAN_SERVER_URL", False, "BATMAN scoring server URL"),
    ("TEMPO_SERVER_URL", False, "TEMPO scoring server URL"),
    ("IMGT_API_URL", False, "Main API gateway URL"),
    ("NCBI_API_KEY", True, "NCBI API key (rate-limited without it)"),
    ("DATABASE_PATH", False, "Path to IMGT SQLite database"),
    ("ESM2_MODEL_NAME", False, "ESM-2 model identifier"),
]

# Service registry: name -> (env_var, default_url, health_path)
SERVICES: Dict[str, Tuple[str, str, str]] = {
    "API Gateway": ("IMGT_API_URL", "http://localhost:8000", "/health"),
    "Tier 1 ML": ("TIER1_SERVER_URL", "http://localhost:8110", "/health"),
    "Structural": ("STRUCTURAL_SERVER_URL", "http://localhost:8120", "/health"),
    "HLA Server": ("HLA_SERVER_URL", "http://localhost:8101", "/health"),
    "TCR Server": ("TCR_SERVER_URL", "http://localhost:8102", "/health"),
    "VDJdb Server": ("VDJDB_SERVER_URL", "http://localhost:8103", "/health"),
    "IEDB Server": ("IEDB_SERVER_URL", "http://localhost:8104", "/health"),
}

# Model paths relative to project root
MODEL_CHECKPOINT_PATH = "data/models/tier1/best_model.pt"
TRAINING_DATA_PATTERNS = [
    "data/processed/class2_pairs.parquet",
    "data/processed/training_pairs.parquet",
]
EMBEDDING_DIRS = [
    "data/features/esm2_peptide",
    "data/features/esm2_tcr",
]

# Health-check timeout in seconds
HEALTH_TIMEOUT = 5.0


# ---------------------------------------------------------------------------
# Result data classes
# ---------------------------------------------------------------------------

@dataclass
class CheckResult:
    """Result of a single check item."""
    name: str
    passed: bool
    warning: bool = False
    message: str = ""
    detail: str = ""

    @property
    def symbol(self) -> str:
        if self.passed:
            return "\u2713"
        if self.warning:
            return "~"
        return "\u2717"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "passed": self.passed,
            "warning": self.warning,
            "message": self.message,
            "detail": self.detail,
        }


@dataclass
class SectionResult:
    """Result of a check section (group of checks)."""
    name: str
    checks: List[CheckResult] = field(default_factory=list)

    @property
    def passed(self) -> int:
        return sum(1 for c in self.checks if c.passed)

    @property
    def warnings(self) -> int:
        return sum(1 for c in self.checks if c.warning and not c.passed)

    @property
    def failures(self) -> int:
        return sum(1 for c in self.checks if not c.passed and not c.warning)

    @property
    def overall_symbol(self) -> str:
        if all(c.passed for c in self.checks):
            return "\u2713"
        if any(not c.passed and not c.warning for c in self.checks):
            return "\u2717"
        return "~"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "overall": self.overall_symbol,
            "passed": self.passed,
            "warnings": self.warnings,
            "failures": self.failures,
            "checks": [c.to_dict() for c in self.checks],
        }


@dataclass
class ValidationReport:
    """Full validation report."""
    sections: List[SectionResult] = field(default_factory=list)

    @property
    def total_checks(self) -> int:
        return sum(len(s.checks) for s in self.sections)

    @property
    def total_passed(self) -> int:
        return sum(s.passed for s in self.sections)

    @property
    def total_warnings(self) -> int:
        return sum(s.warnings for s in self.sections)

    @property
    def total_failures(self) -> int:
        return sum(s.failures for s in self.sections)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "total_checks": self.total_checks,
            "passed": self.total_passed,
            "warnings": self.total_warnings,
            "failures": self.total_failures,
            "sections": [s.to_dict() for s in self.sections],
        }


# ---------------------------------------------------------------------------
# Check functions
# ---------------------------------------------------------------------------

def check_file_system(root: Path, *, fix: bool = False) -> SectionResult:
    """Validate that required files and directories exist."""
    section = SectionResult(name="File System")

    # Pipeline scripts
    missing_scripts: List[str] = []
    for script in PIPELINE_SCRIPTS:
        if not (root / script).exists():
            missing_scripts.append(script)

    if missing_scripts:
        section.checks.append(CheckResult(
            name="pipeline_scripts",
            passed=False,
            message=f"{len(missing_scripts)}/{len(PIPELINE_SCRIPTS)} pipeline scripts missing",
            detail=", ".join(missing_scripts),
        ))
    else:
        section.checks.append(CheckResult(
            name="pipeline_scripts",
            passed=True,
            message=f"All {len(PIPELINE_SCRIPTS)} pipeline scripts present",
        ))

    # Required directories
    missing_dirs: List[str] = []
    created_dirs: List[str] = []
    for d in REQUIRED_DIRS:
        p = root / d
        if not p.is_dir():
            if fix:
                p.mkdir(parents=True, exist_ok=True)
                created_dirs.append(d)
            else:
                missing_dirs.append(d)

    if created_dirs:
        section.checks.append(CheckResult(
            name="directories",
            passed=True,
            message=f"Created {len(created_dirs)} missing directories",
            detail=", ".join(created_dirs),
        ))
    elif missing_dirs:
        section.checks.append(CheckResult(
            name="directories",
            passed=False,
            message=f"{len(missing_dirs)} required directories missing",
            detail=", ".join(missing_dirs),
        ))
    else:
        section.checks.append(CheckResult(
            name="directories",
            passed=True,
            message=f"All {len(REQUIRED_DIRS)} required directories exist",
        ))

    # Python module imports
    failed_imports: List[str] = []
    for mod_name in IMPORTABLE_MODULES:
        try:
            importlib.import_module(mod_name)
        except Exception:
            failed_imports.append(mod_name)

    if failed_imports:
        section.checks.append(CheckResult(
            name="python_imports",
            passed=False,
            warning=True,
            message=f"{len(failed_imports)} modules failed to import",
            detail=", ".join(failed_imports),
        ))
    else:
        section.checks.append(CheckResult(
            name="python_imports",
            passed=True,
            message=f"All {len(IMPORTABLE_MODULES)} Python modules importable",
        ))

    # DVC initialised
    dvc_dir = root / ".dvc"
    section.checks.append(CheckResult(
        name="dvc_initialized",
        passed=dvc_dir.is_dir(),
        warning=not dvc_dir.is_dir(),
        message="DVC initialized" if dvc_dir.is_dir() else ".dvc/ directory not found",
    ))

    return section


def _load_dotenv(root: Path) -> Dict[str, str]:
    """Read .env file if it exists, returning key=value pairs."""
    env_file = root / ".env"
    result: Dict[str, str] = {}
    if not env_file.exists():
        return result
    for line in env_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            key, _, value = line.partition("=")
            result[key.strip()] = value.strip().strip("\"'")
    return result


def check_environment(root: Path) -> SectionResult:
    """Validate environment variables."""
    section = SectionResult(name="Environment")

    dotenv_vals = _load_dotenv(root)
    # Merge: dotenv values are overridden by real env
    merged = {**dotenv_vals}
    merged.update({k: v for k, v in os.environ.items()})

    for var_name, required, description in ENV_VARS:
        value = merged.get(var_name)
        if value:
            section.checks.append(CheckResult(
                name=var_name,
                passed=True,
                message=f"{var_name} = {value}",
            ))
        elif required:
            section.checks.append(CheckResult(
                name=var_name,
                passed=False,
                message=f"{var_name} not set ({description})",
            ))
        else:
            section.checks.append(CheckResult(
                name=var_name,
                passed=False,
                warning=True,
                message=f"{var_name} not set ({description} - optional)",
            ))

    return section


def check_services() -> SectionResult:
    """Check health endpoints for configured services."""
    section = SectionResult(name="Services")

    try:
        import httpx  # noqa: F811
        _httpx_available = True
    except ImportError:
        _httpx_available = False

    for svc_name, (env_var, default_url, health_path) in SERVICES.items():
        url = os.environ.get(env_var)
        if url is None:
            url = default_url
            configured = False
        else:
            configured = True

        if not _httpx_available:
            section.checks.append(CheckResult(
                name=svc_name,
                passed=False,
                warning=True,
                message=f"{svc_name}: httpx not installed, cannot check",
            ))
            continue

        full_url = url.rstrip("/") + health_path
        try:
            resp = httpx.get(full_url, timeout=HEALTH_TIMEOUT)
            if resp.status_code < 400:
                section.checks.append(CheckResult(
                    name=svc_name,
                    passed=True,
                    message=f"{svc_name}: up (HTTP {resp.status_code})",
                    detail=full_url,
                ))
            else:
                section.checks.append(CheckResult(
                    name=svc_name,
                    passed=False,
                    message=f"{svc_name}: unhealthy (HTTP {resp.status_code})",
                    detail=full_url,
                ))
        except Exception:
            status = "not configured" if not configured else "not running"
            section.checks.append(CheckResult(
                name=svc_name,
                passed=False,
                message=f"{svc_name}: {status}",
                detail=full_url,
            ))

    return section


def check_model_readiness(root: Path) -> SectionResult:
    """Check model checkpoint, training data, and embeddings."""
    section = SectionResult(name="Model Readiness")

    # Model checkpoint
    ckpt = root / MODEL_CHECKPOINT_PATH
    section.checks.append(CheckResult(
        name="model_checkpoint",
        passed=ckpt.exists(),
        message="Model checkpoint found" if ckpt.exists()
                else "No model checkpoint found (model not yet trained)",
        detail=str(ckpt),
    ))

    # Training data
    found_data = False
    for pattern in TRAINING_DATA_PATTERNS:
        if (root / pattern).exists():
            found_data = True
            section.checks.append(CheckResult(
                name="training_data",
                passed=True,
                message=f"Training data exists ({Path(pattern).name})",
                detail=str(root / pattern),
            ))
            break
    if not found_data:
        section.checks.append(CheckResult(
            name="training_data",
            passed=False,
            warning=True,
            message="No training data found",
            detail=", ".join(TRAINING_DATA_PATTERNS),
        ))

    # ESM-2 embeddings
    emb_found = any((root / d).is_dir() and any((root / d).iterdir())
                     for d in EMBEDDING_DIRS
                     if (root / d).is_dir())
    if emb_found:
        section.checks.append(CheckResult(
            name="embeddings",
            passed=True,
            message="ESM-2 embeddings generated",
        ))
    else:
        section.checks.append(CheckResult(
            name="embeddings",
            passed=False,
            warning=True,
            message="ESM-2 embeddings not generated",
            detail=", ".join(EMBEDDING_DIRS),
        ))

    return section


def check_compute() -> SectionResult:
    """Detect available compute devices."""
    section = SectionResult(name="Compute")

    # PyTorch availability
    try:
        import torch
        torch_version = torch.__version__
        section.checks.append(CheckResult(
            name="pytorch",
            passed=True,
            message=f"PyTorch {torch_version} available",
        ))
    except ImportError:
        section.checks.append(CheckResult(
            name="pytorch",
            passed=False,
            message="PyTorch not installed",
        ))
        return section

    # Device detection
    if torch.cuda.is_available():
        device_name = torch.cuda.get_device_name(0)
        try:
            vram_mb = torch.cuda.get_device_properties(0).total_mem // (1024 * 1024)
            detail = f"{device_name}, {vram_mb} MB VRAM"
        except Exception:
            detail = device_name
        section.checks.append(CheckResult(
            name="device",
            passed=True,
            message=f"Device: CUDA ({device_name})",
            detail=detail,
        ))
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        section.checks.append(CheckResult(
            name="device",
            passed=True,
            message="Device: MPS (Apple Silicon)",
        ))
    else:
        section.checks.append(CheckResult(
            name="device",
            passed=True,
            warning=False,
            message="Device: CPU only",
            detail="GPU not detected; training will be slow",
        ))

    return section


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def format_text(report: ValidationReport, *, verbose: bool = False) -> str:
    """Format the report as coloured terminal text."""
    lines: List[str] = [
        "TCRpredictor Configuration Check",
        "=" * 35,
        "",
    ]
    for section in report.sections:
        lines.append(f"[{section.overall_symbol}] {section.name}")
        for check in section.checks:
            prefix = f"    {check.symbol} "
            lines.append(f"{prefix}{check.message}")
            if verbose and check.detail:
                lines.append(f"        {check.detail}")
        lines.append("")

    lines.append(
        f"Summary: {report.total_passed}/{report.total_checks} checks passed, "
        f"{report.total_warnings} warnings, {report.total_failures} critical"
    )
    return "\n".join(lines)


def format_json(report: ValidationReport) -> str:
    """Format the report as JSON."""
    return json_mod.dumps(report.to_dict(), indent=2)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run_all_checks(
    root: Optional[Path] = None,
    *,
    fix: bool = False,
    skip_services: bool = False,
) -> ValidationReport:
    """Execute every check section and return a consolidated report."""
    if root is None:
        root = PROJECT_ROOT

    report = ValidationReport()
    report.sections.append(check_file_system(root, fix=fix))
    report.sections.append(check_environment(root))
    if not skip_services:
        report.sections.append(check_services())
    report.sections.append(check_model_readiness(root))
    report.sections.append(check_compute())
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Pre-flight configuration validation for TCRpredictor",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output for each check",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output results as JSON (for CI integration)",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Attempt to fix issues (create missing dirs, etc.)",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Entry-point for the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    report = run_all_checks(fix=args.fix)

    if args.json_output:
        print(format_json(report))
    else:
        print(format_text(report, verbose=args.verbose))

    # Exit code: 0 if no critical failures, 1 otherwise
    return 1 if report.total_failures > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
