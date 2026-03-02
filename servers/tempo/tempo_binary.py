"""TEMPO binary wrapper — calls the official TEMPO CLI for predictions.

Falls back to the Python scorer (scorer.py) when the binary is not available.
"""
from __future__ import annotations

import csv
import io
import logging
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_BINARY_DIR = Path(__file__).resolve().parent.parent.parent / "bin" / "tempo"
_COMMENT_PREFIX = "#"


def _find_binary() -> Optional[Path]:
    """Locate the TEMPO binary for the current platform."""
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "darwin":
        if machine == "arm64":
            candidate = _BINARY_DIR / "TEMPO-macos-arm64-v1.0.0" / "TEMPO"
        else:
            candidate = _BINARY_DIR / "TEMPO-macos-x86_64-v1.0.0" / "TEMPO"
    elif system == "linux":
        candidate = _BINARY_DIR / "TEMPO-linux-x86_64-v1.0.0" / "TEMPO"
    else:
        return None

    if candidate.exists() and os.access(candidate, os.X_OK):
        return candidate

    # Also check PATH
    which = shutil.which("TEMPO")
    if which:
        return Path(which)

    return None


TEMPO_BINARY = _find_binary()


@dataclass
class TempoPrediction:
    """Result from the TEMPO binary for a single TCR."""
    perc_rank: Optional[float]
    problem: str = ""
    # Original input echoed back
    trav: str = ""
    traj: str = ""
    cdr3_tra: str = ""
    trbv: str = ""
    trbj: str = ""
    cdr3_trb: str = ""


def is_available() -> bool:
    """Check if the TEMPO binary is available."""
    return TEMPO_BINARY is not None


def list_epitopes() -> list[str]:
    """Return list of available epitope model IDs."""
    if not TEMPO_BINARY:
        return []
    try:
        result = subprocess.run(
            [str(TEMPO_BINARY), "list_epitopes", "--csv"],
            capture_output=True, text=True, timeout=10,
        )
        epitopes = []
        for line in result.stdout.splitlines():
            if line.startswith("#") or line.startswith("model"):
                continue
            parts = line.split(",")
            if parts:
                epitopes.append(parts[0].strip())
        return [e for e in epitopes if e]
    except Exception as e:
        logger.warning("Failed to list TEMPO epitopes: %s", e)
        return []


def predict(
    tcrs: list[dict[str, str]],
    epitope_id: str,
    chain: str = "AB",
    species: str = "HomoSapiens",
) -> list[TempoPrediction]:
    """Run TEMPO binary prediction on a list of TCRs.

    Args:
        tcrs: List of dicts with keys matching TEMPO CSV columns.
              For chain AB: trav, traj, cdr3_tra, trbv, trbj, cdr3_trb
              For chain A: trav, traj, cdr3_tra
              For chain B: trbv, trbj, cdr3_trb
        epitope_id: TEMPO epitope model ID (e.g., "A0201_GILGFVFTL")
        chain: "AB", "A", or "B"
        species: "HomoSapiens" or "MusMusculus"

    Returns:
        List of TempoPrediction results, one per input TCR.
    """
    if not TEMPO_BINARY:
        raise RuntimeError("TEMPO binary not found")

    if not tcrs:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "input.csv"
        output_path = Path(tmpdir) / "output.csv"

        # Write input CSV
        if chain == "AB":
            fieldnames = ["TRAV", "TRAJ", "cdr3_TRA", "TRBV", "TRBJ", "cdr3_TRB"]
        elif chain == "A":
            fieldnames = ["TRAV", "TRAJ", "cdr3_TRA"]
        else:
            fieldnames = ["TRBV", "TRBJ", "cdr3_TRB"]

        with open(input_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for tcr in tcrs:
                row = {}
                if "TRAV" in fieldnames:
                    row["TRAV"] = tcr.get("trav", tcr.get("TRAV", tcr.get("v_a_gene", "")))
                    row["TRAJ"] = tcr.get("traj", tcr.get("TRAJ", tcr.get("j_a_gene", "")))
                    row["cdr3_TRA"] = tcr.get("cdr3_tra", tcr.get("cdr3_TRA", tcr.get("cdr3_a", "")))
                if "TRBV" in fieldnames:
                    row["TRBV"] = tcr.get("trbv", tcr.get("TRBV", tcr.get("v_b_gene", tcr.get("v_gene", ""))))
                    row["TRBJ"] = tcr.get("trbj", tcr.get("TRBJ", tcr.get("j_b_gene", tcr.get("j_gene", ""))))
                    row["cdr3_TRB"] = tcr.get("cdr3_trb", tcr.get("cdr3_TRB", tcr.get("cdr3_b", tcr.get("cdr3", ""))))
                writer.writerow(row)

        # Run TEMPO
        cmd = [
            str(TEMPO_BINARY), "predict",
            str(input_path), str(output_path),
            epitope_id,
            "--chain", chain,
            "--species", species,
            "--no-qc",
        ]
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=120,
                cwd=str(TEMPO_BINARY.parent),
            )
            if result.returncode != 0:
                logger.error("TEMPO failed: %s", result.stderr)
                return [TempoPrediction(perc_rank=None, problem=f"TEMPO error: {result.stderr[:200]}")]
        except subprocess.TimeoutExpired:
            logger.error("TEMPO timed out")
            return [TempoPrediction(perc_rank=None, problem="timeout")]

        # Parse output CSV
        if not output_path.exists():
            return [TempoPrediction(perc_rank=None, problem="no output file")]

        predictions = []
        with open(output_path, "r") as f:
            # Skip comment lines
            lines = [line for line in f if not line.startswith(_COMMENT_PREFIX)]

        # The TEMPO binary sometimes uses a split-column CSV layout where each
        # logical row (header or data) spans two physical lines.  The first line
        # carries the TCR-gene columns; the second line carries the score columns
        # (perc_rank, problem) at a shifted offset (the first field is empty).
        #
        # Example split-column header:
        #   Line A: "TRAV,TRAJ,cdr3_TRA,TRBV,TRBJ,cdr3_TRB"
        #   Line B: ",perc_rank,problem"
        #
        # Detect this by checking whether "perc_rank" is absent from the first
        # non-comment line but present in the second.
        if len(lines) >= 2:
            first_header_fields = next(csv.reader([lines[0]]))
            second_header_fields = next(csv.reader([lines[1]]))

            if "perc_rank" not in first_header_fields and "perc_rank" in second_header_fields:
                # Build a unified column list: take all fields from the first
                # header line, then append the non-empty fields from the second.
                unified_header = list(first_header_fields)
                score_fields = [f for f in second_header_fields if f.strip()]
                unified_header.extend(score_fields)

                # The offset within the second line where the score values live
                # is determined by the position of the first non-empty field.
                score_offset = next(
                    (i for i, f in enumerate(second_header_fields) if f.strip()), 0
                )

                data_lines = lines[2:]  # Skip the two header lines
                for idx in range(0, len(data_lines) - 1, 2):
                    tcr_fields = next(csv.reader([data_lines[idx]]))
                    score_line_fields = next(csv.reader([data_lines[idx + 1]]))

                    row: dict[str, str] = {}
                    for col_i, col_name in enumerate(first_header_fields):
                        row[col_name] = tcr_fields[col_i] if col_i < len(tcr_fields) else ""

                    for s_i, col_name in enumerate(score_fields):
                        actual_i = score_offset + s_i
                        row[col_name] = (
                            score_line_fields[actual_i]
                            if actual_i < len(score_line_fields)
                            else ""
                        )

                    perc_rank_str = row.get("perc_rank", "").strip()
                    try:
                        perc_rank = float(perc_rank_str) if perc_rank_str else None
                    except ValueError:
                        perc_rank = None

                    predictions.append(TempoPrediction(
                        perc_rank=perc_rank,
                        problem=row.get("problem", ""),
                        trav=row.get("TRAV", ""),
                        traj=row.get("TRAJ", ""),
                        cdr3_tra=row.get("cdr3_TRA", ""),
                        trbv=row.get("TRBV", ""),
                        trbj=row.get("TRBJ", ""),
                        cdr3_trb=row.get("cdr3_TRB", ""),
                    ))
                return predictions

        # Standard single-header CSV (all columns on one header line)
        reader = csv.DictReader(io.StringIO("".join(lines)))
        for row in reader:
            perc_rank_str = row.get("perc_rank", "").strip()
            try:
                perc_rank = float(perc_rank_str) if perc_rank_str else None
            except ValueError:
                perc_rank = None

            predictions.append(TempoPrediction(
                perc_rank=perc_rank,
                problem=row.get("problem", ""),
                trav=row.get("TRAV", ""),
                traj=row.get("TRAJ", ""),
                cdr3_tra=row.get("cdr3_TRA", ""),
                trbv=row.get("TRBV", ""),
                trbj=row.get("TRBJ", ""),
                cdr3_trb=row.get("cdr3_TRB", ""),
            ))

        return predictions


def perc_rank_to_score(perc_rank: float) -> float:
    """Convert TEMPO percentile rank (lower=better) to 0-1 score (higher=better)."""
    return max(0.0, min(1.0, 1.0 - perc_rank / 100.0))
