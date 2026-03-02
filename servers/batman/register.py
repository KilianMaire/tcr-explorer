"""Class II peptide register modeling.

Class II peptides bind in multiple registers (shifted 9-mer cores within
longer peptides). This module identifies viable registers and scores them.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class RegisterPrediction:
    """A single register prediction for a Class II peptide."""
    core: str
    core_start: int
    core_end: int
    pfr_n: str  # N-terminal peptide flanking residue(s)
    pfr_c: str  # C-terminal peptide flanking residue(s)
    score: float = 0.0


def identify_registers(peptide: str, core_length: int = 9) -> list[RegisterPrediction]:
    """Identify all possible 9-mer binding registers within a peptide."""
    if len(peptide) <= core_length:
        return [RegisterPrediction(
            core=peptide[:core_length] if len(peptide) >= core_length else peptide,
            core_start=0,
            core_end=min(len(peptide), core_length),
            pfr_n="",
            pfr_c="",
        )]
    registers = []
    for start in range(len(peptide) - core_length + 1):
        end = start + core_length
        registers.append(RegisterPrediction(
            core=peptide[start:end],
            core_start=start,
            core_end=end,
            pfr_n=peptide[:start],
            pfr_c=peptide[end:],
        ))
    return registers


def score_registers(registers: list[RegisterPrediction]) -> list[RegisterPrediction]:
    """Sort registers by score descending (best first)."""
    return sorted(registers, key=lambda r: r.score, reverse=True)
