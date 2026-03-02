"""Convert IEDB T-cell assay records → BATMAN training DataFrame."""
from __future__ import annotations

from typing import Any

import pandas as pd

# IEDB qualitative_measure → activation integer
_MEASURE_MAP: dict[str, int] = {
    "positive": 2,
    "positive-high": 2,
    "strong positive": 2,
    "intermediate": 1,
    "positive-low": 1,
    "weak positive": 1,
    "negative": 0,
}


def _map_activation(measure: str | None) -> int | None:
    if not measure:
        return None
    return _MEASURE_MAP.get(measure.lower().strip())


def iedb_records_to_batman_df(
    records: list[dict[str, Any]],
    tcr_id: str,
    index_peptide: str,
) -> pd.DataFrame:
    """Map IEDB T-cell assay records to pybatman CSV format.

    Drops rows with unknown activation labels.
    Deduplicates peptides by keeping the highest activation value.
    """
    rows = []
    for rec in records:
        peptide = (rec.get("sequence") or "").upper().strip()
        measure = (rec.get("metadata") or {}).get("qualitative_measure")
        activation = _map_activation(measure)
        if not peptide or activation is None:
            continue
        rows.append({"tcr": tcr_id, "index": index_peptide,
                     "peptide": peptide, "activation": activation})

    if not rows:
        return pd.DataFrame(columns=["tcr", "index", "peptide", "activation"])

    df = pd.DataFrame(rows)
    # Deduplicate: keep highest activation per peptide
    df = df.sort_values("activation", ascending=False).drop_duplicates("peptide")
    return df.reset_index(drop=True)


def merge_user_data(
    iedb_df: pd.DataFrame,
    user_data: list[dict[str, Any]],
    tcr_id: str,
    index_peptide: str,
) -> pd.DataFrame:
    """Overlay user-provided activation data on top of IEDB data.

    User entries take precedence over IEDB entries for the same peptide.
    """
    if not user_data:
        return iedb_df

    user_df = pd.DataFrame([
        {"tcr": tcr_id, "index": index_peptide,
         "peptide": row["peptide"].upper().strip(),
         "activation": int(row["activation"])}
        for row in user_data
        if "peptide" in row and "activation" in row
    ])

    # Remove IEDB rows that user has overridden
    overridden = set(user_df["peptide"].str.upper())
    base = iedb_df[~iedb_df["peptide"].str.upper().isin(overridden)]
    return pd.concat([base, user_df], ignore_index=True)
