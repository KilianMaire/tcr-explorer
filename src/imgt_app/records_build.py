"""Harmonize VDJdb, IEDB, McPAS, and TCR3d raw snapshots into one record index.

Each `harmonize_*` function reads one raw source and returns a DataFrame with
columns `SCHEMA_COLUMNS`, one row per chain. Deposited sequences (nt, full aa,
CDR1/2) are kept verbatim when the source provides them; nothing is fabricated.
A value that cannot be filled is None.
"""
from __future__ import annotations

import json
import zipfile
from pathlib import Path
from typing import Optional

import pandas as pd

from .cdr_enricher import _gene_to_chain
from .input_router import _normalize_gene

SCHEMA_COLUMNS: list[str] = [
    "source",
    "source_record_id",
    "pairing_key",
    "species",
    "chain",
    "cdr3_aa",
    "cdr3_nt",
    "full_aa",
    "full_nt",
    "v_gene",
    "d_gene",
    "j_gene",
    "cdr1_aa",
    "cdr2_aa",
    "epitope_aa",
    "antigen",
    "antigen_organism",
    "mhc_class",
    "mhc_a",
    "mhc_b",
    "pdb_id",
    "reference_pmid",
    "external_url",
    "score",
]

# The design spec's table lists 23 columns and the brief's "22 columns above"
# comment undercounts by one; both omit full_nt even though the brief's own
# IEDB mapping rule requires it (deposited domain nt is not the same field as
# the CDR3 nt). SCHEMA_COLUMNS here has 24 entries: the spec's 23 plus full_nt.

_NA_STRINGS = {"", "na", "nan", "-"}

_SPECIES_MAP = {
    "homosapiens": "human",
    "musmusculus": "mouse",
}


def _na(v: object) -> Optional[str]:
    """Normalize a raw cell into None when it carries no real value."""
    if v is None:
        return None
    if isinstance(v, float) and pd.isna(v):
        return None
    s = str(v).strip()
    if s.lower() in _NA_STRINGS:
        return None
    return s


def normalize_species(raw: str) -> str:
    key = (raw or "").strip().lower().replace(" ", "")
    return _SPECIES_MAP.get(key, key)


# IEDB carries species per chain only as an NCBITaxon ontology IRI (no plain
# species name column exists in tcr_full_v3.csv), e.g.
# "http://purl.obolibrary.org/obo/NCBITaxon_9606". Map the common lab taxa
# explicitly; anything unmapped keeps the taxon id, never guessed.
_TAXON_SPECIES_MAP = {
    "9606": "human",
    "10090": "mouse",
}


def _species_from_organism_iri(raw: object) -> Optional[str]:
    s = _na(raw)
    if s is None:
        return None
    tail = s.rstrip("/").rsplit("/", 1)[-1]
    taxon_id = tail.rsplit("_", 1)[-1] if "_" in tail else tail
    return _TAXON_SPECIES_MAP.get(taxon_id, tail.lower())


def _normalize_gene_or_none(raw: object) -> Optional[str]:
    g = _na(raw)
    if g is None:
        return None
    try:
        return _normalize_gene(g)
    except Exception:
        return g


def _chain_from_gene(*genes: object) -> Optional[str]:
    """Alpha/beta from a gene's TRA/TRB prefix.

    ``cdr_enricher._gene_to_chain`` defaults unmatched genes to "TRB", so we
    only trust it when the raw string itself starts with a TRA/TRB segment
    prefix (V/D/J/C) — otherwise a gamma/delta or garbage gene would silently
    resolve as beta.
    """
    for g in genes:
        s = _na(g)
        if not s:
            continue
        gu = s.strip().upper()
        if gu.startswith("TRA") or gu.startswith("TRB"):
            chain = _gene_to_chain(gu)
            if chain == "TRA":
                return "alpha"
            if chain == "TRB":
                return "beta"
    return None


def _looks_like_pmid(raw: object) -> Optional[str]:
    s = _na(raw)
    if s is None:
        return None
    stripped = s.replace("PMID:", "").replace("PMID", "").strip()
    if stripped.isdigit():
        return stripped
    if stripped.lower().startswith("10.") and "/" in stripped:
        return stripped
    return None


def _empty_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=SCHEMA_COLUMNS)


def _finalize(rows: list[dict]) -> pd.DataFrame:
    if not rows:
        return _empty_frame()
    df = pd.DataFrame(rows)
    for col in SCHEMA_COLUMNS:
        if col not in df.columns:
            df[col] = None
    df = df[SCHEMA_COLUMNS]
    df["cdr3_aa"] = df["cdr3_aa"].astype(str).str.strip().str.upper()
    df = df[df["cdr3_aa"].str.len() > 0].reset_index(drop=True)
    return df


def _blank_record() -> dict:
    return {c: None for c in SCHEMA_COLUMNS}


# ---------------------------------------------------------------------------
# VDJdb
# ---------------------------------------------------------------------------
def harmonize_vdjdb(txt_path: str) -> pd.DataFrame:
    raw = pd.read_csv(txt_path, sep="\t", dtype=str, keep_default_na=False)
    rows: list[dict] = []
    for i, r in raw.iterrows():
        gene = _na(r.get("gene"))
        chain = {"TRA": "alpha", "TRB": "beta"}.get((gene or "").upper())
        if chain is None:
            continue
        cdr3 = _na(r.get("cdr3"))
        if not cdr3:
            continue
        complex_id = _na(r.get("complex.id"))
        if complex_id is None or complex_id == "0":
            pairing_key = f"vdjdb:u{i}"
        else:
            pairing_key = f"vdjdb:c{complex_id}"

        rec = _blank_record()
        rec.update(
            {
                "source": "vdjdb",
                "source_record_id": f"vdjdb:{i}",
                "pairing_key": pairing_key,
                "species": normalize_species(_na(r.get("species")) or ""),
                "chain": chain,
                "cdr3_aa": cdr3,
                "v_gene": _normalize_gene_or_none(r.get("v.segm")),
                "j_gene": _normalize_gene_or_none(r.get("j.segm")),
                "epitope_aa": _na(r.get("antigen.epitope")),
                "antigen": _na(r.get("antigen.gene")),
                "antigen_organism": _na(r.get("antigen.species")),
                "mhc_class": _na(r.get("mhc.class")),
                "mhc_a": _na(r.get("mhc.a")),
                "mhc_b": _na(r.get("mhc.b")),
                "reference_pmid": _looks_like_pmid(r.get("reference.id")),
                "external_url": "https://vdjdb.cdr3.net/search",
                "score": _na(r.get("vdjdb.score")),
            }
        )
        rows.append(rec)
    return _finalize(rows)


# ---------------------------------------------------------------------------
# IEDB
# ---------------------------------------------------------------------------
def _iedb_flatten(path: str) -> pd.DataFrame:
    raw = pd.read_csv(path, header=[0, 1], dtype=str, keep_default_na=False)
    top = raw.columns.get_level_values(0)
    second = raw.columns.get_level_values(1)

    receptor_cols = {second[i]: i for i in range(len(top)) if top[i] in ("Receptor", "Reference", "Epitope")}
    chain1_cols = {second[i]: i for i in range(len(top)) if top[i] == "Chain 1"}
    chain2_cols = {second[i]: i for i in range(len(top)) if top[i] == "Chain 2"}
    return raw, receptor_cols, chain1_cols, chain2_cols


def _iedb_chain_row(
    values: pd.Series,
    cols: dict,
    *,
    receptor_id: str,
    pairing_key: str,
    epitope: Optional[str],
    antigen: Optional[str],
    antigen_organism: Optional[str],
    mhc_a: Optional[str],
    iri: Optional[str],
) -> Optional[dict]:
    def get(name: str) -> Optional[str]:
        idx = cols.get(name)
        if idx is None:
            return None
        return _na(values.iloc[idx])

    v_gene_raw = get("Curated V Gene") or get("Calculated V Gene")
    d_gene_raw = get("Curated D Gene") or get("Calculated D Gene")
    j_gene_raw = get("Curated J Gene") or get("Calculated J Gene")

    chain = _chain_from_gene(v_gene_raw, j_gene_raw)
    if chain is None:
        return None

    cdr3 = get("CDR3 Curated") or get("CDR3 Calculated")
    if not cdr3:
        return None

    rec = _blank_record()
    rec.update(
        {
            "source": "iedb",
            "source_record_id": f"iedb:{receptor_id}",
            "pairing_key": pairing_key,
            "species": _species_from_organism_iri(get("Organism IRI")),
            "chain": chain,
            "cdr3_aa": cdr3,
            "full_aa": get("Protein Sequence"),
            "full_nt": get("Nucleotide Sequence"),
            "v_gene": _normalize_gene_or_none(v_gene_raw),
            "d_gene": _normalize_gene_or_none(d_gene_raw),
            "j_gene": _normalize_gene_or_none(j_gene_raw),
            "cdr1_aa": get("CDR1 Curated") or get("CDR1 Calculated"),
            "cdr2_aa": get("CDR2 Curated") or get("CDR2 Calculated"),
            "epitope_aa": epitope,
            "antigen": antigen,
            "antigen_organism": antigen_organism,
            "mhc_a": mhc_a,
            "external_url": iri or f"https://www.iedb.org/receptor/{receptor_id}",
        }
    )
    return rec


def harmonize_iedb(csv_path: str) -> pd.DataFrame:
    raw, receptor_cols, chain1_cols, chain2_cols = _iedb_flatten(csv_path)
    rows: list[dict] = []
    for _, values in raw.iterrows():

        def rget(name: str) -> Optional[str]:
            idx = receptor_cols.get(name)
            if idx is None:
                return None
            return _na(values.iloc[idx])

        receptor_id = rget("IEDB Receptor ID")
        if not receptor_id:
            continue
        pairing_key = f"iedb:{receptor_id}"
        epitope = rget("Name")
        antigen = rget("Source Molecule")
        antigen_organism = rget("Source Organism")
        mhc_a = rget("MHC Allele Names")
        iri = rget("IEDB IRI")

        for cols in (chain1_cols, chain2_cols):
            rec = _iedb_chain_row(
                values,
                cols,
                receptor_id=receptor_id,
                pairing_key=pairing_key,
                epitope=epitope,
                antigen=antigen,
                antigen_organism=antigen_organism,
                mhc_a=mhc_a,
                iri=iri,
            )
            if rec is not None:
                rows.append(rec)
    return _finalize(rows)


# ---------------------------------------------------------------------------
# McPAS
# ---------------------------------------------------------------------------
def harmonize_mcpas(csv_path: str) -> pd.DataFrame:
    raw = pd.read_csv(csv_path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    rows: list[dict] = []
    for i, r in raw.iterrows():
        pairing_key = f"mcpas:{i}"
        species = normalize_species(_na(r.get("Species")) or "")
        epitope = _na(r.get("Epitope.peptide"))
        antigen = _na(r.get("Antigen.protein"))
        antigen_organism = _na(r.get("Pathology"))
        mhc_a = _na(r.get("MHC"))
        pmid = _looks_like_pmid(r.get("PubMed.ID"))

        alpha_cdr3 = _na(r.get("CDR3.alpha.aa"))
        if alpha_cdr3:
            rec = _blank_record()
            rec.update(
                {
                    "source": "mcpas",
                    "source_record_id": f"mcpas:{i}",
                    "pairing_key": pairing_key,
                    "species": species,
                    "chain": "alpha",
                    "cdr3_aa": alpha_cdr3,
                    "cdr3_nt": _na(r.get("CDR3.alpha.nt")),
                    "v_gene": _normalize_gene_or_none(r.get("TRAV")),
                    "j_gene": _normalize_gene_or_none(r.get("TRAJ")),
                    "epitope_aa": epitope,
                    "antigen": antigen,
                    "antigen_organism": antigen_organism,
                    "mhc_a": mhc_a,
                    "reference_pmid": pmid,
                    "external_url": "https://friedmanlab.weizmann.ac.il/McPAS-TCR/",
                }
            )
            rows.append(rec)

        beta_cdr3 = _na(r.get("CDR3.beta.aa"))
        if beta_cdr3:
            rec = _blank_record()
            rec.update(
                {
                    "source": "mcpas",
                    "source_record_id": f"mcpas:{i}",
                    "pairing_key": pairing_key,
                    "species": species,
                    "chain": "beta",
                    "cdr3_aa": beta_cdr3,
                    "cdr3_nt": _na(r.get("CDR3.beta.nt")),
                    "v_gene": _normalize_gene_or_none(r.get("TRBV")),
                    "d_gene": _normalize_gene_or_none(r.get("TRBD")),
                    "j_gene": _normalize_gene_or_none(r.get("TRBJ")),
                    "epitope_aa": epitope,
                    "antigen": antigen,
                    "antigen_organism": antigen_organism,
                    "mhc_a": mhc_a,
                    "reference_pmid": pmid,
                    "external_url": "https://friedmanlab.weizmann.ac.il/McPAS-TCR/",
                }
            )
            rows.append(rec)
    return _finalize(rows)


# ---------------------------------------------------------------------------
# TCR3d
# ---------------------------------------------------------------------------
def _tcr3d_chain_lookup(chain_tsv: str) -> dict:
    """pdb_id -> {"alpha": (cdr1, cdr2), "beta": (cdr1, cdr2)} best effort fill-in."""
    try:
        raw = pd.read_csv(chain_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}
    lookup: dict = {}
    for _, r in raw.iterrows():
        pdb_id = _na(r.get("pdb_id"))
        tcr_type = (_na(r.get("tcr_type")) or "").lower()
        if not pdb_id or tcr_type not in ("alpha", "beta"):
            continue
        cdr1 = _na(r.get("cdr1_sequences"))
        cdr2 = _na(r.get("cdr2_sequences"))
        lookup.setdefault(pdb_id, {})[tcr_type] = (cdr1, cdr2)
    return lookup


def harmonize_tcr3d(complexes_tsv: str, chain_tsv: str) -> pd.DataFrame:
    raw = pd.read_csv(complexes_tsv, sep="\t", dtype=str, keep_default_na=False)
    chain_lookup = _tcr3d_chain_lookup(chain_tsv)
    rows: list[dict] = []
    for _, r in raw.iterrows():
        pdb_id = _na(r.get("PDB_ID"))
        if not pdb_id:
            continue
        pairing_key = f"tcr3d:{pdb_id}"
        species = normalize_species(_na(r.get("TCR_organism")) or "")
        epitope = _na(r.get("Epitope"))
        mhc_a = _na(r.get("MHC_allele"))
        pmid = _looks_like_pmid(r.get("Pubmed"))
        fallback = chain_lookup.get(pdb_id, {})

        alpha_cdr3 = _na(r.get("CDR3_alpha"))
        if alpha_cdr3:
            cdr1 = _na(r.get("CDR1_alpha")) or (fallback.get("alpha", (None, None))[0])
            cdr2 = _na(r.get("CDR2_alpha")) or (fallback.get("alpha", (None, None))[1])
            rec = _blank_record()
            rec.update(
                {
                    "source": "tcr3d",
                    "source_record_id": f"tcr3d:{pdb_id}:alpha",
                    "pairing_key": pairing_key,
                    "species": species,
                    "chain": "alpha",
                    "cdr3_aa": alpha_cdr3,
                    "v_gene": _normalize_gene_or_none(r.get("TRAV_gene")),
                    "j_gene": _normalize_gene_or_none(r.get("TRAJ_gene")),
                    "cdr1_aa": cdr1,
                    "cdr2_aa": cdr2,
                    "epitope_aa": epitope,
                    "mhc_a": mhc_a,
                    "pdb_id": pdb_id,
                    "reference_pmid": pmid,
                    "external_url": f"https://www.rcsb.org/structure/{pdb_id}",
                }
            )
            rows.append(rec)

        beta_cdr3 = _na(r.get("CDR3_beta"))
        if beta_cdr3:
            cdr1 = _na(r.get("CDR1_beta")) or (fallback.get("beta", (None, None))[0])
            cdr2 = _na(r.get("CDR2_beta")) or (fallback.get("beta", (None, None))[1])
            rec = _blank_record()
            rec.update(
                {
                    "source": "tcr3d",
                    "source_record_id": f"tcr3d:{pdb_id}:beta",
                    "pairing_key": pairing_key,
                    "species": species,
                    "chain": "beta",
                    "cdr3_aa": beta_cdr3,
                    "v_gene": _normalize_gene_or_none(r.get("TRBV_gene")),
                    "j_gene": _normalize_gene_or_none(r.get("TRBJ_gene")),
                    "cdr1_aa": cdr1,
                    "cdr2_aa": cdr2,
                    "epitope_aa": epitope,
                    "mhc_a": mhc_a,
                    "pdb_id": pdb_id,
                    "reference_pmid": pmid,
                    "external_url": f"https://www.rcsb.org/structure/{pdb_id}",
                }
            )
            rows.append(rec)
    return _finalize(rows)


# ---------------------------------------------------------------------------
# build_index
# ---------------------------------------------------------------------------
_DEDUP_KEY = ["source", "chain", "cdr3_aa", "v_gene", "j_gene", "epitope_aa", "mhc_a"]


def _resolve_vdjdb_path(raw_dir: Path) -> Optional[str]:
    zips = list(raw_dir.glob("vdjdb-*.zip"))
    if zips:
        zpath = zips[0]
        with zipfile.ZipFile(zpath) as zf:
            members = [n for n in zf.namelist() if n.endswith("vdjdb.slim.txt")]
            if members:
                extracted = raw_dir / "_vdjdb_slim_extracted.txt"
                extracted.write_bytes(zf.read(members[0]))
                return str(extracted)
    flat = raw_dir / "vdjdb_slim.txt"
    if flat.exists():
        return str(flat)
    return None


def _resolve_iedb_path(raw_dir: Path) -> Optional[str]:
    zips = list(raw_dir.glob("iedb_receptor_full_v3.zip"))
    if zips:
        with zipfile.ZipFile(zips[0]) as zf:
            if "tcr_full_v3.csv" in zf.namelist():
                extracted = raw_dir / "_iedb_tcr_extracted.csv"
                extracted.write_bytes(zf.read("tcr_full_v3.csv"))
                return str(extracted)
    flat = raw_dir / "iedb_tcr.csv"
    if flat.exists():
        return str(flat)
    return None


def build_index(raw_dir: str, out_parquet: str, out_meta: str) -> dict:
    raw_path = Path(raw_dir)
    frames: list[pd.DataFrame] = []
    per_source: dict[str, int] = {}

    vdjdb_path = _resolve_vdjdb_path(raw_path)
    if vdjdb_path:
        df = harmonize_vdjdb(vdjdb_path)
        frames.append(df)
        per_source["vdjdb"] = len(df)

    iedb_path = _resolve_iedb_path(raw_path)
    if iedb_path:
        df = harmonize_iedb(iedb_path)
        frames.append(df)
        per_source["iedb"] = len(df)

    mcpas_path = raw_path / "mcpas.csv"
    if mcpas_path.exists():
        df = harmonize_mcpas(str(mcpas_path))
        frames.append(df)
        per_source["mcpas"] = len(df)

    complexes_path = raw_path / "tcr3d_complexes.tsv"
    chain_path = raw_path / "tcr3d_chain.tsv"
    if complexes_path.exists():
        df = harmonize_tcr3d(str(complexes_path), str(chain_path))
        frames.append(df)
        per_source["tcr3d"] = len(df)

    combined = pd.concat(frames, ignore_index=True) if frames else _empty_frame()
    combined = combined[SCHEMA_COLUMNS]
    combined = combined.drop_duplicates(subset=_DEDUP_KEY).reset_index(drop=True)

    # per_source counts reflect rows actually kept after dedup, not pre-dedup harmonized counts.
    if not combined.empty:
        per_source = combined["source"].value_counts().to_dict()
    else:
        per_source = {}

    out_parquet_path = Path(out_parquet)
    out_parquet_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out_parquet_path, index=False)

    meta = {
        "built_columns": SCHEMA_COLUMNS,
        "rows_total": len(combined),
        "per_source": per_source,
        "snapshots": {
            "vdjdb": "vdjdb-2026-06-03",
            "iedb": "receptor_full_v3",
            "mcpas": "mcpas.csv",
            "tcr3d": "tcr3d_complexes.tsv",
        },
    }
    out_meta_path = Path(out_meta)
    out_meta_path.parent.mkdir(parents=True, exist_ok=True)
    out_meta_path.write_text(json.dumps(meta, indent=2))

    return meta
