"""Streamlit frontend for the IMGT Search Engine.

Run with:
    streamlit run src/imgt_app/frontend.py

Calls the main API at http://127.0.0.1:8000 (override with IMGT_API_URL env var).
"""
from __future__ import annotations

import os
from typing import Any

import httpx
import pandas as pd
import streamlit as st

API_BASE = os.getenv("IMGT_API_URL", "http://127.0.0.1:8000")


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _search(payload: dict[str, Any]) -> dict[str, Any]:
    with httpx.Client(timeout=30.0) as client:
        resp = client.post(f"{API_BASE}/search", json=payload)
        resp.raise_for_status()
        return resp.json()


def _predict_cdr(v_gene: str, species: str) -> dict[str, Any]:
    with httpx.Client(timeout=15.0) as client:
        resp = client.get(
            f"{API_BASE}/predict/cdr",
            params={"v_gene": v_gene, "species": species},
        )
        resp.raise_for_status()
        return resp.json()


# ---------------------------------------------------------------------------
# Data shaping
# ---------------------------------------------------------------------------

def _records_to_df(records: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for rec in records:
        meta: dict[str, Any] = rec.get("metadata") or {}
        rows.append(
            {
                "source": rec.get("source", ""),
                "species": rec.get("species", ""),
                "gene_name": rec.get("gene_name", ""),
                "allele_name": rec.get("allele_name") or "",
                "region": rec.get("region") or "",
                "CDR3": rec.get("sequence", ""),
                "CDR1": meta.get("cdr1_aa") or "",
                "CDR2": meta.get("cdr2_aa") or "",
                "antigen_epitope": (
                    rec.get("antigen_epitope")
                    or meta.get("antigen_epitope")
                    or ""
                ),
                "mhc_a": meta.get("mhc_a") or "",
                "mhc_class": meta.get("mhc_class") or "",
                "j_segm": meta.get("j_segm") or "",
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Page: Search
# ---------------------------------------------------------------------------

def _page_search() -> None:
    st.header("Search")

    # Sidebar filters
    st.sidebar.header("Filters")
    source = st.sidebar.selectbox(
        "Source",
        options=["", "hla", "tcr", "vdjdb", "iedb"],
        format_func=lambda v: v or "(any)",
    )
    species = st.sidebar.selectbox(
        "Species",
        options=["", "human", "mouse", "other"],
        format_func=lambda v: v or "(any)",
    )
    gene_name = st.sidebar.text_input("Gene / V-gene name", placeholder="e.g. TRBV12-3")
    sequence_contains = st.sidebar.text_input(
        "CDR3 / sequence contains", placeholder="e.g. CASSQ"
    )
    antigen_epitope = st.sidebar.text_input(
        "Antigen epitope contains", placeholder="e.g. GILG"
    )
    limit = st.sidebar.slider("Max results", min_value=5, max_value=200, value=50, step=5)
    search_clicked = st.sidebar.button("Search", type="primary")

    if not search_clicked:
        st.info("Set filters in the sidebar and click **Search** to begin.")
        return

    payload: dict[str, Any] = {"limit": limit}
    if source:
        payload["source"] = source
    if species:
        payload["species"] = species
    if gene_name:
        payload["gene_name"] = gene_name.strip()
    if sequence_contains:
        payload["sequence_contains"] = sequence_contains.strip()
    if antigen_epitope:
        payload["antigen_epitope"] = antigen_epitope.strip()

    with st.spinner("Querying API…"):
        try:
            data = _search(payload)
        except httpx.HTTPStatusError as exc:
            st.error(f"API error {exc.response.status_code}: {exc.response.text}")
            return
        except Exception as exc:
            st.error(f"Could not reach API at {API_BASE} — {exc}")
            return

    records: list[dict[str, Any]] = data.get("records", [])
    total: int = data.get("total", len(records))

    st.success(f"**{total}** record(s) found — showing up to {limit}")

    if not records:
        st.info("No records matched your query.")
        return

    df = _records_to_df(records)
    st.dataframe(df, use_container_width=True)

    csv_bytes = df.to_csv(index=False).encode()
    st.download_button(
        label="Download CSV",
        data=csv_bytes,
        file_name="imgt_results.csv",
        mime="text/csv",
    )


# ---------------------------------------------------------------------------
# Page: CDR Predict
# ---------------------------------------------------------------------------

def _page_cdr_predict() -> None:
    st.header("CDR1 / CDR2 Prediction")
    st.caption(
        "Predict CDR1 and CDR2 amino acid sequences for a TCR V gene "
        "using stitchr IMGT germline data."
    )

    col1, col2 = st.columns(2)
    with col1:
        v_gene = st.text_input("V-gene name", placeholder="e.g. TRBV12-3")
    with col2:
        species = st.selectbox("Species", options=["human", "mouse", "other"])

    if st.button("Predict", type="primary") and v_gene:
        with st.spinner("Predicting CDRs…"):
            try:
                result = _predict_cdr(v_gene.strip(), species)
            except httpx.HTTPStatusError as exc:
                st.error(f"API error {exc.response.status_code}: {exc.response.text}")
                return
            except Exception as exc:
                st.error(f"Could not reach API at {API_BASE} — {exc}")
                return

        st.subheader("Result")
        cols = st.columns(4)
        cols[0].metric("V-gene", result.get("v_gene", ""))
        cols[1].metric("Allele", result.get("allele") or "—")
        cols[2].metric("CDR1 (aa)", result.get("cdr1_aa") or "—")
        cols[3].metric("CDR2 (aa)", result.get("cdr2_aa") or "—")

        if result.get("cdr1_nt") or result.get("cdr2_nt"):
            st.subheader("Nucleotide sequences")
            nt_cols = st.columns(2)
            nt_cols[0].text_area("CDR1 (nt)", result.get("cdr1_nt") or "", height=80)
            nt_cols[1].text_area("CDR2 (nt)", result.get("cdr2_nt") or "", height=80)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(
        page_title="IMGT Search Engine",
        page_icon=":dna:",
        layout="wide",
    )
    st.title("IMGT Search Engine")
    st.caption(
        "Search HLA (IMGT/HLA via EBI), TCR (NCBI Entrez), VDJdb, and IEDB "
        "immune receptor / epitope databases."
    )
    st.caption(f"API: `{API_BASE}`")

    page = st.sidebar.radio(
        "Page",
        options=["Search", "CDR Predict"],
        label_visibility="collapsed",
    )

    if page == "Search":
        _page_search()
    else:
        _page_cdr_predict()


if __name__ == "__main__":
    main()
