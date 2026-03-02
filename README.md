# IMGT Search Engine

A four-tier distributed search API for immune receptor and epitope databases:
**HLA** (IMGT/HLA via EBI), **TCR** (NCBI Entrez), **VDJdb**, and **IEDB**.

Features CDR1/CDR2 enrichment via stitchr IMGT germline data, a local SQLite search
index, bulk VDJdb import, natural-language query parsing, and a Streamlit UI.

---

## Requirements

- Python 3.11+
- Internet access (EBI, NCBI, VDJdb, and IEDB public APIs)

---

## Install

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# One-time download of stitchr IMGT germline data (required for CDR1/CDR2 prediction)
stitchrdl -s human
```

---

## Start the services

Each service runs independently. Open five terminals from the project root:

| Terminal | Command | Port |
|----------|---------|------|
| 1 | `uvicorn servers.hla_server:app --port 8101 --reload` | 8101 |
| 2 | `uvicorn servers.tcr_server:app --port 8102 --reload` | 8102 |
| 3 | `uvicorn servers.vdjdb_server:app --port 8103 --reload` | 8103 |
| 4 | `uvicorn servers.iedb_server:app --port 8104 --reload` | 8104 |
| 5 | `PYTHONPATH=src uvicorn imgt_app.api:app --port 8000 --reload` | 8000 |

Start the tool servers (8101–8104) before the main API (8000).

### Health check

```bash
curl http://localhost:8000/health
# {"status":"ok"}
```

---

## `/search` endpoint

**POST** `http://localhost:8000/search`

All fields are optional — omit any you don't need.

```json
{
  "source": "hla | tcr | vdjdb | iedb",
  "species": "human | mouse | other",
  "gene_name": "HLA-A",
  "region": "exon2",
  "sequence_contains": "GILG",
  "antigen_epitope": "GILGFVFTL",
  "limit": 50
}
```

### Examples

**HLA allele sequences (EBI IMGT/HLA)**
```bash
curl -s -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{"source":"hla","species":"human","gene_name":"HLA-A","limit":5}'
```

**TCR V-gene sequences (NCBI Entrez)**
```bash
curl -s -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{"source":"tcr","species":"human","gene_name":"TRBV12-3","limit":5}'
```

**VDJdb CDR3 sequences filtered by antigen epitope**
```bash
curl -s -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{"source":"vdjdb","antigen_epitope":"GILGFVFTL","limit":10}'
```

**IEDB T-cell assay data**
```bash
curl -s -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{"source":"iedb","gene_name":"HLA-A*02:01","sequence_contains":"GILG","limit":10}'
```

### Table output (CSV or Markdown)

```bash
# CSV
curl -s -X POST "http://localhost:8000/search/table?fmt=csv" \
  -H "Content-Type: application/json" \
  -d '{"source":"vdjdb","limit":10}'

# Markdown
curl -s -X POST "http://localhost:8000/search/table?fmt=md" \
  -H "Content-Type: application/json" \
  -d '{"source":"hla","species":"human","limit":5}'
```

---

## CDR1 / CDR2 prediction

**GET** `http://localhost:8000/predict/cdr`

Returns germline CDR1 and CDR2 amino acid (and nucleotide) sequences for a TCR V gene
using stitchr IMGT data. Requires `stitchrdl -s human` to have been run.

```bash
curl "http://localhost:8000/predict/cdr?v_gene=TRBV12-3&species=human"
```

Response:
```json
{
  "v_gene": "TRBV12-3",
  "species": "human",
  "allele": "TRBV12-3*01",
  "cdr1_aa": "SGHDN",
  "cdr2_aa": "FNNNVP",
  "cdr1_nt": "...",
  "cdr2_nt": "..."
}
```

---

## Bulk VDJdb ingest

Upload a VDJdb TSV export to populate the local SQLite index:

```bash
curl -s -X POST http://localhost:8000/ingest/vdjdb \
  -F "file=@/path/to/vdjdb_full.tsv"
```

Download VDJdb data from <https://vdjdb.cdr3.net> → Browse → Export TSV.

After ingesting, search the local index without hitting the upstream API:

```bash
curl -s -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{"source":"vdjdb","antigen_epitope":"NLVPMVATV","limit":20}'
```

---

## Natural-language query

Requires an OpenAI-compatible LLM running locally (e.g. LM Studio). Falls back to
a heuristic keyword parser when no LLM is available.

```bash
curl -s -X POST http://localhost:8000/query/nl \
  -H "Content-Type: application/json" \
  -d '{"query":"find CDR3 sequences for TRBV12 in human against influenza", "limit":10}'
```

---

## Streamlit UI

```bash
streamlit run src/imgt_app/frontend.py
```

Opens at <http://localhost:8501> with two pages:

- **Search** — filter by source/species/gene/CDR3/epitope, results table with CDR3 + CDR1 + CDR2 + antigen_epitope columns, CSV download
- **CDR Predict** — look up germline CDR1/CDR2 for any TCR V gene

---

## Environment variables

All are optional; defaults work for local development.

| Variable | Default | Description |
|----------|---------|-------------|
| `HLA_SERVER_URL` | `http://127.0.0.1:8101` | HLA tool server URL |
| `TCR_SERVER_URL` | `http://127.0.0.1:8102` | TCR tool server URL |
| `VDJDB_SERVER_URL` | `http://127.0.0.1:8103` | VDJdb tool server URL |
| `IEDB_SERVER_URL` | `http://127.0.0.1:8104` | IEDB tool server URL |
| `IMGT_DB_PATH` | `./imgt.db` | SQLite database file path |
| `NCBI_API_KEY` | *(empty)* | NCBI API key — raises rate limit to 10 req/s |
| `VDJDB_DATA_FILE` | `examples/vdjdb_seed.csv` | Local VDJdb seed/export file |
| `VDJDB_API_URL` | VDJdb public API | Set to `""` to disable upstream VDJdb calls |
| `IEDB_API_URL` | IEDB Query API | Override IEDB endpoint |
| `IMGT_API_URL` | `http://127.0.0.1:8000` | Main API URL (used by Streamlit frontend) |
| `LMSTUDIO_BASE_URL` | `http://127.0.0.1:1234/v1` | LM Studio OpenAI-compatible endpoint |
| `LMSTUDIO_MODEL` | `local-model` | LM Studio model ID |

---

## Run tests

```bash
PYTHONPATH=src pytest tests/ -v
```

---

## Architecture

```
                        ┌─────────────────────┐
                        │  Main API  :8000     │
                        │  src/imgt_app/api.py │
                        └──────────┬──────────┘
              ┌───────────┬────────┴────────┬───────────┐
              ▼           ▼                 ▼           ▼
        HLA :8101   TCR :8102        VDJdb :8103  IEDB :8104
        (EBI API)  (NCBI Entrez)   (seed CSV +   (IEDB Query
                                    VDJdb API +    API proxy)
                                    CDR enricher)
```

The main API merges results from the local SQLite index and the matching tool server.
Tool servers are independently startable — they have no shared state or imports from
the main API package.
