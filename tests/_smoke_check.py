"""Quick smoke check for parse_vdjdb_tsv — run with: python3 tests/_smoke_check.py"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

import requests

from imgt_app.file_ingest import parse_vdjdb_tsv

TSV = (
    "cdr3\tv_segm\tj_segm\tspecies\tmhc_a\tmhc_b\tmhc_class"
    "\tantigen_epitope\tantigen_gene\tantigen_species\tscore\n"
    "CASSIRSSYEQYF\tTRBV19\tTRBJ2-7\tHomoSapiens\tHLA-A*02:01\tB2M\tI"
    "\tGILGFVFTL\tM1\tInfluenzaA\t3\n"
    "CASGDSSYEQYF\tTRBV13-2\tTRBJ2-7\tMusMusculus\tH-2Kb\tB2M\tI"
    "\tSIINFEKL\tOVA\tHomoSapiens\t3\n"
    "CAVLDSNYQLIW\tTRAV41\tTRAJ33\tUnknownSp\tHLA-A*02:01\tB2M\tI"
    "\tNLVPMVATV\tpp65\tCMV\t2\n"
)

records = parse_vdjdb_tsv(TSV.encode())
assert len(records) == 3, f"Expected 3 records, got {len(records)}"
assert records[0].species == "human", records[0].species
assert records[1].species == "mouse", records[1].species
assert records[2].species == "other", records[2].species
assert records[0].antigen_epitope == "GILGFVFTL"
assert records[0].gene_name == "TRBV19"
assert records[0].region == "CDR3"
assert records[0].source == "vdjdb"
assert records[0].sequence == "CASSIRSSYEQYF"
print("All smoke checks passed.")

# BATMAN server health
try:
    r = requests.get("http://localhost:8105/health", timeout=3)
    print(f"batman: {r.json()}")
except Exception as e:
    print(f"batman: OFFLINE ({e})")

# Test /predict/activation
try:
    r = requests.post("http://localhost:8000/predict/activation", json={
        "tcr_cdr3": "CASSIRSSYEQYF",
        "tcr_v_gene": "TRBV19",
        "hla_allele": "HLA-A*02:01",
        "candidate_peptides": ["GILGFVFTL", "NLVPMVATV"],
    }, timeout=30)
    data = r.json()
    print(f"/predict/activation: {len(data['results'])} peptides scored")
    for row in data["results"]:
        print(f"  {row['peptide']}: composite={row['composite_score']}")
except Exception as e:
    print(f"/predict/activation: ERROR ({e})")

# TEMPO server health
try:
    r = requests.get("http://localhost:8106/health", timeout=3)
    print(f"tempo: {r.json()}")
except Exception as e:
    print(f"tempo: OFFLINE ({e})")

# TEMPO /tempo/score
try:
    r = requests.post("http://localhost:8106/tempo/score", json={
        "v_gene": "TRAV12-2",
        "j_gene": "TRAJ30",
        "cdr3": "CAVGDDKIIF",
        "chain": "alpha",
        "species": "human",
    }, timeout=10)
    data = r.json()
    print(f"/tempo/score: log_likelihood={data.get('log_likelihood')}")
except Exception as e:
    print(f"/tempo/score: ERROR ({e})")

# Cross-reactivity prediction via main API
try:
    r = requests.post("http://localhost:8000/predict/crossreactivity", json={
        "reference_peptide": "LLWNGPMAV",
        "variant_peptides": ["ALWNGPMAV"],
        "mhc_allele": "HLA-A*02:01",
        "mhc_class": "I",
    }, timeout=10)
    data = r.json()
    print(f"/predict/crossreactivity: {len(data.get('results', []))} variants scored")
except Exception as e:
    print(f"/predict/crossreactivity: ERROR ({e})")
