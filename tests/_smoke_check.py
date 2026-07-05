"""Quick smoke check for parse_vdjdb_tsv — run with: python3 tests/_smoke_check.py"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "src"))

from tcr_explorer.file_ingest import parse_vdjdb_tsv

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
