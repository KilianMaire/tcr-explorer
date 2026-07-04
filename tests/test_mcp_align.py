from tcr_explorer import mcp_server

def test_align_tool_returns_msaresult():
    out = mcp_server.align_tcr_genes(sequences=[{"name":"a","seq":"CASS"},{"name":"b","seq":"CASF"}], seq_type="aa")
    assert out["n_sequences"] == 2 and "records" in out and "engine" in out
