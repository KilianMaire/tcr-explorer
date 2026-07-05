from __future__ import annotations

import dataclasses
import logging
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, PlainTextResponse

from .cdr_enricher import get_cdr1_cdr2
from .config import settings
from .dossier_models import (
    AlignRequest,
    AskRequest,
    AskResponse,
    DossierRequest,
    MSAResult,
    RecordsRequest,
    RecordsResponse,
    SimilarRequest,
    SimilarResponse,
    PairedSimilarRequest,
    PairedSimilarResponse,
    TCRDossier,
)
from .mcp_clients import ToolServerClient
from .models import (
    AssignRequest, AssignResponse,
    CDRPredictResponse, QueryRequest, QueryResponse, QueryUnderstanding, QueryBlock,
    ReconstructRequest, ReconstructResponse,
    SearchRequest, SearchResponse, Species,
)
from .reconstructor import reconstruct_tcr

_logger = logging.getLogger(__name__)


app = FastAPI(title="TCR Explorer", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
hla_client = ToolServerClient(settings.hla_server_url)
mhc_client = ToolServerClient(settings.mhc_server_url)


def _coerce_search_response(value: "SearchResponse | dict") -> SearchResponse:
    """Accept either a SearchResponse or a raw dict (as returned by mocks in tests)."""
    if isinstance(value, dict):
        return SearchResponse(**value)
    return value


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/predict/cdr", response_model=CDRPredictResponse)
def predict_cdr(v_gene: str, species: Species = "human") -> CDRPredictResponse:
    """Return CDR1 and CDR2 amino acid sequences for a TCR V gene using stitchr IMGT data."""
    result = get_cdr1_cdr2(v_gene, species)
    return CDRPredictResponse(v_gene=v_gene, species=species, **result)


@app.post("/search", response_model=SearchResponse)
async def search(req: SearchRequest) -> SearchResponse:
    """Look up MHC allele sequences from the live hla/mhc proxies.

    TCR record search moved to POST /v1/tcr/records, which serves the vendored
    snapshot. Only source hla and source mhc are handled here.
    """
    if req.source == "hla":
        return _coerce_search_response(await hla_client.search(req))
    if req.source == "mhc":
        return _coerce_search_response(await mhc_client.search(req))
    raise HTTPException(
        status_code=400,
        detail="source must be hla or mhc; use POST /v1/tcr/records for TCR records",
    )


@app.post("/reconstruct", response_model=ReconstructResponse)
def reconstruct(req: ReconstructRequest) -> ReconstructResponse:
    """
    Reconstruct a full TCR coding sequence from a CDR3, optionally with V and J.

    V and J are optional: if either is omitted, it is inferred from the CDR3 by
    tallying database records that carry the exact same CDR3 and taking the most
    common pairing (reported as inferred, with the supporting count and
    alternatives). Uses stitchr IMGT germline for V/J regions; CDR3 is
    back-translated with human-optimised codons. Assembly follows IMGT/VDJdb
    CDR3 boundaries: Cys104 … Phe/Trp118 (inclusive).
    """
    from .records import infer_vj_from_cdr3, _gene_base  # local: avoids import cycle

    v_gene = (req.v_gene or "").strip() or None
    j_gene = (req.j_gene or "").strip() or None
    genes_inferred = False
    inference_support: Optional[int] = None
    inference_alternatives: Optional[list[str]] = None

    if not v_gene or not j_gene:
        candidates = infer_vj_from_cdr3(req.cdr3_aa, req.species)
        # constrain to whatever the caller did pin down
        if v_gene:
            candidates = [c for c in candidates if _gene_base(c["v_gene"]) == _gene_base(v_gene)]
        if j_gene:
            candidates = [c for c in candidates if _gene_base(c["j_gene"]) == _gene_base(j_gene)]
        if not candidates:
            return ReconstructResponse(
                v_gene=v_gene or "", j_gene=j_gene or "", cdr3_aa=req.cdr3_aa,
                species=req.species, full_nt=None, full_aa=None, v_region_nt=None,
                cdr3_nt="", j_region_nt=None, v_found=bool(v_gene), j_found=bool(j_gene),
                genes_inferred=False,
                note=(
                    "No database record carries this exact CDR3, so V and J could "
                    "not be inferred. Provide a V gene and a J gene to reconstruct."
                ),
            )
        top = candidates[0]
        v_gene = v_gene or top["v_gene"]
        j_gene = j_gene or top["j_gene"]
        genes_inferred = True
        inference_support = top["count"]
        inference_alternatives = [
            f"{c['v_gene']}/{c['j_gene']} (n={c['count']})" for c in candidates
        ]

    result = reconstruct_tcr(v_gene, j_gene, req.cdr3_aa, req.species)
    result["genes_inferred"] = genes_inferred
    result["inference_support"] = inference_support
    result["inference_alternatives"] = inference_alternatives
    return ReconstructResponse(**result)


def _dossier_markdown(d: TCRDossier) -> str:
    lines = [f"# TCR Dossier ({d.status})", "", d.summary, "",
             f"- chain: {d.chain}", f"- species: {d.species}"]
    if d.genes.get("v"):
        lines.append(f"- V: {d.genes['v'].call} ({d.genes['v'].score_method})")
    if d.known_epitopes:
        lines.append(f"- known epitopes: {d.known_epitopes_total}")
    if d.warnings:
        lines.append("- warnings: " + ", ".join(w.code for w in d.warnings))
    return "\n".join(lines)


@app.post("/v1/tcr/dossier", response_model=TCRDossier)
def tcr_dossier(req: DossierRequest, request: Request):
    # Synchronous by design: build_dossier is fully synchronous and its epitope
    # lookup (dossier_epitopes._run_search) must create its own event loop, which
    # it can only do off the request's running loop. FastAPI runs a sync route in
    # a threadpool, so known_epitopes actually surface (an async route left them
    # permanently empty).
    from .dossier import build_dossier  # local import: avoids circular import
    # (dossier -> dossier_epitopes -> api.search/_IEDB_HITS_CAP)

    d = build_dossier(req)
    if "text/markdown" in request.headers.get("accept", ""):
        return PlainTextResponse(_dossier_markdown(d), media_type="text/markdown")
    return d


@app.post("/v1/tcr/similar", response_model=SimilarResponse)
def tcr_similar(req: SimilarRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: the engine
    # reads a parquet index off disk, so it runs in FastAPI's threadpool rather
    # than blocking the event loop.
    from .similarity import find_similar_tcrs  # local import: avoids import cycles

    neigh, engine, total, warnings = find_similar_tcrs(
        req.cdr3,
        req.v_gene,
        req.j_gene,
        species=req.species,
        top_k=req.top_k,
        min_similarity=req.min_similarity,
    )
    return SimilarResponse(neighbours=neigh, engine=engine, total_candidates=total, warnings=warnings)


@app.post("/v1/tcr/similar_paired", response_model=PairedSimilarResponse)
def tcr_similar_paired(req: PairedSimilarRequest):
    # Synchronous by design, same rationale as /v1/tcr/similar above. Paired
    # scoring is tcrdist only (the paired distance is alpha tcrdist plus beta
    # tcrdist); without the tcrdist extra it returns no neighbours and says so.
    from .similarity import find_similar_paired_tcrs  # local import: avoids import cycles

    neigh, engine, total, warnings = find_similar_paired_tcrs(
        req.cdr3_a,
        req.v_a,
        req.cdr3_b,
        req.v_b,
        species=req.species,
        top_k=req.top_k,
        min_similarity=req.min_similarity,
    )
    return PairedSimilarResponse(neighbours=neigh, engine=engine, total_candidates=total, warnings=warnings)


@app.post("/v1/tcr/records", response_model=RecordsResponse)
def tcr_records(req: RecordsRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: retrieval
    # reads a parquet index off disk via sync helpers that must not run inside
    # a running event loop, so it runs in FastAPI's threadpool.
    from .records import retrieve_records  # local import: avoids import cycles

    return retrieve_records(req)


@app.post("/v1/tcr/align", response_model=MSAResult)
def tcr_align(req: AlignRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: align may
    # resolve germline sets from disk or shell out to clustalo, so it runs in
    # FastAPI's threadpool rather than blocking the event loop.
    from .msa import align  # local import: avoids import cycles

    return align(req)


@app.post("/v1/tcr/assign", response_model=AssignResponse)
def tcr_assign(req: AssignRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: assign may
    # resolve germline sets from disk, so it runs in FastAPI's threadpool rather
    # than blocking the event loop.
    from . import tcr_align  # local import: avoids import cycles

    result = tcr_align.assign(
        req.sequence, species=req.species, chain=req.chain, want_d=req.want_d
    )
    return AssignResponse(**dataclasses.asdict(result))


@app.post("/v1/tcr/query", response_model=QueryResponse)
def tcr_query(req: QueryRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier above: routed
    # tools may hit disk-backed lookups, so it runs in FastAPI's threadpool
    # rather than blocking the event loop.
    from .query_router import route_query  # local import: avoids import cycles

    result = route_query(req.query, species=req.species, force=req.force)
    understood = QueryUnderstanding(
        input=result.input,
        detected_type=result.detected_type,
        species=result.species,
        tools=result.tools,
        note=result.note,
    )
    blocks = [QueryBlock(tool=b.tool, title=b.title, data=b.data) for b in result.blocks]
    return QueryResponse(understood=understood, blocks=blocks, warnings=result.warnings)


@app.post("/v1/tcr/ask", response_model=AskResponse)
def tcr_ask(req: AskRequest):
    # Synchronous by design, same rationale as /v1/tcr/dossier and /v1/tcr/similar
    # above: the routed intent may hit build_dossier/find_similar_tcrs, both of
    # which are themselves synchronous for the same reasons.
    from .ask import answer  # local import: avoids circular import

    return answer(req)


_UI_HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>TCR Explorer</title><meta name="viewport" content="width=device-width, initial-scale=1">
<style>
 body{font-family:system-ui,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;color:#111}
 h1{font-size:1.3rem;text-align:center} input,select,button{font-size:1rem;padding:.4rem}
 #q{width:100%;box-sizing:border-box;margin-bottom:.5rem;font-size:1.1rem}
 .hero{text-align:center;margin-bottom:1.2rem}
 .guidance{color:#555;max-width:640px;margin:0 auto}
 .searchbar{background:#f7f7fb;border:1px solid #e2e2e8;border-radius:10px;padding:1rem}
 .searchrow{display:flex;flex-wrap:wrap;gap:.6rem;align-items:center;justify-content:center}
 .card{border:1px solid #ddd;border-radius:8px;padding:1rem;margin:1rem 0}
 .warn{color:#a15c00} .syn{color:#7a3e00;font-style:italic} table{border-collapse:collapse;width:100%}
 .loading{color:#0b5;font-weight:600} .loading::after{content:'';display:inline-block;width:.7em;height:.7em;margin-left:.4em;border:2px solid #0b5;border-top-color:transparent;border-radius:50%;animation:spin .7s linear infinite;vertical-align:middle} @keyframes spin{to{transform:rotate(360deg)}} button:disabled{opacity:.6}
 td,th{border:1px solid #eee;padding:.3rem;text-align:left;font-size:.9rem} .muted{color:#666}
 h3{margin:.3rem 0}
 .echo{background:#eef7ee;border:1px solid #cfe8cf;border-radius:6px;padding:.4rem .7rem;margin:.6rem 0;font-size:.9rem}
 .section{margin-top:1rem}
 .rec{border:1px solid #e2e2e2;border-radius:6px;padding:.6rem .8rem;margin:.5rem 0;background:#fafafa}
 .rec.neigh{background:#f5f8ff;border-style:dashed}
 .badge{display:inline-block;background:#333;color:#fff;font-size:.75rem;font-weight:700;padding:.1rem .4rem;border-radius:4px;letter-spacing:.03em}
 .kind{color:#666;font-size:.8rem}
 .xspecies{display:inline-block;background:#ffe6cc;color:#8a4b00;font-size:.72rem;font-weight:600;padding:.05rem .4rem;border-radius:4px;margin-left:.4rem}
 .comp{font-family:monospace;font-size:.85rem;background:#eee;padding:.2rem .4rem;border-radius:4px;display:inline-block}
 .partner{margin-top:.4rem;padding-left:.6rem;border-left:3px solid #ccc}
 .db-inference{margin-top:.6rem;padding-top:.4rem;border-top:1px dashed #ccc} .db-inference h4{margin:.2rem 0;font-size:.95rem}
 .chiprow{display:flex;flex-wrap:wrap;gap:.5rem;justify-content:center;margin:.8rem 0}
 .chip{background:#fff;border:1px solid #bbb;border-radius:999px;padding:.3rem .8rem;font-size:.85rem;cursor:pointer}
 .chip:hover{background:#f0f0f0}
 .hiddenform{display:none;margin-top:1.5rem}
 .hiddenform h2{margin:.2rem 0 .4rem;font-size:1.25rem}
 #rc_seq{width:100%;box-sizing:border-box;margin-bottom:.5rem} #rc_v,#rc_j{flex:1;min-width:9rem}
 .onboard pre{white-space:pre-wrap;word-break:break-word;font-family:monospace;font-size:.82rem;background:#eee;padding:.5rem;border-radius:6px}
 .onboard h4{margin:.8rem 0 .2rem}
 .copy-btn{margin-top:.2rem}
</style></head><body>
<div class="hero">
<h1>TCR Explorer</h1>
<p class="guidance muted">Ask anything about a TCR in one box: a gene, a CDR3, a species plus a CDR3, a V+J+CDR3 phrase, a full chain, or a database id.</p>
</div>
<div class="searchbar">
<form id="f">
<input id="q" placeholder="e.g. mouse CASSGGTGEQYF, TRBV20-1, human CASSLGTEAFF TRBJ2-7, vdjdb:12345" value="TRBV20-1">
<div class="searchrow">
<select id="sp"><option value="">auto (detect from text)</option><option value="human">human</option><option value="mouse">mouse</option></select>
<button type="submit">Search</button>
</div>
</form>
</div>
<div class="chiprow" id="chips">
<button type="button" class="chip" data-tool="records">records</button>
<button type="button" class="chip" data-tool="assign">assign</button>
<button type="button" class="chip" data-tool="dossier">dossier</button>
<button type="button" class="chip" data-tool="similar">similar</button>
<button type="button" class="chip" data-tool="paired">paired</button>
<button type="button" class="chip" data-tool="reconstruct">reconstruct</button>
<button type="button" class="chip" data-tool="align">align</button>
</div>
<div id="out"></div>
<div class="searchbar hiddenform" id="rcform">
<h2>Reconstruct or assign a TCR</h2>
<p class="muted">Paste a CDR3, a V region, or a full chain (nucleotide or amino acid). It is aligned against IMGT germline and called at the allele level: V, J, D and constant calls with identity, aligned span and the full tie set; per region identity; the extracted CDR3; and, when enough sequence is present, a reconstructed full chain (mouse constant is oracle-validated, human is not). A bare CDR3 has no framework to call a V allele from, so it is refused and backed instead by a database frequency inference, clearly labeled as weaker. Fill in both V and J below to skip assignment and reconstruct directly from those explicit genes (allele defaults to *01 unless written into the gene name, e.g. TRAV7-4*02).</p>
<form id="rcf">
<input id="rc_seq" placeholder="paste a CDR3, a V region, or a full chain (nucleotide or amino acid), e.g. CASSLGTEAFF">
<div class="searchrow">
<select id="rc_sp"><option>human</option><option>mouse</option></select>
<input id="rc_v" placeholder="V gene (optional; set both V and J to reconstruct directly)">
<input id="rc_j" placeholder="J gene (optional; set both V and J to reconstruct directly)">
<button type="submit">Assign</button>
</div>
</form>
<div id="rc_out"></div>
</div>
<div class="searchbar hiddenform" id="pairform">
<h2>Find paired (alpha and beta) neighbours</h2>
<p class="muted">Give both chains of a receptor. TCR Explorer scores paired references (reconstructed from the index by pairing key) with the paired tcrdist, which is the sum of the alpha and beta single chain distances. These neighbours' epitopes are an inferred, weaker signal, not confirmed specificity. Paired scoring needs the tcrdist extra installed and both query V genes present in the reference table; otherwise it returns nothing with an explanatory warning.</p>
<form id="pf">
<div class="searchrow">
<input id="p_ca" placeholder="alpha CDR3 e.g. CAVNFGGGKLIF">
<input id="p_va" placeholder="alpha V gene e.g. TRAV12-1">
</div>
<div class="searchrow">
<input id="p_cb" placeholder="beta CDR3 e.g. CASSIRSSYEQYF">
<input id="p_vb" placeholder="beta V gene e.g. TRBV19">
</div>
<div class="searchrow">
<select id="p_sp"><option>human</option><option>mouse</option></select>
<button type="submit">Find paired</button>
</div>
</form>
<div id="p_out"></div>
</div>
<div class="searchbar hiddenform" id="alignform">
<h2>Align a gene set</h2>
<p class="muted">Align a germline set (species + chain + segment) or a gene list. V/J/C come from the germline source; D is not available there.</p>
<form id="af">
<select id="a_sp"><option>human</option><option>mouse</option></select>
<input id="a_chain" placeholder="chain e.g. TRB" value="TRB">
<select id="a_seg"><option>V</option><option>D</option><option>J</option><option>C</option></select>
<label><input type="checkbox" id="a_translate"> translate</label>
<button type="submit">Align</button></form>
<div id="a_out"></div>
</div>
<div class="card onboard" id="onboard">
<h3>Ask in plain English, use your own AI assistant</h3>
<p class="muted">TCR Explorer also runs as an MCP server, so your own AI assistant (any MCP compatible client) can call these tools directly while you keep chatting in plain English. Paste either artifact below into your assistant.</p>
<h4>MCP server config</h4>
<pre id="mcpConfigOut"></pre>
<button type="button" class="copy-btn" id="copyConfigBtn">Copy config</button>
<h4>Install prompt</h4>
<pre id="mcpPromptOut"></pre>
<button type="button" class="copy-btn" id="copyPromptBtn">Copy prompt</button>
</div>
<script>
function esc(s){return (s==null?'':String(s)).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));}
function sourceBadge(s){return '<span class="badge">'+esc((s||'').toUpperCase())+'</span>';}
function compStrip(c){if(!c)return '';
 return '<div class="comp">'+esc(c.v_germline_aa||'-')+' | '+esc(c.cdr3_aa||'-')+' | '+esc(c.j_germline_aa||'-')+'</div>';}
function partnerLine(rec,pairs){if(!pairs||!pairs.length)return '';
 for(const p of pairs){
  if(p.pairing_key!==rec.pairing_key)continue;
  const other=(rec.chain==='alpha')?p.beta:p.alpha;
  if(other&&other.cdr3_aa!==rec.cdr3_aa){
   return '<div class="partner">paired '+esc(other.chain)+': <code>'+esc(other.cdr3_aa)+'</code> ('+esc(other.source)+')</div>';
  }
 }
 return '';}
function recordCard(rec,pairs,cls){let h='<div class="rec '+cls+'">';
 h+=sourceBadge(rec.source)+' <span class="muted">'+esc(rec.chain)+' / '+esc(rec.species)+'</span>';
 if(rec.mhc_is_cross_species)h+='<span class="xspecies">human HLA (transgenic)</span>';
 if(rec.external_url)h+=' <a href="'+esc(rec.external_url)+'" target="_blank" rel="noopener">source</a>';
 h+='<br>aa: <code>'+esc(rec.full_aa||rec.cdr3_aa)+'</code>';
 if(rec.full_aa_kind)h+=' <span class="kind">('+esc(rec.full_aa_kind)+')</span>';
 const nt=rec.full_nt||rec.cdr3_nt;
 const ntKind=rec.full_nt?rec.full_nt_kind:rec.cdr3_nt_kind;
 if(nt&&ntKind){h+='<br>nt: <code>'+esc(nt)+'</code> <span class="kind">('+esc(ntKind)+')</span>';}
 const genes=[];
 if(rec.v_gene)genes.push('V '+esc(rec.v_gene));
 if(rec.d_gene)genes.push('D '+esc(rec.d_gene));
 if(rec.j_gene)genes.push('J '+esc(rec.j_gene));
 if(genes.length)h+='<br>'+genes.join(' &middot; ');
 const cdrs=[];
 if(rec.cdr1_aa)cdrs.push('CDR1 '+esc(rec.cdr1_aa));
 if(rec.cdr2_aa)cdrs.push('CDR2 '+esc(rec.cdr2_aa));
 if(rec.cdr3_aa)cdrs.push('CDR3 '+esc(rec.cdr3_aa));
 if(cdrs.length)h+='<br>'+cdrs.join(' &middot; ');
 if(rec.epitope_aa||rec.mhc_a){
  h+='<br>epitope: '+esc(rec.epitope_aa||'unknown');
  if(rec.mhc_class||rec.mhc_a)h+=' &middot; MHC '+esc(rec.mhc_class||'')+' '+esc(rec.mhc_a||'');
 }
 if(rec.composition)h+='<br>'+compStrip(rec.composition);
 h+=partnerLine(rec,pairs);
 return h+'</div>';}
function renderRecords(data){if(!data)return '';
 let h='<div class="section"><h3>Exact records</h3>';
 if(data.exact&&data.exact.length){for(const rec of data.exact)h+=recordCard(rec,data.pairs,'exact');}
 else h+='<p class="muted">No exact records found.</p>';
 h+='</div><div class="section"><h3>Near neighbours</h3>';
 if(data.neighbours&&data.neighbours.length){for(const rec of data.neighbours)h+=recordCard(rec,data.pairs,'neigh');}
 else h+='<p class="muted">No near neighbours found.</p>';
 h+='</div>';
 if(data.warnings&&data.warnings.length){h+='<p class="warn">'+esc(data.warnings.map(w=>w.message||w.code).join('; '))+'</p>';}
 return h;}
function neighTable(ns){if(!ns||!ns.length)return '<p class="muted">No matching known TCRs found in the reference database for this CDR3.</p>';
 let h='<h3>Similar TCRs (inferred, not confirmed specificity)</h3><table><tr><th>CDR3</th><th>V</th><th>sim</th><th>epitope</th></tr>';
 for(const n of ns)h+=`<tr><td>${esc(n.cdr3_b_aa)}</td><td>${esc(n.v_b_gene)}</td><td>${esc(n.similarity)}</td><td>${esc(n.epitope_aa)}</td></tr>`;return h+'</table>';}
function pairedTable(b){
 const ns=b.neighbours||[];
 let h='<div class="card"><h3>Paired neighbours <span class="muted">(engine '+esc(b.engine)+', '+esc(b.total_candidates)+' paired candidates)</span></h3>';
 if(!ns.length){h+='<p class="muted">No paired neighbours returned.</p>';}
 else{h+='<p class="muted">Inferred, not confirmed specificity. Distance is the absolute paired tcrdist (lower is closer).</p>';
  h+='<table><tr><th>alpha CDR3</th><th>alpha V</th><th>beta CDR3</th><th>beta V</th><th>dist</th><th>sim</th><th>epitope</th></tr>';
  for(const n of ns)h+=`<tr><td><code>${esc(n.cdr3_a_aa)}</code></td><td>${esc(n.v_a_gene)}</td><td><code>${esc(n.cdr3_b_aa)}</code></td><td>${esc(n.v_b_gene)}</td><td>${esc(n.distance)}</td><td>${esc(n.similarity)}</td><td>${esc(n.epitope_aa)}</td></tr>`;
  h+='</table>';}
 if(b.warnings&&b.warnings.length)h+='<p class="warn">'+esc(b.warnings.map(w=>w.message||w.code).join('; '))+'</p>';
 return h+'</div>';}
function render(b){let h=`<div class="card"><h3>intent: ${esc(b.intent)} <span class="muted">(source ${esc(b.plan_source)}, llm ${b.llm_used})</span></h3>`;
 if(b.dossier){const d=b.dossier;h+=`<p><b>${esc(d.summary)}</b></p><p>chain: ${esc(d.chain)} &middot; species: ${esc(d.species)} &middot; status: ${esc(d.status)}</p>`;
  if(d.status==='partial' && !(d.genes&&d.genes.v) && !(d.known_epitopes&&d.known_epitopes.length)){h+='<p class="muted">A CDR3 on its own cannot identify V/D/J. Provide the V and J genes (V+J+CDR3), a gene name, or a full V(D)J sequence to get an annotation.</p>';}
  if(d.genes&&d.genes.v){h+=`<p>V: ${esc(d.genes.v.call)} (${esc(d.genes.v.score_method)})</p>`;}
  if(d.regions){h+='<p>';for(const k of ['cdr1','cdr2','cdr3']){if(d.regions[k]&&d.regions[k].aa)h+=`${k}: <code>${esc(d.regions[k].aa)}</code> `;}h+='</p>';}
  if(d.junction&&d.junction.cdr3_nt_is_synthetic){h+=`<p class="syn">cdr3_nt is synthetic (back-translated)</p>`;}
  if(d.known_epitopes&&d.known_epitopes.length){h+='<h3>Known epitopes (retrieved)</h3><table><tr><th>epitope</th><th>MHC</th><th>antigen</th></tr>';
   for(const e of d.known_epitopes)h+=`<tr><td>${esc(e.epitope_sequence)}</td><td>${esc(e.mhc_allele)}</td><td>${esc(e.antigen_name)}</td></tr>`;h+='</table>';}
  if(d.neighbours&&d.neighbours.length){h+=neighTable(d.neighbours);}
  if(d.warnings&&d.warnings.length){h+='<p class="warn">warnings: '+d.warnings.map(w=>esc(w.code)).join(', ')+'</p>';}}
 if(b.neighbours_result){const nr=b.neighbours_result;h+=neighTable(nr.neighbours);if(nr.warnings&&nr.warnings.length){h+='<p class="warn">'+nr.warnings.map(w=>esc(w.message||w.code)).join('; ')+'</p>';}}
 if(b.search_result){h+=`<p class="muted">search returned ${b.search_result.total} records</p>`;}
 return h+'</div>';}
function shade(c){c=c||0;if(c>=0.9)return 'background:#08519c;color:#fff';if(c>=0.7)return 'background:#3182bd;color:#fff';if(c>=0.5)return 'background:#6baed6';if(c>=0.3)return 'background:#bdd7e7';return '';}
function cell(t,s){return '<span style="'+s+'">'+esc(t)+'</span>';}
function renderAlign(b){let h='<div class="card"><h3>engine: '+esc(b.engine)+' <span class="muted">('+esc(b.n_sequences)+' sequences, '+esc(b.mean_pct_identity)+'% identity, view '+esc(b.view)+')</span></h3>';
 h+='<p class="muted">Shading shows per-column conservation (darker is more conserved).</p><pre>';
 const cons=b.conservation||[];
 for(const rec of (b.records||[])){
  const name=esc(rec.name.padEnd(12));
  if(rec.aligned_aa && rec.aligned_nt){
   let aa=name+' aa  ',nt=name+' nt  ';
   for(let i=0;i<rec.aligned_aa.length;i++){aa+=cell(' '+rec.aligned_aa[i]+' ',shade(cons[i]));nt+=cell(rec.aligned_nt.substr(3*i,3),shade(cons[i]));}
   h+=aa+'<br>'+nt+'<br><br>';
  }else{
   const s=rec.aligned_aa||rec.aligned||rec.aligned_nt||'';let row=name+'  ';
   for(let i=0;i<s.length;i++){row+=cell(s[i],shade(cons[i]));}
   h+=row+'<br>';
  }
 }
 h+='</pre>';
 if(b.warnings&&b.warnings.length){h+='<p class="warn">warnings: '+b.warnings.map(w=>esc(w.code)).join(', ')+'</p>';}
 return h+'</div>';}
function seqBlock(label,seq){if(!seq)return '';return '<p class="muted" style="margin:.3rem 0 .1rem">'+esc(label)+' ('+seq.length+' aa)</p><div class="comp" style="word-break:break-all;white-space:pre-wrap">'+esc(seq)+'</div>';}
function renderReconstruct(b){
 if(!b.v_found||!b.j_found){return '<p class="warn">Germline not found'+(!b.v_found?' for V '+esc(b.v_gene):'')+(!b.j_found?' for J '+esc(b.j_gene):'')+'. Check the gene names.</p>';}
 let h='<div class="card"><h3>'+esc(b.v_gene)+' + '+esc(b.cdr3_aa)+' + '+esc(b.j_gene)+' <span class="muted">('+esc(b.species)+', V '+esc(b.v_allele_used||'?')+' / J '+esc(b.j_allele_used||'?')+')</span></h3>';
 if(b.genes_inferred){h+='<p class="echo">V and J inferred from '+esc(b.inference_support)+' database record'+(b.inference_support==1?'':'s')+' carrying this exact CDR3 (most common pairing). This is a frequency inference, not a germline assignment.';
  if(b.inference_alternatives&&b.inference_alternatives.length>1){h+='<br><span class="muted">other pairings: '+esc(b.inference_alternatives.slice(1).join(', '))+'</span>';}
  h+='</p>';}
 h+=seqBlock('variable domain',b.full_aa);
 if(b.full_chain_aa){h+=seqBlock('full chain (variable + constant, reconstructed)',b.full_chain_aa);
  h+='<p class="kind">constant: '+esc(b.constant_source||'')+'</p>';}
 else{h+='<p class="muted">Full chain unavailable: no vendored constant region for '+esc(b.species)+' '+esc((b.v_gene||'').slice(0,3))+'.</p>';}
 h+='<p class="kind">'+esc(b.note||'')+'</p>';
 return h+'</div>';}
function callRow(label,call){
 if(!call)return '';
 return '<p class="kind">'+esc(label)+': '+esc(call.alleles.join(', '))+
  ' <span class="muted">(identity '+(call.identity*100).toFixed(1)+'%, span '+
  esc(call.aligned_span[0])+' to '+esc(call.aligned_span[1])+')</span></p>';}
function renderAssign(b){
 if(!b.chain){return '<p class="warn">No germline alignment found for this sequence'+
  (b.warnings&&b.warnings.length?' ('+b.warnings.map(esc).join(', ')+')':'')+'.</p>';}
 let h='<div class="card"><h3>'+esc(b.chain)+' chain <span class="muted">('+esc(b.species)+
  ', input: '+esc(b.input_kind)+')</span></h3>';
 if(b.v_determinable){h+=callRow('V',b.v_call);}
 else{
  h+='<p class="warn">V not determinable: '+esc(b.v_reason||'')+'</p>';
  if(b.v_db_inference&&b.v_db_inference.length){
   h+='<div class="db-inference"><h4>database frequency inference (weaker signal)</h4><ul>'+
    b.v_db_inference.map(x=>'<li>'+esc(x.chain)+' '+esc(x.v_gene)+' / '+esc(x.j_gene)+
     ' (n='+esc(x.count)+')</li>').join('')+'</ul></div>';}}
 h+=callRow('J',b.j_call);
 if(b.d_call){h+=callRow('D'+(b.d_call.low_confidence?' (low confidence)':''),b.d_call);}
 if(b.constant_call){h+=callRow('constant',b.constant_call);}
 if(b.regions&&Object.keys(b.regions).length){
  h+='<p class="kind">regions: '+Object.entries(b.regions)
   .map(([k,v])=>esc(k)+' '+(v*100).toFixed(1)+'%').join(', ')+'</p>';}
 if(b.cdr3_aa){h+='<p class="kind">CDR3: '+esc(b.cdr3_aa)+'</p>';}
 if(b.reconstruction){h+=renderReconstruct(b.reconstruction);}
 if(b.warnings&&b.warnings.length){h+='<p class="warn">warnings: '+b.warnings.map(esc).join(', ')+'</p>';}
 return h+'</div>';}
function cardWrap(title,inner){return '<div class="card"><h3>'+esc(title)+'</h3>'+inner+'</div>';}
function renderBlock(block){
 const t=block.tool;
 if(t==='records')return cardWrap(block.title,renderRecords(block.data));
 if(t==='assign')return renderAssign(block.data);
 if(t==='dossier')return render({dossier:block.data,intent:'dossier',plan_source:'query_router',llm_used:false});
 if(t==='ask')return render(block.data);
 if(t==='similar')return render({neighbours_result:block.data,intent:'similar',plan_source:'query_router',llm_used:false});
 return cardWrap(block.title,'<pre>'+esc(JSON.stringify(block.data))+'</pre>');
}
function renderQueryResult(b){
 const u=b.understood||{};
 let h='<div class="echo">understood as: <b>'+esc(u.detected_type)+'</b> -&gt; '+esc((u.tools||[]).join(', '));
 if(u.species)h+=' (species: '+esc(u.species)+')';
 h+='</div>';
 if(u.note)h+='<p class="muted">'+esc(u.note)+'</p>';
 for(const block of (b.blocks||[]))h+=renderBlock(block);
 if(b.warnings&&b.warnings.length)h+='<p class="warn">'+esc(b.warnings.join('; '))+'</p>';
 out.innerHTML=h;
}
const f=document.getElementById('f'),out=document.getElementById('out'),chips=document.getElementById('chips');
function speciesOverride(){const v=document.getElementById('sp').value;return v?v:null;}
async function runQuery(q,spOverride,force){
 out.innerHTML='<p class="loading">Searching...</p>';
 try{
  const body={query:q};
  if(spOverride)body.species=spOverride;
  if(force)body.force=force;
  const r=await fetch('/v1/tcr/query',{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify(body)});
  if(!r.ok){out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
  const b=await r.json();renderQueryResult(b);
 }catch(err){out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
}
f.addEventListener('submit',async e=>{e.preventDefault();const btn=f.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Searching...';
 try{
  const q=document.getElementById('q').value,sp=speciesOverride();
  await runQuery(q,sp,null);
 }finally{btn.disabled=false;btn.textContent=t0;}});
function toggleHiddenForm(id){const el=document.getElementById(id);el.style.display=(el.style.display==='none'||!el.style.display)?'block':'none';}
chips.addEventListener('click',e=>{
 const btn=e.target.closest('.chip');if(!btn)return;
 const tool=btn.dataset.tool;
 if(tool==='reconstruct'){toggleHiddenForm('rcform');return;}
 if(tool==='align'){toggleHiddenForm('alignform');return;}
 if(tool==='paired'){toggleHiddenForm('pairform');return;}
 const q=document.getElementById('q').value,sp=speciesOverride();
 runQuery(q,sp,tool);
});
const af=document.getElementById('af'),a_out=document.getElementById('a_out');
af.addEventListener('submit',async e=>{e.preventDefault();const btn=af.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Aligning...';a_out.innerHTML='<p class="loading">Aligning...</p>';
 try{
  const r=await fetch('/v1/tcr/align',{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify({species:document.getElementById('a_sp').value,
    chain:document.getElementById('a_chain').value,
    segment:document.getElementById('a_seg').value,
    translate:document.getElementById('a_translate').checked})});
  if(!r.ok){a_out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
  const b=await r.json();a_out.innerHTML=renderAlign(b);
 }catch(err){a_out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
 finally{btn.disabled=false;btn.textContent=t0;}});
const rcf=document.getElementById('rcf'),rc_out=document.getElementById('rc_out');
rcf.addEventListener('submit',async e=>{e.preventDefault();const btn=rcf.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Working...';rc_out.innerHTML='<p class="loading">Working...</p>';
 try{
  const seq=document.getElementById('rc_seq').value;
  const sp=document.getElementById('rc_sp').value;
  const vg=document.getElementById('rc_v').value;
  const jg=document.getElementById('rc_j').value;
  if(vg&&jg){
   const r=await fetch('/reconstruct',{method:'POST',headers:{'Content-Type':'application/json'},
    body:JSON.stringify({v_gene:vg,j_gene:jg,cdr3_aa:seq,species:sp})});
   if(!r.ok){rc_out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
   const b=await r.json();rc_out.innerHTML=renderReconstruct(b);
  }else{
   const r=await fetch('/v1/tcr/assign',{method:'POST',headers:{'Content-Type':'application/json'},
    body:JSON.stringify({sequence:seq,species:sp,want_d:true})});
   if(!r.ok){rc_out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
   const b=await r.json();rc_out.innerHTML=renderAssign(b);
  }
 }catch(err){rc_out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
 finally{btn.disabled=false;btn.textContent=t0;}});
const pf=document.getElementById('pf'),p_out=document.getElementById('p_out');
pf.addEventListener('submit',async e=>{e.preventDefault();const btn=pf.querySelector('button');const t0=btn.textContent;btn.disabled=true;btn.textContent='Finding...';p_out.innerHTML='<p class="loading">Finding...</p>';
 try{
  const r=await fetch('/v1/tcr/similar_paired',{method:'POST',headers:{'Content-Type':'application/json'},
   body:JSON.stringify({cdr3_a:document.getElementById('p_ca').value,v_a:document.getElementById('p_va').value,
    cdr3_b:document.getElementById('p_cb').value,v_b:document.getElementById('p_vb').value,
    species:document.getElementById('p_sp').value})});
  if(!r.ok){p_out.innerHTML='<p class="warn">Request failed ('+r.status+')</p>';return;}
  const b=await r.json();p_out.innerHTML=pairedTable(b);
 }catch(err){p_out.innerHTML='<p class="warn">Error: '+esc(String(err))+'</p>';}
 finally{btn.disabled=false;btn.textContent=t0;}});
const MCP_CONFIG_TEXT='{"mcpServers":{"tcr-explorer":{"command":"uvx","args":["--from","tcr-explorer[tcrdist]","tcr-explorer-mcp"]}}}';
const MCP_PROMPT_TEXT="Set up the TCR Explorer MCP server so you can answer T cell receptor questions against real immunology databases. Add an MCP server named tcr-explorer that runs `uvx --from 'tcr-explorer[tcrdist]' tcr-explorer-mcp` (if uvx is unavailable, `pip install 'tcr-explorer[tcrdist]'` then run `python -m tcr_explorer.mcp_server`). The tcrdist extra gives authoritative similarity scoring; drop [tcrdist] for a leaner install that falls back to BLOSUM automatically. It exposes these read only tools: retrieve_tcr_records, assign_tcr_alleles, get_tcr_dossier, find_similar_tcrs, find_similar_paired_tcrs, and align_tcr_genes. After adding it, confirm the connection and suggest three example questions I can ask.";
document.getElementById('mcpConfigOut').textContent=MCP_CONFIG_TEXT;
document.getElementById('mcpPromptOut').textContent=MCP_PROMPT_TEXT;
document.getElementById('copyConfigBtn').addEventListener('click',()=>{navigator.clipboard.writeText(MCP_CONFIG_TEXT).catch(()=>{});});
document.getElementById('copyPromptBtn').addEventListener('click',()=>{navigator.clipboard.writeText(MCP_PROMPT_TEXT).catch(()=>{});});
</script></body></html>"""


@app.get("/ui", response_class=HTMLResponse)
def ui() -> HTMLResponse:
    return HTMLResponse(_UI_HTML)
