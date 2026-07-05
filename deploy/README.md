---
title: TCR Explorer
emoji: 🧬
colorFrom: indigo
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
license: mit
---

# TCR Explorer (public instance)

A hosted, public instance of [TCR Explorer](https://github.com/KilianMaire/tcr-explorer),
a federated T cell receptor analysis tool, so researchers can query it from a
chat assistant, a browser, or the REST API without installing anything.

It exposes the same tools three ways:

- **`/mcp`**, the MCP server over HTTP, for MCP capable chat assistants,
- **REST API + `/openapi.json`**, for assistants that call a REST API,
- **`/ui`**, the human web query box.

## What data it serves

This instance serves the **IEDB** database (CC BY 4.0, with attribution) plus
the germline features that need no external records: allele assignment, chain
reconstruction, alignment, and CDR loops. Record retrieval and neighbour search
run against IEDB (about 207,000 receptor chains).

VDJdb, McPAS and TCR3d are **not** served here, their licenses do not permit
public redistribution. To query all sources, install the tool locally:

```bash
pip install "tcr-explorer[tcrdist]"
tcr-explorer-refresh
```

## Connect your chat assistant

### Claude (with a custom connector)

On a plan that supports custom connectors, add a remote MCP server pointing at:

```
https://kmbiotraxion-tcr-explorer.hf.space/mcp
```

Then ask T cell receptor questions in plain language. The assistant sees seven
read only tools: retrieve_tcr_records, assign_tcr_alleles, get_tcr_dossier,
find_similar_tcrs, find_similar_paired_tcrs, align_tcr_genes, and ask_tcr.

### ChatGPT (custom GPT with Actions)

Create a custom GPT, add an Action, and import the OpenAPI schema from:

```
https://kmbiotraxion-tcr-explorer.hf.space/openapi.json
```

The GPT can then call the REST endpoints (retrieve, assign, reconstruct,
similar, align) directly. If your ChatGPT account supports remote MCP, you can
instead point it at the `/mcp` URL above.

### Any MCP capable client

Point it at `https://kmbiotraxion-tcr-explorer.hf.space/mcp` (streamable HTTP,
no authentication).

## Notes

- This is a free tier instance. It sleeps after inactivity, so the first request
  after a while may be slow while it wakes.
- It is unauthenticated and lightly rate limited. Please be gentle.
- The IEDB snapshot is periodic, not live.

## Attribution and license

- IEDB data is used under CC BY 4.0 (Vita R et al., The Immune Epitope Database).
- The IMGT germline reference is used under CC BY 4.0.
- TCR Explorer itself is MIT licensed:
  https://github.com/KilianMaire/tcr-explorer .
- Package: https://pypi.org/project/tcr-explorer/ . DOI:
  https://doi.org/10.5281/zenodo.21204936 .
