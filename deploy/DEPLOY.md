# Deploying the public IEDB + MCP instance to Hugging Face Spaces

This folder deploys a free, public instance of TCR Explorer on Hugging Face
Spaces (Docker SDK). It serves three surfaces from the same functions: the human
web UI (`/ui`), the REST API (`/openapi.json`), and the MCP server over HTTP
(`/mcp`) for chat assistants. It serves IEDB records (CC BY 4.0) plus the
germline features, and serves no other record database (licensing).

## What is in this folder

- `Dockerfile`, installs TCR Explorer from a pinned GitHub commit tarball, copies
  the pre-built IEDB index into `/index`, and runs the FastAPI app on port 7860
  with `TCR_EXPLORER_MCP_HTTP`, `TCR_EXPLORER_DEMO`, and
  `TCR_EXPLORER_RATE_LIMIT_PER_MIN` set.
- `build_iedb_index.py`, builds the IEDB records index locally (the datacenter
  build runners time out downloading from iedb.org, so the index is built here
  and baked into the image instead of fetched at build time).
- `README.md`, the Space card. Its YAML frontmatter (`sdk: docker`,
  `app_port: 7860`) tells Hugging Face how to build and route the Space, and its
  body is the researcher-facing connection guide.

A Space is its own git repository, and Hugging Face expects the `Dockerfile`,
the `README.md`, and the index files at the repository root. So they are copied
to the root of the Space repo, not into a `deploy/` subfolder there.

## Push the deploy

The Space already exists at `KMBioTraxion/TCR_Explorer`. To update it:

```bash
# 1. build the IEDB index locally (needs the package installed, e.g. a venv)
python deploy/build_iedb_index.py /tmp/iedb-index

# 2. clone the Space repo and enable git LFS in it (Hugging Face requires binary
#    files to go through LFS/Xet; the default .gitattributes already LFS-tracks
#    *.parquet, so keep that rule)
git clone https://huggingface.co/spaces/KMBioTraxion/TCR_Explorer /tmp/tcrx-space
cd /tmp/tcrx-space
git lfs install --local

# 3. copy the Space files (Dockerfile, README, and the built index) to its root
cp /path/to/tcr-explorer/deploy/Dockerfile /path/to/tcr-explorer/deploy/README.md .
cp /tmp/iedb-index/records_index.parquet /tmp/iedb-index/records_index.meta.json .

# 4. commit and push (use a Hugging Face write token as the git password; the
#    LFS object uploads during the push)
git add -A
git lfs ls-files            # should list records_index.parquet
git commit -m "TCR Explorer IEDB + MCP instance"
git push
```

Hugging Face rebuilds the image. The build no longer touches iedb.org (it copies
the baked index), so it is faster and reliable, and container cold starts are
fast. LFS files are materialized in the build context, so the `COPY` finds the
real parquet.

Live URLs:

- UI: `https://kmbiotraxion-tcr-explorer.hf.space/ui`
- MCP: `https://kmbiotraxion-tcr-explorer.hf.space/mcp`
- OpenAPI: `https://kmbiotraxion-tcr-explorer.hf.space/openapi.json`

## Updating to newer code

The `Dockerfile` installs from a pinned commit tarball. To ship newer code, bump
the commit SHA in the `archive/<sha>.tar.gz` URL, copy the updated `Dockerfile`
to the Space repo, and push it (Hugging Face rebuilds on Space repo changes, not
on this repo's changes). To refresh the IEDB snapshot, re-run
`build_iedb_index.py`, copy the new parquet and meta into the Space repo, and
push (a plain "Factory rebuild" reuses the same baked index, so it does not
refresh the data).

## Why IEDB only

The record sources have mixed licenses. IMGT (the germline, bundled) and IEDB
are CC BY 4.0 and redistributable with attribution. VDJdb is CC BY-ND (no
derivatives), McPAS has no license (citation only), and TCR3d has no site
license. A public instance that served those would make the operator their
redistributor, which the tool's download-on-first-run design avoids. IEDB is the
one large source that can be served safely, so it is the only one baked in.

## Operational notes

- Free tier: the Space sleeps after inactivity (cold start on the next request).
- Unauthenticated: `/mcp` and the REST API are open, with a light per-IP rate
  limit (`TCR_EXPLORER_RATE_LIMIT_PER_MIN`, 120/min by default here).
- The IEDB snapshot is fixed at image build time, refresh it with a rebuild.
