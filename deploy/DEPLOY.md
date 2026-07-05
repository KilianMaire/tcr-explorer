# Deploying the public germline-only demo to Hugging Face Spaces

This folder holds everything needed to run a free, public, germline-only demo of
TCR Explorer on Hugging Face Spaces (Docker SDK). The demo serves the web UI and
REST API with germline features only. It downloads no record databases, so it
stays within the data licenses (see the note at the end).

## What is in this folder

- `Dockerfile`, installs TCR Explorer from the public GitHub main tarball and
  runs the FastAPI app on port 7860 with `TCR_EXPLORER_DEMO=1` set.
- `README.md`, the Space card. Its YAML frontmatter (`sdk: docker`,
  `app_port: 7860`) is what tells Hugging Face how to build and route the Space.

A Space is its own git repository, and Hugging Face expects `Dockerfile` and
`README.md` at the repository root. So these two files are copied to the root of
the Space repo, not into a `deploy/` subfolder there.

## One time setup

1. Create the Space on https://huggingface.co/new-space
   - Owner, your account.
   - Space name, for example `tcr-explorer`.
   - License, MIT.
   - Space SDK, choose **Docker** (blank template).
   - Visibility, Public.

2. Authenticate git for Hugging Face once (needed to push):
   ```bash
   pip install huggingface_hub
   huggingface-cli login          # paste a write token from https://huggingface.co/settings/tokens
   ```

## Push the demo

From this repository root, with `<user>` replaced by your Hugging Face username:

```bash
# clone the (empty) Space repo next to this repo
git clone https://huggingface.co/spaces/<user>/tcr-explorer /tmp/tcrx-space

# copy the two Space files to its root
cp deploy/Dockerfile deploy/README.md /tmp/tcrx-space/

# commit and push
cd /tmp/tcrx-space
git add Dockerfile README.md
git commit -m "TCR Explorer germline-only demo"
git push
```

Hugging Face builds the image and serves it. The live URL is
`https://<user>-tcr-explorer.hf.space` , and `/` redirects to the `/ui` query box.

## Updating the demo

The `Dockerfile` installs from GitHub `main`, so a rebuild picks up the latest
code. Hugging Face rebuilds when the Space repo changes, not when this repo
changes. To force a rebuild after pushing new commits here, either use the
"Factory rebuild" button in the Space settings, or push any change to the Space
repo (for example bump a comment in its `Dockerfile`).

For a reproducible build, pin the install to a commit: in the `Dockerfile`,
replace `heads/main.tar.gz` with `<commit-sha>.tar.gz`.

## Why germline only

The record sources have mixed licenses. IMGT (the germline, bundled) and IEDB
are CC BY 4.0 and redistributable with attribution. VDJdb is CC BY-ND (no
derivatives), McPAS has no license (citation only), and TCR3d has no site
license. A public instance that served those records would make the operator
their redistributor, which the tool's download-on-first-run design deliberately
avoids. The demo therefore ships no records and runs germline features only. To
add a limited records feature within the licenses, IEDB (CC BY 4.0) is the one
source that could be enabled safely.
