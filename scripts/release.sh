#!/usr/bin/env bash
#
# Release helper for tcr-explorer.
#
#   scripts/release.sh            dry run: tests + build + twine check + clean-install sanity. No upload.
#   scripts/release.sh --test     also upload to TestPyPI  (needs a [testpypi] token in ~/.pypirc)
#   scripts/release.sh --publish  also upload to PyPI and tag the release (needs a [pypi] token)
#
# The version is the single source of truth in pyproject.toml. Bump it there and
# commit before publishing (an upload refuses a dirty tree or an existing version).
# The token is never passed here: twine reads ~/.pypirc (or TWINE_USERNAME/TWINE_PASSWORD).
#
set -euo pipefail

cd "$(dirname "$0")/.."
PY=$([ -x .venv/bin/python ] && echo .venv/bin/python || echo python3)

MODE=dryrun
case "${1:-}" in
  --test)    MODE=test ;;
  --publish) MODE=publish ;;
  "")        MODE=dryrun ;;
  *) echo "usage: scripts/release.sh [--test|--publish]" >&2; exit 2 ;;
esac

VERSION=$("$PY" - <<'PY'
import re, pathlib
t = pathlib.Path("pyproject.toml").read_text()
print(re.search(r'^version\s*=\s*"([^"]+)"', t, re.M).group(1))
PY
)
echo "==> tcr-explorer release: version ${VERSION}, mode ${MODE}"

# Uploads require a clean main and an unused version.
if [ "$MODE" != dryrun ]; then
  branch=$(git rev-parse --abbrev-ref HEAD)
  [ "$branch" = main ] || { echo "refuse: not on main (on ${branch})" >&2; exit 1; }
  [ -z "$(git status --porcelain)" ] || { echo "refuse: working tree not clean; commit first" >&2; exit 1; }
fi
if [ "$MODE" = publish ]; then
  code=$(curl -s -o /dev/null -w '%{http_code}' "https://pypi.org/pypi/tcr-explorer/${VERSION}/json")
  [ "$code" = 404 ] || { echo "refuse: version ${VERSION} already on PyPI (HTTP ${code}). Bump it in pyproject.toml." >&2; exit 1; }
  if git rev-parse "v${VERSION}" >/dev/null 2>&1; then echo "refuse: git tag v${VERSION} already exists" >&2; exit 1; fi
fi

echo "==> ensuring build + twine"
"$PY" -m pip install -q --upgrade build twine

echo "==> running tests"
PYTHONPATH=src:. "$PY" -m pytest tests/ -q

echo "==> building sdist + wheel"
rm -rf dist build src/tcr_explorer.egg-info
"$PY" -m build

echo "==> twine check"
"$PY" -m twine check dist/*

echo "==> clean-install sanity from the built wheel"
TMPV="$(mktemp -d)/venv"
"$PY" -m venv "$TMPV"
"$TMPV/bin/pip" install -q dist/*.whl
( cd "$TMPV" && "$TMPV/bin/python" - "$VERSION" <<'PY'
import sys, importlib.metadata as m
want = sys.argv[1]
got = m.version("tcr-explorer")
assert got == want, f"version mismatch: wheel {got} != pyproject {want}"
from tcr_explorer.cdr_enricher import get_cdr1_cdr2
assert get_cdr1_cdr2("TRBV19", "HUMAN")["allele"] == "TRBV19*01", "bundled germline broken"
print("sanity OK: version matches and bundled germline resolves offline")
PY
)
rm -rf "$(dirname "$TMPV")"

case "$MODE" in
  dryrun)
    echo "==> dry run complete. Artifacts in dist/. Re-run with --test or --publish to upload." ;;
  test)
    echo "==> uploading to TestPyPI"
    "$PY" -m twine upload --repository testpypi dist/*
    echo "TestPyPI: https://test.pypi.org/project/tcr-explorer/${VERSION}/" ;;
  publish)
    echo "==> uploading to PyPI"
    "$PY" -m twine upload dist/*
    git tag -a "v${VERSION}" -m "tcr-explorer ${VERSION}"
    git push origin "v${VERSION}"
    echo "PyPI:   https://pypi.org/project/tcr-explorer/${VERSION}/"
    echo "tagged: v${VERSION}" ;;
esac
