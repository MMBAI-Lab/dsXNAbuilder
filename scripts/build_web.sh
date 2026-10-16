#!/bin/bash
# Build the static web app in web/: self-hosted Pyodide + numpy + NGL, and the
# dsxnabuilder package with its data zipped for the worker.
# Usage: scripts/build_web.sh        (then: python3 scripts/serve.py)
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PYTHON:-.venv/bin/python}
[ -x "$PY" ] || PY=python3

"$PY" scripts/fetch_web_vendor.py

rm -f web/dsxnabuilder.zip
"$PY" - <<'PY'
import pathlib, zipfile
DATA = ["fiber_models.json", "modxna_catalog.json", "aliases.json", "base_pairing.json",
        "residue_index.json", "amber_residues.json"]
with zipfile.ZipFile("web/dsxnabuilder.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(pathlib.Path("dsxnabuilder").rglob("*.py")):
        z.write(p)
    for name in DATA:
        z.write(pathlib.Path("data") / name)
PY
# VERSION names the build: the release ("0.2.0") when built from its tag, otherwise package
# version + source commit (+ ".dirty<time>" for an uncommitted tree), so that every deployment
# has a new version and browsers do not keep serving cached assets under the old one.
# scripts/deploy_danslab.sh stamps it into the deployed copy (asset ?v= and PAGE_VERSION); the
# sources keep "dev".
VER=$("$PY" -c "import dsxnabuilder; print(dsxnabuilder.__version__)")
CLEAN=$([ -z "$(git status --porcelain --untracked-files=no)" ] && echo yes || echo no)
if [ "$CLEAN" = yes ] && [ "$(git describe --exact-match --tags HEAD 2>/dev/null)" = "v$VER" ]; then
  :                                                   # a release, built from its tag
else
  VER="$VER.g$(git rev-parse --short HEAD)"
  [ "$CLEAN" = yes ] || VER="$VER.dirty$(date +%s)"
fi
echo "$VER" > web/VERSION
echo "version $VER"

du -sh web/pyodide web/vendor web/dsxnabuilder.zip
