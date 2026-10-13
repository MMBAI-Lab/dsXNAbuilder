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
VER=$("$PY" -c "import dsxnabuilder; print(dsxnabuilder.__version__)")
echo "$VER" > web/VERSION

# Version-stamp the assets the page loads (GitHub Pages caches for 10 min
# under unchanging names). Idempotent.
sed -i -E "s@(href=\"style\.css)(\?v=[^\"]*)?\"@\1?v=$VER\"@g" web/index.html
sed -i -E "s@(src=\"app\.js)(\?v=[^\"]*)?\"@\1?v=$VER\"@g" web/index.html
sed -i -E "s@(new URL\('worker\.js)(\?v=[^']*)?'@\1?v=$VER'@g" web/app.js

du -sh web/pyodide web/vendor web/dsxnabuilder.zip
