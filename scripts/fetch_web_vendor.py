#!/usr/bin/env python3
# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Download the self-hosted runtime of the web app into web/ (no CDN at run time).

* Pyodide core (GitHub release pyodide-core-<v>.tar.bz2)
* numpy wheel for that Pyodide, from the Pyodide CDN, checked against the
  sha256 recorded in pyodide-lock.json
* NGL viewer (dist/ngl.js from the npm registry tarball, integrity checked)

Idempotent: skips what is already present at the pinned versions.
"""

from __future__ import annotations

import base64
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import urllib.request

PYODIDE = "314.0.7"
NGL = "2.4.0"
WEB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "web"))
CORE_FILES = ["pyodide.js", "pyodide.mjs", "pyodide.asm.mjs", "pyodide.asm.wasm",
              "python_stdlib.zip", "pyodide-lock.json", "package.json"]


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "dsxnabuilder-build"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def pyodide() -> None:
    dest = os.path.join(WEB, "pyodide")
    stamp = os.path.join(dest, ".version")
    if os.path.exists(stamp) and open(stamp).read().strip() == PYODIDE:
        print(f"pyodide {PYODIDE}: present")
    else:
        shutil.rmtree(dest, ignore_errors=True)
        os.makedirs(dest)
        url = (f"https://github.com/pyodide/pyodide/releases/download/{PYODIDE}/"
               f"pyodide-core-{PYODIDE}.tar.bz2")
        print(f"downloading {url}")
        with tarfile.open(fileobj=io.BytesIO(get(url)), mode="r:bz2") as tf:
            for f in CORE_FILES:
                m = tf.getmember(f"pyodide/{f}")
                with open(os.path.join(dest, f), "wb") as fh:
                    fh.write(tf.extractfile(m).read())
        with open(stamp, "w") as fh:
            fh.write(PYODIDE + "\n")
    with open(os.path.join(dest, "pyodide-lock.json")) as fh:
        lock = json.load(fh)
    pkg = lock["packages"]["numpy"]
    for dep in ["numpy"] + list(pkg.get("depends", [])):
        p = lock["packages"][dep]
        path = os.path.join(dest, p["file_name"])
        if os.path.exists(path) and _sha(open(path, "rb").read()) == p["sha256"]:
            continue
        url = f"https://cdn.jsdelivr.net/pyodide/v{PYODIDE}/full/{p['file_name']}"
        print(f"downloading {url}")
        data = get(url)
        if _sha(data) != p["sha256"]:
            sys.exit(f"sha256 mismatch for {p['file_name']}")
        with open(path, "wb") as fh:
            fh.write(data)
    print(f"numpy {pkg['version']}: ok")


def ngl() -> None:
    dest = os.path.join(WEB, "vendor")
    os.makedirs(dest, exist_ok=True)
    stamp = os.path.join(dest, ".ngl-version")
    if os.path.exists(stamp) and open(stamp).read().strip() == NGL:
        print(f"ngl {NGL}: present")
        return
    meta = json.loads(get(f"https://registry.npmjs.org/ngl/{NGL}"))
    tarball = get(meta["dist"]["tarball"])
    algo, b64 = meta["dist"]["integrity"].split("-", 1)
    if base64.b64encode(hashlib.new(algo, tarball).digest()).decode() != b64:
        sys.exit("ngl tarball integrity mismatch")
    with tarfile.open(fileobj=io.BytesIO(tarball), mode="r:gz") as tf:
        data = tf.extractfile(tf.getmember("package/dist/ngl.js")).read()
        lic = tf.extractfile(tf.getmember("package/LICENSE")).read()
    with open(os.path.join(dest, "ngl.js"), "wb") as fh:
        fh.write(data)
    with open(os.path.join(dest, "NGL-LICENSE.txt"), "wb") as fh:
        fh.write(lic)
    with open(stamp, "w") as fh:
        fh.write(NGL + "\n")
    print(f"ngl {NGL}: ok ({len(data) / 1e6:.1f} MB)")


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


if __name__ == "__main__":
    pyodide()
    ngl()
