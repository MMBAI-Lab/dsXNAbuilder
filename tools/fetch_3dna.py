#!/usr/bin/env python3
# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Fetch the 3DNA data files dsXNAbuilder needs (text only, no binaries).

    python tools/fetch_3dna.py                       # pinned GitHub mirror
    python tools/fetch_3dna.py --official x3dna-v2.4-linux-64bit.tar.gz

What is fetched
---------------
* ``fiber/README``            twist/rise table of the 55 fiber models
* ``fiber/strNN/*.pdb|.rpt``  repeating units (Lu & Olson 2003, Table 4)
* ``config/Atomic_X.pdb``     Olson et al. (2001) standard bases, X=A,C,G,T,U

Provenance
----------
The official 3DNA v2.4.8 tarball is distributed through the 3DNA forum
(http://forum.x3dna.org/site-announcements/download-instructions/), whose
download links no longer work (checked 2026-10-11), and is licensed
CC BY-NC 4.0. The data files are therefore read from a third-party mirror of
3DNA pinned to one commit, which is the definitive source of this project;
every file's sha256 is recorded in ``PROVENANCE.json``. ``--official`` is kept
in case an official tarball becomes available: it diffs it against the mirror. The mirror's own LICENSE file cannot relicense 3DNA: the data stay
under the 3DNA terms (see LICENSES/README.md).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tarfile
import urllib.request

MIRROR = "pipepipexiaji/X3DNA"
MIRROR_COMMIT = "1363f200eaa75a7931817029f5f65af15bb0ac0b"   # 2017-09-02
STANDARD_BASES = ["A", "C", "G", "T", "U"]

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(HERE, "..", "external", "3dna")


def wanted(path: str) -> bool:
    if path == "fiber/README":
        return True
    parts = path.split("/")
    if len(parts) == 3 and parts[0] == "fiber" and parts[1].startswith("str") \
            and parts[2].endswith((".pdb", ".rpt", ".dat")):
        return True
    return path in {f"config/Atomic_{b}.pdb" for b in STANDARD_BASES}


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "dsxnabuilder-fetch"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def fetch_mirror(out: str) -> dict:
    tree = json.loads(_get(f"https://api.github.com/repos/{MIRROR}/git/trees/"
                           f"{MIRROR_COMMIT}?recursive=1"))
    paths = sorted(e["path"] for e in tree["tree"]
                   if e["type"] == "blob" and wanted(e["path"]))
    files = {}
    for p in paths:
        data = _get(f"https://raw.githubusercontent.com/{MIRROR}/{MIRROR_COMMIT}/{p}")
        files[p] = _store(out, p, data)
    return {"source": f"github.com/{MIRROR}", "commit": MIRROR_COMMIT,
            "note": "third-party mirror of 3DNA v2.x data files; 3DNA licence applies",
            "files": files}


def fetch_official(tarball: str, out: str) -> dict:
    files = {}
    with tarfile.open(tarball) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            # x3dna-v2.4/fiber/str01/A.pdb -> fiber/str01/A.pdb
            parts = m.name.split("/")
            rel = "/".join(parts[1:]) if len(parts) > 1 else m.name
            if wanted(rel):
                files[rel] = _store(out, rel, tf.extractfile(m).read())
    return {"source": os.path.basename(tarball),
            "sha256_tarball": _sha_file(tarball),
            "note": "official 3DNA distribution (CC BY-NC 4.0)", "files": files}


def _store(out: str, rel: str, data: bytes) -> str:
    dest = os.path.join(out, rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "wb") as fh:
        fh.write(data)
    return hashlib.sha256(data).hexdigest()


def _sha_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--official", help="official 3DNA tarball (x3dna-v2.4*.tar.gz)")
    ap.add_argument("--out", default=DEFAULT_OUT)
    a = ap.parse_args(argv)
    out = os.path.abspath(a.out)
    prov_path = os.path.join(out, "PROVENANCE.json")
    previous = None
    if os.path.exists(prov_path):
        with open(prov_path, encoding="utf-8") as fh:
            previous = json.load(fh)

    prov = fetch_official(a.official, out) if a.official else fetch_mirror(out)
    n_str = len({p.split("/")[1] for p in prov["files"] if p.startswith("fiber/str")})
    print(f"{len(prov['files'])} files, {n_str} fiber models -> {out}")

    if previous and previous.get("source") != prov["source"]:
        old, new = previous["files"], prov["files"]
        changed = sorted(p for p in set(old) & set(new) if old[p] != new[p])
        only_old = sorted(set(old) - set(new))
        only_new = sorted(set(new) - set(old))
        prov["diff_vs_previous"] = {"previous_source": previous["source"],
                                    "changed": changed, "only_previous": only_old,
                                    "only_current": only_new}
        print(f"diff vs {previous['source']}: {len(changed)} changed, "
              f"{len(only_old)} only there, {len(only_new)} only here")

    with open(prov_path, "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=1, sort_keys=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
