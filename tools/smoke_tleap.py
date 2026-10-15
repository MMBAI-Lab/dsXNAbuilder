#!/usr/bin/env python3
"""Offline smoke test of every library residue with real AmberTools.

    python tools/smoke_tleap.py [--amberhome /opt/amber26] [--jobs 8] [--random 12]

For each residue of data/residue_index.json a short duplex is built with that
residue at its terminal position (5'-cap first, 3'-cap last, internal in the
middle) through the same Session the web app uses, and the package's tleap
input is run: it must finish with Errors = 0, without "Created a new atom"
(an atom name the unit does not know), without atoms added by tleap (every
hydrogen must come from the PDB under the unit's name) and with an integer
total charge.
``--random N`` adds N random modified duplexes that also go through the
restrained sander minimization (500 steps) and must end without NaN.

``--prune`` moves the library residues that fail into the "unavailable" list of
residue_index.json, with the reason, and deletes their files: some modXNA
residues build alone but not in a duplex (e.g. a sugar fragment that types O3'
as a hydroxyl, so tleap has no C-O-P angle for the linkage). The curated set
must never fail.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import io
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from dsxnabuilder.catalog import Catalog  # noqa: E402
from dsxnabuilder.fragments import base_families, scaffold_context  # noqa: E402
from dsxnabuilder.reslib import ResidueLibrary  # noqa: E402
from dsxnabuilder.service import Session  # noqa: E402

PAIR = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}


def duplex_for(key: str, res: dict, catalog, families):
    """A short duplex with the residue at its own position (5' end, middle or 3'
    end), on a model with its sugar pucker (dsxnabuilder.fragments.scaffold_context)."""
    bb, sugar, base = res["triplet"]
    cap = res["cap"]
    letter, chem, model = scaffold_context((bb, sugar, base), catalog, families)
    rna = chem == "RNA"
    tok = f"{{{bb}/{sugar}/{base}:{letter}}}"
    flank = "r(GC)" if rna else "d(GC)"
    if cap == "5cap":
        top, seq = tok + flank + flank, letter + "GCGC"
    elif cap == "3cap":
        top, seq = flank + flank + tok, "GCGC" + letter
    else:
        top, seq = flank + tok + flank, "GC" + letter + "GC"
    comp = "".join(PAIR[c] if not (rna and PAIR[c] == "T") else "U" for c in seq)
    bottom = ("r(" if rna else "d(") + comp + ")"
    return {"top": top, "bottom": bottom, "model": model}


def run_package(zbytes: bytes, env: dict, minimize: bool) -> dict:
    d = tempfile.mkdtemp(prefix="smoke-", dir=os.path.join(ROOT, "build"))
    zipfile.ZipFile(io.BytesIO(zbytes)).extractall(d)
    work = os.path.join(d, "dsxnabuilder")
    p = subprocess.run(["tleap", "-f", "leap.in"], cwd=work, env=env, capture_output=True,
                       text=True, timeout=300)
    log = p.stdout
    first = re.search(r"Error!\n([^\n]+)", log)
    out = {"errors": int(re.search(r"Errors = (\d+)", log).group(1)) if "Errors =" in log else -1,
           "first_error": first.group(1).strip() if first else None,
           "new_atom": "Created a new atom" in log,
           # tleap fills atoms missing from the PDB ("Leap added N missing atoms")
           "added": bool(re.search(r"Leap added [1-9]\d* missing atom|Added missing heavy atom",
                                   log))}
    m = re.search(r"Total unperturbed charge:\s+(-?[\d.]+)", log)
    out["charge"] = float(m.group(1)) if m else None
    if minimize and out["errors"] == 0:
        # validation only (the product leaves minimization to the user): 500 steps in vacuum
        with open(os.path.join(work, "min.in"), "w") as fh:
            fh.write("smoke minimization\n &cntrl\n  imin = 1, maxcyc = 500, ncyc = 250,\n"
                     "  ntb = 0, igb = 0, cut = 999.0,\n /\n")
        subprocess.run(["sander", "-O", "-i", "min.in", "-p", "system.parm7", "-c", "system.rst7",
                        "-ref", "system.rst7", "-r", "min.rst7", "-o", "min.out"], cwd=work,
                       env=env, capture_output=True, timeout=1200)
        mo = open(os.path.join(work, "min.out")).read() if os.path.exists(
            os.path.join(work, "min.out")) else ""
        out["min_ok"] = "FINAL RESULTS" in mo and "NaN" not in mo
    subprocess.run(["rm", "-rf", d])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--amberhome", default=os.environ.get("AMBERHOME", "/opt/amber26"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--random", type=int, default=12)
    ap.add_argument("--prune", action="store_true",
                    help="drop the failing library residues from the library (see above)")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    env = dict(os.environ, AMBERHOME=a.amberhome, LC_ALL="C",
               PATH=os.path.join(a.amberhome, "bin") + os.pathsep + os.environ["PATH"])
    os.makedirs(os.path.join(ROOT, "build"), exist_ok=True)
    rl = ResidueLibrary.load_default()
    s = Session(residue_library=rl)

    catalog, families = Catalog.from_json(), base_families()
    jobs = []
    for key, res in sorted(rl.residues.items()):
        jobs.append((key, duplex_for(key, res, catalog, families), False))
    rng = random.Random(a.seed)
    toks_rna = ["mA", "mC", "mG", "mU", "fA", "fC", "fG", "fU", "eA", "eT", "eC^5mC", "+A",
                "+T", "+G", "+C^5mC", "rA", "rG"]
    for i in range(a.random):
        n = rng.randint(6, 14)
        toks, letters = [], []
        for _ in range(n):
            t = rng.choice(toks_rna)
            toks.append(t)
            letters.append(re.findall(r"[ACGTU]", t)[0])
        link = lambda: rng.choice(["", "*1", "*2", "*"])
        top = "".join(t + (link() if j < n - 1 else "") for j, t in enumerate(toks))
        comp = "".join({"A": "U", "U": "A", "T": "A", "G": "C", "C": "G"}[c] for c in letters)
        jobs.append((f"random{i}", {"top": top, "bottom": f"r({comp})", "model": rng.choice([20, 22, 1, 54])}, True))

    def one(job):
        key, req, minimize = job
        out = s.build(dict(req, package=True))
        if not out["ok"]:
            return key, req, {"build_error": [i["code"] for i in out["issues"] if i["level"] == "error"]}
        return key, req, run_package(out["zip"], env, minimize)

    fails, pruned = 0, {}
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for key, req, r in ex.map(one, jobs):
            bad = ("build_error" in r or r["errors"] != 0 or r["new_atom"] or r["added"]
                   or r["charge"] is None or abs(r["charge"] - round(r["charge"])) > 1e-3
                   or r.get("min_ok") is False)
            if bad:
                fails += 1
                print(f"FAIL {key}: {r}  {req}")
                if key in rl.residues:
                    pruned[key] = ("tleap in a duplex: " + (r.get("first_error") or "")
                                   if r.get("first_error") else f"build in a duplex: {r}")
    print(f"{len(jobs)} systems, {fails} failures")
    if a.prune and pruned:
        sys.path.insert(0, os.path.join(ROOT, "tools"))
        from build_residue_library import curated
        bad_curated = [k for k in pruned if tuple(k.split("-")) in set(curated())]
        if bad_curated:
            print(f"curated residues fail, nothing pruned: {bad_curated}")
            return 1
        index_path = os.path.join(ROOT, "data", "residue_index.json")
        with open(index_path, encoding="utf-8") as fh:
            index = json.load(fh)
        for key, why in pruned.items():
            os.remove(os.path.join(ROOT, "web", "lib", "residues", index["residues"][key]["file"]))
            del index["residues"][key]
            index.setdefault("unavailable", {})[key] = why
        with open(index_path, "w", encoding="utf-8") as fh:
            json.dump(index, fh, indent=1, sort_keys=True)
        print(f"pruned {len(pruned)}: {len(index['residues'])} residues left")
        return 1 if fails - len(pruned) else 0
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
