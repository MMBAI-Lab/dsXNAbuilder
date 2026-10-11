#!/usr/bin/env python3
"""
Prebuild the modXNA residue library served with the static web app.

    python tools/build_residue_library.py [--jobs 8] [--modxna external/modxna]
                                          [--cpptraj external/cpptraj/bin]

Runs ``modxna.sh -i <triplet> -m MXR [--5cap|--3cap] --clean`` once per
residue, each in its own scratch directory (modxna.sh writes tmp.* files in
the working directory, so separate directories make the runs independent and
parallel). cpptraj comes first on PATH from the master build
(tools/update_cpptraj.sh): modxna.sh calls plain ``cpptraj`` from PATH and
uses $CPPTRAJ only for its version check. tleap and sander come from
$AMBERHOME (AmberTools 26).

Outputs
-------
web/lib/residues/<KEY>.lib     OFF library, unit named MXR (renamed per job)
web/lib/frcmod.modxna          modXNA parameters used to build them
data/residue_index.json        key -> triplet, cap, charge, heavy atoms, sha256
                               plus modXNA/cpptraj commits and frcmod sha256

KEY = <BACKBONE>-<SUGAR>-<BASE>-<int|5cap|3cap>.

The curated set covers what the web language writes without explicit blocks:
sugars d, r, m (2'-OMe), f (2'-F), e (MOE), + (LNA); bases A C G T U 5mC on
each; phosphodiester and both phosphorothioates; internal, 5'- and 3'-capped.
Canonical capped residues are included because a strand with a modXNA end must
take both ends from modXNA (terminal charges differ from Amber's).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from dsxnabuilder.catalog import CANONICAL_BACKBONE, CANONICAL_SUGAR, CAP_5_OH  # noqa: E402
from dsxnabuilder.offlib import parse_off  # noqa: E402

UNIT = "MXR"
SUGARS = {"DC2": "DNA", "RC3": "RNA", "OME": "RNA", "AF2": "RNA", "MOE": "RNA", "LNA": "RNA"}
BASES = {"DNA": ["DAA", "DCC", "DGG", "DTT", "RUU", "M5C"],
         "RNA": ["RAA", "RCC", "RGG", "RUU", "DTT", "M5C"]}
CHARGE = {"int": -1.0, "5cap": -0.320348, "3cap": -0.679652}


def residue_key(backbone: str, sugar: str, base: str, cap: str) -> str:
    return f"{backbone}-{sugar}-{base}-{cap}"


def curated() -> list:
    out = []
    for sugar, chem in SUGARS.items():
        po = CANONICAL_BACKBONE[chem]
        for base in BASES[chem]:
            canonical = sugar == CANONICAL_SUGAR[chem] and base in (
                ["DAA", "DCC", "DGG", "DTT"] if chem == "DNA" else ["RAA", "RCC", "RGG", "RUU"])
            for bb in (po, "PS1", "PS2"):
                if not (canonical and bb == po):
                    out.append((bb, sugar, base, "int"))
                out.append((bb, sugar, base, "3cap"))
            out.append((CAP_5_OH, sugar, base, "5cap"))
    return out


def build_one(job, modxna_sh: str, env: dict, scratch: str) -> dict:
    bb, sugar, base, cap = job
    key = residue_key(*job)
    work = tempfile.mkdtemp(prefix=key + "-", dir=scratch)
    with open(os.path.join(work, "in.txt"), "w") as fh:
        fh.write(f"{bb} {sugar} {base}\n")
    cmd = ["bash", modxna_sh, "-i", "in.txt", "-m", UNIT, "--clean"]
    if cap == "5cap":
        cmd.append("--5cap")
    elif cap == "3cap":
        cmd.append("--3cap")
    p = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True, timeout=900)
    lib = os.path.join(work, f"{UNIT}.lib")
    if p.returncode != 0 or not os.path.exists(lib):
        return {"key": key, "ok": False, "log": (p.stdout + p.stderr)[-3000:]}
    with open(lib, encoding="utf-8") as fh:
        text = fh.read()
    shutil.rmtree(work, ignore_errors=True)
    return {"key": key, "ok": True, "text": text, "job": job}


def git_head(path: str) -> str:
    try:
        return subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return "?"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--modxna", default=os.path.join(ROOT, "external", "modxna"))
    ap.add_argument("--cpptraj", default=os.path.join(ROOT, "external", "cpptraj", "bin"))
    ap.add_argument("--amberhome", default=os.environ.get("AMBERHOME", "/opt/amber26"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--only", nargs="*",
                    help="build only these keys (BB-SUGAR-BASE-int|5cap|3cap); keys outside "
                         "the curated set are added to the library")
    a = ap.parse_args(argv)

    modxna_sh = os.path.join(a.modxna, "modxna.sh")
    env = dict(os.environ)
    env["AMBERHOME"] = a.amberhome
    env["PATH"] = os.pathsep.join([a.cpptraj, os.path.join(a.amberhome, "bin"), env.get("PATH", "")])
    env["CPPTRAJ"] = os.path.join(a.cpptraj, "cpptraj")
    env["OMP_NUM_THREADS"] = "1"
    # modxna.sh rewrites mol2 lines with awk printf; under a locale with a
    # decimal comma (e.g. es_UY) the 3'-cap O3' line comes out as "-0,0000"
    # and the sugar fragment is lost. Force the C locale.
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    version = subprocess.run([env["CPPTRAJ"], "--version"], capture_output=True, text=True,
                             env=env).stdout.strip()

    out_dir = os.path.join(ROOT, "web", "lib", "residues")
    os.makedirs(out_dir, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="modxna-", dir=os.path.join(ROOT, "build"))
    jobs = curated()
    if a.only:
        jobs = []
        for key in a.only:
            parts = key.split("-")
            if len(parts) != 4 or parts[3] not in CHARGE:
                sys.exit(f"bad key {key!r}: expected BB-SUGAR-BASE-int|5cap|3cap")
            jobs.append(tuple(parts))
    print(f"{len(jobs)} residues, {a.jobs} parallel jobs, {version}")

    index_path = os.path.join(ROOT, "data", "residue_index.json")
    index = {"residues": {}}
    if os.path.exists(index_path) and a.only:
        with open(index_path, encoding="utf-8") as fh:
            index = json.load(fh)
    failures = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(build_one, j, modxna_sh, env, scratch) for j in jobs]
        for n, fut in enumerate(cf.as_completed(futs), 1):
            r = fut.result()
            if not r["ok"]:
                failures.append(r)
                print(f"  FAIL {r['key']}")
                continue
            unit = parse_off(r["text"])[UNIT]
            bb, sugar, base, cap = r["job"]
            q = round(unit.net_charge, 6)
            if abs(q - CHARGE[cap]) > 1e-3:
                failures.append({"key": r["key"], "log": f"net charge {q} != {CHARGE[cap]}"})
                print(f"  FAIL {r['key']} charge {q}")
                continue
            with open(os.path.join(out_dir, r["key"] + ".lib"), "w", encoding="utf-8") as fh:
                fh.write(r["text"])
            index["residues"][r["key"]] = {
                "triplet": [bb, sugar, base], "cap": cap, "file": r["key"] + ".lib",
                "sha256": hashlib.sha256(r["text"].encode()).hexdigest(),
                "charge": q, "heavy": unit.heavy_atoms(), "unit": UNIT,
                "head": unit.atom_names[unit.head] if unit.head is not None else None,
                "tail": unit.atom_names[unit.tail] if unit.tail is not None else None}
            if n % 20 == 0:
                print(f"  {n}/{len(jobs)}")

    frcmod_src = os.path.join(a.modxna, "dat", "frcmod.modxna")
    frcmod_dst = os.path.join(ROOT, "web", "lib", "frcmod.modxna")
    shutil.copyfile(frcmod_src, frcmod_dst)
    with open(frcmod_dst, "rb") as fh:
        frcmod_sha = hashlib.sha256(fh.read()).hexdigest()
    index.update({
        "schema": "dsxnabuilder/residue_index/1",
        "unit": UNIT,
        "modxna_commit": git_head(a.modxna),
        "modxna_version": _modxna_version(modxna_sh),
        "cpptraj": version, "cpptraj_commit": git_head(os.path.dirname(a.cpptraj)),
        "amberhome": a.amberhome, "frcmod_sha256": frcmod_sha,
        "generator": "tools/build_residue_library.py"})
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(index, fh, indent=1, sort_keys=True)
    print(f"{len(index['residues'])} residues in {index_path}; {len(failures)} failures")
    for f in failures:
        print(f"--- {f['key']}\n{f.get('log', '')[-800:]}")
    shutil.rmtree(scratch, ignore_errors=True)
    return 1 if failures else 0


def _modxna_version(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("VERSION="):
                return line.split("=", 1)[1].strip().strip("'\"")
    return "?"


if __name__ == "__main__":
    sys.exit(main())
