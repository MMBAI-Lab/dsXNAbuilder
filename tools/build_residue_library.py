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

The extended set adds every other fragment of the catalogue with canonical
partners, so that each one can be used on its own:
  base      with 2'-deoxyribose/phosphodiester and ribose/phosphodiester,
            internal, 5'- and 3'-capped (protected phosphoramidites excluded);
  sugar     with A C G T U and the phosphodiester of its pucker family,
            internal, 5'- and 3'-capped;
  backbone  with 2'-deoxyribose + dA dC dG dT and ribose + rA rC rG rU,
            internal and 3'-capped (a 5' end has no phosphate).
EXTRA lists specific combinations (the web examples).

Every residue built is then tried on a fiber scaffold (dsxnabuilder.reslib
transplant): those whose stereocentres, bonds or contacts cannot be reproduced
there are left out, like those modxna.sh cannot build (missing parameters);
both are recorded with the reason in residue_index.json ("unavailable").
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from typing import Optional  # noqa: E402

from dsxnabuilder.catalog import (CANONICAL_BACKBONE, CANONICAL_SUGAR, CAP_5_OH,  # noqa: E402
                                  Catalog)
from dsxnabuilder.fibers import load_library  # noqa: E402
from dsxnabuilder.fragments import base_families, test_context  # noqa: E402
from dsxnabuilder.offlib import parse_off  # noqa: E402
from dsxnabuilder.reslib import duplicate_names, tleap_name, transplant  # noqa: E402

UNIT = "MXR"
MODXNA = os.path.join(ROOT, "external", "modxna")
EXTRA = [("PS2", "2AD", "RYC", "int")]       # example "modXNA code" of the web app
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


def p_handedness(mol2_path: str) -> Optional[int]:
    """Sign of V(X, kept O, O5') about P in a backbone fragment that keeps one
    non-bridging oxygen and carries a new group X on P (None otherwise)."""
    import numpy as np
    atoms, bonds, sec = {}, [], None
    for line in open(mol2_path, errors="replace"):
        if line.startswith("@<TRIPOS>"):
            sec = line[9:].strip()
            continue
        f = line.split()
        if sec == "ATOM" and len(f) >= 6:
            atoms[f[0]] = (f[1], np.array([float(x) for x in f[2:5]]))
        elif sec == "BOND" and len(f) >= 3:
            bonds.append((f[1], f[2]))
    pid = [i for i, (n, _) in atoms.items() if n == "P"]
    if not pid:
        return None
    nb = [atoms[b if a == pid[0] else a][0] for a, b in bonds if pid[0] in (a, b)]
    kept = [o for o in ("OP1", "OP2") if o in nb]
    new = [n for n in nb if n not in ("O3'", "O5'", "OP1", "OP2")]
    if "O5'" not in nb or len(kept) != 1 or len(new) != 1:
        return None
    pos = {n: x for n, x in atoms.values()}
    p = pos["P"]
    return int(np.sign(np.dot(pos[new[0]] - p, np.cross(pos[kept[0]] - p, pos["O5'"] - p))))


def ambiguous_backbones(catalog, modxna_dir: str) -> dict:
    """Stereo pairs XX1/XX2 whose fragments have the same configuration at P.
    reslib._phosphorus_rule follows the geometry of the fragment, so the two
    would build the same diastereomer and at least one contradicts its name.
    (Phosphorothioates follow their names instead; see reslib.)"""
    hand = {}
    for f in catalog.by_moiety("backbone"):
        if f.code.startswith("PS"):
            continue
        path = os.path.join(modxna_dir, "dat", "lib_backbone", f.code + ".mol2")
        if os.path.exists(path):
            hand[f.code] = p_handedness(path)
    out = {}
    for code, h in hand.items():
        twin = code[:-1] + {"1": "2", "2": "1"}.get(code[-1], "")
        if code[-1] in "12" and h is not None and hand.get(twin) == h:
            out[code] = (f"{code} and {twin} have the same configuration at P in modXNA: "
                         "the diastereomer cannot be assigned")
    return out


def extended(catalog, review: dict) -> list:
    """Each catalogue fragment outside the curated set, with canonical partners."""
    have = {c for job in curated() for c in job[:3]}
    have |= set(ambiguous_backbones(catalog, MODXNA))
    protected = set(review.get("protected", {}).get("codes", []))
    out = []
    for f in catalog.by_moiety("base"):
        if f.code in have or f.code in protected:
            continue
        for bb, sugar in (("DPO", "DC2"), ("RPO", "RC3")):
            out += [(bb, sugar, f.code, "int"), (bb, sugar, f.code, "3cap"),
                    (CAP_5_OH, sugar, f.code, "5cap")]
    for f in catalog.by_moiety("sugar"):
        if f.code in have:
            continue
        north = f.pucker_domain == "North"
        bb = "RPO" if north else "DPO"
        bases = ["RAA", "RCC", "RGG", "RUU", "DTT"] if north else ["DAA", "DCC", "DGG", "DTT", "RUU"]
        for base in bases:
            out += [(bb, f.code, base, "int"), (bb, f.code, base, "3cap"),
                    (CAP_5_OH, f.code, base, "5cap")]
    for f in catalog.by_moiety("backbone"):
        if f.code in have or f.code == CAP_5_OH:
            continue
        for sugar, bases in (("DC2", ["DAA", "DCC", "DGG", "DTT"]),
                             ("RC3", ["RAA", "RCC", "RGG", "RUU"])):
            for base in bases:
                out += [(f.code, sugar, base, "int"), (f.code, sugar, base, "3cap")]
    return out


def validate_on_scaffold(job, text: str, catalog, families, fibers) -> Optional[str]:
    """None if the residue can be transplanted onto a fiber scaffold, else why not."""
    bb, sugar, base, cap = job
    letter, chem, model = test_context((bb, sugar, base), catalog, families)
    names, xyz, _ = fibers.models[model].nucleotide(3, "I", letter, chem, fibers.standard)
    if cap == "5cap":
        keep = [i for i, n in enumerate(names) if n not in ("P", "OP1", "OP2")]
        names, xyz = [names[i] for i in keep], xyz[keep]
    unit = parse_off(text)[UNIT]
    unit.atom_names = [tleap_name(n) for n in unit.atom_names]
    dup = duplicate_names(unit)
    if dup:       # tleap matches PDB atoms to the unit by name
        return f"duplicate atom names in the modXNA residue: {', '.join(dup)}"
    try:
        _, _, rep = transplant(unit, names, xyz)
    except Exception as exc:          # atoms that cannot be reached from the scaffold
        return f"transplant: {exc}"
    bad = ([f"configuration at {', '.join(rep['stereo'])}"] if rep["stereo"] else []) + \
        rep["distorted"][:2] + rep["contacts"][:2]
    return f"does not fit a fiber scaffold (model #{model}): {'; '.join(bad)}" if bad else None


def revalidate() -> int:
    catalog = Catalog.from_json(os.path.join(ROOT, "data", "modxna_catalog.json"))
    families = base_families()
    fibers = load_library(os.path.join(ROOT, "data", "fiber_models.json"))
    index_path = os.path.join(ROOT, "data", "residue_index.json")
    with open(index_path, encoding="utf-8") as fh:
        index = json.load(fh)
    index.setdefault("unavailable", {})
    out_dir = os.path.join(ROOT, "web", "lib", "residues")
    dropped = 0
    for key, r in sorted(index["residues"].items()):
        path = os.path.join(out_dir, r["file"])
        with open(path, encoding="utf-8") as fh:
            reason = validate_on_scaffold(tuple(r["triplet"]) + (r["cap"],), fh.read(),
                                          catalog, families, fibers)
        if reason:
            if tuple(key.split("-")) in set(curated()):
                print(f"  CURATED FAILS {key}: {reason}")
                return 1
            index["unavailable"][key] = reason
            del index["residues"][key]
            os.remove(path)
            dropped += 1
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(index, fh, indent=1, sort_keys=True)
    print(f"{dropped} dropped; {len(index['residues'])} residues, "
          f"{len(index['unavailable'])} unavailable")
    return 0


def leap_reason(log: str) -> str:
    """The first tleap/modxna.sh error worth showing (missing parameters, ...)."""
    for pat in (r"Could not find (?:angle|bond|torsion) parameter[^\n]*", r"No torsion terms[^\n]*",
                r"Error:[^\n]*"):
        m = re.search(pat, log)
        if m:
            return "modxna.sh: " + " ".join(m.group(0).split())
    return "modxna.sh failed"


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
        leap = os.path.join(work, "leap.log")
        log = open(leap, errors="replace").read() if os.path.exists(leap) else ""
        return {"key": key, "ok": False, "log": (p.stdout + p.stderr)[-3000:],
                "reason": leap_reason(log + p.stdout + p.stderr)}
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
    ap.add_argument("--set", choices=("curated", "extended"), default="extended",
                    help="curated: what the sequence language writes; extended (default): "
                         "plus every catalogue fragment with canonical partners and EXTRA")
    ap.add_argument("--validate-only", action="store_true",
                    help="re-run the scaffold check on the current library (after a change "
                         "to the transplant) and move the residues that fail to 'unavailable'")
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
    if a.validate_only:
        return revalidate()

    out_dir = os.path.join(ROOT, "web", "lib", "residues")
    os.makedirs(out_dir, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="modxna-", dir=os.path.join(ROOT, "build"))
    catalog = Catalog.from_json(os.path.join(ROOT, "data", "modxna_catalog.json"))
    with open(os.path.join(ROOT, "data", "modxna_catalog_review.json"), encoding="utf-8") as fh:
        review = json.load(fh)
    families = base_families()
    fibers = load_library(os.path.join(ROOT, "data", "fiber_models.json"))
    jobs = curated()
    skipped = {}
    if a.set == "extended":
        jobs += extended(catalog, review) + [j for j in EXTRA if j not in jobs]
        skipped = ambiguous_backbones(catalog, a.modxna)
    if a.only:
        jobs = []
        for key in a.only:
            parts = key.split("-")
            if len(parts) != 4 or parts[3] not in CHARGE:
                sys.exit(f"bad key {key!r}: expected BB-SUGAR-BASE-int|5cap|3cap")
            jobs.append(tuple(parts))
    print(f"{len(jobs)} residues, {a.jobs} parallel jobs, {version}")

    index_path = os.path.join(ROOT, "data", "residue_index.json")
    index = {"residues": {}, "unavailable": {}}
    if os.path.exists(index_path) and a.only:
        with open(index_path, encoding="utf-8") as fh:
            index = json.load(fh)
        index.setdefault("unavailable", {})
    else:
        for f in os.listdir(out_dir):          # a full build replaces the library
            if f.endswith(".lib"):
                os.remove(os.path.join(out_dir, f))
    for code, why in skipped.items():
        index["unavailable"][f"{code}-*-*-*"] = why
    failures = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(build_one, j, modxna_sh, env, scratch) for j in jobs]
        for n, fut in enumerate(cf.as_completed(futs), 1):
            r = fut.result()
            if not r["ok"]:
                failures.append(r)
                index["residues"].pop(r["key"], None)
                index["unavailable"][r["key"]] = r["reason"]
                continue
            unit = parse_off(r["text"])[UNIT]
            bb, sugar, base, cap = r["job"]
            q = round(unit.net_charge, 6)
            reason = None
            if abs(q - CHARGE[cap]) > 1e-3:
                reason = f"net charge {q} instead of {CHARGE[cap]}"
            else:
                reason = validate_on_scaffold(r["job"], r["text"], catalog, families, fibers)
            if reason:
                failures.append({"key": r["key"], "log": reason, "reason": reason})
                index["residues"].pop(r["key"], None)
                index["unavailable"][r["key"]] = reason
                continue
            index["unavailable"].pop(r["key"], None)
            with open(os.path.join(out_dir, r["key"] + ".lib"), "w", encoding="utf-8") as fh:
                fh.write(r["text"])
            index["residues"][r["key"]] = {
                "triplet": [bb, sugar, base], "cap": cap, "file": r["key"] + ".lib",
                "sha256": hashlib.sha256(r["text"].encode()).hexdigest(),
                "charge": q, "heavy": unit.heavy_atoms(), "unit": UNIT,
                "head": unit.atom_names[unit.head] if unit.head is not None else None,
                "tail": unit.atom_names[unit.tail] if unit.tail is not None else None}
            if n % 200 == 0:
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
    print(f"{len(index['residues'])} residues in {index_path}; "
          f"{len(index['unavailable'])} unavailable ({len(failures)} in this run)")
    reasons = {}
    for f in failures:
        k = re.sub(r"\d+\.\d+", "#", f["reason"].split(":")[0] + ":" +
                   ":".join(f["reason"].split(":")[1:2]))
        reasons[k] = reasons.get(k, 0) + 1
    for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {v:5d}  {k[:150]}")
    shutil.rmtree(scratch, ignore_errors=True)
    # fragments modxna.sh or the scaffold cannot take are expected in the extended
    # set; the curated set must build completely
    return 1 if any(tuple(f["key"].split("-")) in set(curated()) for f in failures) else 0


def _modxna_version(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("VERSION="):
                return line.split("=", 1)[1].strip().strip("'\"")
    return "?"


if __name__ == "__main__":
    sys.exit(main())
