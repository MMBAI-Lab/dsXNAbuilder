# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Analysis of the use-case MD runs (run from the repository root).

    .venv/bin/python docs/paper/md/analyze_md.py build/md [--out docs/paper/md/results]

For each system: heavy-atom RMSD to the structure built by dsXNAbuilder
(terminal base pairs excluded), base-pair and step parameters (cpptraj
nastruct; modXNA residues mapped to their parent base with resmap), sugar
pseudorotation of every nucleotide (multipucker) and Watson-Crick pair
occupancy. Writes a JSON summary and time series used by Figure 4.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CPPTRAJ = os.path.join(ROOT, "external", "cpptraj", "bin", "cpptraj")
LETTER = {"DAA": "A", "RAA": "A", "DCC": "C", "RCC": "C", "M5C": "C", "DGG": "G", "RGG": "G",
          "DTT": "T", "RUU": "U"}


def read_table(path):
    rows, header = [], None
    for line in open(path):
        if line.startswith("#"):
            header = line[1:].split()
            continue
        if line.strip():
            rows.append(line.split())
    return header, rows


def analyze(sysdir, outdir, stride):
    pkg = os.path.abspath(os.path.join(sysdir, "dsxnabuilder"))
    man = json.load(open(os.path.join(pkg, "manifest.json")))
    chains = man["chains"]
    nA, nB = len(chains["A"]), len(chains["B"])
    n = nA + nB
    resmap = []
    for ch in ("A", "B"):
        for r in chains[ch]:
            if r["resname"] not in ("DA", "DC", "DG", "DT", "A", "C", "G", "U"):
                resmap.append(f"{r['resname']}:{LETTER[r['triplet'][2]]}")
    resmap = sorted(set(resmap))
    # paired columns (exclude overhangs) and terminal pairs for RMSD
    colsA = {r["column"]: k + 1 for k, r in enumerate(chains["A"])}
    colsB = {r["column"]: nA + k + 1 for k, r in enumerate(chains["B"])}
    paired = sorted(set(colsA) & set(colsB))
    inner = paired[1:-1]
    core = sorted([colsA[c] for c in inner] + [colsB[c] for c in inner])
    mask = ":" + ",".join(map(str, core)) + "&!@H="
    # The pairs are known by construction. nastruct's default (geometric search
    # on the first frame) paired two stacked neighbours of the gapmer strand.
    bplist = ",".join(f"{colsA[c]}-{colsB[c]}" for c in paired)
    tag = os.path.basename(os.path.abspath(sysdir))
    o = os.path.abspath(os.path.join(outdir, tag))
    os.makedirs(o, exist_ok=True)
    script = f"""parm {pkg}/system.parm7
reference {pkg}/system.rst7 name built
trajin {pkg}/prod.nc 1 last {stride}
autoimage
rms fit {mask} ref built out {o}/rmsd.dat mass time {10 * stride}
nastruct NA {' '.join('resmap ' + m for m in resmap)} specifiedbp pairs {bplist} naout {o}/na.dat
multipucker P resrange 1-{n} out {o}/pucker.dat range360
run
"""
    with open(os.path.join(o, "analysis.cpptraj"), "w") as fh:
        fh.write(script)
    p = subprocess.run([CPPTRAJ, "-i", os.path.join(o, "analysis.cpptraj")], capture_output=True,
                       text=True, cwd=o)
    if p.returncode != 0:
        sys.exit(p.stdout[-3000:] + p.stderr[-2000:])
    # --- summaries
    _, rows = read_table(f"{o}/rmsd.dat")
    t = np.array([float(r[0]) for r in rows]) / 1000.0      # ns
    rmsd = np.array([float(r[1]) for r in rows])
    h, steps = read_table(f"{o}/BPstep.na.dat")
    hi = {k: i for i, k in enumerate(h)}
    stepvals = {k: np.array([float(r[hi[k]]) for r in steps]) for k in
                ("Shift", "Slide", "Rise", "Tilt", "Roll", "Twist")}
    # nastruct writes one row per pair and frame even when the pair is broken,
    # so occupancy comes from its BP flag and H-bond count, not from row counts.
    hb, bps = read_table(f"{o}/BP.na.dat")
    bi = {k: i for i, k in enumerate(hb)}
    nframes = len(rows)
    npairs_expected = len(paired)
    per_pair = {}
    for r in bps:
        per_pair.setdefault(f"{r[bi['Base1']]}-{r[bi['Base2']]}", []).append(
            (int(r[bi["BP"]]), int(r[bi["HB"]])))
    occ = {k: {"bp": float(np.mean([b for b, _ in v])), "hb2": float(np.mean([h >= 2 for _, h in v]))}
           for k, v in per_pair.items()}
    hp, pk = read_table(f"{o}/pucker.dat")
    pucker = {}
    for j, name in enumerate(hp[1:], start=1):
        vals = np.array([float(r[j]) for r in pk])
        pucker[name] = vals
    # classify residues
    classes = {}
    for ch, off in (("A", 0), ("B", nA)):
        for k, r in enumerate(chains[ch]):
            sug = r["triplet"][1]
            classes[off + k + 1] = {"DC2": "DNA", "RC3": "RNA", "OME": "2'-OMe", "AF2": "2'-F",
                                    "MOE": "MOE", "LNA": "LNA"}.get(sug, sug)
    by_class = {}
    for name, vals in pucker.items():
        resnum = int("".join(c for c in name.split(":")[-1] if c.isdigit()) or 0) if ":" in name else None
        if resnum is None:
            import re
            m = re.search(r"(\d+)", name)
            resnum = int(m.group(1)) if m else None
        cl = classes.get(resnum)
        if cl:
            by_class.setdefault(cl, []).append(vals)
    north = {cl: float(np.mean([(v < 90) | (v > 270) for vs in arr for v in vs]))
             for cl, arr in by_class.items()}
    # conventional pseudorotation quadrants N 315-45, E 45-135, S 135-225, W 225-315
    quadrants = {}
    for cl, arr in by_class.items():
        v = np.concatenate(arr) % 360
        quadrants[cl] = {q: float(np.mean(((v - lo) % 360) < 90))
                         for q, lo in (("N", 315), ("E", 45), ("S", 135), ("W", 225))}
    summary = {
        "system": tag, "frames": nframes, "time_ns": float(t[-1]) if len(t) else 0.0,
        "rmsd_mean": float(rmsd[len(rmsd) // 5:].mean()), "rmsd_sd": float(rmsd[len(rmsd) // 5:].std()),
        "rmsd_max": float(rmsd.max()),
        "paired_columns": npairs_expected, "pairs_tracked": len(per_pair),
        "bp_fraction": float(np.mean([v["bp"] for v in occ.values()])),
        "hb2_fraction": float(np.mean([v["hb2"] for v in occ.values()])),
        "hb2_fraction_min_pair": float(min(v["hb2"] for v in occ.values())),
        "pair_occupancy": occ,
        "steps_mean": {k: float(np.nanmean(v)) for k, v in stepvals.items()},
        "north_fraction_by_class": north,
        "pucker_quadrants_by_class": quadrants,
        "resmap": resmap,
    }
    np.savez(f"{o}/series.npz", t=t, rmsd=rmsd,
             **{f"pucker_{cl}": np.concatenate(arr) for cl, arr in by_class.items()})
    json.dump(summary, open(f"{o}/summary.json", "w"), indent=1)
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mddir")
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "paper", "md", "results"))
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    for s in a.only or ("dickerson", "gapmer", "sirna"):
        d = os.path.join(a.mddir, s)
        if os.path.exists(os.path.join(d, "dsxnabuilder", "prod.nc")):
            print(json.dumps(analyze(d, a.out, a.stride), indent=1))
