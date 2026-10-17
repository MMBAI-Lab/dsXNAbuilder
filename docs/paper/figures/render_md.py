# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Figure 4c renders: final MD frame over the model built by dsXNAbuilder.

    pymol -cq docs/paper/figures/render_md.py -- build/md

The last frame of prod.nc (solvent and ions stripped by cpptraj) is superposed
on the package's structure.pdb, the builder output with the helix axis along z,
using the heavy atoms of the inner base pairs, which is the same core as the
RMSD of analyze_md.py. The built model is drawn as a translucent grey cartoon
under the MD frame.
"""
import json
import os
import subprocess
import sys

import numpy as np
from pymol import cmd

# Run from the repository root (PyMOL does not set __file__ to this script).
D = "docs/paper/figures"
CPPTRAJ = os.path.join("external", "cpptraj", "bin", "cpptraj")
# Same style as render.py, copied: importing render.py would rerun all its renders.
STRAND = {"A": "0x3A6EA5", "B": "0xC8553D"}


def style(obj, highlight_mods=True):
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ladder_mode", 1)
    cmd.set("cartoon_ring_transparency", 0.0)
    cmd.set("cartoon_tube_radius", 0.9)
    for ch, col in STRAND.items():
        cmd.color(col, f"{obj} and chain {ch}")
        cmd.set("cartoon_color", col, f"{obj} and chain {ch}")
    if highlight_mods:
        mods = f"{obj} and (resn X* or resn Y* or resn 3DR)"
        cmd.show("sticks", mods + " and not name P+OP1+OP2+O5'")
        cmd.color("0x2E9E5B", mods + " and elem C")
        cmd.color("red", mods + " and elem O")
        cmd.color("0x3050F8", mods + " and elem N")
        cmd.color("0x90E050", mods + " and elem F")
        cmd.show("spheres", f"{obj} and elem S")
        cmd.set("sphere_scale", 0.45)
        cmd.color("0xE6C200", f"{obj} and elem S")
        cmd.set("stick_radius", 0.18)




def read_pdb(path):
    """Residues in file order: list of (chain, resseq, resname, {atom: xyz})."""
    res, key = [], None
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        k = (line[21], line[22:27], line[17:20].strip())
        if k != key:
            res.append((k[0], k[1].strip(), k[2], {}))
            key = k
        res[-1][3][line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]),
                                                    float(line[46:54])])
    return res


def kabsch(P, Q):
    """R, t such that R @ p + t fits P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def final_frame(pkg, out):
    man = json.load(open(os.path.join(pkg, "manifest.json")))
    nA, nB = len(man["chains"]["A"]), len(man["chains"]["B"])
    tmp = out + ".raw.pdb"
    script = (f"parm {pkg}/system.parm7\ntrajin {pkg}/prod.nc lastframe\nautoimage\n"
              f"strip !:1-{nA + nB}\ntrajout {tmp} pdb\nrun\n")
    p = subprocess.run([CPPTRAJ], input=script, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(p.stdout[-2000:] + p.stderr[-2000:])
    built, md = read_pdb(os.path.join(pkg, "structure.pdb")), read_pdb(tmp)
    os.remove(tmp)
    assert len(built) == len(md) == nA + nB, (len(built), len(md), nA + nB)
    # inner paired residues, as in analyze_md.py
    colsA = {r["column"]: k for k, r in enumerate(man["chains"]["A"])}
    colsB = {r["column"]: nA + k for k, r in enumerate(man["chains"]["B"])}
    paired = sorted(set(colsA) & set(colsB))[1:-1]
    core = [colsA[c] for c in paired] + [colsB[c] for c in paired]
    P, Q = [], []
    for i in core:
        for name, xyz in built[i][3].items():
            if not name.startswith("H") and name in md[i][3]:
                P.append(md[i][3][name])
                Q.append(xyz)
    R, t = kabsch(np.array(P), np.array(Q))
    rmsd = np.sqrt(np.mean(np.sum((np.array(P) @ R.T + t - np.array(Q)) ** 2, axis=1)))
    # residue identifiers of the built model, so both objects share chains/numbers
    n = 0
    with open(out, "w") as fh:
        for (ch, num, rn, _), (_, _, _, atoms) in zip(built, md):
            for name, xyz in atoms.items():
                if name.startswith("H"):
                    continue
                n += 1
                x, y, z = R @ xyz + t
                fh.write(f"ATOM  {n:5d} {name:<4s} {rn:>3s} {ch}{int(num):4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {name[0]:>2s}\n")
        fh.write("END\n")
    return rmsd


def render(name, pkg, out, width=700, height=1100):
    snap = os.path.join(D, "pdb", f"md_{name}.pdb")
    rmsd = final_frame(pkg, snap)
    cmd.reinitialize()
    cmd.load(os.path.join(pkg, "structure.pdb"), "built")
    cmd.load(snap, "md")
    style("md")
    cmd.hide("everything", "built")
    cmd.show("cartoon", "built")
    cmd.color("grey70", "built")
    cmd.set("cartoon_color", "grey70", "built")
    cmd.set("cartoon_transparency", 0.6, "built")
    cmd.hide("everything", "elem H")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_color", "0x404040")
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.2)
    cmd.set("orthoscopic", 1)
    cmd.reset()
    cmd.turn("x", -90)          # helix axis (z) vertical, strand 1 5'->3' upwards
    cmd.zoom("built or md", buffer=1.0, complete=1)
    cmd.png(out, width=width, height=height, dpi=300, ray=1)
    print(f"{name}: final frame vs built, core heavy-atom RMSD {rmsd:.2f} A -> {out}")


mddir = sys.argv[1] if len(sys.argv) > 1 else "build/md"
for name in ("dickerson", "gapmer", "sirna"):
    pkg = os.path.join(mddir, name, "dsxnabuilder")
    if os.path.exists(os.path.join(pkg, "prod.nc")):
        render(name, pkg, os.path.join(D, f"render_md_{name}.png"))
