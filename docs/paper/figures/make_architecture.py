# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Supplementary Figure S2: architecture of dsXNAbuilder.

    .venv/bin/python docs/paper/figures/make_architecture.py

(a) Where each part runs and resides, and how the parts communicate.
(b) The modules of the dsxnabuilder package and the calls between them. Every
arrow is an import found in the source (checked here with ast), except the
dotted one: build uses the residue library through the object it is given.
Writes figureS2.{svg,png,pdf} next to this script.
"""
import ast
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

D = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(D, "..", "..", ".."))
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "savefig.dpi": 300})
MONO, EDGE, GREY = "DejaVu Sans Mono", "#8a8a8a", "#555555"
C_ENTRY, C_SERVICE, C_CORE, C_CHEM, C_OUT, C_BASE, C_DATA = (
    "#f3f3f3", "#eef2f7", "#eef2f7", "#e9f5ee", "#fbeee9", "#ffffff", "#fff8e1")


def imports():
    """module -> package modules it imports (relative imports in dsxnabuilder/)."""
    pkg = os.path.join(ROOT, "dsxnabuilder")
    out = {}
    for f in os.listdir(pkg):
        if f.endswith(".py"):
            tree = ast.parse(open(os.path.join(pkg, f)).read())
            out[f[:-3]] = {n.module for n in ast.walk(tree)
                           if isinstance(n, ast.ImportFrom) and n.level == 1 and n.module}
    return out


def region(ax, x, y, w, h, title, fc, ec=EDGE, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0",
                                fc=fc, ec=ec, lw=0.7, ls=ls))
    ax.text(x + 0.8, y + h - 0.6, title, fontsize=6.8, fontweight="bold", va="top")


def body(ax, x, y, text, size=5.9, family=None):
    ax.text(x, y, text, fontsize=size, va="top", ha="left", linespacing=1.3,
            **({"family": family} if family else {}))


def module(ax, x, y, w, h, name, role, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=0.6",
                                fc=fc, ec=EDGE, lw=0.6))
    ax.text(x + w / 2, y + h - 0.55, name, fontsize=6.3, family=MONO, fontweight="bold",
            ha="center", va="top")
    ax.text(x + w / 2, y + h - 2.05, role, fontsize=5.4, ha="center", va="top", linespacing=1.2,
            color="#333333")


def arrow(ax, p, q, dashed=False, dotted=False, color=GREY, rad=0.0, both=False, lw=0.6):
    ls = (0, (1, 1.2)) if dotted else (0, (3, 2)) if dashed else "-"
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="<|-|>" if both else "-|>", mutation_scale=6,
                                 lw=lw, color=color, ls=ls, shrinkA=0, shrinkB=0,
                                 connectionstyle=f"arc3,rad={rad}"))


def label(ax, x, y, text, ha="center", size=5.5):
    ax.text(x, y, text, fontsize=size, color=GREY, ha=ha, va="center", linespacing=1.2,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))


def panel_a(ax):
    ax.text(0.2, 102.4, "a", fontsize=11, fontweight="bold", va="top")
    # repository and offline tools; command line
    region(ax, 1.5, 80.5, 21.5, 19.0, "Repository (git)", C_DATA, ls=(0, (4, 2)))
    body(ax, 2.3, 97.0,
         "dsxnabuilder/  the Python package\n"
         "data/  JSON tables: fiber models,\n   catalogue, aliases, pairing,\n"
         "   Amber units, residue index\n"
         "web/lib/  244 modXNA .lib +\n   frcmod.modxna\n"
         "web/  page, app.js, worker.js\n"
         "tools/  make data/ and web/lib/\n   from 3DNA, modXNA, CPPTRAJ,\n   AmberTools (offline)\n"
         "tests/, scripts/  checks, build,\n   deployment")
    region(ax, 1.5, 62.5, 21.5, 15.0, "Command line (local)", C_ENTRY)
    body(ax, 2.3, 75.0,
         "CPython ≥ 3.10 + NumPy\n"
         "dsxnabuilder → cli.main\n→ service.Session (+ batch)\n"
         "reads data/ and web/lib/ from\nthe checkout (paths.py)\n"
         "writes a PDB, a zip or a\ndirectory; the Amber package\nruns with AmberTools (tleap)")
    # web server and the two browser contexts
    region(ax, 27.0, 89.5, 45.5, 10.0, "Static web server  (danslab.xyz/dsxnabuilder/)", C_SERVICE)
    body(ax, 27.8, 97.0,
         "index.html, style.css, app.js, worker.js\nvendor/ngl.js  (NGL 2.4.0)\n"
         "VERSION  (package version + commit)")
    body(ax, 49.8, 97.0,
         "pyodide/  (Pyodide 314.0.7: CPython +\n   NumPy in WebAssembly, ~16 MB)\n"
         "dsxnabuilder.zip  (package + 6 tables)\n"
         "lib/residues/*.lib, lib/frcmod.modxna")
    region(ax, 27.0, 62.5, 19.0, 21.0, "Browser: main thread", C_ENTRY)
    body(ax, 27.8, 81.0,
         "app.js\n  input, pairing grid,\n  issues, tabs, high-\n  throughput list,\n"
         "  downloads (Blob),\n  reproducible URL\n"
         "NGL viewer (WebGL)\n"
         "fetches lib/residues/\n  <key>.lib for the\n  residues a build needs\n  (cached per page)")
    region(ax, 53.5, 62.5, 19.0, 21.0, "Browser: Web Worker", C_CHEM)
    body(ax, 54.3, 81.0,
         "worker.js\n  loads Pyodide and\n  NumPy, unpacks\n  dsxnabuilder.zip in\n"
         "  /home/pyodide\n"
         "webapi.call(name, json)\n  → service.Session\n"
         "residue files kept by\n  ResidueLibrary.provide\n"
         "zips returned base64")
    # arrows
    arrow(ax, (23.3, 95.0), (26.7, 95.0))
    label(ax, 25.0, 92.0, "build_web.sh\ndeploy_\ndanslab.sh", size=5.0)
    arrow(ax, (12.2, 80.2), (12.2, 77.8))
    arrow(ax, (36.5, 89.2), (36.5, 83.8))
    label(ax, 41.3, 86.5, "page, scripts;\n.lib on demand", size=5.0)
    arrow(ax, (63.0, 89.2), (63.0, 83.8))
    label(ax, 67.4, 86.5, "Pyodide, zip,\nVERSION", size=5.0)
    arrow(ax, (46.3, 73.6), (53.2, 73.6), both=True)
    ax.text(49.75, 83.0, "calls:\nmodels\nvalidate\nbuild\ncomplement\npieces\nbatch_\nvalidate\n"
            "batch_zip", fontsize=4.9, color="#333333", ha="center", va="top", family=MONO,
            linespacing=1.1)
    ax.text(49.75, 72.4, "postMessage\n(JSON)\n→ {id, name,\n   payload}\n← {id, ok,\n   result}\n"
            "← {status}\n   while\n   loading", fontsize=4.9, color=GREY, ha="center", va="top",
            linespacing=1.15)


def panel_b(ax, imp):
    ax.text(0.2, 59.6, "b", fontsize=11, fontweight="bold", va="top")
    ax.text(3.0, 59.3, "dsxnabuilder package: the main calls between modules "
            "(Supplementary Table S6 lists every import)", fontsize=6.3, va="top", color="#333333")
    boxes = {}

    def put(name, x, y, w, h, role, fc):
        module(ax, x, y, w, h, name, role, fc)
        boxes[name] = (x, y, w, h)

    H = 5.6
    # entry points
    put("webapi", 4.0, 49.5, 16.0, H, "JSON dispatcher of the worker;\ncreates ResidueLibrary", C_ENTRY)
    put("cli", 27.0, 49.5, 16.0, H, "dsxnabuilder models | build;\ncreates ResidueLibrary", C_ENTRY)
    put("analysis", 55.0, 49.5, 17.0, H, "base-pair and step parameters\n(tests and paper figures)", C_BASE)
    # service
    put("service", 9.0, 40.0, 26.0, H, "Session: loads the tables once;\nvalidate · build · batch_validate · batch_zip",
        C_SERVICE)
    put("batch", 44.0, 40.0, 15.0, H, "FASTA-like lists;\nzip, report.tsv", C_SERVICE)
    # pipeline and outputs
    put("build", 3.0, 30.5, 18.0, H, "build_duplex: the steps\nof Supplementary Fig. S1", C_CORE)
    put("pdbout", 27.0, 30.5, 15.0, H, "PDB file with\nREMARK 999 header", C_OUT)
    put("amberio", 48.0, 30.5, 15.0, H, "Amber/modXNA package,\nmanifest, reproduce", C_OUT)
    # the steps
    w, g, y = 8.3, 0.65, 16.0
    row = [("reslib", "modXNA lib,\ntransplant;\nlib/*.lib,\nfrcmod"),
           ("seqlang", "notation →\nResidueSpec;\naliases.json"),
           ("naming", "residue\nnames X01…,\nuser names"),
           ("fibers", "fiber models,\nscrew;\nfiber_models.json"),
           ("nucleotide", "base swap,\nO2′, chi,\npair kinds"),
           ("hydrogens", "H from\ntemplates;\namber_\nresidues.json"),
           ("substituents", "2′-O torsion\nscan"),
           ("fragments", "library tab,\ncatalogue\nnames")]
    for k, (name, role) in enumerate(row):
        fc = C_CHEM if name in ("hydrogens", "substituents", "reslib") else C_CORE
        put(name, 1.0 + k * (w + g), y, w, 8.0, role, fc)
    # foundations
    ax.add_patch(FancyBboxPatch((1.0, 0.8), 71.2, 9.6, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc="#fafafa", ec=EDGE, lw=0.5, ls=(0, (4, 2))))
    ax.text(1.8, 9.8, "foundations, imported by the modules above", fontsize=5.6, color="#333333",
            va="top")
    base = [("geometry", "frames, CEHS,\nKabsch, pucker"), ("structure", "Atom, PDB I/O,\nbase frames"),
            ("graph", "bond-graph\nwalks"), ("offlib", "Amber OFF\n.lib files"),
            ("catalog", "modXNA codes;\nmodxna_catalog"), ("issues", "Issue: code,\nmessage"),
            ("paths", "data/, web/lib/\nlocations")]
    bw, bg = 9.4, 0.7
    for k, (name, role) in enumerate(base):
        put(name, 1.6 + k * (bw + bg), 1.5, bw, 5.6, role, C_BASE)

    def top(n, f=0.5):
        x, y, w_, h = boxes[n]
        return (x + f * w_, y + h + 0.25)

    def bottom(n, f=0.5):
        x, y, w_, h = boxes[n]
        return (x + f * w_, y - 0.25)

    def left(n, f=0.5):
        x, y, w_, h = boxes[n]
        return (x - 0.25, y + f * h)

    def right(n, f=0.5):
        x, y, w_, h = boxes[n]
        return (x + w_ + 0.25, y + f * h)

    drawn = set()

    def edge(a, b, p, q, **kw):
        assert b in imp[a], f"{a} does not import {b}"
        drawn.add((a, b))
        arrow(ax, p, q, **kw)

    edge("webapi", "service", bottom("webapi", 0.6), top("service", 0.25))
    edge("cli", "service", bottom("cli", 0.4), top("service", 0.8))
    edge("cli", "batch", bottom("cli", 0.85), top("batch", 0.3))
    edge("service", "batch", right("service"), left("batch"))
    edge("service", "build", bottom("service", 0.2), top("build", 0.6))
    edge("service", "pdbout", bottom("service", 0.75), top("pdbout", 0.4))
    edge("service", "amberio", bottom("service", 0.95), top("amberio", 0.3))
    edge("amberio", "pdbout", left("amberio"), right("pdbout"))
    edge("pdbout", "build", left("pdbout"), right("build"))
    edge("pdbout", "fragments", bottom("pdbout", 0.8), top("fragments", 0.4))
    for k, name in enumerate(["seqlang", "naming", "fibers", "nucleotide", "hydrogens", "substituents"]):
        edge("build", name, bottom("build", 0.35 + 0.11 * k), top(name))
    # used, not imported: Session hands build the ResidueLibrary that webapi or cli made
    assert "reslib" not in imp["build"]
    arrow(ax, bottom("build", 0.1), top("reslib", 0.6), dotted=True)
    ax.text(0.6, 27.6, "install() on the\nobject Session\npasses", fontsize=4.8, color=GREY,
            ha="left", va="center", style="italic", linespacing=1.1)
    return drawn


def main():
    imp = imports()
    fig = plt.figure(figsize=(7.3, 10.25))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 73)
    ax.set_ylim(0, 102.5)
    ax.axis("off")
    panel_a(ax)
    ax.plot([1.0, 72.0], [61.2, 61.2], color="#cccccc", lw=0.5)
    panel_b(ax, imp)
    for ext in ("svg", "png", "pdf"):
        kw = {"metadata": {"Date": None}} if ext == "svg" else {}
        fig.savefig(os.path.join(D, f"figureS2.{ext}"), facecolor="white", **kw)
    plt.close(fig)
    print(os.path.join(D, "figureS2.png"))


if __name__ == "__main__":
    main()
