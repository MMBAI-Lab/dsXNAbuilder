# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Supplementary Figure S1: flow of the dsXNAbuilder code.

    .venv/bin/python docs/paper/figures/make_flowchart.py

Writes figureS1.{svg,png,pdf} next to this script. Boxes name the functions
of the dsxnabuilder package that run each step; the right column shows the
data files they read and the offline tools that make them.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

D = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "savefig.dpi": 300})

STEP, CHEM, OUT, DATA, ENTRY, ERR = "#eef2f7", "#e9f5ee", "#fbeee9", "#fff8e1", "#f3f3f3", "#fdecea"
EDGE, MONO = "#8a8a8a", "DejaVu Sans Mono"

# main column, right column (data), units of 0.1 inch
XL, XR, W, WR = 2.0, 49.0, 43.0, 22.0
H, GAP = 6.4, 1.6


def box(ax, x, y, w, h, title, module, detail, fc, ec=EDGE, ls="-"):
    """Rounded box: bold title, module in monospace, detail lines."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc=fc, ec=ec, lw=0.6, ls=ls))
    lines = [(title, dict(fontsize=7.3, fontweight="bold"))]
    if module:
        lines.append((module, dict(fontsize=6.2, family=MONO, color="#4a4a4a")))
    if detail:
        lines.append((detail, dict(fontsize=6.3, color="#222222", linespacing=1.25)))
    ax.text(x + 0.9, y + h - 0.75, lines[0][0], va="top", ha="left", **lines[0][1])
    yy = y + h - 2.2
    for text, kw in lines[1:]:
        ax.text(x + 0.9, yy, text, va="top", ha="left", **kw)
        yy -= 1.35 * (text.count("\n") + 1)


def arrow(ax, p, q, dashed=False, label=None, lpos=0.5, color="#555555", rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, lw=0.7, color=color,
                                 ls=(0, (3, 2)) if dashed else "-",
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))
    if label:
        x = p[0] + lpos * (q[0] - p[0])
        y = p[1] + lpos * (q[1] - p[1])
        ax.text(x, y, label, fontsize=6, color=color, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))


def main():
    fig = plt.figure(figsize=(7.3, 10.3))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 73)
    ax.set_ylim(0, 103)
    ax.axis("off")

    # entry points, the high-throughput list and the session
    box(ax, XL, 95.3, 20.5, H, "Web app", "app.js → worker.js → webapi",
        "Pyodide in a Web Worker;\none duplex or up to 50", ENTRY)
    box(ax, XL + 22.5, 95.3, 20.5, H, "Command line", "cli.main",
        "dsxnabuilder build;\none duplex or up to 500", ENTRY)
    box(ax, XL, 85.8, W, H, "Session", "service.Session",
        "loads the catalogue, aliases and fiber models once;\n"
        "validate · build · batch_validate · batch_zip", STEP)
    box(ax, XR, 85.8, WR, 101.7 - 85.8, "High-throughput list", "batch.parse_batch",
        "one record per duplex:\n>name=…; model=…\nstrand 1 (5′→3′)\nstrand 2 (3′→5′)\n"
        "the format is checked first,\nthen one build per record\n\n"
        "web: up to 50 records;\ncommand line: up to 500", ENTRY)
    for x in (XL + 10.25, XL + 32.75):
        arrow(ax, (x, 95.0), (x, 92.5))
    arrow(ax, (XL + 43.3, 98.5), (XR - 0.3, 98.5))
    ax.text(XL + 45.0, 99.3, "list", fontsize=6, color="#555555", ha="center", va="bottom")
    arrow(ax, (XR - 0.3, 89.0), (XL + W + 0.3, 89.0))

    # main pipeline
    steps = [
        ("1  Parse both strands", "seqlang.parse_duplex",
         "tokens → modXNA backbone/sugar/base triplets; strand 2 reversed;\n"
         "a linkage belongs to its 3′ nucleotide; ‘-’ = overhang or abasic site", STEP),
        ("2  Model, orientation and phase", "fibers.FiberLibrary · build.choose_orientation",
         "hybrids built on the strands that match their chemistry; repeat\n"
         "phase fitted to the purine/pyrimidine pattern; sequence checks", STEP),
        ("3  Scaffold nucleotides", "fibers.FiberModel.nucleotide · nucleotide",
         "slot i mod k, moved by the screw operation ⌊i/k⌋ times; base swap\n"
         "in the standard frame; O2′ added or removed; pucker and syn checks", STEP),
        ("4  modXNA residues", "reslib.ResidueLibrary.install → transplant",
         "shared atoms kept; new atoms by local superposition; phosphorus rule\n"
         "(PS1/PS2); ring closures relaxed; stereo, bond and contact checks", CHEM),
        ("5  Chains and residue names", "build._chains · naming.ResidueNamer",
         "chain A = strand 1, chain B = strand 2 (5′→3′); 5′ phosphate removed;\n"
         "modified residues X01, X02, … or names chosen by the user", STEP),
        ("6  Hydrogens", "hydrogens.add_hydrogens",
         "from the unit tleap will match (Amber or modXNA), by superposition\n"
         "anchored on the bonded atom; heavy atoms do not move", CHEM),
        ("7  Flexible 2′-O substituents", "substituents.orient_2o_chains",
         "2′-OMe and MOE torsions scanned out of contacts (2.9 Å; 2.0 Å with H);\n"
         "each chain sees the chains already turned", CHEM),
        ("8  Orientation and final checks", "build._orient_and_centre · _check_result",
         "flipped hybrids turned by Rx(180°), then centred;\n"
         "O3′–P within 1.45–1.75 Å; inter-residue contacts < 2.2 Å", STEP),
    ]
    tops = []
    y = 83.0
    for i, (title, module, detail, fc) in enumerate(steps):
        box(ax, XL, y - H, W, H, title, module, detail, fc)
        tops.append(y)
        y -= H + (6.6 if i == 3 else GAP)
    cx = XL + W / 2
    arrow(ax, (cx, 85.5), (cx, tops[0] + 0.3))
    for i in range(len(steps) - 1):
        if i != 3:
            arrow(ax, (cx, tops[i] - H - 0.3), (cx, tops[i + 1] + 0.3))
    # after the transplant: any error stops the build
    yd = tops[3] - H - 3.3
    ax.add_patch(Polygon([(cx, yd + 2.3), (cx + 5.2, yd), (cx, yd - 2.3), (cx - 5.2, yd)],
                         closed=True, fc="white", ec=EDGE, lw=0.6))
    ax.text(cx, yd, "errors?", fontsize=6.5, ha="center", va="center")
    arrow(ax, (cx, tops[3] - H - 0.3), (cx, yd + 2.4))
    arrow(ax, (cx, yd - 2.4), (cx, tops[4] + 0.3))
    ax.text(cx + 1.2, yd - 2.0, "no", fontsize=6, color="#555555", ha="left", va="center")
    box(ax, XR, yd - 4.0, WR, 5.6, "Build stops", None,
        "issues returned with stable\ncodes (Table S5); no structure", ERR, ec="#c0392b")
    arrow(ax, (cx + 5.3, yd), (XR - 0.3, yd))
    ax.text((cx + 5.3 + XR) / 2, yd + 0.6, "yes", fontsize=6, color="#555555", ha="center",
            va="bottom")

    # data files (right column) and the offline tools that make them
    data = [
        (0, 0, 1.5, "Notation and catalogue", "data/aliases.json",
         "base_pairing.json, modxna_catalog.json\ntools/build_data.py, build_catalog.py"),
        (1, 2, 0.0, "40 fiber models", "data/fiber_models.json",
         "slots, base templates, screw\noperation, puckers, sources;\n"
         "tools/derive_fiber_models.py\n← 3DNA fiber repeats"),
        (3, 3, 0.0, "244 modXNA residues", "web/lib/ + frcmod.modxna",
         "tools/build_residue_library.py\n← modxna.sh, CPPTRAJ"),
        (5, 5, 0.0, "Amber nucleotide units", "data/amber_residues.json",
         "tools/amber_residue_names.py\n← AmberTools DNA/RNA libs"),
    ]
    for first, last, extra, title, module, detail in data:
        top = tops[first] + extra
        bottom = tops[last] - H
        box(ax, XR, bottom, WR, top - bottom, title, module, detail, DATA, ls=(0, (4, 2)))
        for k in range(first, last + 1):
            yk = tops[k] - H / 2
            arrow(ax, (XR - 0.3, yk), (XL + W + 0.3, yk), dashed=True, color="#a08a3c")
    box(ax, XR, tops[6] - H, WR, H, "Offline checks", None,
        "pytest (with CPPTRAJ nastruct);\nsmoke_tleap.py: every residue\n"
        "through tleap/sander; e2e.py", "white", ls=(0, (1, 1.5)))

    # outputs
    yo, ho = 1.0, 6.9
    outs = [(XL, 20.5, "PDB file", "pdbout.structure_pdb",
             "all atoms; REMARK 999 header:\ninput, model, residues, data"),
            (XL + 22.5, 20.5, "Amber/modXNA package", "amberio.package",
             "leap.in, libs/, frcmod.modxna,\nrun.sh, manifest (reproduce)"),
            (XR, WR, "High-throughput zip", "batch.batch_files",
             "pdb/<name>.pdb, report.tsv,\ninput.txt, README.txt")]
    for x, w, title, module, detail in outs:
        box(ax, x, yo, w, ho, title, module, detail, OUT)
    y8 = tops[-1] - H - 0.3
    arrow(ax, (XL + 10.25, y8), (XL + 10.25, yo + ho + 0.3))
    arrow(ax, (XL + 32.75, y8), (XL + 32.75, yo + ho + 0.3), dashed=True)
    ax.text(XL + 33.6, (y8 + yo + ho) / 2, "optional", fontsize=6, color="#555555", va="center")
    ax.add_patch(FancyArrowPatch((XL + W + 0.3, tops[-1] - H / 2), (XR + WR / 2, yo + ho + 0.3),
                                 arrowstyle="-|>", mutation_scale=7, lw=0.7, color="#555555",
                                 ls=(0, (3, 2)), connectionstyle="angle,angleA=0,angleB=90,rad=0",
                                 shrinkA=0, shrinkB=0))
    ax.text(XR + WR / 2 + 0.9, (tops[-1] - H / 2 + yo + ho) / 2, "per record", fontsize=6,
            color="#555555", va="center")

    for ext in ("svg", "png", "pdf"):
        kw = {"metadata": {"Date": None}} if ext == "svg" else {}
        fig.savefig(os.path.join(D, f"figureS1.{ext}"), facecolor="white", **kw)
    plt.close(fig)
    print(os.path.join(D, "figureS1.png"))


if __name__ == "__main__":
    main()
