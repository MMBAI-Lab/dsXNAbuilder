"""Figures 1-3 and the graphical abstract of the dsXNAbuilder paper.

Run from the repository root after gallery_build.py, render.py (PyMOL),
screenshot.py and scripts/bench_browser.py:

    .venv/bin/python docs/paper/figures/make_figures.py

Every number plotted is computed here from the package, its data files and
the cpptraj reference outputs in tests/data/cpptraj.
"""
import json
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

sys.path.insert(0, ".")
from dsxnabuilder.analysis import duplex_parameters
from dsxnabuilder.fibers import load_library
from dsxnabuilder.reslib import ResidueLibrary
from dsxnabuilder.service import Session
from dsxnabuilder.structure import o3p_distances, parse_pdb, read_pdb

D = "docs/paper/figures"
B = "docs/paper/build"
os.makedirs(B, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300})
C_STRAND1, C_STRAND2, C_MOD, C_GREY = "#3A6EA5", "#C8553D", "#2E9E5B", "#6b6b6b"


def panel(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="right")


def show(ax, path):
    ax.imshow(mpimg.imread(path))
    ax.axis("off")


# --------------------------------------------------------------------------- data
session = Session(residue_library=ResidueLibrary.load_default())
lib = session.fibers
raw = json.load(open("data/fiber_models.json"))["models"]
meta = json.load(open(f"{D}/pdb/meta.json"))


def validation_data():
    out = {}
    # (a) helical twist per bp: built models vs 3DNA table (README twist / repeat)
    tw = []
    for mid, m in sorted(lib.models.items()):
        r = raw[str(mid)]
        # 3DNA tabulates the rotation of the whole repeat of k bp, defined modulo
        # 360 (e.g. #3: 0 deg for 11 bp, #14: -90 deg for 6 bp); compare per bp
        # on the equivalent branch (|readme| + 360 n) / k with the sign of the
        # normalized screw
        k, unit = r["k"], r["twist"]
        n = round((k * m.bp_twist - unit) / 360.0)
        table = (unit + 360.0 * n) / k
        tw.append((mid, m.family, table, m.bp_twist, m.bp_rise, m.handedness, n != 0))
    out["twist"] = tw
    # (b) cpptraj comparison
    std = {k: dict(zip(v.names, v.xyz)) for k, v in lib.standard.items()}
    comp = {"pair": [], "step": []}
    for case in ["adna1", "arna20", "bdna4", "cdna7", "ddna14", "hyb22", "sdna49", "zdna15"]:
        atoms = read_pdb(f"tests/data/cpptraj/{case}.pdb")
        pairs, steps, _ = duplex_parameters(atoms, std)
        def table(p):
            return np.array([[float(x) for x in l.split()[3:9]] for l in open(p)
                             if l.strip() and not l.startswith("#")])
        comp["pair"].append((np.array([p.as_tuple() for p in pairs]),
                             table(f"tests/data/cpptraj/{case}.BP.out")))
        comp["step"].append((np.array([s.as_tuple() for s in steps]),
                             table(f"tests/data/cpptraj/{case}.BPstep.out")))
    out["cpptraj"] = comp
    # (c) O3'-P in random builds on every model, canonical and modified
    rng = random.Random(7)
    o3p_can, o3p_mod, nbuild = [], [], 0
    pair = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
    for mid, m in sorted(lib.models.items()):
        c1, c2 = m.strand_chem("I"), m.strand_chem("II")
        for rep in range(5):
            n = rng.randint(10, 24)
            top = "".join(rng.choice("ACGT") for _ in range(n))
            top = top.replace("T", "U") if c1 == "RNA" else top
            bot = "".join(pair[c] for c in top)
            bot = bot.replace("T", "U") if c2 == "RNA" else bot.replace("U", "T")
            res = session.build({"top": top, "bottom": bot, "model": mid})
            if res["ok"]:
                nbuild += 1
                o3p_can += [d for *_, d in o3p_distances(parse_pdb(res["pdb"]))]
    toks = ["mA", "mC", "mG", "mU", "fA", "fC", "fG", "fU", "eA", "eT", "eG", "C^5mC", "dA",
            "dC", "dG", "dT"]
    north = [mid for mid, m in lib.models.items()
             if (m.slots[0].I.pucker or {}).get("domain") == "North"
             and (m.slots[0].II.pucker or {}).get("domain") == "North"]
    for rep in range(60):
        mid = rng.choice(north)
        n = rng.randint(8, 20)
        t = [rng.choice(toks + ["+A", "+G", "+C", "+T"]) for _ in range(n)]
        letters = [x[-1] if "^" not in x else "C" for x in t]
        top = "".join(x + rng.choice(["", "*1", "*2"]) for x in t[:-1]) + t[-1]
        comp_ = "".join({"A": "U", "U": "A", "T": "A", "G": "C", "C": "G"}[c] for c in letters)
        res = session.build({"top": top, "bottom": f"r({comp_})", "model": mid})
        if res["ok"]:
            nbuild += 1
            o3p_mod += [d for *_, d in o3p_distances(parse_pdb(res["pdb"]))]
    out["o3p"] = (np.array(o3p_can), np.array(o3p_mod), nbuild)
    return out


V = validation_data()
json.dump({"twist": V["twist"],
           "cpptraj_max_abs": {k: float(max(np.abs(a - b).max() for a, b in V["cpptraj"][k]))
                               for k in ("pair", "step")},
           "o3p": {"canonical": [float(V["o3p"][0].min()), float(V["o3p"][0].max()), len(V["o3p"][0])],
                   "modified": [float(V["o3p"][1].min()), float(V["o3p"][1].max()), len(V["o3p"][1])],
                   "builds": V["o3p"][2]}},
          open(f"{B}/validation_summary.json", "w"), indent=1)


# --------------------------------------------------------------------------- Figure 1
def figure1():
    fig = plt.figure(figsize=(7.1, 5.2))
    # a: workflow
    ax = fig.add_axes([0.0, 0.83, 1.0, 0.15]); ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.axis("off")
    steps = [("Both strands typed\n5′→3′ over 3′→5′", "#eef2f7"),
             ("Parsed to modXNA\ntriplets; checks", "#eef2f7"),
             ("Fiber model:\nscrew placement", "#eef2f7"),
             ("Base swap, O2′,\nabasic sites", "#eef2f7"),
             ("modXNA transplant\n(shared atoms kept)", "#e9f5ee"),
             ("PDB (+ optional\nAmber/modXNA)", "#fbeee9")]
    w, gap = 14.0, 2.6
    for k, (txt, col) in enumerate(steps):
        x = 1.5 + k * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 1.5), w, 7, boxstyle="round,pad=0.3,rounding_size=1.2",
                                    fc=col, ec="#8a8a8a", lw=0.6))
        ax.text(x + w / 2, 5, txt, ha="center", va="center", fontsize=6.4)
        if k < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.35, 5), (x + w + gap - 0.35, 5),
                                         arrowstyle="-|>", mutation_scale=7, lw=0.8, color="#555"))
    ax.text(0.0, 1.0, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
    # b: input and pairing grid (screenshots of the app)
    img_in = mpimg.imread(f"{D}/screen_input_lna.png")
    img_grid = mpimg.imread(f"{D}/screen_grid_lna.png")
    axb1 = fig.add_axes([0.0, 0.50, 0.40, 0.31]); axb1.imshow(img_in); axb1.axis("off")
    axb2 = fig.add_axes([0.0, 0.06, 0.40, 0.42]); axb2.imshow(img_grid); axb2.axis("off")
    axb1.text(0.0, 1.0, "b", transform=axb1.transAxes, fontsize=11, fontweight="bold", va="top")
    # c: 3D model
    axc = fig.add_axes([0.40, 0.06, 0.24, 0.76]); show(axc, f"{D}/render_gapmer.png")
    axc.text(0.0, 1.0, "c", transform=axc.transAxes, fontsize=11, fontweight="bold", va="top")
    axc.text(0.5, -0.01, "RNA (blue) · LNA/DNA gapmer (red)\nmodXNA residues in sticks, S as spheres",
             transform=axc.transAxes, ha="center", va="top", fontsize=6)
    # d: PDB header
    axd = fig.add_axes([0.65, 0.04, 0.35, 0.78]); axd.axis("off")
    axd.text(0.0, 1.0, "d", transform=axd.transAxes, fontsize=11, fontweight="bold", va="top")
    pdb = open(f"{D}/pdb/gapmer.pdb").read().splitlines()
    rem = [l for l in pdb if l.startswith("REMARK")]
    body = [l[11:] for l in rem]                     # header lines as written (<= 80 col)
    pick = []
    for i, l in enumerate(body):
        if l.startswith(("dsXNAbuilder", "strand", "fiber model", "helical", "modified")) or \
                (body[i - 1].startswith(("strand 2", "fiber model #")) and l.startswith("    ")):
            pick.append(l)
        elif l.lstrip().startswith(("X01", "X02")) or (i and body[i - 1].lstrip().startswith(("X01", "X02"))):
            pick.append(l)
    lines = pick + ["  …"] + [body[-2], body[-1]] + \
        [l[:66] for l in pdb if l.startswith("ATOM")][:3] + ["…"]
    axd.text(0.02, 0.96, "\n".join(lines), family="DejaVu Sans Mono", fontsize=4.3, va="top",
             ha="left", transform=axd.transAxes, linespacing=1.3,
             bbox=dict(boxstyle="round,pad=0.6", fc="#f7f7f7", ec="#bbbbbb", lw=0.5))
    axd.text(0.09, 0.985, "header of the output PDB (REMARK 999 records)", transform=axd.transAxes,
             fontsize=6, va="bottom")
    fig.savefig(f"{D}/figure1.png", bbox_inches="tight")
    fig.savefig(f"{D}/figure1.pdf", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- Figure 2
def figure2():
    names = ["b_dna", "a_dna", "c_dna", "d_dna", "z_dna", "s_dna", "a_rna", "hybrid"]
    fig, axs = plt.subplots(1, 8, figsize=(7.1, 3.0))
    for ax, n in zip(axs, names):
        show(ax, f"{D}/render_{n}.png")
        mid = meta[n]["model"]
        m = lib.models[mid]
        ax.set_title(f"{meta[n]['label'].split(' (')[0]}\n#{mid}", fontsize=7, pad=2)
        ax.text(0.5, -0.02, f"{m.bp_twist:+.1f}°\n{m.bp_rise:.2f} Å", transform=ax.transAxes,
                ha="center", va="top", fontsize=6.3)
    fig.subplots_adjust(wspace=0.02, left=0.01, right=0.99, top=0.86, bottom=0.12)
    fig.savefig(f"{D}/figure2.png", bbox_inches="tight")
    fig.savefig(f"{D}/figure2.pdf", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- Figure 3
def figure3():
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.6))
    # a: twist
    ax = axs[0, 0]
    tw = V["twist"]
    x = np.array([t[2] for t in tw]); y = np.array([t[3] for t in tw])
    left = np.array([t[5] == "left" for t in tw])
    ax.plot([-40, 50], [-40, 50], color="#bbbbbb", lw=0.8, zorder=0)
    ax.scatter(x[~left], y[~left], s=14, color=C_STRAND1, label="right-handed", zorder=2)
    ax.scatter(x[left], y[left], s=14, color=C_STRAND2, label="left-handed (Z-DNA)", zorder=2)
    ax.set_xlabel("helical twist in 3DNA table (°/bp)")
    ax.set_ylabel("helical twist of built duplex (°/bp)")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    wrapped = np.array([t[6] for t in tw])
    ax.scatter(x[wrapped], y[wrapped], s=40, facecolors="none", edgecolors="k", lw=0.6,
               label="repeat rotation taken modulo 360°", zorder=3)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.text(0.98, 0.04, f"40 models\nmax |Δ| = {np.abs(x - y).max():.3f}°",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    panel(ax, "a")
    # b: cpptraj
    ax = axs[0, 1]
    names = ["shift", "slide", "rise", "tilt", "roll", "twist"]
    pnames = ["shear", "stretch", "stagger", "buckle", "propeller", "opening"]
    cols = plt.cm.tab10(np.arange(6))
    allx, ally = [], []
    for kind, labels, marker in (("step", names, "o"), ("pair", pnames, "s")):
        for k in range(6):
            xs = np.concatenate([a[:, k] for a, b in V["cpptraj"][kind]])
            ys = np.concatenate([b[:, k] for a, b in V["cpptraj"][kind]])
            if k >= 3:
                ax.scatter(xs, ys, s=7, marker=marker, color=cols[k], alpha=0.8,
                           label=f"{names[k]} / {pnames[k]}" if kind == "step" else None, lw=0)
            allx.append(xs); ally.append(ys)
    lim = [-60, 45]
    ax.plot(lim, lim, color="#bbbbbb", lw=0.8, zorder=0)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("dsXNAbuilder (°)"); ax.set_ylabel("cpptraj nastruct, 3DNA convention (°)")
    dev = max(float(max(np.abs(a - b).max() for a, b in V["cpptraj"][k])) for k in ("pair", "step"))
    ax.text(0.98, 0.04, f"8 models, all 12 parameters\n(translations not shown)\nmax |Δ| = {dev:.0e} (° or Å)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    ax.text(0.02, 0.98, "angles of steps (○) and pairs (□)", transform=ax.transAxes,
            fontsize=6.5, va="top")
    ax.legend(frameon=False, fontsize=6, loc="center left", handletextpad=0.1)
    panel(ax, "b")
    # c: O3'-P
    ax = axs[1, 0]
    can, mod, nb = V["o3p"]
    bins = np.linspace(1.55, 1.65, 41)
    ax.hist(can, bins=bins, color=C_STRAND1, alpha=0.75, label=f"canonical ({len(can)} bonds)")
    ax.hist(mod, bins=bins, color=C_MOD, alpha=0.75, label=f"modified ({len(mod)} bonds)")
    ax.set_xlabel("O3′–P distance (Å)"); ax.set_ylabel("count")
    ax.set_xlim(1.55, 1.65)
    ax.legend(frameon=False, fontsize=7)
    ax.text(0.98, 0.65, f"{nb} random duplexes\non all 40 models", transform=ax.transAxes,
            ha="right", fontsize=7)
    panel(ax, "c")
    # d: timings
    ax = axs[1, 1]
    bench = json.load(open(f"{B}/bench_browser.json"))
    for kind, col, lab in (("canonical", C_STRAND1, "canonical B-DNA"),
                           ("modified", C_MOD, "2′-OMe/2′-F/PS RNA")):
        ns = sorted(int(k) for k in bench[kind])
        med = [np.median(bench[kind][str(n)]["ms"]) / 1000 for n in ns]
        lo = [min(bench[kind][str(n)]["ms"]) / 1000 for n in ns]
        hi = [max(bench[kind][str(n)]["ms"]) / 1000 for n in ns]
        ax.errorbar(ns, med, yerr=[np.array(med) - lo, np.array(hi) - med], marker="o", ms=3.5,
                    color=col, lw=1, capsize=2, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks([12, 24, 48, 96, 192]); ax.set_xticklabels(["12", "24", "48", "96", "192"])
    ax.set_xlabel("duplex length (bp)"); ax.set_ylabel("time in the browser (s)")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.text(0.98, 0.04, f"first load of the app: {bench['load_s']:.1f} s\n(local server)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5)
    panel(ax, "d")
    fig.tight_layout()
    fig.savefig(f"{D}/figure3.png", bbox_inches="tight")
    fig.savefig(f"{D}/figure3.pdf", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- graphical abstract
def graphical_abstract():
    fig = plt.figure(figsize=(8.0, 3.2))
    ax0 = fig.add_axes([0.01, 0.15, 0.27, 0.7]); ax0.axis("off")
    ax0.text(0.0, 0.92, "Type both strands", fontsize=11, fontweight="bold", va="top")
    ax0.text(0.0, 0.68, "5′-r(GCAAUGUACUGAAUGC)-3′\n3′-+C*+G*+T*dT*dA*dC*dA*dT*\n   dG*dA*dC*dT*dT*+A*+C*+G-5′",
             family="DejaVu Sans Mono", fontsize=7.5, va="top",
             bbox=dict(boxstyle="round,pad=0.5", fc="#f4f6f9", ec="#aab4c3", lw=0.6))
    ax0.text(0.0, 0.12, "2′-OMe · 2′-F · MOE · LNA\nphosphorothioates · 5mC · abasic",
             fontsize=8, color=C_MOD, va="top")
    for k, (n, lab) in enumerate((("b_dna", "B-DNA"), ("z_dna", "Z-DNA"), ("a_rna", "A-RNA"),
                                  ("gapmer", "LNA gapmer : RNA"))):
        ax = fig.add_axes([0.31 + k * 0.125, 0.08, 0.12, 0.8])
        show(ax, f"{D}/render_{n}.png")
        ax.set_title(lab, fontsize=8)
    axr = fig.add_axes([0.82, 0.15, 0.17, 0.7]); axr.axis("off")
    axr.text(0.0, 0.92, "40 fiber models", fontsize=10, fontweight="bold", va="top")
    axr.text(0.0, 0.72, "exact screw symmetry\n244 modXNA residues\nPDB for any program\nAmber/modXNA topology\nruns in the browser",
             fontsize=8, va="top", linespacing=1.5)
    fig.savefig(f"{D}/graphical_abstract.png", dpi=300)
    plt.close(fig)


figure1(); figure2(); figure3(); graphical_abstract()
print(json.dumps(json.load(open(f"{B}/validation_summary.json"))["cpptraj_max_abs"]),
      json.load(open(f"{B}/validation_summary.json"))["o3p"])
