"""Figures 1-3 and the graphical abstract of the dsXNAbuilder paper.

Run from the repository root after gallery_build.py, render.py (PyMOL),
screenshot.py and scripts/bench_browser.py:

    .venv/bin/python docs/paper/figures/make_figures.py

Every number plotted is computed here from the package, its data files and
the cpptraj reference outputs in tests/data/cpptraj.
"""
import json
import math
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
B = "docs/paper/data"
os.makedirs(B, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300,
                     # editable SVG: text stays <text>, not glyph paths
                     "svg.fonttype": "none"})
C_STRAND1, C_STRAND2, C_MOD, C_GREY = "#3A6EA5", "#C8553D", "#2E9E5B", "#6b6b6b"


def panel(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="right")


def save(fig, name, **kw):
    """PNG and PDF for the manuscript, SVG for editing (fixed metadata, no date)."""
    fig.savefig(f"{D}/{name}.png", **kw)
    fig.savefig(f"{D}/{name}.pdf", **kw)
    fig.savefig(f"{D}/{name}.svg", metadata={"Date": None}, **kw)


def show(ax, path):
    ax.imshow(mpimg.imread(path))
    ax.axis("off")


# --------------------------------------------------------------------------- data
RESLIB = ResidueLibrary.load_default()
session = Session(residue_library=RESLIB)
lib = session.fibers
raw = json.load(open("data/fiber_models.json"))["models"]
meta = json.load(open(f"{D}/pdb/meta.json"))


PAIR = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}


def canonical_strands(m, n, rng):
    """Random Watson-Crick strands with the chemistry of the model strands."""
    top = "".join(rng.choice("ACGT") for _ in range(n))
    top = top.replace("T", "U") if m.strand_chem("I") == "RNA" else top
    bot = "".join(PAIR[c] for c in top)
    bot = bot.replace("T", "U") if m.strand_chem("II") == "RNA" else bot.replace("U", "T")
    return top, bot


def built_twist(mid, m, rng):
    """Helical twist per bp of a built duplex: C1' azimuth of strand 1 around
    the helix axis (z), averaged over a whole number of repeats so that the
    variation within a multi-bp repeat cancels."""
    top, bot = canonical_strands(m, m.k * max(1, 24 // m.k) + 1, rng)
    res = session.build({"top": top, "bottom": bot, "model": mid})
    c1 = [a.xyz for a in parse_pdb(res["pdb"]) if a.name == "C1'" and a.chain == "A"]
    az = np.degrees(np.arctan2([x[1] for x in c1], [x[0] for x in c1]))
    return float(np.mean((np.diff(az) + 180.0) % 360.0 - 180.0))


BONDS = {"C2'-F2' (2'-F)": ("AF2", "C2'", "F2'"), "O2'-CM2 (2'-OMe)": ("OME", "O2'", "CM2"),
         "C2'-O2' (LNA)": ("LNA", "C2'", "O2'"), "O2'-C6' (LNA)": ("LNA", "O2'", "C6'"),
         "C6'-C4' (LNA)": ("LNA", "C6'", "C4'"), "C5-C8 (5mC)": ("M5C", "C5", "C8"),
         "P-S (PS1)": ("PS1", "P", "S1"), "P-S (PS2)": ("PS2", "P", "S2")}


def mod_bond_lengths(res, out):
    """Bonds made by transplanted atoms in a built duplex, keyed as BONDS. The
    modXNA triplet of each X.. residue comes from the PDB header (REMARK 999)."""
    trip = {}
    for line in res["pdb"].splitlines():
        if line.startswith("REMARK 999") and "=" in line:
            f = line[10:].split()
            if f and f[0].startswith(("X", "Y")) and len(f[0]) == 3:
                for tok in f[1:]:
                    if tok.count("/") == 2 or tok.count("-") == 2:
                        trip[f[0]] = tok.replace("/", "-").split("-")
                        break
    by = {}
    for a in parse_pdb(res["pdb"]):
        by.setdefault((a.chain, a.resseq, a.resname), {})[a.name] = np.asarray(a.xyz)
    for (_, _, rn), at in by.items():
        t = trip.get(rn)
        if not t:
            continue
        for k, (frag, x, y) in BONDS.items():
            if frag in t and x in at and y in at:
                out.setdefault(k, []).append(float(np.linalg.norm(at[x] - at[y])))


def library_bond_lengths():
    """The same bonds in the prebuilt modXNA residues (reference values)."""
    out = {}
    for key in RESLIB.residues:
        c = RESLIB.unit(key).coords()
        t = key.split("-")
        for k, (frag, x, y) in BONDS.items():
            if frag in t and x in c and y in c:
                out.setdefault(k, []).append(float(np.linalg.norm(c[x] - c[y])))
    return out


def validation_data():
    out = {}
    # (a) helical twist per bp: built duplexes vs 3DNA table (README twist / repeat)
    rng_a = random.Random(3)
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
        tw.append((mid, m.family, table, built_twist(mid, m, rng_a), m.bp_rise, m.handedness,
                   n != 0))
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
    o3p_can, o3p_mod, nbuild, bonds = [], [], 0, {}
    clashes = {"canonical": 0, "modified": 0, "canonical_models": set()}
    for mid, m in sorted(lib.models.items()):
        for rep in range(5):
            top, bot = canonical_strands(m, rng.randint(10, 24), rng)
            res = session.build({"top": top, "bottom": bot, "model": mid})
            if res["ok"]:
                nbuild += 1
                o3p_can += [d for *_, d in o3p_distances(parse_pdb(res["pdb"]))]
                if res["summary"]["clashes"]:
                    clashes["canonical"] += 1
                    clashes["canonical_models"].add(mid)
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
            mod_bond_lengths(res, bonds)
            clashes["modified"] += res["summary"]["clashes"] > 0
    out["o3p"] = (np.array(o3p_can), np.array(o3p_mod), nbuild)
    clashes["canonical_models"] = sorted(clashes["canonical_models"])
    out["builds_with_contacts"] = clashes            # contacts < 2.2 A between residues
    ref = library_bond_lengths()
    out["mod_bonds"] = {k: {"n": len(v), "built": [min(v), max(v)],
                            "library": [min(ref[k]), max(ref[k])],
                            "max_dev": max(abs(x - np.mean(ref[k])) for x in v)}
                        for k, v in sorted(bonds.items())}
    return out


V = validation_data()
json.dump({"twist": V["twist"],
           "twist_max_abs_diff": max(abs((b - t + 180.0) % 360.0 - 180.0)
                                     for _, _, t, b, *_ in V["twist"]),
           "cpptraj_max_abs": {k: float(max(np.abs(a - b).max() for a, b in V["cpptraj"][k]))
                               for k in ("pair", "step")},
           "o3p": {"canonical": [float(V["o3p"][0].min()), float(V["o3p"][0].max()), len(V["o3p"][0])],
                   "modified": [float(V["o3p"][1].min()), float(V["o3p"][1].max()), len(V["o3p"][1])],
                   "builds": V["o3p"][2]},
           "mod_bonds": V["mod_bonds"],
           "builds_with_contacts": V["builds_with_contacts"]},
          open(f"{B}/validation_summary.json", "w"), indent=1)


# --------------------------------------------------------------------------- Figure 1
def figure1():
    W, H = 7.1, 7.6
    fig = plt.figure(figsize=(W, H))
    # a: workflow
    ax = fig.add_axes([0.0, 0.89, 1.0, 0.10]); ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.axis("off")
    # the complete flow of the code is Supplementary Figure S1 (make_flowchart.py)
    steps = [("Both strands\ntyped 5′→3′\nover 3′→5′", "#eef2f7"),
             ("Parsed to\nmodXNA triplets;\nchecks", "#eef2f7"),
             ("Fiber model:\nscrew\nplacement", "#eef2f7"),
             ("Base swap,\nO2′, abasic\nsites", "#eef2f7"),
             ("modXNA\ntransplant\n(shared atoms)", "#e9f5ee"),
             ("Hydrogens;\n2′-O chains\nout of contacts", "#e9f5ee"),
             ("PDB (+ optional\nAmber/modXNA)", "#fbeee9")]
    w, gap = 12.0, 2.3
    for k, (txt, col) in enumerate(steps):
        x = 1.5 + k * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 1.5), w, 7, boxstyle="round,pad=0.3,rounding_size=1.2",
                                    fc=col, ec="#8a8a8a", lw=0.6))
        ax.text(x + w / 2, 5, txt, ha="center", va="center", fontsize=6.5, linespacing=1.15)
        if k < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.35, 5), (x + w + gap - 0.35, 5),
                                         arrowstyle="-|>", mutation_scale=7, lw=0.8, color="#555"))
    ax.text(0.0, 1.0, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
    # left column: b (input and pairing grid, screenshots of the app) over d (PDB header);
    # each image at its own aspect ratio so that no gap opens between them
    img_in = mpimg.imread(f"{D}/screen_input_lna.png")
    img_grid = mpimg.imread(f"{D}/screen_grid_lna.png")
    wl, top = 0.60, 0.875
    h_in = 0.40 * W * img_in.shape[0] / img_in.shape[1] / H
    h_grid = wl * W * img_grid.shape[0] / img_grid.shape[1] / H
    axb1 = fig.add_axes([0.0, top - h_in, 0.40, h_in]); axb1.imshow(img_in); axb1.axis("off")
    y_grid = top - h_in - 0.008 - h_grid
    axb2 = fig.add_axes([0.0, y_grid, wl, h_grid]); axb2.imshow(img_grid); axb2.axis("off")
    axb1.text(0.0, 1.0, "b", transform=axb1.transAxes, fontsize=11, fontweight="bold", va="top")
    # d: PDB header
    axd = fig.add_axes([0.0, 0.0, wl, y_grid - 0.02]); axd.axis("off")
    axd.text(0.0, 1.0, "d", transform=axd.transAxes, fontsize=11, fontweight="bold", va="top")
    # the header as the current version writes it (the render in c is of the same duplex)
    g = meta["gapmer"]
    pdb = session.build({"top": g["top"], "bottom": g["bottom"], "model": g["model"]})["pdb"]
    pdb = pdb.splitlines()
    rem = [l for l in pdb if l.startswith("REMARK")]
    body = [l[11:] for l in rem]                     # header lines as written (<= 80 col)
    pick = []
    for i, l in enumerate(body):
        if l.startswith(("dsXNAbuilder", "strand", "fiber model", "helical", "modified")) or \
                (body[i - 1].startswith(("strand 2", "fiber model #")) and l.startswith("    ")):
            pick.append(l)
        elif l.lstrip().startswith(("X01", "X02")) or (i and body[i - 1].lstrip().startswith(("X01", "X02"))):
            pick.append(l)
    tail = [i for i, l in enumerate(body) if l.startswith("modXNA ")][0]
    lines = pick + ["  …"] + body[tail:] + \
        [l[:54] for l in pdb if l.startswith("ATOM")][:4] + ["…"]
    axd.text(0.06, 0.93, "header of the output PDB (REMARK 999 records)", transform=axd.transAxes,
             fontsize=7, va="bottom")
    axd.text(0.012, 0.88, "\n".join(lines), family="DejaVu Sans Mono", fontsize=7, va="top",
             ha="left", transform=axd.transAxes, linespacing=1.2,
             bbox=dict(boxstyle="round,pad=0.5", fc="#f7f7f7", ec="#bbbbbb", lw=0.5))
    # c: 3D model, right column
    wc = 0.37
    hc = wc * W * 1100 / 700 / H
    yc = 0.5 * (top - hc) + 0.03
    axc = fig.add_axes([1.0 - wc, yc, wc, hc]); show(axc, f"{D}/render_gapmer.png")
    axc.text(0.0, 1.0, "c", transform=axc.transAxes, fontsize=11, fontweight="bold", va="top")
    axc.text(0.5, -0.01, "RNA (blue) · LNA/DNA gapmer (red)\nmodXNA residues in sticks, S as spheres",
             transform=axc.transAxes, ha="center", va="top", fontsize=7)
    save(fig, "figure1", bbox_inches="tight")
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
        ax.text(0.5, -0.02, f"{m.bp_twist:+.1f}°".replace("-", "−") + f"\n{m.bp_rise:.2f} Å",
                transform=ax.transAxes,
                ha="center", va="top", fontsize=7)
    fig.subplots_adjust(wspace=0.02, left=0.01, right=0.99, top=0.86, bottom=0.12)
    save(fig, "figure2", bbox_inches="tight")
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
               label="repeat rotation mod 360°", zorder=3)
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
    ax.set_xlabel("dsXNAbuilder (°)"); ax.set_ylabel("CPPTRAJ nastruct, 3DNA convention (°)")
    dev = max(float(max(np.abs(a - b).max() for a, b in V["cpptraj"][k])) for k in ("pair", "step"))
    ax.text(0.98, 0.04, f"8 models, all 12 parameters\n(translations not shown)\nmax |Δ| = {dev / 10 ** math.floor(math.log10(dev)):.1f} × 10$^{{{math.floor(math.log10(dev))}}}$ (° or Å)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    ax.text(0.02, 0.98, "angles of steps (○) and pairs (□)", transform=ax.transAxes,
            fontsize=7, va="top")
    ax.legend(frameon=False, fontsize=7, loc="center left", handletextpad=0.1)
    panel(ax, "b")
    # c: O3'-P
    ax = axs[1, 0]
    can, mod, nb = V["o3p"]
    bins = np.linspace(1.55, 1.65, 41)
    ax.hist(can, bins=bins, color=C_STRAND1, alpha=0.75,
            label=f"canonical: 200 duplexes, 40 models ({len(can)} bonds)")
    ax.hist(mod, bins=bins, color=C_MOD, alpha=0.75,
            label=f"modified: 60 duplexes, 8 C3′-endo models ({len(mod)} bonds)")
    ax.set_xlabel("O3′–P distance (Å)"); ax.set_ylabel("count")
    ax.set_xlim(1.55, 1.65)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.text(0.98, 0.70, f"{nb} random duplexes,\nno relaxation", transform=ax.transAxes,
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
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
    panel(ax, "d")
    fig.tight_layout()
    save(fig, "figure3", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- graphical abstract
def graphical_abstract():
    """NAR: 5:2 landscape, at least 127 x 50 mm, sans-serif text of 12-16 pt, an
    image not used in the figures (render_ga.py), TIF/EPS/PDF named graphical_abstract."""
    from PIL import Image
    img = mpimg.imread(f"{D}/render_ga.png")[..., :3]
    ink = np.where(img.min(axis=2) < 0.97)                # crop the white margin
    img = img[ink[0].min():ink[0].max() + 1, ink[1].min():ink[1].max() + 1]
    fig = plt.figure(figsize=(10.0, 4.0), facecolor="white")
    ax0 = fig.add_axes([0.02, 0.05, 0.33, 0.9]); ax0.axis("off")
    ax0.text(0.0, 0.97, "dsXNAbuilder", fontsize=16, fontweight="bold", va="top")
    ax0.text(0.0, 0.80, "Type both strands", fontsize=12, va="top")
    ax0.text(0.0, 0.66, "5′-r(GCAAUGUACUGAAUGCAUGC)-3′\n3′-eC*eG*eT*eT*eA*dC…eC*eG-5′",
             family="DejaVu Sans Mono", fontsize=12, va="top", linespacing=1.4,
             bbox=dict(boxstyle="round,pad=0.45", fc="#f4f6f9", ec="#aab4c3", lw=0.8))
    ax0.text(0.0, 0.30, "2′-OMe · 2′-F · MOE · LNA\nphosphorothioates · 5mC\nabasic sites · mismatches",
             fontsize=12, color=C_MOD, va="top", linespacing=1.35)
    axm = fig.add_axes([0.355, 0.08, 0.42, 0.84]); axm.imshow(img); axm.axis("off")
    axm.text(0.5, 1.0, "MOE gapmer : RNA on a DNA:RNA hybrid fiber model", fontsize=12,
             ha="center", va="bottom", transform=axm.transAxes)
    axr = fig.add_axes([0.80, 0.05, 0.19, 0.9]); axr.axis("off")
    axr.text(0.0, 0.88, "40 fiber models\nfrom 3DNA\n\nmodXNA chemistry\non either strand\n\n"
             "PDB for any program\n\nruns in the browser", fontsize=12, va="top", linespacing=1.3)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{D}/graphical_abstract.{ext}", dpi=300, facecolor="white",
                    **({"metadata": {"Date": None}} if ext == "svg" else {}))
    plt.close(fig)
    im = Image.open(f"{D}/graphical_abstract.png").convert("RGB")      # no alpha channel
    im.save(f"{D}/graphical_abstract.png", dpi=(300, 300))
    im.save(f"{D}/graphical_abstract.tif", dpi=(300, 300), compression="tiff_lzw")

figure1(); figure2(); figure3(); graphical_abstract()
print(json.dumps(json.load(open(f"{B}/validation_summary.json"))["cpptraj_max_abs"]),
      json.load(open(f"{B}/validation_summary.json"))["o3p"])
