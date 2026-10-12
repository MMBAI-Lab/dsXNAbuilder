"""Figure 4: the MD use case (run after docs/paper/md/analyze_md.py)."""
import json, os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

R = "docs/paper/md/results"
D = "docs/paper/figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300,
                     "svg.fonttype": "none"})        # editable SVG text
SYS = [("dickerson", "Dickerson–Drew B-DNA", "#3A6EA5"),
       ("gapmer", "LNA/PS gapmer : RNA", "#C8553D"),
       ("sirna", "2′-OMe/2′-F/PS siRNA", "#2E9E5B")]
LS = {"dickerson": "-", "gapmer": "--", "sirna": ":"}
CLASS_COL = {"DNA": "#3A6EA5", "RNA": "#8c564b", "LNA": "#C8553D", "2'-OMe": "#2E9E5B",
             "2'-F": "#9467bd", "MOE": "#e377c2"}


def main():
    fig = plt.figure(figsize=(7.1, 5.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25], hspace=0.38, wspace=0.32)
    ax = fig.add_subplot(gs[0, :2])
    summ = {}
    for name, label, col in SYS:
        f = f"{R}/{name}/series.npz"
        if not os.path.exists(f):
            continue
        s = np.load(f)
        summ[name] = json.load(open(f"{R}/{name}/summary.json"))
        ax.plot(s["t"], s["rmsd"], lw=0.4, color=col, alpha=0.35)
        w = 10                                        # 1-ns running mean (100-ps frames)
        ax.plot(s["t"][w - 1:], np.convolve(s["rmsd"], np.ones(w) / w, "valid"), lw=1.1,
                color=col, label=label)
    ax.set_xlabel("time (ns)"); ax.set_ylabel("RMSD to built model (Å)")
    ax.set_ylim(0, 6.2)
    ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.5)
    ax.text(-0.06, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold")
    # b: pseudorotation phase by sugar class
    axb = fig.add_subplot(gs[0, 2])
    bins = np.linspace(0, 360, 37)
    for name, _, _ in SYS:
        f = f"{R}/{name}/series.npz"
        if not os.path.exists(f):
            continue
        s = np.load(f)
        for k in s.files:
            if k.startswith("pucker_"):
                cl = k[7:]
                h, _ = np.histogram(s[k], bins=bins, density=True)
                axb.plot(0.5 * (bins[1:] + bins[:-1]), h, color=CLASS_COL.get(cl, "k"), lw=0.9,
                         ls=LS[name], label=f"{cl} ({name})")
    axb.axvspan(0, 90, color="#eeeeee", zorder=0); axb.axvspan(270, 360, color="#eeeeee", zorder=0)
    axb.set_xlabel("pseudorotation phase P (°)"); axb.set_ylabel("density")
    axb.set_xticks([0, 90, 180, 270, 360])
    axb.text(45, 1.0, "N", transform=axb.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)
    axb.text(180, 1.0, "S", transform=axb.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)
    axb.legend(frameon=False, fontsize=5.5, loc="upper right", handlelength=2.6)
    axb.text(-0.2, 1.02, "b", transform=axb.transAxes, fontsize=11, fontweight="bold")
    # c: snapshots
    for k, (name, label, col) in enumerate(SYS):
        img = f"{D}/render_md_{name}.png"
        axc = fig.add_subplot(gs[1, k])
        if os.path.exists(img):
            axc.imshow(mpimg.imread(img))
        axc.axis("off")
        axc.set_title(label, fontsize=7.5)
        if k == 0:
            axc.text(-0.2, 1.02, "c", transform=axc.transAxes, fontsize=11, fontweight="bold")
    fig.savefig(f"{D}/figure4.png", bbox_inches="tight"); fig.savefig(f"{D}/figure4.pdf", bbox_inches="tight")
    fig.savefig(f"{D}/figure4.svg", bbox_inches="tight", metadata={"Date": None})
    json.dump(summ, open(f"{R}/all_summaries.json", "w"), indent=1)


if __name__ == "__main__":
    main()
