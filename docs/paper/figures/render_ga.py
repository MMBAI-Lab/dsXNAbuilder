"""Render for the graphical abstract (an image not used in any figure).

    .venv/bin/python docs/paper/figures/render_ga.py      # builds the PDB
    pymol -cq docs/paper/figures/render_ga.py             # renders it

A 5-10-5 MOE/DNA phosphorothioate gapmer bound to its RNA target on the
DNA:RNA hybrid model #22, drawn with the helix axis horizontal.
"""
import os
import sys

D = "docs/paper/figures"
PDB = f"{D}/pdb/ga_moe_gapmer.pdb"
REQ = {"top": "r(GCAAUGUACUGAAUGCAUGC)",
       "bottom": "eC*eG*eT*eT*eA*dC*dA*dT*dG*dA*dC*dT*dT*dA*dC*eG*eT*eA*eC*eG", "model": 22}

try:
    from pymol import cmd
except ImportError:                                   # build step (plain Python)
    sys.path.insert(0, ".")
    from dsxnabuilder.reslib import ResidueLibrary
    from dsxnabuilder.service import Session
    res = Session(residue_library=ResidueLibrary.load_default()).build(REQ)
    assert res["ok"] and res["summary"]["clashes"] == 0, res["issues"]
    open(PDB, "w").write(res["pdb"])
    print(PDB, res["summary"]["n_modified"], "modified nucleotides")
    sys.exit(0)

STRAND = {"A": "0x3A6EA5", "B": "0xC8553D"}
cmd.reinitialize()
cmd.load(PDB, "ga")
cmd.hide("everything", "ga")
cmd.show("cartoon", "ga")
cmd.set("cartoon_ring_mode", 3)
cmd.set("cartoon_ladder_mode", 1)
cmd.set("cartoon_tube_radius", 0.9)
for ch, col in STRAND.items():
    cmd.set("cartoon_color", col, f"ga and chain {ch}")
mods = "ga and resn X* and not name P+OP1+OP2+O5'"
cmd.show("sticks", mods)
cmd.color("0x2E9E5B", mods + " and elem C")
cmd.color("red", mods + " and elem O")
cmd.color("0x3050F8", mods + " and elem N")
cmd.show("spheres", "ga and elem S")
cmd.set("sphere_scale", 0.45)
cmd.color("0xE6C200", "ga and elem S")
cmd.set("stick_radius", 0.2)
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_trace_color", "0x404040")
cmd.set("ray_shadows", 0)
cmd.set("specular", 0.2)
cmd.set("orthoscopic", 1)
cmd.reset()
cmd.turn("x", -90)          # helix axis vertical...
cmd.turn("z", 90)           # ...then horizontal
cmd.turn("x", 25)           # tilt to open the grooves
cmd.zoom("ga", buffer=0.5, complete=1)
cmd.png(f"{D}/render_ga.png", width=2400, height=1050, dpi=300, ray=1)
