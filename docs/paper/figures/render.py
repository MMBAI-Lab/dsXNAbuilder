# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""PyMOL renders for the figures:  pymol -cq docs/paper/figures/render.py"""
import os
from pymol import cmd

D = "docs/paper/figures"
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


def render(name, out, width, height, turn_y=0):
    cmd.reinitialize()
    cmd.load(f"{D}/pdb/{name}.pdb", name)
    style(name)
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
    if turn_y:
        cmd.turn("y", turn_y)
    cmd.zoom(name, buffer=1.0, complete=1)
    cmd.png(out, width=width, height=height, dpi=300, ray=1)


for name in ("b_dna", "a_dna", "c_dna", "d_dna", "z_dna", "s_dna", "a_rna", "hybrid"):
    render(name, f"{D}/render_{name}.png", 520, 1000)
for name in ("gapmer", "sirna", "abasic"):
    render(name, f"{D}/render_{name}.png", 700, 1100, turn_y=0)
