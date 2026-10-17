# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Hydrogens placed from the Amber and modXNA residue templates."""

import os

import numpy as np
import pytest

from dsxnabuilder.hydrogens import Template, add_hydrogens, amber_templates
from dsxnabuilder.structure import parse_pdb, residues

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HAS_LIB = os.path.isdir(os.path.join(ROOT, "web", "lib", "residues"))


def bonded_h(ats):
    """(H, heavy atom it is bonded to, distance) for every H of a residue."""
    heavy = [a for a in ats if a.element != "H"]
    out = []
    for h in (a for a in ats if a.element == "H"):
        d = [float(np.linalg.norm(h.xyz - a.xyz)) for a in heavy]
        k = int(np.argmin(d))
        out.append((h, heavy[k], d[k]))
    return out


def test_heavy_atoms_do_not_move_and_bonds_keep_template_length(lib):
    tpl = amber_templates()["DA"]
    names, xyz, _ = lib.models[4].nucleotide(3, "I", "A", "DNA", lib.standard)
    out_names, out_xyz = add_hydrogens(names, xyz, tpl)
    out = dict(zip(out_names, out_xyz))
    for n, p in zip(names, xyz):
        assert np.array_equal(out[n], p)
    assert set(out_names) == set(tpl.names)
    for h in tpl.hydrogens:
        parent = tpl.nb[h][0]
        assert abs(np.linalg.norm(out[h] - out[parent])
                   - np.linalg.norm(tpl.xyz[h] - tpl.xyz[parent])) < 1e-6, h
    # each hydrogen follows its heavy atom in the output
    for h in tpl.hydrogens:
        assert out_names[out_names.index(h) - 1] in [tpl.nb[h][0]] + tpl.hydrogens


@pytest.mark.parametrize("top,bottom,model", [
    ("CGCGAATTCGCG", "GCGCTTAAGCGC", 4),
    ("GGCGCUUGCGUC", "CCGCGAACGCAG", 20),
    ("CGCGCG", "GCGCGC", 15),
    ("r(GACUAGCUGA)", "d(CTGATCGACT)", 22),
    ("CGCGA-TTCGCG", "GCGCTTAAGCGC", 4),
])
def test_canonical_hydrogens(session, top, bottom, model):
    out = session.build({"top": top, "bottom": bottom, "model": model})
    assert out["ok"]
    atoms = parse_pdb(out["pdb"])
    assert out["summary"]["n_atoms"] == len(atoms)
    for ats in residues(atoms).values():
        for h, x, d in bonded_h(ats):
            assert 0.94 < d < 1.11, (h.resname, h.name, x.name, d)
    hs = np.array([a.xyz for a in atoms if a.element == "H"])
    dd = np.linalg.norm(hs[:, None] - hs[None], axis=2) + 9.0 * np.eye(len(hs))
    assert dd.min() > 1.5


@pytest.mark.skipif(not HAS_LIB, reason="residue library not built")
@pytest.mark.parametrize("top,bottom,model", [
    ("mG*mA*fC*mU*fA*mG*fC*mU*mG*mA", "CUGAUCGACU", 20),
    ("r(GACUAGCUGA)", "eC*eT*eG*dA*dT*dC*dG*eA*eC*eT", 22),     # MOE H6' hit H4' at 1.4 A
    ("r(GCAAUGUACUGAAUGC)", "+C*+G*+T*dT*dA*dC*dA*dT*dG*dA*dC*dT*dT*+A*+C*+G", 22),
    ("CGC^5mCGAA*2TTCGCG", "GCGCTT*1AAGCGC", 4),
])
def test_modxna_hydrogens(top, bottom, model):
    """modXNA residues keep the hydrogens of their unit; the 2'-O chains are
    oriented with them, so the methyl and methylene hydrogens stay clear of
    the sugar."""
    from dsxnabuilder.reslib import ResidueLibrary
    from dsxnabuilder.service import Session
    rl = ResidueLibrary.load_default()
    out = Session(residue_library=rl).build({"top": top, "bottom": bottom, "model": model})
    assert out["ok"], out["issues"]
    atoms = parse_pdb(out["pdb"])
    for ats in residues(atoms).values():
        for h, x, d in bonded_h(ats):
            assert 0.94 < d < 1.12, (h.resname, h.name, x.name, d)
    hs = np.array([a.xyz for a in atoms if a.element == "H"])
    dd = np.linalg.norm(hs[:, None] - hs[None], axis=2) + 9.0 * np.eye(len(hs))
    assert dd.min() > 1.5


def test_abasic_c1_hydrogens(session):
    out = session.build({"top": "CGCGA-TTCGCG", "bottom": "GCGCTTAAGCGC", "model": 4})
    ats = [a for a in parse_pdb(out["pdb"]) if a.resname == "3DR"]
    c = {a.name: a.xyz for a in ats}
    for h in ("H1'1", "H1'2"):
        assert 1.05 < np.linalg.norm(c[h] - c["C1'"]) < 1.12
    v1, v2 = c["H1'1"] - c["C1'"], c["H1'2"] - c["C1'"]
    ang = np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
    assert 100 < ang < 120
    # H1'1 takes the place of the base: on the side of C5' (beta face), H1'2 opposite
    ring = np.array([c[n] for n in ("C1'", "C2'", "C3'", "C4'", "O4'")])
    centre = ring.mean(0)
    normal = np.linalg.svd(ring - centre)[2][2]
    side = lambda n: np.dot(c[n] - centre, normal)
    assert side("H1'1") * side("C5'") > 0 > side("H1'2") * side("C5'")


def test_template_from_amber_units_is_complete():
    t = amber_templates()
    for unit in ("DA", "DC5", "DG3", "DTN", "A", "C5", "G3", "UN"):
        assert isinstance(t[unit], Template) and t[unit].hydrogens
        for h in t[unit].hydrogens:
            assert len(t[unit].nb[h]) == 1
