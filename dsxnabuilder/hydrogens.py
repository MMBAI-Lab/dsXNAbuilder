"""
dsxnabuilder.hydrogens
======================

Hydrogen atoms of a built duplex.

Every residue takes the hydrogens of the template tleap will match it to: the
Amber unit for canonical nucleotides (DNA.OL15 / RNA.lib, chosen by chain
position: DA5, DA, DA3, ...) and the modXNA unit for modified ones. The names
are therefore the ones tleap expects and it keeps the coordinates as written.

Each hydrogen is placed like the new atoms of the modXNA transplant: the
template coordinates of the heavy atoms within two bonds of it (three when
fewer than three are present) are superposed on the built ones and the
hydrogen follows. Heavy atoms do not move; bond lengths and angles, and the
torsions of methyl and hydroxyl hydrogens, come from the template.

Abasic sites take the sugar-phosphate hydrogens of the canonical unit of their
chemistry plus a second hydrogen on C1' along the former glycosidic bond. They
are named as in the PDB chemical component 3DR, where C1' carries H1'1 (on the
side of the missing base) and H1'2, and C4' carries H4'1.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import kabsch
from .graph import within
from .paths import DATA

AMBER_JSON = os.path.join(DATA, "amber_residues.json")
ABASIC_NAMES = {"H1'": "H1'2", "H4'": "H4'1"}     # Amber -> PDB component 3DR


class Template:
    """Atom names, elements, coordinates and bonds of a residue template."""

    def __init__(self, names: Sequence[str], elements: Sequence[int], xyz, bonds):
        self.names = list(names)
        self.elements = list(elements)
        self.xyz = {n: np.asarray(p, float) for n, p in zip(self.names, np.asarray(xyz, float))}
        self.nb: Dict[str, List[str]] = {n: [] for n in self.names}
        for i, j in bonds:
            a, b = self.names[i], self.names[j]
            self.nb[a].append(b)
            self.nb[b].append(a)
        self.hydrogens = [n for n, e in zip(self.names, self.elements) if e == 1]
        self.heavy = {n for n, e in zip(self.names, self.elements) if e != 1}

    @classmethod
    def from_off(cls, unit) -> "Template":
        """Template of a parsed OFF unit (offlib.OffUnit)."""
        return cls(unit.atom_names, unit.elements, unit.xyz, unit.bonds)


_AMBER: Optional[Dict[str, Template]] = None


def amber_templates(path: str = AMBER_JSON) -> Dict[str, Template]:
    """Amber units of DNA.OL15.lib and RNA.lib by name (read once)."""
    global _AMBER
    if _AMBER is None:
        with open(path, encoding="utf-8") as fh:
            units = json.load(fh)["units"]
        _AMBER = {k: Template(u["atoms"], u["elements"], u["xyz"], u["bonds"])
                  for k, u in units.items() if "atoms" in u}
    return _AMBER


def amber_unit(letter: str, chem: str, cap: str) -> str:
    """Amber unit tleap uses for a canonical nucleotide at this chain position."""
    name = ("D" + letter) if chem == "DNA" else letter
    return name + {"5cap": "5", "3cap": "3"}.get(cap, "")


def add_hydrogens(names: Sequence[str], xyz: np.ndarray, tpl: Template,
                  abasic: bool = False) -> Tuple[List[str], np.ndarray]:
    """Heavy atoms ``names``/``xyz`` with the template's hydrogens added.

    Each hydrogen is written right after the heavy atom it is bonded to.
    Hydrogens whose heavy atom is absent (the base of an abasic site) are
    skipped.
    """
    built = {n: np.asarray(p, float) for n, p in zip(names, xyz)}
    present = set(built) & tpl.heavy
    on: Dict[str, List[Tuple[str, np.ndarray]]] = {}
    for parent in tpl.names:
        hs = [h for h in tpl.hydrogens if tpl.nb[h][0] == parent]
        if not hs or parent not in built:
            continue
        # the hydrogens of one atom share its reference atoms: one fit for all
        for h, p in zip(hs, _fit(tpl, hs, present, built)):
            on.setdefault(parent, []).append((ABASIC_NAMES.get(h, h) if abasic else h, p))
    if abasic and "C1'" in built:
        on.setdefault("C1'", []).insert(0, ("H1'1", _abasic_h1(tpl, present, built)))
    out_names: List[str] = []
    out_xyz: List[np.ndarray] = []
    for n in names:
        out_names.append(n)
        out_xyz.append(built[n])
        for hn, p in on.get(n, []):
            out_names.append(hn)
            out_xyz.append(p)
    return out_names, np.array(out_xyz)


def _fit(tpl: Template, atoms: Sequence[str], present: set,
         built: Dict[str, np.ndarray]) -> List[np.ndarray]:
    """Place ``atoms``, all bonded to the same atom, by superposing the
    reference atoms around it. The fit is anchored on the bonded atom, so the
    bond keeps the template length exactly (a plain fit spreads the misfit of
    the sugar over the anchor and stretches C-H bonds by up to 0.1 A)."""
    refs = _refs(tpl, atoms[0], present, 2)
    if len(refs) < 3:
        refs = _refs(tpl, atoms[0], present, 3)
    if len(refs) < 3:
        raise ValueError(f"cannot place {atoms[0]}: too few heavy atoms around it")
    # refs[0] is the bonded atom; its overwhelming weight makes the fit a pure
    # rotation about it
    weights = [1e6] + [1.0] * (len(refs) - 1)
    R, _, _ = kabsch(np.array([tpl.xyz[r] for r in refs]), np.array([built[r] for r in refs]),
                     weights)
    a0, a1 = tpl.xyz[refs[0]], built[refs[0]]
    return [a1 + R @ (tpl.xyz[a] - a0) for a in atoms]


def _abasic_h1(tpl: Template, present: set, built: Dict[str, np.ndarray]) -> np.ndarray:
    """Second C1' hydrogen of an abasic site: along the template's glycosidic
    bond, at the template's C1'-H1' length."""
    gly = "N9" if "N9" in tpl.heavy else "N1"
    n = _fit(tpl, [gly], present, built)[0]
    c1 = built["C1'"]
    d = float(np.linalg.norm(tpl.xyz["H1'"] - tpl.xyz["C1'"]))
    return c1 + d * (n - c1) / np.linalg.norm(n - c1)


def _refs(tpl: Template, atom: str, present: set, depth: int) -> List[str]:
    """Present heavy atoms within ``depth`` bonds of ``atom``, through heavy atoms."""
    return [a for a in within(tpl.nb, atom, depth, through=tpl.heavy) if a in present]
