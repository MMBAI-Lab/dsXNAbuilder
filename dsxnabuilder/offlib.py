"""
dsxnabuilder.offlib
===================

Minimal reader/renamer for Amber OFF (``.lib``) object files.

Only what dsXNAbuilder needs: per unit, the atom names, types, charges and
elements, the coordinates (``positions``), the bonds (``connectivity``) and
the head/tail atoms (``connect``). Renaming a unit rewrites the unit name in
the ``!!index`` block, in every ``!entry.<name>.`` key and in the residue
table, which is what tleap matches against PDB residue names.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass
class OffUnit:
    name: str
    atom_names: List[str]
    atom_types: List[str]
    charges: List[float]
    elements: List[int]
    xyz: np.ndarray
    bonds: List[Tuple[int, int]]          # 0-based indices
    head: Optional[int] = None            # 0-based, None if absent
    tail: Optional[int] = None
    residue_names: List[str] = field(default_factory=list)

    @property
    def net_charge(self) -> float:
        return float(sum(self.charges))

    def coords(self) -> Dict[str, np.ndarray]:
        return {n: self.xyz[i] for i, n in enumerate(self.atom_names)}

    def neighbours(self) -> Dict[str, List[str]]:
        nb: Dict[str, List[str]] = {n: [] for n in self.atom_names}
        for i, j in self.bonds:
            a, b = self.atom_names[i], self.atom_names[j]
            nb[a].append(b)
            nb[b].append(a)
        return nb

    def heavy_atoms(self) -> List[str]:
        return [n for n, e in zip(self.atom_names, self.elements) if e != 1]


_KEY = re.compile(r"^!entry\.(?P<unit>[^.]+)\.unit\.(?P<field>\S+)\s")


def _split(text: str) -> Dict[str, Dict[str, List[str]]]:
    """{unit: {field: [data lines]}}."""
    out: Dict[str, Dict[str, List[str]]] = {}
    cur = None
    for line in text.splitlines():
        if line.startswith("!"):
            m = _KEY.match(line + " ")
            cur = (m.group("unit"), m.group("field")) if m else None
            if cur:
                out.setdefault(cur[0], {})[cur[1]] = []
            continue
        if cur:
            out[cur[0]][cur[1]].append(line)
    return out


def _tokens(line: str) -> List[str]:
    return re.findall(r'"[^"]*"|\S+', line)


def parse_off(text: str) -> Dict[str, OffUnit]:
    units: Dict[str, OffUnit] = {}
    for name, f in _split(text).items():
        atoms = [_tokens(l) for l in f.get("atoms", []) if l.strip()]
        names = [a[0].strip('"') for a in atoms]
        types = [a[1].strip('"') for a in atoms]
        elmnt = [int(a[6]) for a in atoms]
        chg = [float(a[7]) for a in atoms]
        pos = np.array([[float(x) for x in l.split()] for l in f.get("positions", [])
                        if l.strip()]).reshape(-1, 3)
        bonds = []
        for l in f.get("connectivity", []):
            t = l.split()
            if len(t) >= 2:
                bonds.append((int(t[0]) - 1, int(t[1]) - 1))
        conn = [int(l.split()[0]) for l in f.get("connect", []) if l.strip()]
        head = conn[0] - 1 if len(conn) > 0 and conn[0] > 0 else None
        tail = conn[1] - 1 if len(conn) > 1 and conn[1] > 0 else None
        res = [_tokens(l)[0].strip('"') for l in f.get("residues", []) if l.strip()]
        units[name] = OffUnit(name, names, types, chg, elmnt, pos, bonds,
                              head, tail, res)
    return units


def rename_unit(text: str, old: str, new: str) -> str:
    """Rename unit ``old`` to ``new`` (keys, index and residue table)."""
    text = text.replace(f'!entry.{old}.unit.', f'!entry.{new}.unit.')
    lines = text.splitlines(keepends=True)
    out = []
    section = None
    for line in lines:
        if line.startswith("!"):
            section = line
        elif line.strip() == f'"{old}"' and section and section.startswith("!!index"):
            line = line.replace(f'"{old}"', f'"{new}"')
        elif section and section.startswith(f"!entry.{new}.unit.residues ") \
                and line.lstrip().startswith(f'"{old}"'):
            line = line.replace(f'"{old}"', f'"{new}"', 1)
        elif section and section.startswith(f"!entry.{new}.unit.name "):
            line = line.replace(f'"{old}"', f'"{new}"')
        out.append(line)
    return "".join(out)


def signed_volume(c: Dict[str, np.ndarray], centre: str, a: str, b: str, d: str) -> float:
    """(a-centre) x (b-centre) . (d-centre): sign encodes the handedness."""
    o = c[centre]
    return float(np.dot(np.cross(c[a] - o, c[b] - o), c[d] - o))
