"""
dsxnabuilder.reslib
===================

Prebuilt modXNA residues and the *local transplant* that installs them on the
fiber scaffold.

Transplant
----------
The scaffold nucleotide sits exactly in its fiber frame. A modXNA residue
(built and relaxed in isolation by modxna.sh) shares most heavy atoms with it
by name. The transplant keeps the scaffold coordinates of every shared atom,
so the helical frame and the O3'-P continuity stay exact and the sugar keeps
the model's pucker; only the atoms unique to the modXNA residue (2'-F, the
2'-O-methyl or MOE chain, the LNA bridge, 5-methyl, ...) are placed, one by
one in breadth-first order from the shared core, by superposing the library
coordinates of their already placed neighbours within two bonds (three if
fewer than three are available) onto the placed coordinates. Bond lengths,
angles and the torsions around the new atoms therefore come from the modXNA
residue, and the strain left where a new ring closes onto the scaffold (LNA on
a pucker other than its own) is measured and reported.

Phosphorothioates are placed by an explicit rule rather than by fitting. Both
modXNA fragments PS1 and PS2 put the sulfur on the same pro-chiral position of
their own coordinates and differ only in which oxygen they keep (PS1 keeps
OP2, PS2 keeps OP1), so a fit would give the same diastereomer for both. Here
the sulfur takes the place of the oxygen that the fragment does *not* keep:
PS1 puts S on the scaffold's OP1, PS2 on its OP2, with OP1/OP2 in the Amber
convention (the scaffold is relabelled by chirality in
tools/derive_fiber_models.py). The absolute Rp/Sp labels are not asserted.
"""

from __future__ import annotations

import json
import os
from collections import deque
from typing import Dict, List, Optional, Tuple

import numpy as np

from .catalog import Issue
from .geometry import kabsch
from .offlib import OffUnit, parse_off, rename_unit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
STRAIN_TOL = 0.15     # A, bond length change where a new atom closes a ring
# Scaffold atoms that always keep their fiber coordinates: furanose ring,
# backbone and glycosidic atoms (they define the helical frame and O3'-P).
CORE_ATOMS = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'",
              "N9", "N1"}


def residue_key(backbone: str, sugar: str, base: str, cap: str) -> str:
    return f"{backbone}-{sugar}-{base}-{cap}"


class ResidueLibrary:
    def __init__(self, index: dict, root: Optional[str] = None,
                 frcmod_path: Optional[str] = None):
        self.index = index
        self.residues: Dict[str, dict] = index.get("residues", {})
        self.unit_name = index.get("unit", "MXR")
        self.root = root
        self.frcmod_path = frcmod_path
        self._texts: Dict[str, str] = {}
        self._units: Dict[str, OffUnit] = {}
        self._frcmod: Optional[str] = None

    @classmethod
    def load_default(cls) -> "ResidueLibrary":
        with open(os.path.join(ROOT, "data", "residue_index.json"), encoding="utf-8") as fh:
            index = json.load(fh)
        root = os.path.join(ROOT, "web", "lib", "residues")
        frc = os.path.join(ROOT, "web", "lib", "frcmod.modxna")
        return cls(index, root if os.path.isdir(root) else None,
                   frc if os.path.exists(frc) else None)

    # ------------------------------------------------------------------
    def key_for(self, spec) -> str:
        return residue_key(spec.backbone, spec.sugar, spec.base, spec.cap)

    def has(self, spec) -> bool:
        return self.key_for(spec) in self.residues

    def provide(self, libs: Optional[Dict[str, str]] = None,
                frcmod: Optional[str] = None) -> None:
        """Library texts fetched by the browser (key -> OFF text)."""
        for k, v in (libs or {}).items():
            self._texts[k] = v
            self._units.pop(k, None)
        if frcmod:
            self._frcmod = frcmod

    def text(self, key: str) -> str:
        if key not in self._texts:
            if self.root is None:
                raise KeyError(f"residue library file for {key} was not provided")
            with open(os.path.join(self.root, self.residues[key]["file"]), encoding="utf-8") as fh:
                self._texts[key] = fh.read()
        return self._texts[key]

    def available(self, key: str) -> bool:
        return key in self._texts or self.root is not None

    def unit(self, key: str) -> OffUnit:
        if key not in self._units:
            self._units[key] = parse_off(self.text(key))[self.unit_name]
        return self._units[key]

    def frcmod(self) -> str:
        if self._frcmod is None:
            if self.frcmod_path is None:
                raise KeyError("frcmod.modxna was not provided")
            with open(self.frcmod_path, encoding="utf-8") as fh:
                self._frcmod = fh.read()
        return self._frcmod

    # ------------------------------------------------------------------
    def transplant(self, spec, names: List[str], xyz: np.ndarray
                   ) -> Tuple[List[str], np.ndarray, dict]:
        key = self.key_for(spec)
        unit = self.unit(key)
        new_names, new_xyz, report = transplant(unit, names, xyz)
        issues = []
        if report["strained"]:
            issues.append(Issue("warning", spec.column, "TRANSPLANT_STRAIN",
                                f"{spec.token}: modified atoms close onto the scaffold with "
                                f"strained bonds ({', '.join(report['strained'][:3])}); "
                                "minimize before MD",
                                {"token": spec.token, "bonds": ", ".join(report["strained"][:3])}))
        return new_names, new_xyz, {"key": key, "placed": report["placed"], "issues": issues}

    def package_libs(self, result, provided: Optional[Dict[str, str]] = None
                     ) -> Tuple[Dict[str, str], str]:
        """OFF files renamed to the residue names used in the PDB."""
        self.provide(provided)
        out: Dict[str, str] = {}
        for r in result.modified():
            fname = f"{r.resname}.lib"
            if fname in out:
                continue
            out[fname] = rename_unit(self.text(self.key_for(r.spec)), self.unit_name, r.resname)
        return out, self.frcmod()


# --------------------------------------------------------------------------
def transplant(unit: OffUnit, names: List[str], xyz: np.ndarray
               ) -> Tuple[List[str], np.ndarray, dict]:
    """Install ``unit`` (heavy atoms) on a scaffold nucleotide (names, xyz)."""
    sc = {n: p for n, p in zip(names, xyz)}
    lib = unit.coords()
    heavy = unit.heavy_atoms()
    hset = set(heavy)
    nb = {a: [b for b in unit.neighbours()[a] if b in hset] for a in heavy}
    placed: Dict[str, np.ndarray] = {}
    special = {}
    # phosphorothioate rule (see module docstring)
    if "P" in sc:
        if "S1" in hset and "OP2" in hset and "OP1" in sc:
            special["S1"] = "OP1"
        if "S2" in hset and "OP1" in hset and "OP2" in sc:
            special["S2"] = "OP2"
    for a in heavy:
        if a in sc and a not in special:
            placed[a] = np.array(sc[a], float)
    # Ring closure (LNA/cEt bridge C4'-C6'-O2'): a shared substituent bonded to a
    # new atom that also bonds to the scaffold core must move with the bridge,
    # otherwise the ribose O2' position leaves the bridge 0.6 A too long.
    released = set()
    for a in heavy:
        if a in placed:
            continue
        shared_nb = [b for b in nb[a] if b in placed]
        if len(shared_nb) >= 2:
            for b in shared_nb:
                if b not in CORE_ATOMS:
                    placed.pop(b, None)
                    released.add(b)
    for s_name, o_name in special.items():
        v = sc[o_name] - sc["P"]
        d = float(np.linalg.norm(lib[s_name] - lib["P"]))
        placed[s_name] = sc["P"] + d * v / np.linalg.norm(v)
    new = [a for a in heavy if a not in placed and a not in special]
    order: List[str] = []
    pending = set(new)
    while pending:
        progress = False
        # atoms with the most placed neighbours first
        cands = sorted(pending, key=lambda a: -sum(b in placed for b in nb[a]))
        for a in cands:
            if not any(b in placed for b in nb[a]):
                continue
            refs = _refs(a, nb, placed, 2)
            if len(refs) < 3:
                refs = _refs(a, nb, placed, 3)
            if len(refs) < 3:
                continue
            R, t, _ = kabsch(np.array([lib[r] for r in refs]),
                             np.array([placed[r] for r in refs]))
            placed[a] = R @ lib[a] + t
            order.append(a)
            pending.discard(a)
            progress = True
            break
        if not progress:
            raise ValueError(f"cannot place atoms {sorted(pending)} of {unit.name}")
    if released:
        _refine(list(new), placed, lib, nb)
    strained = []
    newset = set(new) | set(special)
    for i, j in unit.bonds:
        a, b = unit.atom_names[i], unit.atom_names[j]
        if a not in hset or b not in hset or not ({a, b} & newset):
            continue
        d0 = float(np.linalg.norm(lib[a] - lib[b]))
        d1 = float(np.linalg.norm(placed[a] - placed[b]))
        if abs(d1 - d0) > STRAIN_TOL:
            strained.append(f"{a}-{b} {d1:.2f}/{d0:.2f} A")
    out_names = list(heavy)
    out_xyz = np.array([placed[a] for a in out_names])
    return out_names, out_xyz, {"placed": order + list(special), "strained": strained}


def _refs(atom: str, nb: Dict[str, List[str]], placed: Dict[str, np.ndarray],
          depth: int) -> List[str]:
    """Placed atoms within ``depth`` bonds of ``atom`` (graph distance)."""
    seen = {atom: 0}
    q = deque([atom])
    while q:
        a = q.popleft()
        if seen[a] == depth:
            continue
        for b in nb[a]:
            if b not in seen:
                seen[b] = seen[a] + 1
                q.append(b)
    return [a for a, d in sorted(seen.items(), key=lambda kv: kv[1])
            if d > 0 and a in placed]


def _refine(movable: List[str], placed: Dict[str, np.ndarray], lib: Dict[str, np.ndarray],
            nb: Dict[str, List[str]], iters: int = 400) -> None:
    """Relax ring-closing new atoms to the library 1-2 and 1-3 distances.

    Only the new atoms move; the scaffold stays fixed. Used when new atoms
    close a ring onto the scaffold (LNA/cEt bridge), so that the misfit between
    the fiber sugar and the bicyclic library sugar is spread over bond angles
    instead of being left in one bond. Plain gradient descent: a handful of
    atoms, no external dependency (runs in Pyodide).
    """
    mov = set(movable)
    pairs = []
    for a in movable:
        for b in nb[a]:
            pairs.append((a, b, 4.0))
            for c in nb[b]:
                if c != a:
                    pairs.append((a, c, 0.3))
    target = {(a, b): float(np.linalg.norm(lib[a] - lib[b])) for a, b, _ in pairs}
    for _ in range(iters):
        grad = {a: np.zeros(3) for a in movable}
        for a, b, w in pairs:
            r = placed[a] - placed[b]
            d = float(np.linalg.norm(r)) or 1e-9
            g = 2.0 * w * (d - target[(a, b)]) * r / d
            grad[a] += g
            if b in mov:
                grad[b] -= g
        for a in movable:
            placed[a] = placed[a] - 0.05 * grad[a]
