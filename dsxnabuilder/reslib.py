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
Every other backbone modified at phosphorus (methylphosphonate, phosphoryl
guanidine, ...) keeps the handedness its new group has in the modXNA residue
(see _phosphorus_rule): unlike PS1/PS2, those pairs are mirror images at P.

Checks
------
A transplant is refused rather than returned wrong: every carbon or sulfur
stereocentre must keep the handedness it has in the modXNA residue (a shared
substituent such as the O2' of arabinose, which the scaffold puts on the
ribose face, is placed from the residue instead; ring atoms cannot be), every
bond must keep its length within GEOMETRY_TOL and no two atoms three or more
bonds apart may come closer than CONTACT_MIN. modXNA atom names with a leading
digit are read as tleap reads them back from a PDB ("1H15" -> "H151").
"""

from __future__ import annotations

import json
import os
import re
from collections import deque
from typing import Dict, List, Optional, Tuple

import numpy as np

from .catalog import Issue
from .geometry import kabsch
from .offlib import OffUnit, parse_off, rename_unit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
STRAIN_TOL = 0.15     # A, bond length change where a new atom closes a ring
GEOMETRY_TOL = 0.30   # A, bond length change that makes a transplant unusable
CONTACT_MIN = 2.0     # A, non-bonded atoms of one residue (1-4 and beyond)
CHIRAL_MIN = 0.25     # |u1 . (u2 x u3)| of a tetrahedral centre is ~0.5-0.6; planar ~0
# Scaffold atoms that always keep their fiber coordinates: furanose ring,
# backbone and glycosidic atoms (they define the helical frame and O3'-P).
CORE_ATOMS = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'",
              "N9", "N1"}


def tleap_name(name: str) -> str:
    """Atom name as tleap reads it from a PDB: leading digits move to the end
    (old PDB convention, "1H15" -> "H151"). modXNA names some hydrogens with a
    leading digit; written as such, tleap would not match them to the unit."""
    m = re.match(r"^(\d+)(\D.*)$", name)
    return m.group(2) + m.group(1) if m else name


def duplicate_names(unit: OffUnit) -> List[str]:
    """Atom names that occur twice in a unit (as tleap reads them)."""
    seen, out = set(), []
    for n in unit.atom_names:
        if n in seen and n not in out:
            out.append(n)
        seen.add(n)
    return out


def _tleap_names(unit: OffUnit) -> OffUnit:
    unit.atom_names = [tleap_name(n) for n in unit.atom_names]
    return unit


def _rename_atoms_off(text: str) -> str:
    """Apply tleap_name to the quoted atom names of an OFF file."""
    return re.sub(r'"(\d+\D[^"\s]*)"', lambda m: f'"{tleap_name(m.group(1))}"', text)


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
        """The modXNA unit, with atom names as tleap reads them back from the PDB."""
        if key not in self._units:
            self._units[key] = _tleap_names(parse_off(self.text(key))[self.unit_name])
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
        bad = ([f"configuration at {', '.join(report['stereo'])}"] if report["stereo"] else []) + \
            report["distorted"][:3] + report["contacts"][:3]
        if bad:
            issues.append(Issue("error", spec.column, "TRANSPLANT_FAILED",
                                f"{spec.token}: the {'/'.join(spec.triplet())} residue does not fit "
                                f"the fiber scaffold ({'; '.join(bad)})",
                                {"token": spec.token, "triplet": "/".join(spec.triplet()),
                                 "detail": "; ".join(bad)}, None))
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
            # atoms first: a residue name the user chose may start with a digit
            out[fname] = rename_unit(_rename_atoms_off(self.text(self.key_for(r.spec))),
                                     self.unit_name, r.resname)
        return out, self.frcmod()


# --------------------------------------------------------------------------
def transplant(unit: OffUnit, names: List[str], xyz: np.ndarray
               ) -> Tuple[List[str], np.ndarray, dict]:
    """Install ``unit`` (heavy atoms) on a scaffold nucleotide (names, xyz).

    The report lists the atoms placed, strained ring closures (a warning) and
    what makes the result unusable: stereocentres whose configuration differs
    from the modXNA residue (``stereo``), bonds of the residue stretched or
    compressed beyond GEOMETRY_TOL (``distorted``) and non-bonded atoms closer
    than CONTACT_MIN (``contacts``), as when the scaffold ring cannot carry the
    sugar (L- or alpha-sugars, unlocked rings) or base atoms share names with a
    different ring.
    """
    sc = {n: p for n, p in zip(names, xyz)}
    lib = unit.coords()
    heavy = unit.heavy_atoms()
    hset = set(heavy)
    unb = unit.neighbours()
    nb = {a: [b for b in unb[a] if b in hset] for a in heavy}
    placed: Dict[str, np.ndarray] = {}
    elements = dict(zip(unit.atom_names, unit.elements))
    special = _phosphorus_rule(nb, hset, sc, lib, elements)
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
    order = _place(new, placed, lib, nb, unit.name)
    if released:
        _refine(list(new), placed, lib, nb)
    # A shared substituent kept from the scaffold can sit on the wrong face for
    # this sugar (O2' of arabinose: the scaffold puts it where ribose has it).
    # Such atoms are placed from the residue instead; ring atoms cannot be.
    # stereocentres: carbon and sulfur (a three-coordinate N or O inverts freely)
    centres = {n for n, e in zip(unit.atom_names, unit.elements) if e in (6, 16) and n in hset}
    stereo = _stereo_mismatches(lib, placed, nb, centres)
    rings = _ring_atoms(nb)
    redo = sorted({a for c in stereo for a in [c] + nb[c]
                   if a in sc and a not in CORE_ATOMS and a not in special and a not in rings})
    if redo:
        for a in redo:
            placed.pop(a)
        order += _place(redo, placed, lib, nb, unit.name)
        new += redo
        stereo = _stereo_mismatches(lib, placed, nb, centres)
    strained, distorted = [], []
    newset = set(new) | set(special)
    for i, j in unit.bonds:
        a, b = unit.atom_names[i], unit.atom_names[j]
        if a not in hset or b not in hset:
            continue
        d0 = float(np.linalg.norm(lib[a] - lib[b]))
        d1 = float(np.linalg.norm(placed[a] - placed[b]))
        if abs(d1 - d0) > GEOMETRY_TOL:
            distorted.append(f"{a}-{b} {d1:.2f}/{d0:.2f} A")
        elif {a, b} & newset and abs(d1 - d0) > STRAIN_TOL:
            strained.append(f"{a}-{b} {d1:.2f}/{d0:.2f} A")
    out_names = list(heavy)
    out_xyz = np.array([placed[a] for a in out_names])
    if not np.all(np.isfinite(out_xyz)):        # a ring closure that would not relax
        distorted.insert(0, "ring closure did not converge")
    return out_names, out_xyz, {"placed": order + list(special), "strained": strained,
                                "stereo": stereo, "distorted": distorted,
                                "contacts": _close_contacts(out_names, out_xyz, nb)}


def _phosphorus_rule(nb: Dict[str, List[str]], hset: set, sc: Dict[str, np.ndarray],
                     lib: Dict[str, np.ndarray], elements: Dict[str, int]) -> Dict[str, str]:
    """Atoms on P -> scaffold oxygen whose position each one takes.

    The scaffold's non-bridging oxygens follow the Amber convention,
    V(OP1, OP2, O5') < 0 about P. A backbone modified at phosphorus keeps one
    of them and puts a new group X on the other position; X takes the scaffold
    position with the handedness it has in the modXNA residue (OP1 when
    V(X, kept O, O5') < 0, otherwise OP2) and the kept oxygen the other one.
    The oxygen names of the fragments cannot be used for this: most pairs (MP,
    PR, IB, MO, PC, PN, WX, ZV) name them opposite to the Amber convention.

    Phosphorothioates are the exception. modXNA's PS1 and PS2 have the same
    geometry at P, so their names decide: PS1 keeps OP2 and gets S on OP1, PS2
    the reverse (the absolute Rp/Sp is not asserted). PSS: S1 on OP1, S2 on OP2.
    Other patterns keep the geometry of the residue.
    """
    if "P" not in sc or "P" not in hset:
        return {}
    missing = [o for o in ("OP1", "OP2") if o not in hset and o in sc]
    # bridging atoms (bonded to the sugar, e.g. N5' of a phosphoramidate) stay put
    new = sorted(b for b in nb["P"] if b not in sc and "C5'" not in nb[b])
    if missing == ["OP1", "OP2"] and new == ["S1", "S2"]:
        return {"S1": "OP1", "S2": "OP2"}
    if len(missing) != 1 or len(new) != 1:
        return {}
    x, kept = new[0], "OP2" if missing[0] == "OP1" else "OP1"
    if elements.get(x) == 16 or "O5'" not in lib:
        return {x: missing[0]}
    p = lib["P"]
    v = float(np.dot(lib[x] - p, np.cross(lib[kept] - p, lib["O5'"] - p)))
    slot = "OP1" if v < 0 else "OP2"
    if slot == missing[0]:
        return {x: slot}
    return {x: slot, kept: missing[0]}


def _place(atoms: List[str], placed: Dict[str, np.ndarray], lib: Dict[str, np.ndarray],
           nb: Dict[str, List[str]], unit_name: str) -> List[str]:
    """Place ``atoms`` one by one, most-anchored first, by local superposition."""
    order: List[str] = []
    pending = set(atoms)
    while pending:
        progress = False
        # atoms with the most placed neighbours first; ties by name, so that the
        # result does not depend on set order (string hashing differs per process)
        cands = sorted(pending, key=lambda a: (-sum(b in placed for b in nb[a]), a))
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
            raise ValueError(f"cannot place atoms {sorted(pending)} of {unit_name}")
    return order


def _ring_atoms(nb: Dict[str, List[str]]) -> set:
    """Atoms on a ring: two of their neighbours stay connected without them."""
    out = set()
    for a, ns in nb.items():
        if len(ns) < 2:
            continue
        seen, q = {ns[0]}, deque([ns[0]])
        while q:
            x = q.popleft()
            for y in nb.get(x, []):
                if y != a and y not in seen:
                    seen.add(y)
                    q.append(y)
        if any(n in seen for n in ns[1:]):
            out.add(a)
    return out


def _chirality(pos: Dict[str, np.ndarray], c: str, ns: List[str]) -> float:
    u = [(pos[n] - pos[c]) / np.linalg.norm(pos[n] - pos[c]) for n in ns[:3]]
    return float(np.dot(u[0], np.cross(u[1], u[2])))


def _stereo_mismatches(lib: Dict[str, np.ndarray], placed: Dict[str, np.ndarray],
                       nb: Dict[str, List[str]], centres: set) -> List[str]:
    """sp3 ``centres`` (three or more heavy neighbours, not planar in the residue)
    whose handedness differs from the modXNA residue, or that came out flat.
    Phosphorus is not among them: its configuration follows _phosphorus_rule."""
    out = []
    for c, ns in nb.items():
        if c not in centres or len(ns) < 3:
            continue
        ns = sorted(ns)
        v0 = _chirality(lib, c, ns)
        if abs(v0) < CHIRAL_MIN:
            continue
        v1 = _chirality(placed, c, ns)
        if v0 * v1 <= 0 or abs(v1) < CHIRAL_MIN / 2:
            out.append(c)
    return sorted(out)


def _close_contacts(names: List[str], xyz: np.ndarray, nb: Dict[str, List[str]]) -> List[str]:
    """Atom pairs of the residue three or more bonds apart closer than CONTACT_MIN."""
    idx = {n: i for i, n in enumerate(names)}
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    out = []
    for i, j in zip(*np.nonzero(np.triu(d < CONTACT_MIN, 1))):
        a, b = names[i], names[j]
        near = {a}
        for _ in range(2):          # atoms within two bonds of a
            near |= {y for x in near for y in nb.get(x, [])}
        if b not in near:
            out.append(f"{a}...{b} {d[i, j]:.2f} A")
    return out


def _refs(atom: str, nb: Dict[str, List[str]], placed: Dict[str, np.ndarray],
          depth: int) -> List[str]:
    """Placed atoms within ``depth`` bonds of ``atom`` (graph distance)."""
    seen = {atom: 0}
    q = deque([atom])
    while q:
        a = q.popleft()
        if seen[a] == depth:
            continue
        for b in sorted(nb[a]):
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
