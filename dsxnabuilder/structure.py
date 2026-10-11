"""
dsxnabuilder.structure
======================

Atoms, minimal PDB input/output and the standard base reference frame.

The standard reference frame is that of Olson et al. (2001) J Mol Biol
313:229-237 (doi:10.1006/jmbi.2001.4987): x towards the major groove along the
pseudo-dyad, y along the long axis of the pair, z = x cross y. A base frame is
obtained by superposing the ring atoms of a standard base (which by
construction sits in that frame) onto the observed base: the fitted rotation
and translation *are* the frame. Exocyclic atoms are excluded so that the frame
is robust to substituted base analogues (5mC, s4T, Br5U, ...).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import Frame, kabsch

PURINE_RING = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
PYRIMIDINE_RING = ["N1", "C2", "N3", "C4", "C5", "C6"]

# Polar ring/exocyclic atoms that can take part in base-pair hydrogen bonds.
# Used to detect pairing partners in fiber models (Watson-Crick or Hoogsteen).
HBOND_ATOMS = ("N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6")


@dataclass
class Atom:
    name: str
    resname: str
    resseq: int
    chain: str
    xyz: np.ndarray
    element: str = ""

    def copy(self, **kw) -> "Atom":
        d = dict(name=self.name, resname=self.resname, resseq=self.resseq,
                 chain=self.chain, xyz=np.array(self.xyz, float),
                 element=self.element)
        d.update(kw)
        return Atom(**d)


def element_of(name: str) -> str:
    """Element from an Amber/PDB nucleic-acid atom name (no metals expected)."""
    n = name.lstrip("0123456789")
    if n[:2].upper() in ("BR", "CL"):
        return n[:2].capitalize()
    return n[:1].upper()


def pdb_atom_line(serial: int, a: Atom, record: str = "ATOM") -> str:
    # PDB column 13-16: names shorter than 4 characters start in column 14
    # unless the element symbol has two letters.
    nm = a.name if len(a.name) >= 4 else f" {a.name:<3s}"
    el = a.element or element_of(a.name)
    return (f"{record:<6s}{serial % 100000:5d} {nm:<4s} {a.resname:>3s} "
            f"{a.chain:1s}{a.resseq % 10000:4d}    "
            f"{a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}"
            f"  1.00  0.00          {el:>2s}")


def read_pdb(path: str) -> List[Atom]:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return parse_pdb(fh.read())


def parse_pdb(text: str) -> List[Atom]:
    atoms: List[Atom] = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        atoms.append(Atom(
            name=name, resname=line[17:21].strip(), resseq=int(line[22:26]),
            chain=(line[21:22].strip() or "A"),
            xyz=np.array([float(line[30:38]), float(line[38:46]),
                          float(line[46:54])]),
            element=line[76:78].strip() or element_of(name)))
    return atoms


def format_pdb(atoms: Sequence[Atom], remarks: Iterable[str] = ()) -> str:
    out = [f"REMARK 999 {r}" for r in remarks]
    serial = 0
    prev = None
    for a in atoms:
        if prev is not None and a.chain != prev:
            out.append("TER")
        serial += 1
        out.append(pdb_atom_line(serial, a))
        prev = a.chain
    out += ["TER", "END", ""]
    return "\n".join(out)


def write_pdb(atoms: Sequence[Atom], path: str, remarks: Iterable[str] = ()) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(format_pdb(atoms, remarks))


def residues(atoms: Sequence[Atom]) -> Dict[Tuple[str, int], List[Atom]]:
    """Group atoms by (chain, resseq), preserving order."""
    out: Dict[Tuple[str, int], List[Atom]] = {}
    for a in atoms:
        out.setdefault((a.chain, a.resseq), []).append(a)
    return out


def base_class_of(names: Iterable[str]) -> str:
    return "purine" if "N9" in set(names) else "pyrimidine"


def ring_atoms(names: Iterable[str]) -> List[str]:
    return PURINE_RING if base_class_of(names) == "purine" else PYRIMIDINE_RING


def standard_base_frame(coords: Dict[str, np.ndarray],
                        reference: Dict[str, np.ndarray]) -> Tuple[Frame, float]:
    """Standard reference frame of an observed base.

    ``reference`` holds the standard base in its own frame (Atomic_X). The ring
    atoms common to both are superposed; returns (frame, ring-fit rmsd).
    """
    ring = [a for a in ring_atoms(reference) if a in coords and a in reference]
    if len(ring) < 3:
        raise ValueError(f"not enough ring atoms to define a base frame: {sorted(coords)}")
    R, t, rmsd = kabsch(np.array([reference[a] for a in ring]),
                        np.array([coords[a] for a in ring]))
    return Frame(t, R), rmsd


def o3p_distances(atoms: Sequence[Atom]) -> List[Tuple[str, int, float]]:
    """O3'(i)-P(i+1) distances along each chain: [(chain, i, d)]."""
    by_res: Dict[Tuple[str, int], Dict[str, np.ndarray]] = {}
    for a in atoms:
        by_res.setdefault((a.chain, a.resseq), {})[a.name] = a.xyz
    out = []
    for (chain, i), res in sorted(by_res.items()):
        nxt = by_res.get((chain, i + 1))
        if nxt is None or "O3'" not in res or "P" not in nxt:
            continue
        out.append((chain, i, float(np.linalg.norm(res["O3'"] - nxt["P"]))))
    return out
