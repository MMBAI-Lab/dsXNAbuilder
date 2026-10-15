"""
dsxnabuilder.nucleotide
=======================

Nucleotide-level chemistry shared by the offline fiber derivation and the
runtime builder: atom-name conventions, sugar-phosphate vs base partition,
glycosidic torsion, base swapping in the standard frame and ribose O2'
placement.

All coordinates handled here for templates are *local*: expressed in the
standard reference frame of the residue's own base (Olson et al. 2001). In
that frame a standard base (3DNA ``Atomic_X.pdb``) needs no transformation,
which is what makes base swapping a pure substitution of atoms.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import dihedral

# Sugar-phosphate atoms (Amber names). Everything else in a nucleotide is base.
BACKBONE_ATOMS = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
                  "C2'", "O2'", "C1'")
PHOSPHATE_ATOMS = ("P", "OP1", "OP2")

# 3DNA/PDB-v2 -> Amber atom names.
NAME_MAP = {"O1P": "OP1", "O2P": "OP2", "C5M": "C7", "O5*": "O5'", "C5*": "C5'",
            "C4*": "C4'", "O4*": "O4'", "C3*": "C3'", "O3*": "O3'", "C2*": "C2'",
            "O2*": "O2'", "C1*": "C1'"}

WATSON_CRICK = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
WC_PAIRS = {frozenset("AT"), frozenset("AU"), frozenset("GC")}
WOBBLE_PAIRS = {frozenset("GU"), frozenset("GT")}
PURINES = {"A", "G", "I", "X"}         # I, X: inosine and xanthosine of 3DNA fiber repeats
# standard base whose ring atoms define the reference frame of each base letter
RING_LETTER = {"A": "A", "G": "G", "I": "G", "X": "G", "C": "C", "T": "T", "U": "U"}

# Ideal bond lengths (A) for template QC; Parkinson et al. (1996) Acta Cryst
# D52:57-64 / Gelbin et al. (1996) JACS 118:519-529 averages, rounded.
IDEAL_BONDS = {
    ("P", "OP1"): 1.485, ("P", "OP2"): 1.485, ("P", "O5'"): 1.593,
    ("O5'", "C5'"): 1.440, ("C5'", "C4'"): 1.510, ("C4'", "O4'"): 1.446,
    ("C4'", "C3'"): 1.524, ("C3'", "O3'"): 1.423, ("C3'", "C2'"): 1.525,
    ("C2'", "C1'"): 1.521, ("C1'", "O4'"): 1.420, ("C2'", "O2'"): 1.413,
}
GLYCOSIDIC_BOND = 1.47                 # C1'-N9/N1, same sources


def normalize_name(name: str) -> str:
    """Amber name of a 3DNA/PDB-v2 atom name."""
    return NAME_MAP.get(name, name)


def chi(c: Dict[str, np.ndarray]) -> Optional[float]:
    """Glycosidic torsion O4'-C1'-N9-C4 (purines) / O4'-C1'-N1-C2 (pyrimidines)."""
    if "N9" in c and "C4" in c:
        keys = ("O4'", "C1'", "N9", "C4")
    elif "N1" in c and "C2" in c:
        keys = ("O4'", "C1'", "N1", "C2")
    else:
        return None
    if not all(k in c for k in keys):
        return None
    return dihedral(*(c[k] for k in keys))


def is_syn(chi_deg: Optional[float]) -> bool:
    """syn when chi is in (-45, 135).

    The IUPAC-IUB (1983) split at +-90 deg puts high-anti nucleotides (chi
    around -60 to -90, e.g. the C-DNA fiber model #8 at -84.5) on the syn
    side. Here syn means the base sits over the sugar, which is what matters
    for sequence constraints (Z-DNA syn G at chi ~ +60), so the boundary is
    moved to keep the whole high-anti range in anti.
    """
    return chi_deg is not None and -45.0 < chi_deg < 135.0


def split_residue(names: Sequence[str]) -> Tuple[List[int], List[int]]:
    """Indices of (sugar-phosphate, base) atoms."""
    bb = [i for i, n in enumerate(names) if n in BACKBONE_ATOMS]
    base = [i for i, n in enumerate(names) if n not in BACKBONE_ATOMS]
    return bb, base


def base_swap(names: Sequence[str], xyz: np.ndarray,
              std_names: Sequence[str], std_xyz: np.ndarray
              ) -> Tuple[List[str], np.ndarray]:
    """Replace the base of a local-frame nucleotide by a standard base.

    Sugar-phosphate atoms (including the template's own C1') are kept; the
    standard base atoms (``Atomic_X`` minus its C1') are inserted. This is the
    same operation as 3DNA ``mutate_bases``: the standard frame is built so
    that C1' and the glycosidic N of the four bases nearly coincide, so the
    glycosidic bond survives the swap (checked in tests/test_base_swap.py).
    """
    bb, _ = split_residue(names)
    keep_names = [names[i] for i in bb]
    keep_xyz = [xyz[i] for i in bb]
    for n, p in zip(std_names, std_xyz):
        if n == "C1'":
            continue
        keep_names.append(n)
        keep_xyz.append(p)
    return _order(keep_names, np.array(keep_xyz))


# Output order close to Amber's residue templates: phosphate, sugar, base.
_ORDER = {n: i for i, n in enumerate(
    ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C1'"])}
_TAIL = {n: i for i, n in enumerate(["C3'", "C2'", "O2'", "O3'"])}


def _order(names: List[str], xyz: np.ndarray) -> Tuple[List[str], np.ndarray]:
    def key(i):
        n = names[i]
        if n in _ORDER:
            return (0, _ORDER[n])
        if n in _TAIL:
            return (2, _TAIL[n])
        return (1, i)
    idx = sorted(range(len(names)), key=key)
    return [names[i] for i in idx], xyz[idx]


def tetrahedral_positions(centre: np.ndarray, a: np.ndarray, b: np.ndarray,
                          bond: float) -> Tuple[np.ndarray, np.ndarray]:
    """The two remaining sp3 positions on ``centre`` bonded to ``a`` and ``b``."""
    u = (a - centre) / np.linalg.norm(a - centre)
    v = (b - centre) / np.linalg.norm(b - centre)
    bis = -(u + v)
    bis /= np.linalg.norm(bis)
    perp = np.cross(u, v)
    perp /= np.linalg.norm(perp)
    # Ideal tetrahedral half-angle between the two new bonds: 109.47/2.
    half = math.radians(109.47 / 2.0)
    p1 = centre + bond * (math.cos(half) * bis + math.sin(half) * perp)
    p2 = centre + bond * (math.cos(half) * bis - math.sin(half) * perp)
    return p1, p2


def place_o2prime(c: Dict[str, np.ndarray],
                  bond: float = IDEAL_BONDS[("C2'", "O2'")]) -> np.ndarray:
    """O2' of a ribose built on a (deoxy)ribose sugar.

    In ribonucleosides O2' lies on the alpha face of the furanose, opposite to
    the base (beta face) and cis to O3'. Of the two free tetrahedral positions
    on C2' we take the one on the opposite side of the C1'-C2'-C3' plane from
    the glycosidic nitrogen. In 2'-deoxyribose that position is H2''.
    """
    gly = c["N9"] if "N9" in c else c["N1"]
    p1, p2 = tetrahedral_positions(c["C2'"], c["C1'"], c["C3'"], bond)
    normal = np.cross(c["C1'"] - c["C2'"], c["C3'"] - c["C2'"])
    side_n = np.dot(gly - c["C2'"], normal)
    return p1 if np.dot(p1 - c["C2'"], normal) * side_n < 0 else p2


def to_rna(names: List[str], xyz: np.ndarray) -> Tuple[List[str], np.ndarray]:
    """The nucleotide with O2' added (ribose), unless it has one already."""
    if "O2'" in names:
        return names, xyz
    c = {n: xyz[i] for i, n in enumerate(names)}
    o2 = place_o2prime(c)
    return _order(names + ["O2'"], np.vstack([xyz, o2]))


def to_dna(names: List[str], xyz: np.ndarray) -> Tuple[List[str], np.ndarray]:
    """The nucleotide without O2' (2'-deoxyribose)."""
    keep = [i for i, n in enumerate(names) if n != "O2'"]
    return [names[i] for i in keep], xyz[keep]


def pair_kind(a: str, b: str) -> str:
    """'wc', 'wobble' or 'mismatch' for two canonical letters."""
    s = frozenset((a, b))
    if s in WC_PAIRS:
        return "wc"
    if s in WOBBLE_PAIRS:
        return "wobble"
    return "mismatch"
