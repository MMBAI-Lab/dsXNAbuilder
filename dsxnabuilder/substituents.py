"""
dsxnabuilder.substituents
=========================

Torsions of flexible 2'-O substituents (2'-OMe, MOE) in the built duplex.

The transplant copies the substituent from the modXNA residue, minimized as
an isolated nucleotide, so in a duplex the MOE chain can run into the next
nucleotide (C6' 1.2 A from its C5'). The library conformation is kept when
it makes no contact shorter than CHAIN_CONTACT; otherwise the chain
torsions are scanned in CHAIN_STEP steps, two bonds at a time from O2'
outwards, and the contact-free conformation closest to the library one is
taken. Rotations are rigid about each bond, so bond lengths and angles stay
those of the library. The LNA bridge closes a ring and is not touched.

Hydrogens are already placed and turn with their chain atoms; contacts
involving a hydrogen count below CHAIN_H_CONTACT, which also staggers the
terminal methyl against the sugar. A bare 2'-OH is left as placed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .graph import far_from
from .structure import element_of

CHAIN_CONTACT = 2.9        # A, heavy-atom distance a 2'-O substituent should keep
CHAIN_H_CONTACT = 2.0      # A, the same for pairs involving a hydrogen
CHAIN_STEP = 30.0          # degrees, torsion scan step
ENV_RADIUS = 10.0          # A, atoms of other residues considered around O2'
# Atoms of the same residue count from this many bonds away: 1-5 contacts such
# as CM2...C4' at 2.9 A belong to the minimized modXNA residue.
OWN_MIN_BONDS = 5
SCAN_ROUNDS = 2            # passes over the rotatable bonds
TOL = 1e-9                 # penalties this close are equal


@dataclass
class _Chain:
    """The substituent of one residue and what it must keep away from."""

    atoms: List[str]                            # from O2' outwards (breadth first)
    bonds: List[Tuple[str, str, List[str]]]     # rotatable bond u-v and the atoms it moves
    start: Dict[str, np.ndarray]                # residue coordinates before the scan
    env: np.ndarray                             # atoms around the chain
    cut: np.ndarray                             # chain atom x env atom: closest allowed (0: any)
    pairs: List[Tuple[int, int, float]]         # chain atoms OWN_MIN_BONDS apart: i, j, limit

    def conformation(self, angles: Sequence[float]) -> Dict[str, np.ndarray]:
        """Residue coordinates with each rotatable bond turned by its angle."""
        pos = dict(self.start)
        for (u, v, moving), ang in zip(self.bonds, angles):
            if ang:
                pts = np.array([pos[m] for m in moving])
                pts = rotate_about(pts, pos[v], pos[v] - pos[u], ang)
                for m, p in zip(moving, pts):
                    pos[m] = p
        return pos

    def penalty(self, pos: Dict[str, np.ndarray]) -> float:
        """Sum of squared shortfalls below the contact limits (0: no contact)."""
        c = np.array([pos[a] for a in self.atoms])
        d = np.linalg.norm(c[:, None, :] - self.env[None, :, :], axis=2)
        p = float(np.sum(np.clip(self.cut - d, 0.0, None) ** 2))
        for i, j, lim in self.pairs:
            p += max(0.0, lim - float(np.linalg.norm(c[i] - c[j]))) ** 2
        return p


def orient_2o_chains(residues: Sequence, residue_library) -> None:
    """Re-orient the 2'-O substituents of the modXNA residues in place.

    ``residues`` are build.BuiltResidue (names, xyz, flags); the torsions
    applied are recorded in ``flags["chain_torsions"]``.
    """
    allxyz = np.vstack([r.xyz for r in residues])
    first = np.cumsum([0] + [len(r.names) for r in residues])
    owner = np.concatenate([[k] * len(r.names) for k, r in enumerate(residues)])
    is_h = np.array([element_of(n) == "H" for r in residues for n in r.names])
    for k, r in enumerate(residues):
        if not r.flags.get("modxna") or "O2'" not in r.names:
            continue
        ch = _chain_of(r, residue_library.unit(r.flags["key"]), allxyz[owner != k],
                       is_h[owner != k])
        if ch is None:
            continue
        best = _scan(ch)
        if best is None:
            continue
        idx = {n: i for i, n in enumerate(r.names)}
        pos = ch.conformation(best)
        for a in ch.atoms:
            r.xyz[idx[a]] = pos[a]
        # the residues scanned next see this chain where it is now, so that two
        # neighbouring chains do not both turn into the same space
        allxyz[first[k]:first[k + 1]] = r.xyz
        r.flags["chain_torsions"] = [round(x, 1) for x in best]


def _chain_of(r, unit, other_xyz: np.ndarray, other_h: np.ndarray) -> Optional[_Chain]:
    """The substituent on O2' of residue ``r`` with its contact limits, or None
    (no substituent, only hydrogens, or a ring closure)."""
    present = set(r.names)
    unb = unit.neighbours()
    nb = {a: [b for b in unb.get(a, []) if b in present] for a in r.names}
    tree = _substituent_tree(nb)
    if not tree:
        return None
    parent, chain = tree
    if all(element_of(a) == "H" for a in chain):
        return None
    idx = {n: i for i, n in enumerate(r.names)}
    desc = {a: [b for b in chain if _ancestor(parent, b, a)] for a in ["O2'"] + chain}
    bonds = [(parent[a] if a != "O2'" else "C2'", a, desc[a]) for a in ["O2'"] + chain
             if desc[a]]
    # the environment: other residues' atoms near O2', and own atoms not in the chain
    o2 = r.xyz[idx["O2'"]]
    near = np.linalg.norm(other_xyz - o2, axis=1) < ENV_RADIUS
    own = [n for n in r.names if n not in chain]
    env = np.vstack([other_xyz[near], r.xyz[[idx[n] for n in own]]])
    env_h = np.concatenate([other_h[near], [element_of(n) == "H" for n in own]])
    far = {a: far_from(nb, [a], OWN_MIN_BONDS) for a in chain}
    counted = np.array([[True] * int(near.sum()) + [b in far[a] for b in own]
                        for a in chain])
    chain_h = np.array([element_of(a) == "H" for a in chain])
    cut = np.where(chain_h[:, None] | env_h[None, :], CHAIN_H_CONTACT, CHAIN_CONTACT)
    cut = np.where(counted, cut, 0.0)
    pairs = [(i, j, CHAIN_H_CONTACT if chain_h[i] or chain_h[j] else CHAIN_CONTACT)
             for i, a in enumerate(chain) for j, b in enumerate(chain)
             if i < j and b in far[a]]
    start = {a: r.xyz[idx[a]].copy() for a in r.names}
    return _Chain(chain, bonds, start, env, cut, pairs)


def _turn(angles: Sequence[float]) -> float:
    """How far a set of torsion changes is from the library conformation."""
    return sum(min(x % 360.0, 360.0 - x % 360.0) for x in angles)


def _scan(ch: _Chain) -> Optional[List[float]]:
    """Torsion changes for a contact-free chain closest to the library one
    (or with the fewest contacts); None when the library conformation is free.

    Bonds are scanned two at a time, on a CHAIN_STEP grid, in SCAN_ROUNDS
    passes from O2' outwards.
    """
    best = [0.0] * len(ch.bonds)
    best_p = ch.penalty(ch.conformation(best))
    if best_p == 0.0:
        return None
    grid = [CHAIN_STEP * i for i in range(int(round(360.0 / CHAIN_STEP)))]
    for _ in range(SCAN_ROUNDS):
        for s in range(0, len(ch.bonds), 2):
            pair = range(s, min(s + 2, len(ch.bonds)))
            cands = [[]]
            for _j in pair:
                cands = [c + [g] for c in cands for g in grid]
            # closest to the library first (stable, so ties keep the grid order):
            # once a contact-free candidate is reached no later one can improve,
            # which ends most scans after a few trials (a staggered methyl)
            cands.sort(key=_turn)
            for c in cands:
                trial = list(best)
                for j, g in zip(pair, c):
                    trial[j] = g
                if best_p <= TOL and _turn(trial) >= _turn(best):
                    break
                p = ch.penalty(ch.conformation(trial))
                if p < best_p - TOL or (abs(p - best_p) <= TOL and _turn(trial) < _turn(best)):
                    best, best_p = trial, p
    return best


def _substituent_tree(nb: Dict[str, List[str]]):
    """Atoms hanging from O2' (not through C2') as a tree: (parent, atoms in BFS
    order), or None when there are none or they close a ring (LNA/cEt)."""
    parent = {"O2'": None}
    order: List[str] = []
    queue = ["O2'"]
    while queue:
        a = queue.pop(0)
        for b in nb.get(a, []):
            if a == "O2'" and b == "C2'":
                continue
            if b in parent:
                if parent[a] != b:
                    return None
                continue
            parent[b] = a
            order.append(b)
            queue.append(b)
    return (parent, order) if order else None


def _ancestor(parent: Dict[str, Optional[str]], atom: str, anc: str) -> bool:
    """Whether ``anc`` lies on the path from ``atom`` back to the root."""
    a = parent.get(atom)
    while a is not None:
        if a == anc:
            return True
        a = parent.get(a)
    return False


def rotate_about(points: np.ndarray, origin: np.ndarray, axis: np.ndarray,
                 deg: float) -> np.ndarray:
    """Rotate points by ``deg`` about the axis through ``origin`` (Rodrigues)."""
    k = axis / np.linalg.norm(axis)
    t = math.radians(deg)
    p = points - origin
    return (origin + p * math.cos(t) + np.cross(k, p) * math.sin(t)
            + np.outer(p @ k, k) * (1.0 - math.cos(t)))
