# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.graph
==================

Walks on the bond graph of a residue, given as atom name -> bonded atom names.
"""

from __future__ import annotations

from collections import deque
from typing import Collection, Dict, List, Optional, Sequence, Set

Graph = Dict[str, List[str]]


def within(nb: Graph, start: str, depth: int,
           through: Optional[Collection[str]] = None) -> List[str]:
    """Atoms within ``depth`` bonds of ``start``, nearest first.

    ``start`` itself is not included. Neighbours are visited in name order, so
    the result does not depend on the order of the bonds. With ``through``,
    paths only run over those atoms (the heavy atoms of a residue, say).
    """
    dist = {start: 0}
    queue = deque([start])
    while queue:
        a = queue.popleft()
        if dist[a] == depth:
            continue
        for b in sorted(nb.get(a, [])):
            if b not in dist and (through is None or b in through):
                dist[b] = dist[a] + 1
                queue.append(b)
    return [a for a, d in sorted(dist.items(), key=lambda kv: kv[1]) if d > 0]


def far_from(nb: Graph, sources: Sequence[str], n: int) -> Set[str]:
    """Atoms at least ``n`` bonds away from every atom in ``sources``."""
    near = set(sources)
    for s in sources:
        near.update(within(nb, s, n - 1))
    return {a for a in nb if a not in near}


def ring_atoms(nb: Graph) -> Set[str]:
    """Atoms on a ring: two of their neighbours stay connected without them."""
    out = set()
    for a, ns in nb.items():
        reached: Set[str] = set()
        for s in ns:
            if s in reached:            # joined to an earlier neighbour around a
                out.add(a)
                break
            seen, queue = {s}, deque([s])
            while queue:
                x = queue.popleft()
                for y in nb.get(x, []):
                    if y != a and y not in seen:
                        seen.add(y)
                        queue.append(y)
            reached |= seen
    return out
