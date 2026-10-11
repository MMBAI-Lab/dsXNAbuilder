"""
dsxnabuilder.analysis
=====================

Base-pair and step parameters of a built duplex, in the 3DNA (CEHS)
convention of :mod:`dsxnabuilder.geometry`.

Base frames are fitted on the ring atoms of the Olson et al. (2001) standard
bases; base pairs are taken from the duplex layout written by the builder
(chain A residue k pairs with chain B residue n-1-k, both chains 5'->3' and
of equal length); the pair frames are oriented along the 5'->3' direction of
chain A, so flipped pairs (Z-DNA) are reported as 3DNA does. Validated
against cpptraj ``nastruct`` in tests/test_cpptraj_parameters.py.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from .geometry import Frame, PairParams, StepParams, base_frames_to_pair, frames_to_step_params
from .structure import Atom, residues, standard_base_frame

RING_LETTER = {"A": "A", "G": "G", "C": "C", "T": "T", "U": "U", "I": "G"}


def _letter(resname: str) -> str:
    return resname.strip()[-1].upper()


def duplex_parameters(atoms: Sequence[Atom], standard: Dict[str, Dict[str, np.ndarray]]
                      ) -> Tuple[List[PairParams], List[StepParams], List[Frame]]:
    """(pair parameters, step parameters, pair frames) along chain A.

    ``standard`` maps A/C/G/T/U to the standard base coordinates (name -> xyz).
    Residues whose name does not end in a base letter (modXNA residues named
    X01, ...) are not supported here.
    """
    by_chain: Dict[str, List[List[Atom]]] = {}
    for (chain, _), ats in residues(atoms).items():
        by_chain.setdefault(chain, []).append(ats)
    a, b = by_chain["A"], by_chain["B"]
    if len(a) != len(b):
        raise ValueError("duplex_parameters needs two chains of equal length")
    n = len(a)

    def frame(ats: List[Atom]) -> Frame:
        c = {x.name: x.xyz for x in ats}
        return standard_base_frame(c, standard[RING_LETTER[_letter(ats[0].resname)]])[0]

    c1 = [next(x.xyz for x in ats if x.name == "C1'") for ats in a]
    direction = c1[-1] - c1[0]
    pairs, frames = [], []
    for k in range(n):
        f, pp = base_frames_to_pair(frame(a[k]), frame(b[n - 1 - k]), direction=direction)
        frames.append(f)
        pairs.append(pp)
    steps = [frames_to_step_params(frames[k], frames[k + 1]) for k in range(n - 1)]
    return pairs, steps, frames
