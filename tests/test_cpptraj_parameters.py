# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""geometry.py reproduces 3DNA base-pair and step parameters.

Reference values come from cpptraj ``nastruct`` (3DNA conventions) on
duplexes built on eight fiber models (tests/data/cpptraj/README.md). cpptraj
writes four decimals, so agreement is required within 5e-4 (A or deg).
"""

import os

import numpy as np
import pytest

from dsxnabuilder.analysis import duplex_parameters
from dsxnabuilder.structure import read_pdb

DATA = os.path.join(os.path.dirname(__file__), "data", "cpptraj")
CASES = ["adna1", "arna20", "bdna4", "cdna7", "ddna14", "hyb22", "sdna49", "zdna15"]
TOL = 5e-4


def table(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                rows.append([float(x) for x in line.split()[3:9]])
    return np.array(rows)


@pytest.fixture(scope="module")
def standard(lib):
    return {k: dict(zip(v.names, v.xyz)) for k, v in lib.standard.items()}


@pytest.mark.parametrize("case", CASES)
def test_matches_cpptraj_nastruct(case, standard):
    atoms = read_pdb(os.path.join(DATA, f"{case}.pdb"))
    pairs, steps, _ = duplex_parameters(atoms, standard)
    ref_bp = table(os.path.join(DATA, f"{case}.BP.out"))
    ref_st = table(os.path.join(DATA, f"{case}.BPstep.out"))
    got_bp = np.array([p.as_tuple() for p in pairs])
    got_st = np.array([s.as_tuple() for s in steps])
    assert got_bp.shape == ref_bp.shape and got_st.shape == ref_st.shape
    assert np.abs(got_bp - ref_bp).max() < TOL, np.abs(got_bp - ref_bp).max(0)
    assert np.abs(got_st - ref_st).max() < TOL, np.abs(got_st - ref_st).max(0)


def test_z_dna_reported_as_3dna_does(standard):
    """Left-handed Z-DNA: positive rise, negative twist (pairs reoriented)."""
    _, steps, _ = duplex_parameters(read_pdb(os.path.join(DATA, "zdna15.pdb")), standard)
    assert all(s.rise > 3.0 for s in steps)
    twists = sorted({round(s.twist) for s in steps})
    assert twists == [-51, -9]
