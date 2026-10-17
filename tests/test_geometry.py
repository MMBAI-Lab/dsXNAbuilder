# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Geometric formalism.

The decisive test is not the round trip (any self-consistent forward/inverse
pair passes it) but the *fiber symmetry constraint*: an A or B duplex has a
two-fold axis perpendicular to the helix axis through every base pair, which
forces tilt and shift of every step to vanish. This tells apart decomposition
conventions that otherwise look equally valid.
"""

import math

import numpy as np
import pytest

from dsxnabuilder.geometry import (Frame, PairParams, StepParams, apply_step,
                                   base_frames_to_pair, fiber_frames,
                                   frames_to_step_params, kabsch,
                                   pair_to_base_frames, pseudorotation, rot_x,
                                   rot_z, step_rotation, PUCKER_ENVELOPES)

FIBER_CASES = {
    # name: (twist, rise, x_displacement, inclination)
    "B-DNA": (36.0, 3.38, 0.0, 2.0),
    "A-DNA": (32.7, 2.56, -4.17, 19.1),
    "A-RNA": (32.7, 2.81, -4.40, 15.8),
}


def _random_proper_rotation(rng):
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 2] *= -1
    return Q


def test_step_roundtrip_exact():
    rng = np.random.default_rng(20260101)
    worst = 0.0
    for _ in range(2000):
        sp = StepParams(shift=rng.uniform(-3, 3), slide=rng.uniform(-3, 3),
                        rise=rng.uniform(2, 5), tilt=rng.uniform(-20, 20),
                        roll=rng.uniform(-25, 25), twist=rng.uniform(10, 60))
        f0 = Frame(rng.uniform(-5, 5, 3), _random_proper_rotation(rng))
        back = frames_to_step_params(f0, apply_step(f0, sp))
        worst = max(worst, float(np.abs(np.array(sp.as_tuple())
                                        - np.array(back.as_tuple())).max()))
    # measured (CEHS/3DNA convention): 3.3e-11 over 5000 cases with twist in (-170, 170)
    assert worst < 1e-10, f"round trip not exact: {worst}"


def test_pair_roundtrip_exact():
    rng = np.random.default_rng(7)
    for _ in range(500):
        pp = PairParams(shear=rng.uniform(-1, 1), stretch=rng.uniform(-1, 1),
                        stagger=rng.uniform(-1, 1), buckle=rng.uniform(-20, 20),
                        propeller=rng.uniform(-25, 25), opening=rng.uniform(-15, 15))
        bp = Frame(rng.uniform(-5, 5, 3), _random_proper_rotation(rng))
        f1, f2 = pair_to_base_frames(bp, pp)
        bp_back, pp_back = base_frames_to_pair(f1, f2)
        assert np.allclose(bp.origin, bp_back.origin, atol=1e-9)
        assert np.allclose(bp.triad, bp_back.triad, atol=1e-9)
        assert np.allclose(np.array(pp.as_tuple()),
                           np.array(pp_back.as_tuple()), atol=1e-9)


@pytest.mark.parametrize("name", sorted(FIBER_CASES))
def test_fiber_dyad_symmetry(name):
    """tilt = 0 and shift = 0 for every uniform helix with a two-fold axis."""
    twist, rise, dx, inc = FIBER_CASES[name]
    frames = fiber_frames(24, twist, rise, dx, inc)
    for i in range(len(frames) - 1):
        sp = frames_to_step_params(frames[i], frames[i + 1])
        # measured (CEHS/3DNA convention): max|tilt| = 1.8e-15 deg, max|shift| = 4.2e-15 A
        assert abs(sp.tilt) < 1e-12, f"{name}: tilt = {sp.tilt}"
        assert abs(sp.shift) < 1e-12, f"{name}: shift = {sp.shift}"


@pytest.mark.parametrize("name", sorted(FIBER_CASES))
def test_fiber_is_exact_screw(name):
    """The frames are the exact orbit of the screw operator S(Omega, h)."""
    twist, rise, dx, inc = FIBER_CASES[name]
    frames = fiber_frames(12, twist, rise, dx, inc)
    S = rot_z(twist)
    for i in range(len(frames) - 1):
        assert np.allclose(frames[i + 1].triad, S @ frames[i].triad, atol=1e-12)
        assert np.allclose(frames[i + 1].origin,
                           S @ frames[i].origin + np.array([0, 0, rise]),
                           atol=1e-12)


def test_fiber_step_params_are_uniform():
    """All steps of a fiber have the same parameters."""
    frames = fiber_frames(12, 32.7, 2.81, -4.40, 15.8)
    vals = np.array([frames_to_step_params(frames[i], frames[i + 1]).as_tuple()
                     for i in range(len(frames) - 1)])
    assert vals.std(axis=0).max() < 1e-9


def test_a_form_step_parameters_are_physical():
    """Fiber A-RNA has a large negative slide and a positive roll.

    That is the known signature of the A form, and a check that the sign
    convention is not inverted.
    """
    sp = frames_to_step_params(*fiber_frames(2, 32.7, 2.81, -4.40, 15.8))
    assert -2.5 < sp.slide < -1.0, sp.slide
    assert 5.0 < sp.roll < 15.0, sp.roll
    assert 3.0 < sp.rise < 3.6, sp.rise
    assert 29.0 < sp.twist < 33.0, sp.twist


def test_b_form_step_parameters_are_physical():
    sp = frames_to_step_params(*fiber_frames(2, 36.0, 3.38, 0.0, 2.0))
    assert abs(sp.slide) < 0.5, sp.slide
    assert abs(sp.roll) < 4.0, sp.roll
    assert 3.2 < sp.rise < 3.5, sp.rise
    assert 34.0 < sp.twist < 37.0, sp.twist


def test_complementary_strand_frame_is_proper():
    """The Rx(180) turn of strand II does not invert handedness."""
    f1, f2 = pair_to_base_frames(Frame([0, 0, 0], rot_z(17.0)),
                                 PairParams(propeller=-13.7, buckle=3.0))
    assert np.isclose(np.linalg.det(f1.triad), 1.0)
    assert np.isclose(np.linalg.det(f2.triad), 1.0)


def test_kabsch_never_inverts():
    rng = np.random.default_rng(3)
    mob = rng.normal(size=(8, 3))
    ref = -mob                      # mirror image
    R, t, rmsd = kabsch(mob, ref)
    assert np.isclose(np.linalg.det(R), 1.0), "Kabsch returned a reflection"


@pytest.mark.parametrize("model_id, envelopes", [
    (1, {"C3'-endo"}), (4, {"C2'-endo"}), (20, {"C3'-endo"}),
    (15, {"C3'-endo", "C2'-endo"}),          # Z-DNA: G C3'-endo (P ~ 28), C C2'-endo
])
def test_pseudorotation_of_fiber_sugars(lib, model_id, envelopes):
    """Domains and envelopes as 3DNA and cpptraj name them: P in [0, 36) is
    C3'-endo, [36, 72) C4'-exo, ..., [144, 180) C2'-endo. The pucker stored
    with each model slot is the one of its coordinates."""
    m = lib.get(model_id)
    got = set()
    for bp, slot in enumerate(m.slots):
        names, xyz, _ = m.nucleotide(bp, "I", slot.I.native, m.strand_chem("I"), lib.standard)
        pr = pseudorotation(dict(zip(names, xyz)))
        assert pr["envelope"] == PUCKER_ENVELOPES[int(pr["P"] // 36)]
        assert pr["domain"] == ("North" if pr["P"] < 90 or pr["P"] > 270 else "South")
        assert (pr["envelope"], pr["domain"]) == (slot.I.pucker["envelope"],
                                                  slot.I.pucker["domain"])
        got.add(pr["envelope"])
    assert got == envelopes


def test_pseudorotation_handles_4_thio():
    coords = {"C1'": [0, 0, 0], "C2'": [1.5, 0, 0], "C3'": [2.2, 1.3, 0.3],
              "C4'": [1.2, 2.4, -0.1], "S4'": [0.0, 1.5, 0.4]}
    out = pseudorotation(coords)
    assert out is not None and out["ring_heteroatom"] == "S4'"
    assert pseudorotation({"C1'": [0, 0, 0]}) is None


# --------------------------------------------------------------------------
# 3DNA (CEHS) convention: signs, hinge middle frame, flipped pairs
# --------------------------------------------------------------------------
def test_sign_convention_matches_3dna():
    """tilt > 0 turns about +x and roll > 0 about +y (El Hassan & Calladine
    1995; 3DNA); buckle > 0 turns base I about +x relative to base II."""
    s = math.sin(math.radians(10.0))
    z = apply_step(Frame(), StepParams(tilt=10.0, twist=0.0, rise=0.0)).triad[:, 2]
    assert np.allclose(z, [0.0, -s, math.cos(math.radians(10.0))])
    z = apply_step(Frame(), StepParams(roll=10.0, twist=0.0, rise=0.0)).triad[:, 2]
    assert np.allclose(z, [s, 0.0, math.cos(math.radians(10.0))])
    f1, f2 = pair_to_base_frames(Frame(), PairParams(buckle=10.0))
    rel = (f2.triad @ rot_x(180.0)).T @ f1.triad
    assert np.allclose(rel[:, 2], [0.0, -s, math.cos(math.radians(10.0))])


def test_middle_step_frame_is_hinge_bisector():
    """The middle frame bisects the z axes and contains the hinge axis."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        sp = StepParams(shift=rng.uniform(-2, 2), slide=rng.uniform(-2, 2),
                        rise=rng.uniform(2, 4), tilt=rng.uniform(-15, 15),
                        roll=rng.uniform(-20, 20), twist=rng.uniform(-60, 60))
        f0 = Frame(rng.uniform(-5, 5, 3), _random_proper_rotation(rng))
        f1 = apply_step(f0, sp)
        _, M = step_rotation(sp)
        T = f0.triad @ M
        zsum = f0.triad[:, 2] + f1.triad[:, 2]
        assert np.allclose(T[:, 2], zsum / np.linalg.norm(zsum), atol=1e-12)
        hinge = np.cross(f0.triad[:, 2], f1.triad[:, 2])
        assert abs(np.dot(hinge, T[:, 2])) < 1e-12
        # the step translation is expressed in that frame
        d = T.T @ (f1.origin - f0.origin)
        assert np.allclose(d, [sp.shift, sp.slide, sp.rise], atol=1e-10)


def test_flipped_pair_is_reoriented_along_strand_direction():
    """Flipped pairs (Z-DNA): the pair frame turns 180 degrees about y so that
    its z follows the 5'->3' direction of strand I; shear, stagger, buckle and
    opening change sign."""
    pp = PairParams(shear=0.4, stretch=-0.1, stagger=0.2, buckle=5.0,
                    propeller=-12.0, opening=3.0)
    bp = Frame([1.0, 2.0, 3.0], rot_z(40.0))
    f1, f2 = pair_to_base_frames(bp, pp)
    same, p_same = base_frames_to_pair(f1, f2, direction=[0, 0, 1])
    assert np.allclose(same.triad, bp.triad) and np.isclose(p_same.shear, 0.4)
    flip, p_flip = base_frames_to_pair(f1, f2, direction=[0, 0, -1])
    assert np.isclose(np.dot(flip.triad[:, 2], [0, 0, -1]), 1.0)
    got = np.array(p_flip.as_tuple())
    assert np.allclose(got, [-0.4, -0.1, -0.2, -5.0, -12.0, -3.0])
