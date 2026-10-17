# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.geometry
=====================

Geometric formalism: reference frames, base-pair and step parameters.

Conventions
-----------
The *standard reference frame* of Olson et al. (2001) J Mol Biol 313:229-237
(doi:10.1006/jmbi.2001.4987) and the base-pair / step parameter algebra of
El Hassan & Calladine (1995) J Mol Biol 251:648-664
(doi:10.1006/jmbi.1995.0462) as implemented in 3DNA (Lu & Olson 2003 NAR
31:5108-5121), i.e. the CEHS scheme. All rotations are active, in
right-handed systems; a triad ``T`` holds the (x, y, z) unit vectors of the
frame as columns in global coordinates, so that

    r_global = o + T . r_local

A ``Frame`` is the pair (origin o, triad T).

Step parameters
---------------
For (shift Dx, slide Dy, rise Dz, tilt tau, roll rho, twist Omega):

    Gamma = sqrt(tau^2 + rho^2)          bending angle (roll-tilt)
    phi   = atan2(tau, rho)              phase of the bending

    R_step = Rz(Omega/2 - phi) . Ry(Gamma) . Rz(Omega/2 + phi)
    M      = Rz(Omega/2 - phi) . Ry(Gamma/2) . Rz(phi)

    T_{i+1} = T_i . R_step
    T_mst   = T_i . M                    (middle-step triad)
    o_{i+1} = o_i + T_mst . (Dx, Dy, Dz)^T

M is the CEHS "hinge" middle frame: both frames are bent by -/+ Gamma/2 about
the hinge axis (the roll-tilt axis, which lies in the xy plane of T_mst at
angle phi from its y axis) so that their z axes coincide, and T_mst is their
average. Positive tilt is a rotation about +x of T_mst and positive roll about
+y. For Gamma = 0 the step is the screw operator S(Omega, Dz) of a fiber
model.

The inverse decomposes T_i^T T_{i+1} = Rz(alpha) Ry(Gamma) Rz(beta) (ZYZ),
with Omega = alpha + beta and phi = (beta - alpha)/2 on the branch
Omega in (-180, 180] (see ``_unwrap_twist_phase``), and rebuilds M from
(Omega, Gamma, phi) to express o_{i+1} - o_i.

Base-pair parameters
--------------------
Same algebra between the two bases of a pair. With base II "unflipped" (its
triad times Rx(180), so that both bases point the same way):

    Gamma_bp = sqrt(kappa^2 + omega^2),   phi_bp = atan2(kappa, omega)
    R_pair   = Rz(sigma/2 - phi_bp) . Ry(Gamma_bp) . Rz(sigma/2 + phi_bp)
             = T_II,unflipped^T . T_I
    T_bp     = T_II,unflipped . Rz(sigma/2 - phi_bp) . Ry(Gamma_bp/2) . Rz(phi_bp)
    o_bp     = (o_I + o_II) / 2,    (Sx, Sy, Sz) = T_bp^T (o_I - o_II)

for (shear Sx, stretch Sy, stagger Sz, buckle kappa, propeller omega,
opening sigma).

Orientation of the pair frame
-----------------------------
3DNA reports parameters with the base-pair z axis along the 5'->3'
direction of strand I. When the bases are flipped relative to that direction
(Z-DNA), ``base_frames_to_pair(..., direction=...)`` turns the pair frame by
180 deg about its y axis, which changes the sign of shear, stagger, buckle
and opening; step parameters between such frames then give the familiar
positive rise and negative twist of Z-DNA.

Verification
------------
Pair and step parameters of duplexes built on fiber models agree with cpptraj
``nastruct`` (3DNA conventions) within its four-decimal output
(tests/test_cpptraj_parameters.py). A uniform helix with a dyad axis through
each base pair has tau = 0 and Dx = 0, and the forward and inverse operations
are exact inverses (tests/test_geometry.py).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

__all__ = [
    "rot_x", "rot_y", "rot_z", "Frame", "StepParams", "PairParams",
    "step_rotation", "apply_step", "pair_rotation", "pair_to_base_frames",
    "frames_to_step_params", "base_frames_to_pair", "fiber_frames", "kabsch",
    "dihedral", "signed_volume", "pseudorotation", "PUCKER_ENVELOPES",
]

DEG = math.pi / 180.0


# --------------------------------------------------------------------------
# Elementary rotations (angles in degrees)
# --------------------------------------------------------------------------
def rot_x(a_deg: float) -> np.ndarray:
    c, s = math.cos(a_deg * DEG), math.sin(a_deg * DEG)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)


def rot_y(a_deg: float) -> np.ndarray:
    c, s = math.cos(a_deg * DEG), math.sin(a_deg * DEG)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)


def rot_z(a_deg: float) -> np.ndarray:
    c, s = math.cos(a_deg * DEG), math.sin(a_deg * DEG)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


# --------------------------------------------------------------------------
# Frames and parameters
# --------------------------------------------------------------------------
@dataclass
class Frame:
    """Right-handed orthonormal reference frame (origin, triad)."""

    origin: np.ndarray = field(default_factory=lambda: np.zeros(3))
    triad: np.ndarray = field(default_factory=lambda: np.eye(3))

    def __post_init__(self) -> None:
        self.origin = np.asarray(self.origin, dtype=float).reshape(3)
        self.triad = np.asarray(self.triad, dtype=float).reshape(3, 3)

    def to_global(self, local: np.ndarray) -> np.ndarray:
        """Local coordinates, (N, 3) or (3,), in the global system."""
        local = np.asarray(local, dtype=float)
        return (self.triad @ local.reshape(-1, 3).T).T + self.origin

    def to_local(self, glob: np.ndarray) -> np.ndarray:
        """Global coordinates, (N, 3) or (3,), in the frame."""
        glob = np.asarray(glob, dtype=float)
        return (self.triad.T @ (glob.reshape(-1, 3) - self.origin).T).T


@dataclass
class StepParams:
    """Step parameters between consecutive base pairs (angstrom, degrees)."""

    shift: float = 0.0   # Dx
    slide: float = 0.0   # Dy
    rise: float = 0.0    # Dz
    tilt: float = 0.0    # tau
    roll: float = 0.0    # rho
    twist: float = 0.0   # Omega

    def as_tuple(self) -> Tuple[float, ...]:
        return (self.shift, self.slide, self.rise, self.tilt, self.roll, self.twist)


@dataclass
class PairParams:
    """Parameters between the two bases of a pair (angstrom, degrees)."""

    shear: float = 0.0      # Sx
    stretch: float = 0.0    # Sy
    stagger: float = 0.0    # Sz
    buckle: float = 0.0     # kappa
    propeller: float = 0.0  # omega
    opening: float = 0.0    # sigma

    def as_tuple(self) -> Tuple[float, ...]:
        return (self.shear, self.stretch, self.stagger,
                self.buckle, self.propeller, self.opening)


# --------------------------------------------------------------------------
# Step: forward and inverse
# --------------------------------------------------------------------------
def _bend(a: float, b: float) -> Tuple[float, float]:
    """(Gamma, phi) of a bend from its two components: (tilt, roll) for a
    step, (buckle, propeller) for a pair. Degrees."""
    gamma = math.hypot(a, b)
    phi = math.degrees(math.atan2(a, b)) if gamma > 1e-12 else 0.0
    return gamma, phi


def _bend_components(gamma: float, phi: float) -> Tuple[float, float]:
    """Inverse of :func:`_bend`: (tilt, roll) or (buckle, propeller)."""
    return gamma * math.sin(phi * DEG), gamma * math.cos(phi * DEG)


def _cehs(omega: float, gamma: float, phi: float) -> Tuple[np.ndarray, np.ndarray]:
    """(R, M): R = Rz(w/2 - phi) Ry(g) Rz(w/2 + phi), M = Rz(w/2 - phi) Ry(g/2) Rz(phi).

    M is the CEHS hinge middle frame relative to the first frame (module
    docstring); angles in degrees.
    """
    R = rot_z(omega / 2.0 - phi) @ rot_y(gamma) @ rot_z(omega / 2.0 + phi)
    M = rot_z(omega / 2.0 - phi) @ rot_y(gamma / 2.0) @ rot_z(phi)
    return R, M


def step_rotation(sp: StepParams) -> Tuple[np.ndarray, np.ndarray]:
    """Return (R_step, M) for a step: T_{i+1} = T_i R_step, T_mst = T_i M."""
    gamma, phi = _bend(sp.tilt, sp.roll)
    return _cehs(sp.twist, gamma, phi)


def apply_step(frame: Frame, sp: StepParams) -> Frame:
    """Propagate a base-pair frame by one step."""
    R_step, M = step_rotation(sp)
    T_mst = frame.triad @ M
    origin = frame.origin + T_mst @ np.array([sp.shift, sp.slide, sp.rise])
    return Frame(origin, frame.triad @ R_step)


def _zyz_from_matrix(R: np.ndarray) -> Tuple[float, float, float]:
    """(alpha, Gamma, beta) with R = Rz(alpha) Ry(Gamma) Rz(beta), in degrees."""
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, R[2, 2]))))
    if abs(math.sin(gamma * DEG)) < 1e-10:
        # no bend: the whole rotation is about z
        alpha = 0.0
        if R[2, 2] > 0:
            beta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
        else:
            beta = math.degrees(math.atan2(-R[1, 0], -R[0, 0]))
    else:
        alpha = math.degrees(math.atan2(R[1, 2], R[0, 2]))
        beta = math.degrees(math.atan2(R[2, 1], -R[2, 0]))
    return alpha, gamma, beta


def _unwrap_twist_phase(alpha: float, beta: float) -> Tuple[float, float]:
    """Resolve the degeneracy (Gamma, phi, Omega) ~ (Gamma, phi - 180, Omega + 360).

    The ZYZ decomposition returns alpha and beta in (-180, 180], so
    Omega = alpha + beta can come out 360 deg away from the physical twist.
    Shifting alpha by -360 k (the same rotation) brings Omega into
    (-180, 180] and shifts phi = (beta - alpha)/2 by +180 k, which keeps
    R = Rz(Omega/2 - phi) Ry(Gamma) Rz(Omega/2 + phi) exact and flips the
    signs of tilt and roll only where they were spurious. This branch is also
    the one on which M (the hinge middle frame) is the bisecting frame.
    """
    omega_raw = alpha + beta
    k = math.floor((omega_raw + 180.0) / 360.0)
    return omega_raw - 360.0 * k, (beta - alpha) / 2.0 + 180.0 * k


def _decompose(R: np.ndarray) -> Tuple[float, float, float]:
    """(Omega, Gamma, phi) with R = Rz(Omega/2 - phi) Ry(Gamma) Rz(Omega/2 + phi)."""
    alpha, gamma, beta = _zyz_from_matrix(R)
    omega, phi = _unwrap_twist_phase(alpha, beta)
    return omega, gamma, phi


def frames_to_step_params(f_i: Frame, f_j: Frame) -> StepParams:
    """Exact inverse of :func:`apply_step` (3DNA/CEHS step parameters)."""
    omega, gamma, phi = _decompose(f_i.triad.T @ f_j.triad)
    _, M = _cehs(omega, gamma, phi)
    d = (f_i.triad @ M).T @ (f_j.origin - f_i.origin)
    tilt, roll = _bend_components(gamma, phi)
    return StepParams(shift=float(d[0]), slide=float(d[1]), rise=float(d[2]),
                      tilt=float(tilt), roll=float(roll), twist=float(omega))


# --------------------------------------------------------------------------
# Base pair: forward and inverse
# --------------------------------------------------------------------------
def pair_rotation(pp: PairParams) -> Tuple[np.ndarray, np.ndarray]:
    """Return (R_pair, M): R_pair = T_IIu^T T_I and T_bp = T_IIu M (T_IIu = unflipped base II)."""
    gamma, phi = _bend(pp.buckle, pp.propeller)
    return _cehs(pp.opening, gamma, phi)


def pair_to_base_frames(bp: Frame, pp: PairParams) -> Tuple[Frame, Frame]:
    """Frames of the two bases from the base-pair frame and parameters."""
    R_pair, M = pair_rotation(pp)
    T_IIu = bp.triad @ M.T
    disp = bp.triad @ (np.array([pp.shear, pp.stretch, pp.stagger]) / 2.0)
    f1 = Frame(bp.origin + disp, T_IIu @ R_pair)
    f2 = Frame(bp.origin - disp, T_IIu @ rot_x(180.0))
    return f1, f2


def base_frames_to_pair(f1: Frame, f2: Frame,
                        direction: Sequence[float] | None = None) -> Tuple[Frame, PairParams]:
    """Inverse of :func:`pair_to_base_frames`: base-pair frame and parameters.

    ``f2`` is the strand II base frame as fitted (z antiparallel to base I);
    the Rx(180) flip is undone first, as ``nastruct``/3DNA do. If
    ``direction`` (the 5'->3' direction of strand I) is given and the pair's
    z axis points against it, the pair frame is turned by 180 deg about its y
    axis (3DNA convention for flipped pairs, e.g. Z-DNA).
    """
    T2u = f2.triad @ rot_x(180.0)
    omega, gamma, phi = _decompose(T2u.T @ f1.triad)
    _, M = _cehs(omega, gamma, phi)
    T_bp = T2u @ M
    o_bp = (f1.origin + f2.origin) / 2.0
    d = T_bp.T @ (f1.origin - f2.origin)
    shear, stretch, stagger = float(d[0]), float(d[1]), float(d[2])
    buckle, propeller = (float(x) for x in _bend_components(gamma, phi))
    opening = float(omega)
    if direction is not None and float(np.dot(T_bp[:, 2], np.asarray(direction, float))) < 0:
        T_bp = T_bp @ rot_y(180.0)
        shear, stagger, buckle, opening = -shear, -stagger, -buckle, -opening
    pp = PairParams(shear=shear, stretch=stretch, stagger=stagger, buckle=buckle,
                    propeller=propeller, opening=opening)
    return Frame(o_bp, T_bp), pp


# --------------------------------------------------------------------------
# Fiber models: exact screw symmetry
# --------------------------------------------------------------------------
def fiber_frames(n_bp: int, twist: float, rise: float,
                 x_displacement: float = 0.0,
                 inclination: float = 0.0,
                 tip: float = 0.0) -> List[Frame]:
    """Base-pair frames of a fiber model with screw symmetry.

    The helix axis is the global z axis. Frame i is

        T_i = Rz(i.Omega) . Ry(tip) . Rx(inclination)
        o_i = Rz(i.Omega) . (dx, 0, 0)^T + (0, 0, i.h)^T

    With inclination = tip = 0 and dx = 0 this is the orbit of the identity
    frame under the screw operator S(Omega, h).
    """
    base_T = rot_y(tip) @ rot_x(inclination)
    frames = []
    for i in range(n_bp):
        Rz = rot_z(i * twist)
        o = Rz @ np.array([x_displacement, 0.0, 0.0]) + np.array([0.0, 0.0, i * rise])
        frames.append(Frame(o, Rz @ base_T))
    return frames


# --------------------------------------------------------------------------
# Superposition and torsions
# --------------------------------------------------------------------------
def kabsch(mobile: np.ndarray, ref: np.ndarray,
           weights: Sequence[float] | None = None) -> Tuple[np.ndarray, np.ndarray, float]:
    """Least-squares superposition: (R, t, rmsd) minimising
    sum_a w_a ||R m_a + t - r_a||^2.

    R is a proper rotation (det R = +1), so the fit never inverts chirality.
    """
    m = np.asarray(mobile, float).reshape(-1, 3)
    r = np.asarray(ref, float).reshape(-1, 3)
    if m.shape != r.shape or m.shape[0] < 3:
        raise ValueError(f"need at least 3 matching points, got {m.shape} and {r.shape}")
    w = np.ones(len(m)) if weights is None else np.asarray(weights, float)
    w = w / w.sum()
    mc = (w[:, None] * m).sum(0)
    rc = (w[:, None] * r).sum(0)
    H = (w[:, None] * (m - mc)).T @ (r - rc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = rc - R @ mc
    diff = (R @ m.T).T + t - r
    rmsd = float(np.sqrt((w * (diff ** 2).sum(1)).sum()))
    return R, t, rmsd


def dihedral(p0, p1, p2, p3) -> float:
    """Torsion angle p0-p1-p2-p3 in degrees, in (-180, 180]."""
    p0, p1, p2, p3 = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def signed_volume(c: Dict[str, np.ndarray], centre: str, a: str, b: str, d: str) -> float:
    """(a - centre) x (b - centre) . (d - centre): its sign is the handedness."""
    o = c[centre]
    return float(np.dot(np.cross(c[a] - o, c[b] - o), c[d] - o))


# --------------------------------------------------------------------------
# Sugar pseudorotation (Altona & Sundaralingam 1972)
# --------------------------------------------------------------------------
# Names of the ten 36-degree sectors of P starting at 0, as 3DNA and cpptraj
# label them: [0, 36) is C3'-endo, [144, 180) C2'-endo.
PUCKER_ENVELOPES = ["C3'-endo", "C4'-exo", "O4'-endo", "C1'-exo", "C2'-endo",
                    "C3'-exo", "C4'-endo", "O4'-exo", "C1'-endo", "C2'-exo"]


def pseudorotation(coords: Dict[str, Iterable[float]]) -> dict | None:
    """Phase P, amplitude nu_max and North/South domain of a furanose ring.

    Rings with O4' or S4' (4'-thio sugars) are accepted. Returns None when
    there is no furanose ring (an acyclic sugar, a morpholino).
    """
    ring_o = "O4'" if "O4'" in coords else ("S4'" if "S4'" in coords else None)
    if ring_o is None or not all(k in coords for k in ("C1'", "C2'", "C3'", "C4'")):
        return None
    c = {k: np.asarray(v, float) for k, v in coords.items()}
    C1, C2, C3, C4, O4 = c["C1'"], c["C2'"], c["C3'"], c["C4'"], c[ring_o]
    nu = [dihedral(C4, O4, C1, C2), dihedral(O4, C1, C2, C3), dihedral(C1, C2, C3, C4),
          dihedral(C2, C3, C4, O4), dihedral(C3, C4, O4, C1)]
    num = (nu[4] + nu[1]) - (nu[3] + nu[0])
    den = 2.0 * nu[2] * (math.sin(36 * DEG) + math.sin(72 * DEG))
    P = math.degrees(math.atan2(num, den)) % 360.0
    numax = nu[2] / math.cos(P * DEG) if abs(math.cos(P * DEG)) > 1e-6 else float("nan")
    domain = "North" if (P < 90.0 or P > 270.0) else "South"
    return {"nu": [round(x, 2) for x in nu], "P": round(P, 1),
            "nu_max": round(numax, 1), "ring_heteroatom": ring_o,
            "envelope": PUCKER_ENVELOPES[int(P // 36) % 10],
            "domain": domain}
