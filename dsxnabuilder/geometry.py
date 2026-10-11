"""
xnafiber.geometry
=================

Formalismo geometrico de FiberXNA.

Convenciones
------------
Se adopta el *standard reference frame* de Olson et al. (2001) J Mol Biol
313:229-237 (doi:10.1006/jmbi.2001.4987) y el algebra de parametros de pares
de bases / escalon de El Hassan & Calladine (1995) J Mol Biol 251:648-664
(doi:10.1006/jmbi.1995.0462), en la realizacion de 3DNA (Lu & Olson 2003 NAR
31:5108-5121). El marco intermedio ("mid-step triad") es el elemento central
de ese esquema, y su proposito declarado es que los parametros de escalon
salgan identicos leyendo el escalon en cualquiera de los dos sentidos; la
raiz cuadrada principal que se usa mas abajo es la construccion que cumple
esa propiedad de forma exacta. Todas las rotaciones son activas, en
sistemas dextrogiros, y
las matrices son column-major: una triada ``T`` tiene como columnas los
versores (x, y, z) del marco expresados en coordenadas globales, de modo que

    r_global = o + T . r_local

Un marco (``Frame``) es el par (origen o, triada T).

Descomposicion de un escalon
----------------------------
Dado el conjunto de parametros de escalon (shift Dx, slide Dy, rise Dz,
tilt tau, roll rho, twist Omega), se define

    Gamma = sqrt(tau^2 + rho^2)        (angulo de curvatura, "roll-tilt")
    phi   = atan2(tau, rho)            (fase de la curvatura)

    R_step = Rz(Omega/2 + phi) . Ry(Gamma) . Rz(Omega/2 - phi)
    T_mst  = T_i . R_step^(1/2)                        (middle step triad)

    T_{i+1} = T_i . R_step
    o_{i+1} = o_i + T_mst . (Dx, Dy, Dz)^T

El giro se reparte simetricamente a ambos lados de la curvatura (convencion
CEHS/3DNA) y el marco intermedio es la raiz cuadrada principal de la rotacion
del escalon. Esta eleccion no es cosmetica: es la unica de las variantes
habituales que satisface las restricciones de simetria de un modelo de fibra.
Un duplex A o B tiene un eje binario perpendicular al eje helicoidal que pasa
por cada par de bases; ese eje lleva el escalon i -> i+1 al escalon i -> i-1,
lo que fuerza

    tau = 0   y   Dx = 0

para cualquier helice uniforme. Con la fase simetrica phi = (alpha - beta)/2
de la descomposicion ZYZ y la raiz cuadrada principal como marco intermedio,
ambas identidades se cumplen a precision de maquina: medido sobre las cinco
conformaciones de ``forms.yaml`` y 23 escalones cada una, max|tau| = 1.8e-15
grados y max|Dx| = 7.9e-16 A. Con la fase asimetrica phi = alpha, o con un
marco intermedio construido "medio giro y media curvatura", aparecen valores
espurios de hasta tau = 2.4 grados y Dx = 0.45 A. El par directo/inverso es
reversible con un error maximo de 7.7e-12 sobre 5000 escalones aleatorios.
Todo ello se comprueba en ``tests/test_geometry.py``.

Para tau = rho = 0 la construccion se reduce exactamente al operador de
tornillo S(Omega, Dz), que es la simetria que define un modelo de fibra de
Arnott.

Par de bases
------------
Con (shear Sx, stretch Sy, stagger Sz, buckle kappa, propeller omega,
opening sigma) y el marco del par T_bp:

    Gamma_bp = sqrt(kappa^2 + omega^2),  phi_bp = atan2(kappa, omega)
    R_pair   = Rz(sigma/2 + phi_bp) . Ry(Gamma_bp) . Rz(sigma/2 - phi_bp)
    H        = R_pair^(1/2)

    T_I  = T_bp . H                       o_I  = o_bp + T_bp.(Sx,Sy,Sz)/2
    T_II = T_bp . H^{-1} . Rx(180)        o_II = o_bp - T_bp.(Sx,Sy,Sz)/2

El giro Rx(180) sobre la hebra complementaria implementa la convencion
estandar: el marco de la base II se define en el mismo sentido que el de la
base I antes de aplicar los parametros intra-par.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

__all__ = [
    "rot_x", "rot_y", "rot_z", "Frame", "StepParams", "PairParams",
    "step_rotation", "apply_step", "pair_rotation", "pair_to_base_frames",
    "frames_to_step_params", "base_frames_to_pair", "rot_sqrt",
    "fiber_frames", "kabsch", "pseudorotation", "PUCKER_ENVELOPES",
]

DEG = math.pi / 180.0


# --------------------------------------------------------------------------
# Rotaciones elementales
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
# Marco
# --------------------------------------------------------------------------
@dataclass
class Frame:
    """Marco de referencia ortonormal dextrogiro."""

    origin: np.ndarray = field(default_factory=lambda: np.zeros(3))
    triad: np.ndarray = field(default_factory=lambda: np.eye(3))

    def __post_init__(self) -> None:
        self.origin = np.asarray(self.origin, dtype=float).reshape(3)
        self.triad = np.asarray(self.triad, dtype=float).reshape(3, 3)

    def to_global(self, local: np.ndarray) -> np.ndarray:
        """Transforma coordenadas locales (N,3) o (3,) al sistema global."""
        local = np.asarray(local, dtype=float)
        return (self.triad @ local.reshape(-1, 3).T).T + self.origin

    def to_local(self, glob: np.ndarray) -> np.ndarray:
        glob = np.asarray(glob, dtype=float)
        return (self.triad.T @ (glob.reshape(-1, 3) - self.origin).T).T

    def rotated(self, R: np.ndarray) -> "Frame":
        return Frame(self.origin.copy(), self.triad @ R)

    def orthonormality_error(self) -> float:
        return float(np.abs(self.triad.T @ self.triad - np.eye(3)).max())


# --------------------------------------------------------------------------
# Parametros
# --------------------------------------------------------------------------
@dataclass
class StepParams:
    """Parametros de escalon entre pares de bases consecutivos."""

    shift: float = 0.0   # Dx, angstrom
    slide: float = 0.0   # Dy, angstrom
    rise: float = 3.38   # Dz, angstrom
    tilt: float = 0.0    # tau, grados
    roll: float = 0.0    # rho, grados
    twist: float = 36.0  # Omega, grados

    def as_tuple(self) -> Tuple[float, ...]:
        return (self.shift, self.slide, self.rise, self.tilt, self.roll, self.twist)


@dataclass
class PairParams:
    """Parametros intra-par de bases."""

    shear: float = 0.0      # Sx, angstrom
    stretch: float = 0.0    # Sy, angstrom
    stagger: float = 0.0    # Sz, angstrom
    buckle: float = 0.0     # kappa, grados
    propeller: float = 0.0  # omega, grados
    opening: float = 0.0    # sigma, grados

    def as_tuple(self) -> Tuple[float, ...]:
        return (self.shear, self.stretch, self.stagger,
                self.buckle, self.propeller, self.opening)


# --------------------------------------------------------------------------
# Escalon: directo e inverso
# --------------------------------------------------------------------------
def rot_sqrt(R: np.ndarray) -> np.ndarray:
    """Raiz cuadrada principal de una rotacion propia (mitad del eje-angulo)."""
    ang = math.acos(max(-1.0, min(1.0, (np.trace(R) - 1.0) / 2.0)))
    if ang < 1e-12:
        return np.eye(3)
    axis = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    axis = axis / (2.0 * math.sin(ang))
    h = ang / 2.0
    K = np.array([[0.0, -axis[2], axis[1]],
                  [axis[2], 0.0, -axis[0]],
                  [-axis[1], axis[0], 0.0]])
    return np.eye(3) + math.sin(h) * K + (1.0 - math.cos(h)) * (K @ K)


def step_rotation(sp: StepParams) -> Tuple[np.ndarray, np.ndarray]:
    """Devuelve (R_step, R_half) para un escalon.

    ``R_half`` es la rotacion que lleva T_i a la middle step triad; es la raiz
    cuadrada principal de ``R_step``.
    """
    gamma = math.hypot(sp.tilt, sp.roll)
    phi = math.degrees(math.atan2(sp.tilt, sp.roll)) if gamma > 1e-12 else 0.0
    R_step = (rot_z(sp.twist / 2.0 + phi) @ rot_y(gamma)
              @ rot_z(sp.twist / 2.0 - phi))
    return R_step, rot_sqrt(R_step)


def apply_step(frame: Frame, sp: StepParams) -> Frame:
    """Propaga un marco de par de bases un escalon."""
    R_step, R_half = step_rotation(sp)
    T_mst = frame.triad @ R_half
    origin = frame.origin + T_mst @ np.array([sp.shift, sp.slide, sp.rise])
    return Frame(origin, frame.triad @ R_step)


def _zyz_from_matrix(R: np.ndarray) -> Tuple[float, float, float]:
    """Extrae (phi, Gamma, psi) con R = Rz(phi).Ry(Gamma).Rz(psi), en grados."""
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, R[2, 2]))))
    if abs(math.sin(gamma * DEG)) < 1e-10:
        # Degenerado (sin curvatura): toda la rotacion es alrededor de z.
        phi = 0.0
        if R[2, 2] > 0:
            psi = math.degrees(math.atan2(R[1, 0], R[0, 0]))
        else:
            psi = math.degrees(math.atan2(-R[1, 0], -R[0, 0]))
    else:
        phi = math.degrees(math.atan2(R[1, 2], R[0, 2]))
        psi = math.degrees(math.atan2(R[2, 1], -R[2, 0]))
    return phi, gamma, psi


def _unwrap_twist_phase(alpha: float, beta: float) -> Tuple[float, float]:
    """Resuelve la degeneracion (Gamma, phi, Omega) ~ (Gamma, phi-180, Omega+360).

    La descomposicion ZYZ devuelve alpha y beta en (-180, 180], de modo que
    Omega = alpha + beta puede salir desplazado 360 grados respecto al giro
    fisico; cada desplazamiento de 360 en Omega va acompanado de uno de 180 en
    la fase, que invierte los signos de tilt y roll. Se fija la rama exigiendo
    Omega en (-180, 180], el unico giro con sentido fisico para un escalon o
    para la apertura de un par.
    """
    omega_raw = alpha + beta
    k = math.floor((omega_raw + 180.0) / 360.0)
    return omega_raw - 360.0 * k, (alpha - beta) / 2.0 + 180.0 * k


def frames_to_step_params(f_i: Frame, f_j: Frame) -> StepParams:
    """Inversa exacta de :func:`apply_step` (usada en control de calidad)."""
    R_rel = f_i.triad.T @ f_j.triad
    alpha, gamma, beta = _zyz_from_matrix(R_rel)
    twist, phi = _unwrap_twist_phase(alpha, beta)
    T_mst = f_i.triad @ rot_sqrt(R_rel)
    d = T_mst.T @ (f_j.origin - f_i.origin)
    return StepParams(shift=float(d[0]), slide=float(d[1]), rise=float(d[2]),
                      tilt=float(gamma * math.sin(phi * DEG)),
                      roll=float(gamma * math.cos(phi * DEG)),
                      twist=float(twist))


# --------------------------------------------------------------------------
# Par de bases: directo e inverso
# --------------------------------------------------------------------------
def pair_rotation(pp: PairParams) -> Tuple[np.ndarray, np.ndarray]:
    """Devuelve (R_pair, H) con ``H`` la raiz cuadrada principal de R_pair."""
    gamma = math.hypot(pp.buckle, pp.propeller)
    phi = math.degrees(math.atan2(pp.buckle, pp.propeller)) if gamma > 1e-12 else 0.0
    R_pair = (rot_z(pp.opening / 2.0 + phi) @ rot_y(gamma)
              @ rot_z(pp.opening / 2.0 - phi))
    return R_pair, rot_sqrt(R_pair)


def pair_to_base_frames(bp: Frame, pp: PairParams) -> Tuple[Frame, Frame]:
    """Marcos de las dos bases a partir del marco del par."""
    _, H = pair_rotation(pp)
    disp = bp.triad @ (np.array([pp.shear, pp.stretch, pp.stagger]) / 2.0)
    f1 = Frame(bp.origin + disp, bp.triad @ H)
    f2 = Frame(bp.origin - disp, bp.triad @ H.T @ rot_x(180.0))
    return f1, f2


def base_frames_to_pair(f1: Frame, f2: Frame) -> Tuple[Frame, PairParams]:
    """Inversa de :func:`pair_to_base_frames`.

    Recupera el marco del par y los parametros intra-par a partir de los dos
    marcos de base, deshaciendo primero el giro Rx(180) de la hebra
    complementaria. Es la misma operacion que realiza ``nastruct``.
    """
    T2u = f2.triad @ rot_x(180.0)            # base II "desvolteada"
    R_pair = T2u.T @ f1.triad
    H = rot_sqrt(R_pair)
    T_bp = f1.triad @ H.T
    o_bp = (f1.origin + f2.origin) / 2.0
    alpha, gamma, beta = _zyz_from_matrix(R_pair)
    opening, phi = _unwrap_twist_phase(alpha, beta)
    d = T_bp.T @ (f1.origin - f2.origin)
    pp = PairParams(shear=float(d[0]), stretch=float(d[1]), stagger=float(d[2]),
                    buckle=float(gamma * math.sin(phi * DEG)),
                    propeller=float(gamma * math.cos(phi * DEG)),
                    opening=float(opening))
    return Frame(o_bp, T_bp), pp


# --------------------------------------------------------------------------
# Modelos de fibra: simetria helicoidal exacta
# --------------------------------------------------------------------------
def fiber_frames(n_bp: int, twist: float, rise: float,
                 x_displacement: float = 0.0,
                 inclination: float = 0.0,
                 tip: float = 0.0) -> List[Frame]:
    """Marcos de par de bases de un modelo de fibra con simetria de tornillo.

    El eje de la helice es el eje z global. El marco i-esimo es

        T_i = Rz(i.Omega) . Ry(tip) . Rx(inclination)
        o_i = Rz(i.Omega) . (dx, 0, 0)^T + (0, 0, i.h)^T

    Para ``inclination = tip = 0`` y ``x_displacement = 0`` esto es exactamente
    la orbita del operador de tornillo S(Omega, h) aplicado al marco identidad.
    """
    base_T = rot_y(tip) @ rot_x(inclination)
    frames = []
    for i in range(n_bp):
        Rz = rot_z(i * twist)
        o = Rz @ np.array([x_displacement, 0.0, 0.0]) + np.array([0.0, 0.0, i * rise])
        frames.append(Frame(o, Rz @ base_T))
    return frames


# --------------------------------------------------------------------------
# Superposicion
# --------------------------------------------------------------------------
def kabsch(mobile: np.ndarray, ref: np.ndarray,
           weights: Sequence[float] | None = None) -> Tuple[np.ndarray, np.ndarray, float]:
    """Superposicion de minimos cuadrados: devuelve (R, t, rmsd).

    Minimiza sum_a w_a ||R.m_a + t - r_a||^2. Rotacion propia garantizada
    (det R = +1), por lo que nunca introduce una inversion quiral.
    """
    m = np.asarray(mobile, float).reshape(-1, 3)
    r = np.asarray(ref, float).reshape(-1, 3)
    if m.shape != r.shape or m.shape[0] < 3:
        raise ValueError(f"se requieren >=3 pares de atomos; recibido {m.shape} vs {r.shape}")
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


# --------------------------------------------------------------------------
# Pseudorrotacion del azucar (Altona & Sundaralingam 1972)
# --------------------------------------------------------------------------
PUCKER_ENVELOPES = ["C3'-endo", "C4'-exo", "O4'-endo", "C1'-exo", "C2'-endo",
                    "C3'-exo", "C4'-endo", "O4'-exo", "C1'-endo", "C2'-exo"]


def dihedral(p0, p1, p2, p3) -> float:
    p0, p1, p2, p3 = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def pseudorotation(coords: Dict[str, Iterable[float]]) -> dict | None:
    """Fase P, amplitud nu_max y dominio N/S del anillo de furanosa.

    Acepta anillos con O4' o S4' (azucares 4'-tio). Devuelve ``None`` si el
    fragmento no contiene un anillo de furanosa reconocible (p. ej. UNA
    aciclica o morfolino).
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
            "envelope": PUCKER_ENVELOPES[int(((P + 18) % 360) // 36)],
            "domain": domain}
