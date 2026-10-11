#!/usr/bin/env python3
"""
Derive ``data/fiber_models.json`` from the 3DNA fiber repeating units.

    python tools/derive_fiber_models.py [--threedna external/3dna] [--out data/fiber_models.json]

Inputs are the text files fetched by ``tools/fetch_3dna.py``. For every
double-stranded model the script

1. reads the repeating unit of strand I (``S1.pdb`` or the per-base files
   ``A/C/G/T.pdb`` of 3DNA's generic models) and of strand II (``S2.pdb``, or
   the dyad image (x,y,z)->(x,-y,-z) of the complementary per-base file);
2. finds the screw S(Omega, h) that maps the unit onto the next one along
   strand I 5'->3'. The README table gives |Omega| and |h| but neither the
   direction of the strands nor, for long units, the handedness, so the four
   sign combinations are tried and the one that closes the backbone
   (O3'(last)-P(first of next unit) = 1.55-1.65 A) is kept; it must be unique;
3. normalizes the model so that strand I runs 5'->3' along +z (Rx(180) when
   3DNA wrote it the other way round);
4. relabels OP1/OP2 by chirality to the Amber convention (the 3DNA files are
   not consistent: some models, and even residues within one model, have
   them swapped), counting the swaps;
5. pairs strand I residues with strand II residues (or their screw images) by
   base-base hydrogen-bond contacts. Neither the closest C1'-C1' distance nor
   the S1->S2 dyad identifies the partner in Z-DNA or D-DNA;
6. computes base standard frames (ring fit onto Olson 2001 ``Atomic_X``),
   nucleotide templates in those frames, chi/syn, sugar pucker, pair and step
   parameters, and runs a QC gate (bond lengths, backbone closure inside and
   across repeats, clashes, ring fit, pairing). Models that fail are kept in
   the JSON with ``status: excluded`` and the reasons, so the UI can show them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)

from dsxnabuilder.geometry import (Frame, base_frames_to_pair,  # noqa: E402
                                   frames_to_step_params, pseudorotation,
                                   rot_x, rot_z)
from dsxnabuilder.nucleotide import (BACKBONE_ATOMS, GLYCOSIDIC_BOND,  # noqa: E402
                                     IDEAL_BONDS, WATSON_CRICK, chi, is_syn,
                                     normalize_name)
from dsxnabuilder.offlib import signed_volume  # noqa: E402
from dsxnabuilder.structure import (HBOND_ATOMS, parse_pdb, residues,  # noqa: E402
                                    standard_base_frame)

SCHEMA = "dsxnabuilder/fiber_models/1"
O3P_OK = (1.55, 1.65)
BOND_TOL = 0.08
CLASH = 2.5
HB_MAX = 3.5
# Amber DNA.OL15 / RNA.OL3 libraries: (OP1-P)x(OP2-P).(O5'-P) = -2.435 A^3.
AMBER_OP_SIGN = -1.0

REFERENCES = {
    "Lu2003": "Lu XJ, Olson WK (2003) 3DNA: a software package for the analysis, "
              "rebuilding and visualization of three-dimensional nucleic acid "
              "structures. Nucleic Acids Res 31:5108-5121. doi:10.1093/nar/gkg680",
    "Lu2008": "Lu XJ, Olson WK (2008) 3DNA: a versatile, integrated software system "
              "for the analysis, rebuilding and visualization of three-dimensional "
              "nucleic-acid structures. Nat Protoc 3:1213-1227. doi:10.1038/nprot.2008.104",
    "Olson2001": "Olson WK et al. (2001) A standard reference frame for the "
                 "description of nucleic acid base-pair geometry. J Mol Biol "
                 "313:229-237. doi:10.1006/jmbi.2001.4987",
    "Arnott1999": "Arnott S (1999) Polynucleotide secondary structures: an "
                  "historical perspective. In: Neidle S (ed) Oxford Handbook of "
                  "Nucleic Acid Structure, pp 1-38. Oxford University Press",
    "Chandrasekaran1990": "Chandrasekaran R, Arnott S (1990) The structures of DNA "
                          "and RNA helices in oriented fibers. In: Saenger W (ed) "
                          "Landolt-Bornstein, Group VII, vol 1b, pp 31-170. Springer",
    "Alexeev1987": "Alexeev DG, Lipanov AA, Skuratovskii IYa (1987) The structure "
                   "of poly(dA).poly(dT) as revealed by an X-ray fibre diffraction. "
                   "J Biomol Struct Dyn 4:989-1011",
    "vanDam2000": "van Dam L, Levitt MH (2000) BII nucleotides in the B and C forms "
                  "of natural-sequence polymeric DNA: a new model for the C form of "
                  "DNA. J Mol Biol 304:541-561",
    "Premilat2001": "Premilat S, Albiser G (2001) A new D-DNA form of "
                    "poly(dA-dT).poly(dA-dT): an A-DNA type structure with reversed "
                    "Hoogsteen pairing. Eur Biophys J 30:404-410",
}

# Family labels (from the 3DNA fiber/README table). 30-36 and 42 are triplexes,
# a quadruplex and a single strand: not duplexes, never built.
FAMILY = {}
for _ids, _fam in (((1, 2, 3, 54), "A-DNA"), ((4, 5, 6, 46, 55), "B-DNA"),
                   ((7, 8, 9, 10, 11, 47, 53), "C-DNA"), ((12, 13, 14), "D-DNA"),
                   ((48,), "D(A)-DNA"), ((52,), "D(B)-DNA"), ((15, 16), "Z-DNA"),
                   ((17,), "L-DNA"), ((18, 19, 37, 38, 39, 40, 41, 43), "B'-DNA"),
                   ((44, 45), "B-DNA (poly dA:dT salts)"), ((49, 50), "S-DNA"),
                   ((51,), "B*-DNA"), ((20,), "A-RNA"), ((21,), "A'-RNA"),
                   ((22, 23, 24, 25), "DNA:RNA hybrid"), ((26, 27), "poly(X)"),
                   ((28, 29), "poly(s2U)")):
    for _i in _ids:
        FAMILY[_i] = _fam
NOT_DUPLEX = {30, 31, 32, 33, 34, 35, 36, 42}
GENERIC = {1, 4, 7, 46, 47, 53, 54, 55}      # README: "generic sequence"
DEPRECATED = {53}                            # README: "(depreciated)"
FIXED_STRAND2 = {44, 45}                     # README note 4


def references_for(mid: int) -> List[str]:
    refs = ["Lu2003", "Lu2008"]
    if mid <= 41:
        refs.append("Arnott1999")
    elif mid <= 43:
        refs.append("Chandrasekaran1990")
    elif mid <= 45:
        refs.append("Alexeev1987")
    elif mid <= 47:
        refs.append("vanDam2000")
    else:
        refs.append("Premilat2001")
    return refs


# --------------------------------------------------------------------------
@dataclass
class Res:
    resname: str
    names: List[str]
    xyz: np.ndarray

    def c(self) -> Dict[str, np.ndarray]:
        return {n: self.xyz[i] for i, n in enumerate(self.names)}

    def moved(self, R: np.ndarray, t=None) -> "Res":
        t = np.zeros(3) if t is None else t
        return Res(self.resname, list(self.names), (R @ self.xyz.T).T + t)

    @property
    def letter(self) -> str:
        return self.resname.strip().upper()[:1]

    @property
    def chem(self) -> str:
        return "RNA" if "O2'" in self.names else "DNA"


DYAD = rot_x(180.0)
RING_LETTER = {"A": "A", "G": "G", "I": "G", "X": "G", "C": "C", "T": "T", "U": "U"}


def read_res(path: str) -> List[Res]:
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    out = []
    for _, ats in residues(parse_pdb(text)).items():
        names = [normalize_name(a.name) for a in ats]
        out.append(Res(ats[0].resname, names, np.array([a.xyz for a in ats])))
    return out


def parse_readme(path: str) -> Dict[int, Tuple[float, float, str]]:
    table: Dict[int, Tuple[float, float, str]] = {}
    last = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    inside = False
    for line in lines:
        if line.startswith("-----"):
            inside = not inside
            continue
        if not inside:
            continue
        m = re.match(r"^\s*(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(.*\S)\s*$", line)
        if m:
            last = int(m.group(1))
            table[last] = (float(m.group(2)), float(m.group(3)), m.group(4))
        elif last is not None and line.strip():
            tw, h, d = table[last]
            table[last] = (tw, h, d + " " + line.strip())
    return table


def screw_apply(res: Res, tw: float, h: float, n: int = 1) -> Res:
    return res.moved(rot_z(n * tw), np.array([0.0, 0.0, n * h]))


def link(a: Res, b: Res) -> float:
    ca, cb = a.c(), b.c()
    if "O3'" not in ca or "P" not in cb:
        return float("nan")
    return float(np.linalg.norm(ca["O3'"] - cb["P"]))


def relabel_op(res: Res) -> int:
    c = res.c()
    if not all(k in c for k in ("P", "OP1", "OP2", "O5'")):
        return 0
    if signed_volume(c, "P", "OP1", "OP2", "O5'") * AMBER_OP_SIGN > 0:
        return 0
    i, j = res.names.index("OP1"), res.names.index("OP2")
    res.names[i], res.names[j] = "OP2", "OP1"
    return 1


def base_frame(res: Res, std: Dict[str, Dict[str, np.ndarray]]) -> Tuple[Frame, float]:
    ref = std[RING_LETTER.get(res.letter, "A")]
    return standard_base_frame(res.c(), ref)


def hbond_contacts(a: Res, b: Res) -> List[Tuple[str, str, float]]:
    ca, cb = a.c(), b.c()
    out = []
    for x in HBOND_ATOMS:
        if x not in ca:
            continue
        for y in HBOND_ATOMS:
            if y in cb:
                d = float(np.linalg.norm(ca[x] - cb[y]))
                if d <= HB_MAX:
                    out.append((x, y, round(d, 3)))
    return out


def pair_type(a: Res, b: Res, contacts) -> str:
    pur, pyr = (a, b) if a.letter in "AGIX" and b.letter not in "AGIX" else \
               (b, a) if b.letter in "AGIX" and a.letter not in "AGIX" else (None, None)
    if pur is None:
        return "other"
    pur_atoms = {x if pur is a else y for x, y, _ in contacts}
    if "N1" in pur_atoms:
        return "wc"
    if "N7" in pur_atoms:
        return "hoogsteen"
    return "other"


def nonstandard(res: Res, std: Dict[str, Dict[str, np.ndarray]]) -> bool:
    if res.letter not in ("A", "C", "G", "T", "U"):
        return True
    expected = set(BACKBONE_ATOMS) | (set(std[res.letter]) - {"C1'"})
    return any(n not in expected for n in res.names)


def missing_atoms(res: Res, std) -> List[str]:
    need = [n for n in BACKBONE_ATOMS if n != "O2'"]
    if res.letter in std and not nonstandard(res, std):
        need += [n for n in std[res.letter] if n != "C1'"]
    elif res.letter in RING_LETTER:
        need += ["N9" if RING_LETTER[res.letter] in "AG" else "N1"]
    return [n for n in need if n not in res.names]


def bond_problems(res: Res) -> List[str]:
    c = res.c()
    out = []
    pairs = dict(IDEAL_BONDS)
    gly = "N9" if "N9" in c else "N1"
    pairs[("C1'", gly)] = GLYCOSIDIC_BOND
    for (a, b), d0 in pairs.items():
        if a in c and b in c:
            d = float(np.linalg.norm(c[a] - c[b]))
            if abs(d - d0) > BOND_TOL:
                out.append(f"{res.resname}:{a}-{b}={d:.2f}")
    return out


def frame_json(f: Frame) -> dict:
    return {"o": [round(float(x), 6) for x in f.origin],
            "R": [[round(float(x), 8) for x in row] for row in f.triad]}


def tmpl_json(f: Frame, res: Res) -> dict:
    local = f.to_local(res.xyz)
    return {"names": list(res.names),
            "xyz": [[round(float(v), 5) for v in p] for p in local]}


def sugar_info(res: Res) -> Optional[dict]:
    p = pseudorotation(res.c())
    if p is None:
        return None
    return {"P": p["P"], "domain": p["domain"], "envelope": p["envelope"]}


# --------------------------------------------------------------------------
def derive_model(mid: int, mdir: str, readme, std, sha: Dict[str, str]) -> dict:
    tw0, h0, desc = readme[mid]
    out = {"id": mid, "name": f"{FAMILY.get(mid, 'fiber')} #{mid}",
           "family": FAMILY.get(mid, "other"), "description": desc,
           "readme_twist": tw0, "readme_rise": h0,
           "references": references_for(mid), "issues": [], "status": "ok"}
    files = sorted(os.listdir(mdir))
    out["source_files"] = {f: sha.get(f"fiber/str{mid:02d}/{f}", "") for f in files}

    if mid in NOT_DUPLEX:
        out.update(status="excluded", issues=["not a double helix (triplex, "
                                              "quadruplex or single strand)"])
        return out

    per_base = {L: read_res(os.path.join(mdir, f"{L}.pdb"))[0]
                for L in "ACGTU" if f"{L}.pdb" in files}
    if "S1.pdb" in files:
        kind = "repeat"
        U1 = read_res(os.path.join(mdir, "S1.pdb"))
        U2 = read_res(os.path.join(mdir, "S2.pdb"))
        per_base = {}
    elif mid in FIXED_STRAND2:
        # README note 4: these per-base files are strand I (A) and strand II
        # (T) as they sit in the duplex, not a dyad-related generic pair.
        kind = "repeat"
        U1, U2 = [per_base["A"]], [per_base["T"]]
        per_base = {}
    elif per_base:
        kind = "per-base"
        ref = "A" if "A" in per_base else sorted(per_base)[0]
        comp = WATSON_CRICK[ref]
        if comp not in per_base:
            out.update(status="excluded", issues=[f"no {comp} file to pair with {ref}"])
            return out
        U1 = [per_base[ref]]
        U2 = [per_base[comp].moved(DYAD)]
    else:
        out.update(status="excluded", issues=["unrecognized file layout"])
        return out
    out["kind"] = kind
    k = len(U1)
    out["k"] = k
    issues = out["issues"]

    # -- 2. screw sign from backbone closure --------------------------------
    cands = []
    for st in (1, -1):
        for sh in (1, -1):
            tw, h = st * tw0, sh * h0
            if any(abs(((tw - c[0]) + 180) % 360 - 180) < 1e-9 and abs(h - c[1]) < 1e-9
                   for c in cands):
                continue
            d = link(U1[-1], screw_apply(U1[0], tw, h))
            cands.append((tw, h, d))
    good = [c for c in cands if O3P_OK[0] - 0.05 <= c[2] <= O3P_OK[1] + 0.05]
    out["screw_candidates"] = [[c[0], c[1], round(c[2], 3)] for c in cands]
    if len(good) != 1:
        best = min(cands, key=lambda c: abs(c[2] - 1.6))
        out.update(status="excluded")
        issues.append(f"backbone does not close uniquely across repeats "
                      f"({len(good)} screw sign(s) fit; best O3'-P {best[2]:.2f} A)")
        return out
    tw, h, _ = good[0]

    # strand II order: its 5'->3' direction maps a unit by S^-1
    def closes_II(u):
        inter = link(u[-1], screw_apply(u[0], tw, h, -1))
        intra = [link(u[j], u[j + 1]) for j in range(len(u) - 1)]
        return inter, intra
    inter2, intra2 = closes_II(U2)
    if not (O3P_OK[0] - 0.05 <= inter2 <= O3P_OK[1] + 0.05):
        r2 = list(reversed(U2))
        i2, ii2 = closes_II(r2)
        if O3P_OK[0] - 0.05 <= i2 <= O3P_OK[1] + 0.05:
            U2, inter2, intra2 = r2, i2, ii2
    # -- 3. normalize to strand I 5'->3' along +z ----------------------------
    flipped = h < 0
    if flipped:
        U1 = [r.moved(DYAD) for r in U1]
        U2 = [r.moved(DYAD) for r in U2]
        per_base = {L: r.moved(DYAD) for L, r in per_base.items()}
        tw, h = -tw, -h
    out["flipped_to_plus_z"] = flipped
    out["twist"], out["rise"] = tw, h

    # -- 4. OP1/OP2 by chirality ---------------------------------------------
    swaps = sum(relabel_op(r) for r in U1 + U2 + list(per_base.values()))
    out["op_relabelled"] = swaps

    # -- 5. pairing ------------------------------------------------------------
    pairs = []
    used = set()
    for j, a in enumerate(U1):
        fa, _ = base_frame(a, std)
        best = None
        for p, b0 in enumerate(U2):
            for m in range(-2, 3):
                if (p, m) in used:
                    continue
                b = screw_apply(b0, tw, h, m)
                fb, _ = base_frame(b, std)
                dz = abs(fa.origin[2] - fb.origin[2])
                if dz > 2.5 or np.linalg.norm(fa.origin - fb.origin) > 8.0:
                    continue
                cts = hbond_contacts(a, b)
                if not cts:
                    continue
                score = (len(cts), -sum(d for *_, d in cts))
                if best is None or score > best[0]:
                    best = (score, p, m, b, cts)
        if best is None:
            out.update(status="excluded")
            issues.append(f"no hydrogen-bonded partner for strand I residue {j + 1} "
                          f"({a.resname})")
            return out
        used.add((best[1], best[2]))
        pairs.append((a, best[3], best[4], best[1], best[2]))

    # -- 6. frames, templates, parameters, QC ---------------------------------
    slots, bp_frames = [], []
    seq_policy = "any" if (mid in GENERIC or set("ACGT") <= set(per_base)) else "swap"
    for j, (a, b, cts, p, m) in enumerate(pairs):
        fa, ra = base_frame(a, std)
        fb, rb = base_frame(b, std)
        ptype = pair_type(a, b, cts)
        # 3DNA orientation: pair z along strand I 5'->3' (+z after normalization)
        bpf, pp = base_frames_to_pair(fa, fb, direction=(0.0, 0.0, 1.0))
        bp_frames.append(bpf)
        slot = {"pairing": ptype,
                "contacts": [list(c) for c in cts],
                "pair_params": [round(x, 3) for x in pp.as_tuple()],
                "I": strand_json(a, fa, ra, std, per_base, None),
                "II": strand_json(b, fb, rb, std, per_base, (tw, h, m))}
        if ptype != "wc":
            seq_policy = "pattern"
        if slot["I"]["syn"] or slot["II"]["syn"]:
            seq_policy = "pattern"
        slots.append(slot)
    if k > 1 and seq_policy != "any":
        seq_policy = "pattern"
    homo = all(s["pairing"] == "other" for s in slots)
    if homo:
        out.update(status="excluded")
        issues.append("non-canonical homo base pairs (e.g. X:X, s2U:s2U); "
                      "not buildable with A/C/G/T/U")
    out["slots"] = slots
    out["sequence_policy"] = seq_policy

    steps = []
    for j in range(k):
        f1 = bp_frames[j]
        if j + 1 < k:
            f2 = bp_frames[j + 1]
        else:
            R = rot_z(tw)
            f2 = Frame(R @ bp_frames[0].origin + np.array([0, 0, h]),
                       R @ bp_frames[0].triad)
        sp = frames_to_step_params(f1, f2)
        steps.append([round(x, 3) for x in sp.as_tuple()])
    out["step_params"] = steps
    out["local_step_twist_mean"] = round(float(np.mean([s[5] for s in steps])), 3)
    hel = helical_twist(U1, tw, h)
    out["bp_twist"] = round(hel, 3)
    out["bp_rise"] = round(h / k, 4)
    out["handedness"] = "right" if hel > 0 else "left"

    qc = model_qc(U1, U2, pairs, tw, h, std)
    out["qc"] = qc
    issues.extend(qc.pop("failures"))
    if issues and out["status"] == "ok":
        out["status"] = "excluded"
    if mid in DEPRECATED and out["status"] == "ok":
        out["status"] = "excluded"
        issues.append("marked deprecated in 3DNA (use #7 or #47)")
    readme_bp = abs(tw0) / k if k else 0
    out["meta"] = {"readme_twist_per_bp": round(readme_bp, 3)}
    return out


def helical_twist(U1: List[Res], tw: float, h: float) -> float:
    """Mean helical twist per residue of strand I (rise > 0 after normalization).

    Measured from the azimuth of C1' around the helix axis between consecutive
    residues, wrapped to (-180, 180]. Unlike the local step twist between
    base-pair frames this is the helical twist of the fiber literature (e.g.
    32.7 for A-DNA, where inclined pairs give a smaller local twist) and its
    sign is the handedness: Z-DNA pairs are flipped relative to B-DNA, so the
    local step parameters of a left-handed helix can come out positive.
    """
    seq = U1 + [screw_apply(U1[0], tw, h)]
    az = [math.degrees(math.atan2(r.c()["C1'"][1], r.c()["C1'"][0])) for r in seq]
    d = [((az[i + 1] - az[i] + 180.0) % 360.0) - 180.0 for i in range(len(seq) - 1)]
    return float(np.mean(d))


def strand_json(res: Res, f: Frame, rmsd: float, std, per_base, screw) -> dict:
    ns = nonstandard(res, std)
    d = {"frame": frame_json(f), "native": res.letter, "resname": res.resname,
         "chem": res.chem, "nonstandard": ns, "ring_rmsd": round(rmsd, 4)}
    ch = chi(res.c())
    d["chi"] = round(ch, 1) if ch is not None else None
    d["syn"] = is_syn(ch)
    d["pucker"] = sugar_info(res)
    templates = {}
    if per_base:
        for L, r in per_base.items():
            rr = r
            if screw is not None:
                # strand II per-base residue: dyad image (already flipped with
                # the rest), then the same screw power as the paired residue
                tw, h, m = screw
                rr = screw_apply(_strand2_from_file(r), tw, h, m)
            templates[L] = tmpl_json(f, rr)
    elif not ns:
        templates[res.letter] = tmpl_json(f, res)
    else:
        templates["_native"] = tmpl_json(f, res)
    d["templates"] = templates
    return d


def _strand2_from_file(r: Res) -> Res:
    # per_base residues are stored in the (possibly flipped) strand I frame;
    # the strand II copy is the dyad image, which commutes with the flip.
    return r.moved(DYAD)


def model_qc(U1, U2, pairs, tw, h, std) -> dict:
    fails: List[str] = []
    qc: Dict[str, object] = {}
    # backbone closure, strand I (5'->3' = +S) and strand II (5'->3' = -S)
    o3p = []
    for u, s in ((U1, 1), (U2, -1)):
        for j in range(len(u) - 1):
            o3p.append(link(u[j], u[j + 1]))
        o3p.append(link(u[-1], screw_apply(u[0], tw, h, s)))
    qc["o3p"] = [round(x, 3) for x in o3p]
    bad = [x for x in o3p if not (O3P_OK[0] <= x <= O3P_OK[1])]
    if bad:
        fails.append(f"O3'-P outside {O3P_OK[0]}-{O3P_OK[1]} A: "
                     + ", ".join(f"{x:.2f}" for x in bad))
    # bonds and atoms
    bonds, miss = [], []
    for r in U1 + U2:
        bonds += bond_problems(r)
        m = missing_atoms(r, std)
        if m:
            miss.append(f"{r.resname}: {','.join(m)}")
    if bonds:
        fails.append("bond lengths off by >0.08 A: " + "; ".join(bonds[:6]))
    if miss:
        fails.append("missing atoms: " + "; ".join(miss))
    # ring fit
    rms = [base_frame(r, std)[1] for r in U1 + U2]
    qc["max_ring_rmsd"] = round(max(rms), 4)
    # clashes on three consecutive repeats of the pairs
    atoms = []
    for n in (-1, 0, 1):
        for j, (a, b, *_rest) in enumerate(pairs):
            for strand, r in (("I", a), ("II", b)):
                rr = screw_apply(r, tw, h, n)
                for nm, p in zip(rr.names, rr.xyz):
                    atoms.append((strand, n, j, nm, p))
    k = len(pairs)
    clashes = []
    xyz = np.array([x[4] for x in atoms])
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if d[i, j] >= CLASH:
                continue
            a, b = atoms[i], atoms[j]
            if (a[0], a[1], a[2]) == (b[0], b[1], b[2]):
                continue
            if _bonded_across(a, b, k):
                continue
            clashes.append(f"{a[0]}{a[1] * k + a[2]}:{a[3]}-{b[0]}{b[1] * k + b[2]}:{b[3]}="
                           f"{d[i, j]:.2f}")
    qc["clashes"] = clashes[:10]
    if clashes:
        fails.append(f"{len(clashes)} heavy-atom contacts < {CLASH} A")
    qc["failures"] = fails
    return qc


_LINK_O3 = {"O3'", "C3'"}
_LINK_P = {"P", "OP1", "OP2", "O5'"}


def _bonded_across(a, b, k) -> bool:
    if a[0] != b[0]:
        return False
    ia, ib = a[1] * k + a[2], b[1] * k + b[2]
    if abs(ia - ib) != 1:
        return False
    # strand I: residue i -> i+1 is 5'->3'; strand II: i+1 -> i
    if a[0] == "I":
        first, second = (a, b) if ia < ib else (b, a)
    else:
        first, second = (a, b) if ia > ib else (b, a)
    return first[3] in _LINK_O3 and second[3] in _LINK_P


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="derive data/fiber_models.json")
    ap.add_argument("--threedna", default=os.path.join(ROOT, "external", "3dna"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "fiber_models.json"))
    a = ap.parse_args(argv)

    with open(os.path.join(a.threedna, "PROVENANCE.json"), encoding="utf-8") as fh:
        prov = json.load(fh)
    sha = prov["files"]
    readme = parse_readme(os.path.join(a.threedna, "fiber", "README"))
    std_res = {}
    for L in "ACGTU":
        with open(os.path.join(a.threedna, "config", f"Atomic_{L}.pdb"), encoding="utf-8") as fh:
            ats = parse_pdb(fh.read())
        std_res[L] = {normalize_name(x.name): x.xyz for x in ats}

    models = {}
    for mid in sorted(readme):
        mdir = os.path.join(a.threedna, "fiber", f"str{mid:02d}")
        if not os.path.isdir(mdir):
            continue
        try:
            models[str(mid)] = derive_model(mid, mdir, readme, std_res, sha)
        except Exception as exc:  # keep going; report the model as excluded
            models[str(mid)] = {"id": mid, "name": f"#{mid}", "family": FAMILY.get(mid, "other"),
                                "description": readme[mid][2], "status": "excluded",
                                "issues": [f"derivation error: {exc!r}"]}

    payload = {
        "schema": SCHEMA,
        "provenance": {"source": prov.get("source"), "commit": prov.get("commit"),
                       "note": prov.get("note"),
                       "readme_sha256": sha.get("fiber/README"),
                       "generator": "tools/derive_fiber_models.py"},
        "references": REFERENCES,
        "standard_bases": {L: {"names": list(v), "xyz": [[round(float(x), 4) for x in p]
                                                          for p in v.values()]}
                           for L, v in std_res.items()},
        "models": models,
    }
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    ok = [m for m in models.values() if m["status"] == "ok"]
    print(f"{len(models)} models, {len(ok)} ok -> {a.out} "
          f"({os.path.getsize(a.out) / 1e6:.2f} MB)")
    for m in models.values():
        flag = "ok " if m["status"] == "ok" else "-- "
        extra = "" if m["status"] == "ok" else "  " + "; ".join(m.get("issues", []))[:150]
        print(f"  {flag}{m['id']:>2} {m.get('family', ''):<22} k={m.get('k', '-')!s:<2} "
              f"tw/bp={m.get('bp_twist', '-')!s:<8} {m.get('kind', ''):<8} "
              f"{m.get('sequence_policy', '')}{extra}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
