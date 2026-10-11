"""The derived fiber models reproduce the 3DNA repeats and are physically sane."""

import math
import os

import numpy as np
import pytest

from conftest import THREEDNA, needs_3dna
from dsxnabuilder.geometry import rot_x, rot_z
from dsxnabuilder.offlib import signed_volume
from dsxnabuilder.structure import parse_pdb, residues

EXPECTED_OK = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21,
               22, 23, 24, 25, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 49, 50, 51,
               52, 54, 55}
LEFT_HANDED = {15, 16}
GENERIC = {1, 4, 7, 46, 47, 54, 55}


def test_expected_models_buildable(lib):
    assert set(lib.models) == EXPECTED_OK


def test_excluded_models_have_reasons(lib):
    for mid, m in lib.excluded.items():
        assert m["issues"], f"model {mid} excluded without a reason"


def test_helical_twist_matches_3dna_table(raw_models):
    """Per-bp helical twist equals the README twist of the unit / k (mod 360)."""
    for key, m in raw_models["models"].items():
        if m["status"] != "ok":
            continue
        per_unit = m["bp_twist"] * m["k"]
        diff = ((per_unit - m["twist"] + 180) % 360) - 180
        assert abs(diff) < 6e-4 * m["k"], key   # bp_twist is stored to 3 decimals
        assert abs(abs(m["twist"]) - abs(m["readme_twist"])) < 1e-9 or \
            abs((abs(m["twist"]) + abs(m["readme_twist"])) % 360) < 1e-9, key
        assert abs(m["rise"] - abs(m["readme_rise"])) < 1e-9, key
        assert m["rise"] > 0


def test_handedness(lib):
    for mid, m in lib.models.items():
        assert (m.handedness == "left") == (mid in LEFT_HANDED), mid


def test_backbone_closes_everywhere(lib, catalog, aliases):
    """O3'(i)-P(i+1) in 1.55-1.65 A along both strands of a 3-repeat duplex."""
    from dsxnabuilder.structure import o3p_distances
    from dsxnabuilder.build import build_duplex
    for mid, m in lib.models.items():
        chem1, chem2 = m.strand_chem("I"), m.strand_chem("II")
        canon = {"I": "G", "X": "G"}
        top = "".join(canon.get(s.I.native, s.I.native) for s in m.slots) * 3
        top = top.replace("U", "T") if chem1 == "DNA" else top.replace("T", "U")
        comp = {"A": "T" if chem2 == "DNA" else "U", "T": "A", "U": "A", "G": "C", "C": "G"}
        bottom = "".join(comp[c] for c in top)
        res = build_duplex(top, bottom, mid, lib, catalog, aliases)
        d = [x for *_, x in res.o3p]
        assert len(d) == 2 * (len(top) - 1)
        assert min(d) >= 1.55 and max(d) <= 1.65, (mid, min(d), max(d))


def test_op1_op2_follow_amber_chirality(lib):
    for mid, m in lib.models.items():
        for s in m.slots:
            for ss in (s.I, s.II):
                for t in ss.templates.values():
                    c = dict(zip(t.names, t.xyz))
                    if all(k in c for k in ("P", "OP1", "OP2", "O5'")):
                        assert signed_volume(c, "P", "OP1", "OP2", "O5'") < 0, mid


def test_watson_crick_c1_distance(lib):
    for mid, m in lib.models.items():
        for j, s in enumerate(m.slots):
            if s.pairing != "wc":
                continue
            _, f1 = m.slot_strand(j, "I")
            _, f2 = m.slot_strand(j, "II")
            a = next(iter(s.I.templates.values()))
            b = next(iter(s.II.templates.values()))
            c1a = f1.to_global(a.xyz[a.names.index("C1'")])[0]
            c1b = f2.to_global(b.xyz[b.names.index("C1'")])[0]
            assert 10.0 < np.linalg.norm(c1a - c1b) < 11.2, (mid, j)


def test_generic_models_dyad_symmetry(raw_models):
    """Single-bp generic repeats have a dyad: local tilt and shift ~ 0."""
    for mid in GENERIC:
        steps = raw_models["models"][str(mid)]["step_params"]
        shift, tilt = steps[0][0], steps[0][3]
        # 3DNA coordinates have 3 decimals; exact symmetry is tested in
        # test_geometry.py::test_fiber_dyad_symmetry on ideal frames
        assert abs(shift) < 0.02 and abs(tilt) < 0.1, (mid, shift, tilt)


@needs_3dna
@pytest.mark.parametrize("mid", sorted(GENERIC))
def test_reproduces_3dna_generic_coordinates(lib, mid):
    """Strand I residue of a generic model, rebuilt for every letter, equals the
    3DNA per-base file after the same normalization (<1e-3 A)."""
    m = lib.models[mid]
    flip = m.meta.get("flipped_to_plus_z")
    raw = __import__("json").load(open(os.path.join(
        os.path.dirname(__file__), "..", "data", "fiber_models.json")))["models"][str(mid)]
    flip = raw["flipped_to_plus_z"]
    for letter in "ACGT":
        path = os.path.join(THREEDNA, "fiber", f"str{mid:02d}", f"{letter}.pdb")
        with open(path) as fh:
            atoms = parse_pdb(fh.read())
        ref = np.array([a.xyz for a in atoms])
        if flip:
            ref = (rot_x(180.0) @ ref.T).T
        names = [a.name.replace("O1P", "OP1").replace("O2P", "OP2").replace("C5M", "C7")
                 for a in atoms]
        t = m.slots[0].I.templates[letter]
        _, frame = m.slot_strand(0, "I")
        got = dict(zip(t.names, frame.to_global(t.xyz)))
        for n, p in zip(names, ref):
            if n in ("OP1", "OP2"):
                # may be relabelled by chirality; compare the pair as a set
                q = min(np.linalg.norm(got["OP1"] - p), np.linalg.norm(got["OP2"] - p))
            else:
                q = np.linalg.norm(got[n] - p)
            assert q < 1e-3, (mid, letter, n, q)
