"""modXNA residue library and the local transplant."""

import io
import json
import os
import zipfile

import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
pytestmark = pytest.mark.skipif(
    not os.path.isdir(os.path.join(ROOT, "web", "lib", "residues")),
    reason="residue library not built (tools/build_residue_library.py)")

CAP_CHARGE = {"int": -1.0, "5cap": -0.320348, "3cap": -0.679652}
with open(os.path.join(ROOT, "data", "amber_residues.json")) as fh:
    AMBER = json.load(fh)["units"]


@pytest.fixture(scope="module")
def rl():
    from dsxnabuilder.reslib import ResidueLibrary
    return ResidueLibrary.load_default()


@pytest.fixture(scope="module")
def msession(rl):
    from dsxnabuilder.service import Session
    return Session(residue_library=rl)


def test_index_complete_and_charged(rl):
    assert len(rl.residues) == 244
    for key, r in rl.residues.items():
        assert abs(r["charge"] - CAP_CHARGE[r["cap"]]) < 1e-4, key
        assert key == "-".join(r["triplet"] + [r["cap"]])
        assert (r["head"] is None) == (r["cap"] == "5cap")
        assert (r["tail"] is None) == (r["cap"] == "3cap")


def scaffold(lib, model=20, letter="A", chem="RNA"):
    m = lib.models[model]
    names, xyz, _ = m.nucleotide(3, "I", letter, chem, lib.standard)
    return names, xyz


def place(rl, lib, key, letter, chem="RNA", model=20):
    from dsxnabuilder.reslib import transplant
    names, xyz = scaffold(lib, model, letter, chem)
    out_names, out_xyz, rep = transplant(rl.unit(key), names, xyz)
    return dict(zip(names, xyz)), dict(zip(out_names, out_xyz)), rep


@pytest.mark.parametrize("key,letter", [("PS1-OME-RAA-int", "A"), ("PS2-AF2-RCC-int", "C"),
                                        ("RPO-MOE-DTT-int", "T"), ("PS1-RC3-M5C-int", "C")])
def test_shared_atoms_do_not_move(rl, lib, key, letter):
    sc, out, rep = place(rl, lib, key, letter)
    for n, p in out.items():
        if n in sc and n not in rep["placed"]:
            assert np.allclose(p, sc[n]), (key, n)
    assert not rep["strained"]


def test_phosphorothioate_stereo_rule(rl, lib):
    from dsxnabuilder.offlib import signed_volume
    sc, ps1, _ = place(rl, lib, "PS1-RC3-RAA-int", "A")
    _, ps2, _ = place(rl, lib, "PS2-RC3-RAA-int", "A")
    for out, s_name, o_name in ((ps1, "S1", "OP1"), (ps2, "S2", "OP2")):
        v = out[s_name] - out["P"]
        u = sc[o_name] - sc["P"]
        cos = np.dot(v, u) / np.linalg.norm(v) / np.linalg.norm(u)
        assert cos > 0.9999, s_name
        assert 1.85 < np.linalg.norm(v) < 2.10
    # (S, retained O, O5') has opposite handedness: two diastereomers
    v1 = signed_volume(ps1, "P", "S1", "OP2", "O5'")
    v2 = signed_volume(ps2, "P", "S2", "OP1", "O5'")
    assert v1 < 0 < v2


def test_substituent_geometry(rl, lib):
    from dsxnabuilder.nucleotide import place_o2prime
    sc, f, _ = place(rl, lib, "RPO-AF2-RCC-int", "C")
    assert 1.36 < np.linalg.norm(f["F2'"] - f["C2'"]) < 1.42
    # 2'-F-ribo: fluorine on the alpha face, where ribose O2' sits
    ref = place_o2prime({k: v for k, v in sc.items() if k != "O2'"})
    assert np.linalg.norm(f["F2'"] - ref) < 0.3
    _, m, _ = place(rl, lib, "RPO-OME-RGG-int", "G")
    assert 1.38 < np.linalg.norm(m["CM2"] - m["O2'"]) < 1.46
    _, e, _ = place(rl, lib, "RPO-MOE-RAA-int", "A")
    for a, b in (("O2'", "C6'"), ("C6'", "C7'"), ("C7'", "O7'"), ("O7'", "C8'")):
        assert 1.38 < np.linalg.norm(e[a] - e[b]) < 1.56, (a, b)
    _, c5, _ = place(rl, lib, "DPO-DC2-M5C-int", "C", "DNA", 4)
    assert 1.46 < np.linalg.norm(c5["C8"] - c5["C5"]) < 1.54
    ring = np.array([c5[n] for n in ("N1", "C2", "N3", "C4", "C5", "C6")])
    centre = ring.mean(0)
    normal = np.linalg.svd(ring - centre)[2][-1]
    assert abs(np.dot(c5["C8"] - centre, normal)) < 0.1


def test_lna_bridge_closes(rl, lib):
    for model in (20, 54):
        chem = "RNA" if model == 20 else "DNA"
        _, l, rep = place(rl, lib, "PS1-LNA-RCC-int", "C", chem, model)
        for a, b in (("C2'", "O2'"), ("O2'", "C6'"), ("C6'", "C4'")):
            assert 1.38 < np.linalg.norm(l[a] - l[b]) < 1.58, (model, a, b)


@pytest.mark.parametrize("top,bottom,model", [
    ("mG*mA*fC*mU*fA*mG*fC*mU*mG*mA", "CUGAUCGACU", 20),
    ("r(GACUAGCUGA)", "+C*+T*+G*dA*dT*dC*dG*+A*+C*+T", 22),
    ("eG*eA*C^5mC*T*A*G*C*eT*eA*eG", "CTGATCGATC", 54),
    ("CGC^5mCGAA*2TTCGCG", "GCGCTT*1AAGCGC", 4),
])
def test_modified_duplexes(msession, rl, top, bottom, model):
    out = msession.build({"top": top, "bottom": bottom, "model": model, "package": True})
    assert out["ok"], out["issues"]
    sm = out["summary"]
    assert 1.55 <= sm["o3p"][0] and sm["o3p"][1] <= 1.65
    assert sm["clashes"] == 0
    z = zipfile.ZipFile(io.BytesIO(out["zip"]))
    names = z.namelist()
    assert "dsxnabuilder/frcmod.modxna" in names
    leap = z.read("dsxnabuilder/leap.in").decode()
    from dsxnabuilder.offlib import parse_off
    from dsxnabuilder.structure import parse_pdb, residues
    units = {}
    for resname in set(sm["name_map"].values()):
        text = z.read(f"dsxnabuilder/libs/{resname}.lib").decode()
        assert f"loadoff libs/{resname}.lib" in leap
        units[resname] = parse_off(text)[resname]
    # every residue: atoms written == heavy atoms of its unit; integer strand charge
    pdb = parse_pdb(z.read("dsxnabuilder/structure.pdb").decode())
    by_chain = {}
    for (chain, seq), ats in residues(pdb).items():
        by_chain.setdefault(chain, []).append(ats)
    for chain, rs in by_chain.items():
        q = 0.0
        for k, ats in enumerate(rs):
            rn = ats[0].resname
            names_ = {a.name for a in ats}
            if rn in units:
                assert names_ == set(units[rn].heavy_atoms()), (rn,)
                q += units[rn].net_charge
            else:
                unit = rn + ("5" if k == 0 else "3" if k == len(rs) - 1 else "")
                assert names_ == set(AMBER[unit]["heavy"]), unit
                q += AMBER[unit]["charge"]
        assert abs(q - round(q)) < 1e-3, (chain, q)
