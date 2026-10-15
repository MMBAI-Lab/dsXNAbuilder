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
    import sys
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    from build_residue_library import curated
    # the web app serves the curated set: what the language writes without explicit blocks
    # (tools/build_residue_library.py --set extended can add every catalogue fragment)
    assert set(rl.residues) == {"-".join(j) for j in curated()}
    assert all(isinstance(v, str) and v for v in rl.index.get("unavailable", {}).values())
    assert not set(rl.residues) & set(rl.index.get("unavailable", {}))
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
    from dsxnabuilder.geometry import signed_volume
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
    ("GCGAeUAGCG", "r(CGCUAUCGC)", 20),          # library MOE chain hit C5' of i+1
    ("GCGAeCAGCG", "r(CGCUGUCGC)", 22),
    ("eC*eG*eT*eT*eA*dC*dA*dT*dG*dA*dC*dT*dT*dA*dC*eG*eT*eA*eC*eG",
     "r(GCAAUGUACUGAAUGCAUGC)", 22),             # 5-10-5 MOE gapmer
])
def test_moe_chain_avoids_contacts(msession, rl, top, bottom, model):
    """2'-O chains are rotated out of contacts; bonds and angles stay the library's."""
    from dsxnabuilder.structure import parse_pdb, residues
    out = msession.build({"top": top, "bottom": bottom, "model": model})
    assert out["ok"], out["issues"]
    assert out["summary"]["clashes"] == 0
    assert "CLASHES" not in [i["code"] for i in out["issues"]]
    ref = rl.unit("RPO-MOE-RAA-int").coords()
    chain = ["O2'", "C6'", "C7'", "O7'", "C8'"]
    for ats in residues(parse_pdb(out["pdb"])).values():
        c = {a.name: np.asarray(a.xyz) for a in ats}
        if "C7'" not in c:
            continue
        # the transplant keeps bonds within 0.03 A of the library (O2' is a scaffold
        # atom); the torsion scan is rigid and must not add to that
        for a, b in zip(chain, chain[1:]):
            assert abs(np.linalg.norm(c[a] - c[b]) - np.linalg.norm(ref[a] - ref[b])) < 0.03
        for a, b, d in zip(chain, chain[1:], chain[2:]):   # bond angles via 1-3 distances
            assert abs(np.linalg.norm(c[a] - c[d]) - np.linalg.norm(ref[a] - ref[d])) < 0.06


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
    # every residue: atoms written == atoms of its unit (hydrogens included);
    # integer strand charge
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
                assert names_ == set(units[rn].atom_names), (rn,)
                q += units[rn].net_charge
            else:
                unit = rn + ("5" if k == 0 else "3" if k == len(rs) - 1 else "")
                assert names_ == set(AMBER[unit]["atoms"]), unit
                q += AMBER[unit]["charge"]
        assert abs(q - round(q)) < 1e-3, (chain, q)


def _c2_handedness(xyz):
    """Sign of the C2' centre (C1', C3', O2'), as in reslib._chirality."""
    o = xyz["C2'"]
    u = [(xyz[n] - o) / np.linalg.norm(xyz[n] - o) for n in ("C1'", "C3'", "O2'")]
    return np.sign(np.dot(u[0], np.cross(u[1], u[2])))


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "web", "lib", "residues",
                                                   "RPO-A5O-RAA-int.lib")),
                    reason="arabinose is only in the extended library (--set extended)")
def test_arabinose_o2_on_its_own_face(msession, rl):
    """The scaffold puts O2' on the ribose face; arabinose has it on the other one, so the
    transplant must place it from the residue, not keep the scaffold's."""
    from dsxnabuilder.structure import parse_pdb, residues
    lib = rl.unit("RPO-A5O-RAA-int").coords()
    rib = rl.unit("RPO-OME-RAA-int").coords()
    assert _c2_handedness(lib) != _c2_handedness(rib)
    out = msession.build({"top": "r(GC){RPO/A5O/RAA:A}r(GC)", "bottom": "r(GCUGC)",
                          "model": 20})
    assert out["ok"], out["issues"]
    ats = residues(parse_pdb(out["pdb"]))[("A", 3)]
    built = {a.name: np.asarray(a.xyz) for a in ats}
    assert _c2_handedness(built) == _c2_handedness(lib)


def test_transplant_refuses_what_the_scaffold_cannot_carry(rl, lib):
    """LNA on a C2'-endo scaffold: the bridge would run through O3'."""
    from dsxnabuilder.reslib import transplant
    names, xyz = scaffold(lib, 4, "C", "RNA")
    _, _, rep = transplant(rl.unit("PS1-LNA-RCC-int"), names, xyz)
    assert rep["contacts"]
    _, _, rep = transplant(rl.unit("PS1-LNA-RCC-int"), *scaffold(lib, 20, "C", "RNA"))
    assert not (rep["contacts"] or rep["distorted"] or rep["stereo"])


def test_transplant_is_independent_of_hash_seed():
    """Ties in the placement order are broken by name, not by set order."""
    import subprocess
    import sys
    code = ("import sys, hashlib; sys.path.insert(0, %r)\n"
            "from dsxnabuilder.service import Session\n"
            "from dsxnabuilder.reslib import ResidueLibrary\n"
            "s = Session(residue_library=ResidueLibrary.load_default())\n"
            "o = s.build({'top': 'mG*eA*+C^5mC*fU*mG', 'bottom': 'CUGGC', 'model': 20})\n"
            "print(o['ok'], hashlib.sha256(o['pdb'].encode()).hexdigest())" % ROOT)
    outs = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONHASHSEED=str(seed))).stdout
            for seed in (1, 2, 3)}
    assert len(outs) == 1 and next(iter(outs)).startswith("True"), outs


def test_atom_names_as_tleap_reads_them(rl):
    from dsxnabuilder.reslib import _rename_atoms_off, duplicate_names, tleap_name
    assert tleap_name("1H15") == "H151" and tleap_name("C1'") == "C1'"
    assert _rename_atoms_off('"1H15" "HC" 0\n"MXR"') == '"H151" "HC" 0\n"MXR"'
    for key in rl.residues:
        u = rl.unit(key)
        assert not duplicate_names(u), key
        assert not any(n[0].isdigit() for n in u.atom_names), key


def test_modxna_code_example(msession):
    out = msession.build({"top": "CG{PS2/OME/M5C:C}GAA*2TTC-CG--",
                          "bottom": "GCGCTT*1AAGCGCAA", "model": 4})
    assert out["ok"], out["issues"]
    codes = {i["code"] for i in out["issues"]}
    assert {"ABASIC", "OVERHANG", "PUCKER_FORCED"} <= codes and "CLASHES" not in codes
    assert out["summary"]["named"][0]["key"] == "PS2/OME/M5C/int"
    assert "5-methylcytosine" in out["summary"]["named"][0]["description"]


def test_phosphorus_follows_the_modxna_geometry(rl, lib):
    """A group on P keeps the handedness it has in the modXNA residue (WX2 came out
    inverted when placed by oxygen name: its OP1/OP2 are named opposite to Amber's);
    phosphorothioates follow their names, since modXNA's PS1 and PS2 coincide."""
    from dsxnabuilder.reslib import transplant

    def hand(c, x, kept):
        p = c["P"]
        return np.sign(np.dot(c[x] - p, np.cross(c[kept] - p, c["O5'"] - p)))
    names, xyz = scaffold(lib, 4, "C", "DNA")
    if "WX2-DC2-DCC-int" in rl.residues:          # extended library only
        u = rl.unit("WX2-DC2-DCC-int")
        out_names, out_xyz, _ = transplant(u, names, xyz)
        built = dict(zip(out_names, out_xyz))
        assert hand(built, "N2", "OP2") == hand(u.coords(), "N2", "OP2")
    for key, s in (("PS1-DC2-DCC-int", "S1"), ("PS2-DC2-DCC-int", "S2")):
        out_names, out_xyz, _ = transplant(rl.unit(key), names, xyz)
        built = dict(zip(out_names, out_xyz))
        target = "OP1" if s == "S1" else "OP2"       # by name, on the Amber scaffold
        assert np.linalg.norm(built[s] - built["P"] - 2.02 * (dict(zip(names, xyz))[target] - built["P"])
                              / np.linalg.norm(dict(zip(names, xyz))[target] - built["P"])) < 0.05


CHAIN_ATOMS = {"CM2", "HM'1", "HM'2", "HM'3",                     # 2'-OMe
               "C6'", "C7'", "O7'", "C8'", "H6'", "H6''", "H7'", "H7''",
               "H8'1", "H8'2", "H8'3"}                            # 2'-MOE


def test_neighbouring_2o_chains_do_not_overlap(msession):
    """Each 2'-O chain is turned against the chains already turned around it,
    so chains on both strands of a fully modified duplex keep apart."""
    from dsxnabuilder.structure import element_of, parse_pdb
    from dsxnabuilder.substituents import CHAIN_CONTACT, CHAIN_H_CONTACT
    out = msession.build({"top": "eA*mA*eU*mC*eU*eG*eC*mG*eC*eG*eC*eA*eG",
                          "bottom": "eU*eU*mA*eG*mA*eC*mG*eC*mG*eC*mG*eU*mC", "model": 20})
    assert out["ok"], out["issues"]
    atoms = [a for a in parse_pdb(out["pdb"]) if a.name in CHAIN_ATOMS]
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            if (a.chain, a.resseq) == (b.chain, b.resseq):
                continue
            hydrogen = "H" in (element_of(a.name), element_of(b.name))
            limit = CHAIN_H_CONTACT if hydrogen else CHAIN_CONTACT
            d = float(np.linalg.norm(a.xyz - b.xyz))
            assert d >= limit - 1e-3, (a.chain, a.resseq, a.name, b.chain, b.resseq, b.name, d)
