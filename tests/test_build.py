"""Duplex construction, packaging and Amber preflight (no AmberTools needed)."""

import io
import json
import os
import zipfile

import numpy as np
import pytest

from dsxnabuilder.build import BuildError, build_duplex
from dsxnabuilder.nucleotide import PHOSPHATE_ATOMS

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
with open(os.path.join(ROOT, "data", "amber_residues.json")) as fh:
    AMBER = json.load(fh)["units"]


def B(top, bottom, model, lib, catalog, aliases):
    return build_duplex(top, bottom, model, lib, catalog, aliases)


def amber_unit(r, first, last):
    base = r.resname
    return base + ("5" if first else "3" if last else "")


def preflight(res):
    """Every atom written, hydrogens included, must exist in the Amber unit
    tleap will match, and none may be missing (tleap would add it)."""
    by_chain = {}
    for r in res.residues:
        by_chain.setdefault(r.chain, []).append(r)
    for chain, rs in by_chain.items():
        for k, r in enumerate(rs):
            unit = amber_unit(r, k == 0, k == len(rs) - 1)
            allowed = set(AMBER[unit]["atoms"])
            extra = set(r.names) - allowed
            assert not extra, (chain, r.resseq, unit, extra)
            missing = allowed - set(r.names)
            assert not missing, (chain, r.resseq, unit, missing)


def test_dickerson_b_dna(lib, catalog, aliases):
    res = B("CGCGAATTCGCG", "GCGCTTAAGCGC", 4, lib, catalog, aliases)
    assert len(res.residues) == 24 and res.n_clashes == 0
    assert all(1.55 <= d <= 1.65 for *_, d in res.o3p)
    preflight(res)
    first = res.residues[0]
    assert first.spec.cap == "5cap" and not set(first.names) & set(PHOSPHATE_ATOMS)
    # integer charge with Amber terminal units
    charge = 0.0
    for chain in "AB":
        rs = [r for r in res.residues if r.chain == chain]
        for k, r in enumerate(rs):
            charge += AMBER[amber_unit(r, k == 0, k == len(rs) - 1)]["charge"]
    assert abs(charge - (-22.0)) < 1e-3


@pytest.mark.parametrize("model", [1, 4, 7, 46, 47, 54, 55, 12, 18, 51])
def test_random_sequences_on_dna_models(lib, catalog, aliases, model):
    rng = np.random.default_rng(model)
    seq = "".join(rng.choice(list("ACGT"), 16))
    comp = "".join({"A": "T", "T": "A", "G": "C", "C": "G"}[c] for c in seq)
    res = B(seq, comp, model, lib, catalog, aliases)
    assert all(1.55 <= d <= 1.65 for *_, d in res.o3p)
    if lib.models[model].sequence_policy == "any":
        assert res.n_clashes == 0, (model, res.n_clashes)
    else:
        # sequence-specific fibers (e.g. poly dA:dT) can be tight for G:C;
        # the contacts must then be reported
        assert res.n_clashes == 0 or any(i.code == "CLASHES" for i in res.issues)
    preflight(res)


def test_rna_and_hybrid_orientation(lib, catalog, aliases):
    res = B("GGCGCUUGCGUC", "CCGCGAACGCAG", 20, lib, catalog, aliases)
    assert all(r.resname in ("A", "C", "G", "U") for r in res.residues)
    preflight(res)
    # model 23 has DNA on strand I: an RNA strand 1 is flipped onto strand II
    res = B("r(GACU)", "d(CTGA)", 23, lib, catalog, aliases)
    assert res.flipped
    a = [r for r in res.residues if r.chain == "A"]
    assert all("O2'" in r.names for r in a)
    z = [r.xyz[r.names.index("C1'")][2] for r in a]
    assert z == sorted(z), "chain A must run 5'->3' along +z after the flip"
    preflight(res)


def test_rna_on_b_dna_gets_o2prime_and_warning(lib, catalog, aliases):
    res = B("r(GACU)", "r(AGUC)", 4, lib, catalog, aliases)
    assert all("O2'" in r.names for r in res.residues)
    assert any(i.code == "PUCKER_FORCED" for i in res.issues)
    preflight(res)


def test_z_dna_phase_and_handedness(lib, catalog, aliases):
    res = B("CGCGCG", "GCGCGC", 15, lib, catalog, aliases)
    assert res.phase == 1
    assert not [i for i in res.issues if i.code in ("SEQ_PATTERN", "SYN_PYRIMIDINE")]
    res = B("CCGCGC", "GGCGCG", 15, lib, catalog, aliases)
    assert any(i.code in ("SEQ_PATTERN", "SYN_PYRIMIDINE") for i in res.issues)


def test_mismatch_and_overhang_issues(lib, catalog, aliases):
    res = B("--ACGTGCA", "TTTGCACAT", 4, lib, catalog, aliases)
    codes = {i.code for i in res.issues}
    assert "OVERHANG" in codes and "MISMATCH" in codes
    assert len([r for r in res.residues if r.chain == "A"]) == 7
    assert len([r for r in res.residues if r.chain == "B"]) == 9
    assert all(1.55 <= d <= 1.65 for *_, d in res.o3p)


def test_modified_residue_needs_library(lib, catalog, aliases):
    with pytest.raises(BuildError) as e:
        B("mA*mCG", "UGC", 20, lib, catalog, aliases)
    assert any(i.code == "NOT_IN_LIBRARY" for i in e.value.issues)


def test_lna_in_south_model_is_an_error(lib, catalog, aliases):
    with pytest.raises(BuildError) as e:
        B("+ACG", "TGC", 4, lib, catalog, aliases)
    assert any(i.code == "PUCKER_LOCKED" for i in e.value.issues)


def test_deterministic(lib, catalog, aliases):
    a = B("ACGTTGCA", "TGCAACGT", 4, lib, catalog, aliases)
    b = B("ACGTTGCA", "TGCAACGT", 4, lib, catalog, aliases)
    assert all(np.array_equal(x.xyz, y.xyz) for x, y in zip(a.residues, b.residues))


def test_default_output_is_a_self_describing_pdb(session):
    out = session.build({"top": "5'-CGCGAATTCGCG-3'", "bottom": "3'-GCGCTTAAGCGC-5'",
                         "model": 4})
    assert out["ok"] and "zip" not in out
    lines = out["pdb"].splitlines()
    atoms = [l for l in lines if l.startswith("ATOM")]
    # 486 heavy atoms + 272 hydrogens: the 758 atoms of the Amber topology
    assert len(atoms) == 758 and sum(l[76:78].strip() == "H" for l in atoms) == 272
    assert out["summary"]["n_atoms"] == 758 and out["summary"]["n_heavy"] == 486
    assert any("CGCGAATTCGCG" in l for l in lines if l.startswith("REMARK"))
    assert all(len(l) <= 80 for l in lines)


def test_service_package(session):
    out = session.build({"top": "5'-CGCGAATTCGCG-3'", "bottom": "3'-GCGCTTAAGCGC-5'",
                         "model": 4, "package": True, "amber": {"water": "opc"}})
    assert out["ok"]
    z = zipfile.ZipFile(io.BytesIO(out["zip"]))
    names = set(z.namelist())
    for f in ("structure.pdb", "leap.in", "run.sh", "manifest.json", "README.txt"):
        assert f"dsxnabuilder/{f}" in names
    assert "dsxnabuilder/min.in" not in names            # simulation is up to the user
    assert z.read("dsxnabuilder/structure.pdb").decode() == out["pdb"]
    man = json.loads(z.read("dsxnabuilder/manifest.json"))
    assert man["model"]["id"] == 4 and man["fiber_data"]["fiber_models_sha256"]
    leap = z.read("dsxnabuilder/leap.in").decode()
    assert "leaprc.DNA.OL15" in leap and "solvateOct" in leap and "OPCBOX" in leap
    assert (z.getinfo("dsxnabuilder/run.sh").external_attr >> 16) & 0o111


def test_service_validate_reports_columns(session):
    v = session.validate({"top": "ACGGG", "bottom": "TGCTA", "model": 4})
    assert v["ok"]
    assert [c["pair"] for c in v["columns"]] == ["wc", "wc", "wc", "wobble", "mismatch"]
    v = session.validate({"top": "ACG", "bottom": "TG", "model": 4})
    assert not v["ok"] and v["issues"][0]["code"] == "LENGTH_MISMATCH"


def test_abasic_site(session):
    out = session.build({"top": "CGCGA-TTCGCG", "bottom": "GCGCTTAAGCGC", "model": 4})
    assert out["ok"]
    assert any(i["code"] == "ABASIC" for i in out["issues"])
    res = [l for l in out["pdb"].splitlines() if l.startswith("ATOM") and l[17:20] == "3DR"]
    names = {l[12:16].strip() for l in res if l[76:78].strip() != "H"}
    assert names == {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
    # hydrogens named as in the PDB component 3DR (two on C1', H4'1)
    hyd = {l[12:16].strip() for l in res if l[76:78].strip() == "H"}
    assert hyd == {"H5'", "H5''", "H4'1", "H3'", "H2'", "H2''", "H1'1", "H1'2"}
    assert 1.55 <= out["summary"]["o3p"][0] and out["summary"]["o3p"][1] <= 1.65
    assert out["summary"]["clashes"] == 0
    remarks = [l for l in out["pdb"].splitlines() if l.startswith("REMARK")]
    assert any("abasic sites" in l for l in remarks) and any("3DR:" in l for l in remarks)
    pkg = session.build({"top": "CGCGA-TTCGCG", "bottom": "GCGCTTAAGCGC", "model": 4,
                         "package": True})
    assert "zip" not in pkg and any(i["code"] == "ABASIC_NO_PARAMETERS" for i in pkg["issues"])


def test_abasic_rna_keeps_o2prime(session):
    out = session.build({"top": "r(GGCG-UUGCGUC)", "bottom": "r(CCGCGAACGCAG)", "model": 20})
    assert out["ok"]
    abasic = [l for l in out["pdb"].splitlines()
              if l.startswith("ATOM") and l[17:20].strip() == "X01"]
    assert {l[12:16].strip() for l in abasic} >= {"C1'", "C2'", "O2'"}
    assert not {"N9", "N1"} & {l[12:16].strip() for l in abasic}


@pytest.mark.skipif(not os.path.isdir(os.path.join(ROOT, "web", "lib", "residues")),
                    reason="residue library not built (tools/build_residue_library.py)")
def test_residue_names_chosen_by_the_user(session):
    from dsxnabuilder.reslib import ResidueLibrary
    from dsxnabuilder.service import Session
    s = Session(residue_library=ResidueLibrary.load_default())
    req = {"top": "mA*mC*fG", "bottom": "UGC", "model": 20,
           "names": {"PS1/OME/RCC/int": "omc", "PS1/AF2/RGG/3cap": "X01"}}
    out = s.build(dict(req, package=True))
    assert out["ok"], out["issues"]
    names = {n["key"]: n["name"] for n in out["summary"]["named"]}
    # user names are applied (upper case); automatic ones skip the names taken
    assert names == {"5PO/OME/RAA/5cap": "X02", "PS1/OME/RCC/int": "OMC",
                     "PS1/AF2/RGG/3cap": "X01"}
    resnames = {l[17:20].strip() for l in out["pdb"].splitlines() if l.startswith("ATOM")}
    assert {"OMC", "X01", "X02"} <= resnames
    z = zipfile.ZipFile(io.BytesIO(out["zip"]))
    assert "dsxnabuilder/libs/OMC.lib" in z.namelist()
    man = json.loads(z.read("dsxnabuilder/manifest.json"))
    assert '--name "PS1/OME/RCC/int=OMC"' in man["reproduce"]
    for bad in ({"PS1/OME/RCC/int": "DA"}, {"PS1/OME/RCC/int": "ABCD"},
                {"PS1/OME/RCC/int": "QQ1", "PS1/AF2/RGG/3cap": "QQ1"}):
        out = s.build(dict(req, names=bad))
        assert not out["ok"] and out["issues"][0]["code"] == "RESIDUE_NAME", bad
