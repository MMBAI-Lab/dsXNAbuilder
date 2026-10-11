"""Base swapping in the standard frame keeps a chemically sane glycosidic bond.

In the generic models 3DNA provides every letter on the same backbone, so a
standard base placed in the A slot frame can be compared with 3DNA's own
nucleotide for that letter.
"""

import numpy as np
import pytest

from dsxnabuilder.nucleotide import base_swap, chi, place_o2prime, to_rna

GENERIC = [1, 4, 7, 46, 47, 54, 55]


@pytest.mark.parametrize("mid", GENERIC)
def test_swap_matches_3dna_letter(lib, mid):
    m = lib.models[mid]
    ss = m.slots[0].I
    a = ss.templates["A"]
    for letter in "CGT":
        std = lib.standard[letter]
        names, xyz = base_swap(a.names, a.xyz, std.names, std.xyz)
        got = dict(zip(names, xyz))
        ref = dict(zip(ss.templates[letter].names, ss.templates[letter].xyz))
        base_atoms = [n for n in ref if n in got and n not in
                      ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")]
        rmsd = np.sqrt(np.mean([np.sum((got[n] - ref[n]) ** 2) for n in base_atoms]))
        # The fiber sources used their own base geometries (e.g. the C of van
        # Dam & Levitt #46/47 is 0.18 A RMSD from Olson's Atomic_C), so this
        # bounds the base-geometry difference; placement itself is checked by
        # the glycosidic bond and chi below.
        assert rmsd < 0.25, (mid, letter, rmsd)
        gly = "N9" if "N9" in got else "N1"
        bond = np.linalg.norm(got["C1'"] - got[gly])
        assert 1.42 < bond < 1.53, (mid, letter, bond)
        assert abs(((chi(got) - chi(ref) + 180) % 360) - 180) < 10.0


def test_swap_on_every_slot_keeps_glycosidic_bond(lib):
    for mid, m in lib.models.items():
        for s in m.slots:
            for ss in (s.I, s.II):
                t = next(iter(ss.templates.values()))
                for letter in "ACGTU":
                    std = lib.standard[letter]
                    names, xyz = base_swap(t.names, t.xyz, std.names, std.xyz)
                    c = dict(zip(names, xyz))
                    gly = "N9" if "N9" in c else "N1"
                    d = np.linalg.norm(c["C1'"] - c[gly])
                    assert 1.38 < d < 1.56, (mid, letter, d)


def test_o2prime_matches_rna_fiber(lib):
    """O2' built on the A-RNA sugar lands on the 3DNA O2' (same face)."""
    for mid in (20, 21):
        t = lib.models[mid].slots[0].I.templates
        t = next(iter(t.values()))
        c = dict(zip(t.names, t.xyz))
        assert "O2'" in c
        built = place_o2prime({k: v for k, v in c.items() if k != "O2'"})
        assert np.linalg.norm(built - c["O2'"]) < 0.25, (mid, np.linalg.norm(built - c["O2'"]))


def test_to_rna_adds_one_atom(lib):
    t = lib.models[4].slots[0].I.templates["A"]
    names, xyz = to_rna(t.names, t.xyz)
    assert names.count("O2'") == 1 and len(names) == len(t.names) + 1
    c = dict(zip(names, xyz))
    assert 1.38 < np.linalg.norm(c["O2'"] - c["C2'"]) < 1.45
