"""The pieces of the modXNA library tab: only what is prebuilt, in the short notation."""

import itertools
import json
import os

import pytest

from dsxnabuilder.fragments import pieces
from dsxnabuilder.seqlang import parse_duplex

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.fixture(scope="module")
def table():
    return pieces()


def test_pieces(table):
    assert [s["prefix"] for s in table["sugars"]] == ["d", "r", "m", "f", "e", "+"]
    assert [b["write"] for b in table["bases"]] == ["A", "C", "G", "T", "U", "C^5mC"]
    assert [lk["write"] for lk in table["linkages"]] == ["", "*1", "*2"]
    assert table["n_residues"] == 244


def test_every_combination_offered_is_prebuilt_and_parses(table, catalog, aliases):
    """Each sugar x base x linkage the tab writes parses to the residue it names, and that
    residue is in the library (or is a canonical Amber nucleotide) at every position."""
    with open(os.path.join(ROOT, "data", "residue_index.json")) as fh:
        prebuilt = set(json.load(fh)["residues"])
    for sug, base, link in itertools.product(table["sugars"], table["bases"], table["linkages"]):
        nt = sug["prefix"] + base["write"]
        # middle of strand 1: the linkage goes before the nucleotide
        dup = parse_duplex(f"dG{link['write']}{nt}dG", "CGC", catalog, aliases, ("DNA", "DNA"))
        r = dup.top.residues[1]
        assert r.sugar == sug["code"], (nt, r.sugar)
        assert r.backbone in link["codes"], (nt, link, r.backbone)
        if not r.is_canonical:
            for cap, bb in (("int", r.backbone), ("3cap", r.backbone), ("5cap", "5PO")):
                assert f"{bb}-{r.sugar}-{r.base}-{cap}" in prebuilt, (nt, link, cap)
        # strand 2 is typed 3'->5': the linkage goes after the nucleotide
        dup = parse_duplex("CGC", f"dG{nt}{link['write']}dG", catalog, aliases, ("DNA", "DNA"))
        assert dup.bottom.residues[1].backbone in link["codes"], (nt, link)
