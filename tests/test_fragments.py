"""The modXNA fragment list of the library tab."""

import json
import os

import pytest

from dsxnabuilder.fragments import fragment_table
from dsxnabuilder.seqlang import parse_duplex

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.fixture(scope="module")
def table():
    return fragment_table()


def by_code(table):
    return {c: (s, e) for s in table["sections"] for e in s["entries"] for c in e["codes"]}


def test_every_fragment_of_the_release_is_listed_once(table):
    with open(os.path.join(ROOT, "data", "modxna_catalog.json")) as fh:
        release = {(f["moiety"], f["code"]) for f in json.load(fh)["fragments"]}
    listed = [(s["moiety"], c) for s in table["sections"] for e in s["entries"] for c in e["codes"]]
    assert len(listed) == len(set(listed)) == len(release) == table["n_fragments"]
    assert set(listed) == release
    assert [s["title"] for s in table["sections"]][:2] == ["Backbone", "Sugar / Linker"]
    assert [s["parent"] for s in table["sections"] if s["moiety"] == "base"] == \
        ["A", "C", "G", "T", "U", "other"]
    assert len(table["names_source"]["sha256"]) == 64


def test_available_fragments_are_those_of_the_residue_library(table):
    with open(os.path.join(ROOT, "data", "residue_index.json")) as fh:
        used = {c for r in json.load(fh)["residues"].values() for c in r["triplet"]}
    avail = {c for c, (_, e) in by_code(table).items() if e["available"]}
    assert avail == used
    # fragments with a symbol in the language show it; the others are written as triplets
    for c in ("DC2", "RC3", "OME", "AF2", "MOE", "LNA", "PS1", "PS2", "M5C", "DAA", "RUU"):
        assert by_code(table)[c][1]["write"], c
    for c, (_, e) in by_code(table).items():
        assert e["available"] or e["reason"], c
    assert set(table["prebuilt"]) == set(json.load(open(os.path.join(
        ROOT, "data", "residue_index.json")))["residues"])


def test_notation_shown_is_accepted_by_the_parser(table, catalog, aliases):
    """Each prefix, linkage and base symbol the tab shows writes that fragment."""
    codes = by_code(table)
    for code in ("DC2", "RC3", "OME", "AF2", "MOE", "LNA"):
        for sym in codes[code][1]["write"].split(" or "):
            dup = parse_duplex(f"{sym}A{sym}C", "UG", catalog, aliases, ("RNA", "RNA"))
            assert {r.sugar for r in dup.top.residues} == {code}, (code, sym)
    for code in ("PS1", "PS2"):
        for sym in codes[code][1]["write"].split(" or "):
            dup = parse_duplex(f"rA{sym}rC", "UG", catalog, aliases, ("RNA", "RNA"))
            assert dup.top.residues[1].backbone == code, (code, sym)
    dup = parse_duplex(codes["M5C"][1]["write"] + "G", "GC", catalog, aliases, ("DNA", "DNA"))
    assert dup.top.residues[0].base == "M5C"


def test_names_and_flags(table):
    codes = by_code(table)
    assert codes["OME"][1]["name"].endswith("O-methyl")
    # backbones named only by their label on the page take the name from the notes
    assert codes["PS1"][1]["name"].startswith("Thiophosphate")
    assert any(f["kind"] == "stereo" and "not confirmed" in f["text"] for f in codes["PS1"][1]["flags"])
    assert any(f["kind"] == "name" and "thymine" in f["text"] for f in codes["DTT"][1]["flags"])
    assert any(f["kind"] == "protected" for f in codes["FWC"][1]["flags"])
