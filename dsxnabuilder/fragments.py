"""
dsxnabuilder.fragments
======================

The pieces of the web app's "modXNA library" tab, and helpers for the
offline residue-library tools.

``pieces()`` lists what the sequence language writes and the prebuilt
residue library has in every combination: the sugars (prefix letters), the
bases (letters and C^5mC) and the linkages (phosphodiester, *1, *2), with
their notation. Every sugar x base x linkage combination is prebuilt,
internal and at either end, so the tab offers only what can be built.

The other helpers (catalogue names and base families, the scaffold used to
try a residue) serve tools/build_residue_library.py, which can also build
the extended set of every catalogue fragment with canonical partners.
"""

from __future__ import annotations

import json
import os
import re
from typing import Dict, List, Optional

from .catalog import CANONICAL

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
# the page gives some backbones only as a label ("PS1 (S)") and their name in the notes
_LABEL = re.compile(r"^[A-Za-z0-9]{3}\s*(\((S|R)\))?$")

# Names of the pieces as shown in the tab; prefixes in the order they are offered.
SUGARS = [("d", "DC2", "2′-deoxyribose (DNA)"), ("r", "RC3", "ribose (RNA)"),
          ("m", "OME", "2′-O-methylribose (2′-OMe)"),
          ("f", "AF2", "2′-deoxy-2′-fluororibose (2′-F)"),
          ("e", "MOE", "2′-O-methoxyethylribose (2′-MOE)"),
          ("+", "LNA", "LNA: ribose with a 2′-O,4′-C-methylene bridge")]
BASES = {"A": "adenine", "C": "cytosine", "G": "guanine", "T": "thymine", "U": "uracil"}
LINKAGES = [("", ["DPO", "RPO"], "phosphodiester"),
            ("*1", ["PS1"], "phosphorothioate PS1: sulfur on OP1, keeps OP2"),
            ("*2", ["PS2"], "phosphorothioate PS2: sulfur on OP2, keeps OP1")]


def _load(name: str, data_dir: str) -> dict:
    with open(os.path.join(data_dir, name), encoding="utf-8") as fh:
        return json.load(fh)


def pieces(data_dir: Optional[str] = None) -> dict:
    """Sugars, bases and linkages of the prebuilt library, with their notation."""
    d = data_dir or DATA
    aliases = _load("aliases.json", d)
    pairing = _load("base_pairing.json", d)["pairing"]
    index = _load("residue_index.json", d)
    used = {m: set() for m in ("backbone", "sugar", "base")}
    for r in index["residues"].values():
        for moiety, code in zip(("backbone", "sugar", "base"), r["triplet"]):
            used[moiety].add(code)
    other = {}                                   # code -> other prefixes that write it
    for prefix, code in aliases.get("sugar", {}).items():
        other.setdefault(code, []).append(prefix)
    sugars = [{"prefix": p, "code": c, "name": n,
               "also": [x for x in other.get(c, []) if x != p]}
              for p, c, n in SUGARS if c in used["sugar"] and aliases["sugar"].get(p) == c]
    bases = [{"write": letter, "name": name} for letter, name in BASES.items()]
    for alias, code in aliases.get("base", {}).items():
        if code in used["base"]:
            bases.append({"write": f"{pairing[code]}^{alias}", "name": _BASE_NAMES.get(code, code)})
    linkages = [{"write": w, "codes": cs, "name": n} for w, cs, n in LINKAGES
                if all(c in used["backbone"] for c in cs if c not in ("DPO", "RPO"))]
    return {"sugars": sugars, "bases": bases, "linkages": linkages,
            "n_residues": len(index["residues"]),
            "modxna": {"version": index.get("modxna_version"),
                       "commit": index.get("modxna_commit")}}


_BASE_NAMES = {"M5C": "5-methylcytosine"}


def notation(aliases: dict, pairing: Dict[str, str]) -> Dict[str, str]:
    """How to write each fragment that has a symbol in the sequence language
    (a base alias goes after the letter it pairs as, e.g. C^5mC)."""
    out: Dict[str, List[str]] = {}
    for prefix, code in aliases.get("sugar", {}).items():
        out.setdefault(code, []).append(prefix)
    unverified = set(aliases.get("unverified", []))
    for sym, code in aliases.get("linkage", {}).items():
        if sym not in unverified:
            out.setdefault(code, []).append(sym)
    for alias, code in aliases.get("base", {}).items():
        out.setdefault(code, []).append(pairing.get(code, "") + "^" + alias)
    for (_, letter), code in CANONICAL.items():
        out.setdefault(code, []).append(letter)
    write = {code: " or ".join(dict.fromkeys(syms)) for code, syms in out.items()}
    write.update({"DPO": "default linkage", "RPO": "default linkage", "5PO": "5′ end (automatic)"})
    return write


def _display(name: str, notes: List[str]):
    """(name, remaining notes): a bare label ("PS1 (S)") gives way to the first note."""
    if _LABEL.match(name) and notes:
        return notes[0], notes[1:]
    return name, notes


def fragment_names(data_dir: Optional[str] = None) -> Dict[str, str]:
    """Fragment code -> name shown in the catalogue."""
    names = _load("modxna_names.json", data_dir or DATA)
    out = {}
    for s in names["sections"]:
        for e in s["entries"]:
            for c in e["codes"]:
                out[c] = _display(e["name"], e["notes"])[0]
    return out


def base_families(data_dir: Optional[str] = None) -> Dict[str, str]:
    """Base code -> parent base of its catalogue section (A, C, G, T, U or other)."""
    names = _load("modxna_names.json", data_dir or DATA)
    return {c: s["parent"] for s in names["sections"] if s["moiety"] == "base"
            for e in s["entries"] for c in e["codes"]}


CANONICAL_LETTER = {code: letter for (_, letter), code in CANONICAL.items()}


def test_context(triplet, catalog, families: Dict[str, str]):
    """(letter, chem, model) to try a residue on a fiber scaffold: the pairing
    letter of its base (canonical, tabulated, or the catalogue family; purine
    or pyrimidine for "other" bases), the strand chemistry its sugar implies
    (ribo-type when 2'-substituted) and a model with its sugar pucker (A-RNA
    #20 for North, B-DNA #4 otherwise)."""
    bb, sugar, base = triplet
    letter = CANONICAL_LETTER.get(base) or {"M5C": "C"}.get(base) or families.get(base)
    if letter in (None, "other"):
        letter = "A" if catalog[base].base_class == "purine" else "C"
    chem = "RNA" if set(catalog[sugar].atom_names) & {"O2'", "F2'"} else "DNA"
    model = 20 if catalog[sugar].pucker_domain == "North" else 4
    return letter, chem, model
