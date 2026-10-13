"""
dsxnabuilder.fragments
======================

The modXNA fragment list of the web app's "modXNA library" tab.

Every fragment of the modXNA release the residue library was built from
(data/modxna_catalog.json) is listed by section, as in the modXNA catalogue
page: backbones, sugars/linkers and base modifications grouped by parent base.
Names and notes come from that page (data/modxna_names.json, written by
tools/fetch_catalog_names.py); known errors and duplicates of the catalogue
are flagged from data/modxna_catalog_review.json.

A fragment is *available* when the prebuilt residue library uses it, at least
with canonical partners; fragments with a symbol in the sequence language also
get the notation that writes it (from the alias table). The exact residues
prebuilt (``prebuilt``, keys BACKBONE-SUGAR-BASE-cap) and those that could not
be built or do not fit a fiber scaffold (``unavailable``, with the reason) let
the residue composer of the web app say whether a given combination exists.
"""

from __future__ import annotations

import collections
import json
import os
import re
from typing import Dict, List, Optional

from .catalog import CANONICAL

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
# the page gives some backbones only as a label ("PS1 (S)") and their name in the notes
_LABEL = re.compile(r"^[A-Za-z0-9]{3}\s*(\((S|R)\))?$")
_STEREO_NOTE = ("modXNA labels it {label}. In dsXNAbuilder PS1 puts the sulfur on OP1 and PS2 "
                "on OP2 (Amber convention); their absolute Rp/Sp configuration in a duplex is "
                "not confirmed yet.")


def _load(name: str, data_dir: str) -> dict:
    with open(os.path.join(data_dir, name), encoding="utf-8") as fh:
        return json.load(fh)


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


def fragment_table(data_dir: Optional[str] = None) -> dict:
    d = data_dir or DATA
    names = _load("modxna_names.json", d)
    review = _load("modxna_catalog_review.json", d)
    catalog = _load("modxna_catalog.json", d)
    index = _load("residue_index.json", d)
    write = notation(_load("aliases.json", d), _load("base_pairing.json", d)["pairing"])
    in_release = {(f["moiety"], f["code"]) for f in catalog["fragments"]}
    used = {m: set() for m in ("backbone", "sugar", "base")}
    for r in index["residues"].values():
        for moiety, code in zip(("backbone", "sugar", "base"), r["triplet"]):
            used[moiety].add(code)
    protected = set(review.get("protected", {}).get("codes", []))
    unavailable = index.get("unavailable", {})
    why: Dict[str, collections.Counter] = {}
    for key, reason in unavailable.items():
        for code in key.split("-")[:3]:
            why.setdefault(code, collections.Counter())[_short_reason(reason)] += 1
    sections = []
    for s in names["sections"]:
        entries = []
        for e in s["entries"]:
            label = e["name"] if _LABEL.match(e["name"]) and e["notes"] else None
            name, notes = _display(e["name"], list(e["notes"]))
            codes = e["codes"]
            flags = []
            for c in codes:
                if c in review["entries"]:
                    flags.append(dict(review["entries"][c], code=c))
                if c in protected:
                    flags.append({"kind": "protected", "code": c,
                                  "text": review["protected"]["text"]})
                if c in ("PS1", "PS2"):
                    flags.append({"kind": "stereo", "code": c,
                                  "text": _STEREO_NOTE.format(label=label or c)})
            available = all(c in used[s["moiety"]] for c in codes)
            reason = None
            if not available:
                if any(c in protected for c in codes):
                    reason = "protected phosphoramidite: not prebuilt"
                else:
                    counts = sum((why.get(c, collections.Counter()) for c in codes),
                                 collections.Counter())
                    reason = counts.most_common(1)[0][0] if counts else "not prebuilt"
            entries.append({
                "codes": codes, "name": name, "notes": notes,
                "available": available, "reason": reason,
                "write": write.get(codes[0]) if available else None,
                "in_release": all((s["moiety"], c) in in_release for c in codes),
                "flags": flags})
        sections.append({"title": s["title"], "moiety": s["moiety"], "parent": s["parent"],
                         "entries": entries})
    return {"sections": sections, "names_source": names["source"],
            "prebuilt": sorted(index["residues"]), "unavailable": unavailable,
            "review_source": review["source"],
            "modxna": {"version": index.get("modxna_version"),
                       "commit": index.get("modxna_commit")},
            "n_fragments": len(in_release)}


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


def _short_reason(reason: str) -> str:
    """One reason per kind: drop the distances and the model of the scaffold check."""
    reason = re.sub(r" \(model #\d+\)", "", reason)
    reason = re.sub(r"\s+-?\d+\.\d+(/-?\d+\.\d+)? A", "", reason)
    return reason
