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

A fragment is *available* when the prebuilt residue library uses it. Every
combination of the available backbones, sugars and bases is prebuilt, so each
available fragment also gets the notation that writes it (from the alias
table); the rest are listed for reference.
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
    sections = []
    for s in names["sections"]:
        entries = []
        for e in s["entries"]:
            name, notes = e["name"], list(e["notes"])
            if _LABEL.match(name) and notes:
                label, name = name, notes.pop(0)
            else:
                label = None
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
            entries.append({
                "codes": codes, "name": name, "notes": notes,
                "available": available,
                "write": write.get(codes[0]) if available else None,
                "in_release": all((s["moiety"], c) in in_release for c in codes),
                "flags": flags})
        sections.append({"title": s["title"], "moiety": s["moiety"], "parent": s["parent"],
                         "entries": entries})
    return {"sections": sections, "names_source": names["source"],
            "review_source": review["source"],
            "modxna": {"version": index.get("modxna_version"),
                       "commit": index.get("modxna_commit")},
            "n_fragments": len(in_release)}
