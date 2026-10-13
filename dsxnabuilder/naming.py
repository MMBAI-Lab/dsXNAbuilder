"""
dsxnabuilder.naming
===================

Residue names for the output PDB and the Amber libraries.

Canonical residues keep their Amber names (DA, DC, DG, DT; A, C, G, U); tleap
turns the first and last residue of each chain into the 5'/3' variants
through ``addPdbResMap``. Every other residue is a modXNA residue and gets a
deterministic three-character name (X01..X99, Y01.., ...) in order of first
appearance; the map goes into the manifest.

The key includes the terminal position: the 5'-capped, internal and 3'-capped
versions of the same triplet are different OFF units (different atoms and
charges), so they must not share a name.

The user may name any of these residues (``custom``: key -> name). A name is
one to three letters or digits, unique, and none of the names tleap or the
PDB already give a meaning to (Amber nucleotides, 3DR, water); automatic names
skip the ones taken.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from .seqlang import ResidueSpec

_NAME = re.compile(r"^[A-Z0-9]{1,3}$")
RESERVED = {n + s for n in ("DA", "DC", "DG", "DT", "A", "C", "G", "U")
            for s in ("", "3", "5", "N")} | {"3DR", "OHE", "WAT", "HOH", "DOD", "NA", "CL",
                                              "K", "MG", "ZN", "CA", "UNK"}


def check_names(custom: Dict[str, str]) -> List[str]:
    """Problems with user-chosen residue names (empty list when all are valid)."""
    out = []
    seen: Dict[str, str] = {}
    for key, name in custom.items():
        if not _NAME.match(name):
            out.append(f"{name!r} ({key}): use one to three letters or digits")
        elif name in RESERVED:
            out.append(f"{name} ({key}) is already a standard residue name")
        elif name in seen:
            out.append(f"{name} is given to both {seen[name]} and {key}")
        seen.setdefault(name, key)
    return out


class ResidueNamer:
    PREFIXES = "XYZQWVUT"

    def __init__(self, custom: Optional[Dict[str, str]] = None) -> None:
        self.custom = {k: v.strip().upper() for k, v in (custom or {}).items() if v and v.strip()}
        self._map: Dict[str, str] = {}
        self._n = 0

    def name_for(self, spec: ResidueSpec, modxna: bool) -> str:
        """``modxna`` forces a modXNA unit even for a canonical triplet (used
        when a strand end must come from modXNA to keep integer charges).

        A 2'-deoxy abasic site with a phosphodiester is the PDB chemical
        component 3DR (1',2'-dideoxyribofuranose-5'-phosphate); other abasic
        sites get an X name described in the PDB header."""
        if spec.is_abasic and spec.sugar == "DC2" and spec.backbone == "DPO":
            return "3DR"
        if not modxna and spec.is_canonical:
            return ("D" + spec.letter) if spec.chem == "DNA" else spec.letter
        key = spec.key()
        if key not in self._map:
            if key in self.custom:
                self._map[key] = self.custom[key]
                return self._map[key]
            taken = set(self.custom.values())
            while True:
                p, k = divmod(self._n, 99)
                if p >= len(self.PREFIXES):
                    raise RuntimeError("more than 792 distinct modified residues")
                name = f"{self.PREFIXES[p]}{k + 1:02d}"
                self._n += 1
                if name not in taken:
                    break
            self._map[key] = name
        return self._map[key]

    @property
    def used_custom(self) -> Dict[str, str]:
        """User-chosen names that were applied (keys present in this duplex)."""
        return {k: v for k, v in self.custom.items() if self._map.get(k) == v}

    @property
    def mapping(self) -> Dict[str, str]:
        return dict(self._map)
