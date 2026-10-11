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
"""

from __future__ import annotations

from typing import Dict

from .seqlang import ResidueSpec


class ResidueNamer:
    PREFIXES = "XYZQWVUT"

    def __init__(self) -> None:
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
            p, k = divmod(self._n, 99)
            if p >= len(self.PREFIXES):
                raise RuntimeError("more than 792 distinct modified residues")
            self._map[key] = f"{self.PREFIXES[p]}{k + 1:02d}"
            self._n += 1
        return self._map[key]

    @property
    def mapping(self) -> Dict[str, str]:
        return dict(self._map)
