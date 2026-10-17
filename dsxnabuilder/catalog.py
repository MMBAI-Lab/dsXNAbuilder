# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.catalog
====================

The modXNA fragment catalogue and the canonical fragment codes.

Each modXNA fragment is a mol2 file whose title line encodes its connection
points:

    modXNA[-5'|-3']-fragment:<NAME>:HEAD01:@<atom>:HEAD01STRIP:@<atom,...>
                             [:TAIL01:@<atom>:TAIL01STRIP:@<atom,...>]
                             [:HEAD03:@<atom>:HEAD03STRIP:@<atom,...>]

  backbone : HEAD01 = @P    (bonds to O3' of the preceding residue)
             TAIL01 = @O5'  (bonds to C5' of the sugar)
  sugar    : HEAD01 = @C5'  (bonds to TAIL01 of the backbone)
             TAIL01 = @O3'  (bonds to P of the next residue)
             HEAD03 = @C1'  (bonds to the base)
  base     : HEAD01 = @N9 (purines) / @N1 (pyrimidines)

The ``*STRIP`` atoms are the caps modXNA removes when it joins the fragments.
data/modxna_catalog.json (tools/build_catalog.py) keeps, per fragment, what
the builder needs: its atom names, the pucker of a sugar and the class of a
base.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .paths import DATA

# (strand chemistry, letter) -> modXNA base code
CANONICAL: Dict[Tuple[str, str], str] = {
    ("DNA", "A"): "DAA", ("DNA", "C"): "DCC", ("DNA", "G"): "DGG", ("DNA", "T"): "DTT",
    ("RNA", "A"): "RAA", ("RNA", "C"): "RCC", ("RNA", "G"): "RGG", ("RNA", "U"): "RUU",
}
CANONICAL_SUGAR = {"DNA": "DC2", "RNA": "RC3"}      # 2'-deoxyribose / ribose
CANONICAL_BACKBONE = {"DNA": "DPO", "RNA": "RPO"}   # phosphodiester
CAP_5_OH = "5PO"

# Names written in the PDB header and shown in the web app.
SUGAR_NAMES = {"DC2": "2'-deoxyribose", "RC3": "ribose", "OME": "2'-O-methylribose",
               "AF2": "2'-deoxy-2'-fluororibose", "MOE": "2'-O-(2-methoxyethyl)ribose",
               "LNA": "LNA (2'-O,4'-C-methylene) ribose", "CET": "cEt",
               "DS4": "4'-thio-2'-deoxyribose", "RS4": "4'-thioribose"}
BASE_NAMES = {"DAA": "adenine", "RAA": "adenine", "DCC": "cytosine", "RCC": "cytosine",
              "DGG": "guanine", "RGG": "guanine", "DTT": "thymine", "RUU": "uracil",
              "M5C": "5-methylcytosine"}
BACKBONE_NAMES = {"DPO": "phosphodiester", "RPO": "phosphodiester", "5PO": "5'-OH",
                  "PS1": "phosphorothioate PS1 (S on OP1 position)",
                  "PS2": "phosphorothioate PS2 (S on OP2 position)"}


def canonical_base(chem: str, letter: str) -> str:
    """modXNA code of a canonical base for a strand chemistry and a letter.

    The D/R prefix of the modXNA codes is the context (deoxy or ribo) in which
    the charges were derived, not a different base. Thymine exists only as
    DTT and uracil only as RUU, so T on a ribo-type sugar (LNA-T, 2'-OMe-T)
    resolves to DTT and U on 2'-deoxyribose to RUU, as in real antisense
    oligonucleotides.
    """
    code = CANONICAL.get((chem, letter))
    if code is not None:
        return code
    other = "DNA" if chem == "RNA" else "RNA"
    code = CANONICAL.get((other, letter))
    if code is None:
        raise KeyError(f"not a canonical base letter: {letter!r}")
    return code


def default_chem(chem: str) -> str:
    """Strand chemistry for plain letters on a model strand: DNA or RNA, and
    DNA for anything else (a model strand of mixed or unknown chemistry)."""
    return chem if chem in ("DNA", "RNA") else "DNA"


@dataclass
class Fragment:
    code: str
    moiety: str                          # "backbone" | "sugar" | "base"
    atom_names: List[str] = field(default_factory=list)
    sugar_geom: Optional[dict] = None    # pseudorotation of a sugar (geometry.pseudorotation)
    base_class: Optional[str] = None     # "purine" | "pyrimidine" | "other"

    @property
    def pucker_domain(self) -> Optional[str]:
        """North or South for a sugar with a furanose ring, else None."""
        return (self.sugar_geom or {}).get("domain")

    @property
    def is_ribo_type(self) -> bool:
        """A sugar substituted at 2' (O2' or F2'): takes the RNA base set."""
        return bool(set(self.atom_names) & {"O2'", "F2'"})


class Catalog:
    """The modXNA fragments by code."""

    def __init__(self, fragments: Sequence[Fragment]):
        self.fragments = list(fragments)
        self._by_code: Dict[str, Fragment] = {f.code: f for f in self.fragments}

    @classmethod
    def from_json(cls, path: Optional[str] = None) -> "Catalog":
        """Read data/modxna_catalog.json (or ``path``)."""
        with open(path or os.path.join(DATA, "modxna_catalog.json"), encoding="utf-8") as fh:
            blob = json.load(fh)
        return cls([Fragment(code=row["code"], moiety=row["moiety"],
                             atom_names=row.get("atom_names", []),
                             sugar_geom=row.get("sugar_geom"),
                             base_class=row.get("base_class"))
                    for row in blob["fragments"]])

    def __getitem__(self, code: str) -> Fragment:
        try:
            return self._by_code[code]
        except KeyError:
            raise KeyError(f"unknown modXNA fragment: {code!r}") from None

    def __contains__(self, code: str) -> bool:
        return code in self._by_code
