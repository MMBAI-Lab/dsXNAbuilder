"""
dsxnabuilder.pdbout
===================

The default output: one self-describing PDB file.

All atoms, hydrogens included, with Amber-style residue and atom names
(DA/DC/DG/DT, A/C/G/U; OP1/OP2; thymine C7, H71-H73), chain A = strand 1 and
chain B = strand 2, both 5'->3'. Modified nucleotides keep the modXNA atom
names (hydrogens included) and get a
three-character residue name (X01, X02, ...); REMARK 999 records spell out
what each one is, together with the input, the fiber model and the data
provenance, so the file can be used and reproduced with any program.
"""

from __future__ import annotations

import textwrap
from typing import Dict, List, Optional

from . import __version__
from .build import BuildResult
from .structure import format_pdb

SUGAR_NAMES = {"DC2": "2'-deoxyribose", "RC3": "ribose", "OME": "2'-O-methylribose",
               "AF2": "2'-deoxy-2'-fluororibose", "MOE": "2'-O-(2-methoxyethyl)ribose",
               "LNA": "LNA (2'-O,4'-C-methylene) ribose", "CET": "cEt", "DS4": "4'-thio-2'-deoxyribose",
               "RS4": "4'-thioribose"}
BASE_NAMES = {"DAA": "adenine", "RAA": "adenine", "DCC": "cytosine", "RCC": "cytosine",
              "DGG": "guanine", "RGG": "guanine", "DTT": "thymine", "RUU": "uracil",
              "M5C": "5-methylcytosine"}
BACKBONE_NAMES = {"DPO": "phosphodiester", "RPO": "phosphodiester", "5PO": "5'-OH",
                  "PS1": "phosphorothioate PS1 (S on OP1 position)",
                  "PS2": "phosphorothioate PS2 (S on OP2 position)"}
CAP_NAMES = {"int": "internal", "5cap": "5'-terminal", "3cap": "3'-terminal"}


_CATALOGUE_NAMES: Optional[Dict[str, str]] = None


def catalogue_name(code: str) -> str:
    """Name of a fragment in the modXNA catalogue (the code when it has none)."""
    global _CATALOGUE_NAMES
    if _CATALOGUE_NAMES is None:
        try:
            from .fragments import fragment_names
            _CATALOGUE_NAMES = fragment_names()
        except OSError:             # data/modxna_names.json not shipped
            _CATALOGUE_NAMES = {}
    return _CATALOGUE_NAMES.get(code, code)


def describe(triplet, cap: str) -> str:
    bb, sugar, base = triplet
    return (f"{BASE_NAMES.get(base) or catalogue_name(base)}, "
            f"{SUGAR_NAMES.get(sugar) or catalogue_name(sugar)}, "
            f"{BACKBONE_NAMES.get(bb) or catalogue_name(bb)}; {CAP_NAMES.get(cap, cap)}")


def remarks(res: BuildResult, fiber_provenance: Optional[dict] = None,
            residue_index: Optional[dict] = None) -> List[str]:
    m = res.model
    out = [f"dsXNAbuilder {__version__}",
           f"strand 1 (chain A) 5'-{res.duplex.top.raw.strip()}-3'",
           f"strand 2 (chain B) 3'-{res.duplex.bottom.raw.strip()}-5'",
           f"fiber model #{m.id} {m.family}: {m.description}"[:200],
           f"helical twist {m.bp_twist} deg/bp, rise {m.bp_rise} A/bp, {m.handedness}-handed",
           "3DNA fiber models: Lu & Olson (2003) NAR 31:5108; Lu & Olson (2008) Nat Protoc 3:1213"]
    if fiber_provenance and fiber_provenance.get("fiber_models_sha256"):
        out.append(f"fiber data sha256 {fiber_provenance['fiber_models_sha256'][:16]}")
    mods: Dict[str, object] = {}
    abasic: Dict[str, object] = {}
    for r in res.residues:
        if r.flags.get("abasic") and r.resname not in abasic:
            abasic[r.resname] = r
        elif r.flags.get("modxna") and r.resname not in mods:
            mods[r.resname] = r
    if abasic:
        out.append("abasic sites (sugar-phosphate without base):")
        for name, r in sorted(abasic.items()):
            sugar = ("1',2'-dideoxyribose (tetrahydrofuran-type)" if r.spec.chem == "DNA"
                     else "1'-deoxyribose")
            ccd = " = PDB chemical component 3DR" if name == "3DR" else ""
            out.append(f"  {name}: {sugar}, {BACKBONE_NAMES.get(r.spec.backbone, r.spec.backbone)}"
                       f"{ccd}")
    if mods:
        out.append("modified nucleotides (modXNA; Love et al. (2024) JCTC 20:9354):")
        for name, r in sorted(mods.items()):
            out.append(f"  {name} = {'/'.join(r.spec.triplet())} {r.spec.cap}: "
                       f"{describe(r.spec.triplet(), r.spec.cap)}")
        if res.custom_names:
            out.append("residue names chosen by the user: " + ", ".join(
                f"{v}={k}" for k, v in sorted(res.custom_names.items())))
        if residue_index:
            out.append(f"modXNA {residue_index.get('modxna_version', '?')} "
                       f"@{str(residue_index.get('modxna_commit', '?'))[:7]}")
    out.append("hydrogens placed on the built heavy atoms from the residue templates "
               "(Amber DNA.OL15/RNA.lib for canonical nucleotides, modXNA for modified "
               "ones) under the names tleap expects")
    return out


def wrap(lines: List[str], width: int = 69) -> List[str]:
    """Keep REMARK records within the 80 columns of the PDB format."""
    out: List[str] = []
    for line in lines:
        indent = len(line) - len(line.lstrip())
        out += textwrap.wrap(line, width, subsequent_indent=" " * (indent + 4),
                             break_long_words=True, break_on_hyphens=False) or [""]
    return out


def structure_pdb(res: BuildResult, fiber_provenance: Optional[dict] = None,
                  residue_index: Optional[dict] = None) -> str:
    return format_pdb(res.atoms(), wrap(remarks(res, fiber_provenance, residue_index)))
