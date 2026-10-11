"""
dsxnabuilder.pdbout
===================

The default output: one self-describing PDB file.

Heavy atoms only, with Amber-style residue and atom names (DA/DC/DG/DT,
A/C/G/U; OP1/OP2; thymine C7), chain A = strand 1 and chain B = strand 2,
both 5'->3'. Modified nucleotides keep the modXNA atom names and get a
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


def describe(triplet, cap: str) -> str:
    bb, sugar, base = triplet
    return (f"{BASE_NAMES.get(base, base)}, {SUGAR_NAMES.get(sugar, sugar)}, "
            f"{BACKBONE_NAMES.get(bb, bb)}; {CAP_NAMES.get(cap, cap)}")


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
    for r in res.residues:
        if r.flags.get("modxna") and r.resname not in mods:
            mods[r.resname] = r
    if mods:
        out.append("modified nucleotides (modXNA; Love et al. (2024) JCTC 20:9354):")
        for name, r in sorted(mods.items()):
            out.append(f"  {name} = {'/'.join(r.spec.triplet())} {r.spec.cap}: "
                       f"{describe(r.spec.triplet(), r.spec.cap)}")
        if residue_index:
            out.append(f"modXNA {residue_index.get('modxna_version', '?')} "
                       f"@{str(residue_index.get('modxna_commit', '?'))[:7]}")
    out.append("heavy atoms only (add hydrogens with your modeling or MD program)")
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
