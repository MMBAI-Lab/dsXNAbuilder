"""
dsxnabuilder.amberio
====================

The downloadable package: a self-contained directory that rebuilds the Amber
topology outside the web app.

    structure.pdb   heavy atoms, Amber residue/atom names (tleap adds H)
    leap.in         tleap input (force field, modXNA libraries, solvation)
    libs/*.lib      modXNA residues used, renamed to their PDB residue names
    frcmod.modxna   modXNA parameters (only when modified residues are present)
    min.in          restrained sander minimization
    run.sh          tleap (+ sander) driver
    manifest.json   inputs, model provenance, residue map, issues, file hashes
    README.txt

The tleap sources mirror what modxna.sh uses to parametrize its residues
(leaprc.protein.ff19SB, leaprc.DNA.OL15, leaprc.RNA.OL3, leaprc.gaff,
frcmod.modxna), so the modXNA units see the same atom types they were built
with. The DNA force field can be switched among the leaprc files shipped with
AmberTools.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from . import __version__
from .build import BuildResult
from .structure import format_pdb

DNA_FF = {"OL15": "leaprc.DNA.OL15", "OL21": "leaprc.DNA.OL21",
          "OL24": "leaprc.DNA.OL24", "bsc1": "leaprc.DNA.bsc1"}
RNA_FF = {"OL3": "leaprc.RNA.OL3"}
WATER = {"none": None,
         "opc": ("leaprc.water.opc", "OPCBOX"),
         "tip3p": ("leaprc.water.tip3p", "TIP3PBOX"),
         "spce": ("leaprc.water.spce", "SPCBOX")}


@dataclass
class AmberOptions:
    dna_ff: str = "OL15"
    rna_ff: str = "OL3"
    water: str = "none"          # none (implicit GB) | opc | tip3p | spce
    box: str = "octahedral"      # octahedral | cubic
    buffer: float = 10.0
    neutralize: bool = True
    maxcyc: int = 2000
    ncyc: int = 500
    restraint_wt: float = 5.0
    igb: int = 8

    def validate(self) -> None:
        if self.dna_ff not in DNA_FF:
            raise ValueError(f"unknown DNA force field {self.dna_ff!r}")
        if self.rna_ff not in RNA_FF:
            raise ValueError(f"unknown RNA force field {self.rna_ff!r}")
        if self.water not in WATER:
            raise ValueError(f"unknown water model {self.water!r}")
        if self.box not in ("octahedral", "cubic"):
            raise ValueError(f"unknown box {self.box!r}")
        if not (5.0 <= self.buffer <= 30.0):
            raise ValueError("buffer must be within 5-30 A")


def _remarks(res: BuildResult) -> List[str]:
    m = res.model
    return [f"dsXNAbuilder {__version__}",
            f"fiber model #{m.id} {m.family}: {m.description}"[:200],
            f"helical twist {m.bp_twist} deg/bp, rise {m.bp_rise} A/bp, {m.handedness}-handed",
            "3DNA fiber data: Lu & Olson (2003) NAR 31:5108; Lu & Olson (2008) Nat Protoc 3:1213",
            "heavy atoms only; tleap adds hydrogens"]


def structure_pdb(res: BuildResult) -> str:
    return format_pdb(res.atoms(), _remarks(res))


def leap_script(res: BuildResult, opts: AmberOptions, libs: List[str]) -> str:
    has_mod = bool(libs)
    lines = ["# dsXNAbuilder: Amber topology", "# run: tleap -f leap.in"]
    if has_mod:
        lines.append("source leaprc.protein.ff19SB")
    lines += [f"source {DNA_FF[opts.dna_ff]}", f"source {RNA_FF[opts.rna_ff]}"]
    if has_mod:
        lines += ["source leaprc.gaff", "loadamberparams frcmod.modxna"]
        lines += [f"loadoff libs/{name}" for name in libs]
    water = WATER[opts.water]
    if water:
        lines.append(f"source {water[0]}")
    else:
        lines.append("set default PBradii mbondi3" if opts.igb == 8
                     else "set default PBradii mbondi2")
    lines += ["", "mol = loadpdb structure.pdb", "check mol", "charge mol"]
    if water:
        solv = "solvateOct" if opts.box == "octahedral" else "solvateBox"
        if opts.neutralize:
            lines += ["addIons mol Na+ 0", "addIons mol Cl- 0"]
        lines.append(f"{solv} mol {water[1]} {opts.buffer:.1f}")
    lines += ["saveamberparm mol system.parm7 system.rst7", "savepdb mol system.pdb",
              "quit", ""]
    return "\n".join(lines)


def restraint_mask(res: BuildResult) -> str:
    """Heavy atoms of the unmodified residues (Amber residue numbering)."""
    canon, k = [], 0
    for r in res.residues:
        k += 1
        if not r.flags.get("modxna"):
            canon.append(k)
    if not canon:
        return ""
    return f":{_ranges(canon)}&!@H="


def min_input(res: BuildResult, opts: AmberOptions) -> str:
    implicit = WATER[opts.water] is None
    body = ["Restrained minimization of the dsXNAbuilder duplex", " &cntrl",
            f"  imin = 1, maxcyc = {opts.maxcyc}, ncyc = {opts.ncyc},"]
    if implicit:
        body += [f"  ntb = 0, igb = {opts.igb}, cut = 999.0,"]
    else:
        body += ["  ntb = 1, cut = 9.0,"]
    mask = restraint_mask(res)
    if mask:
        body += ["  ntr = 1,", f"  restraint_wt = {opts.restraint_wt},",
                 f"  restraintmask = '{mask}',"]
    body += [" /", ""]
    return "\n".join(body)


RUN_SH = """#!/bin/bash
# dsXNAbuilder package: Amber topology (tleap) and restrained minimization (sander).
# Usage: ./run.sh            topology + minimization
#        ./run.sh --no-min   topology only
set -euo pipefail
cd "$(dirname "$0")"
if ! command -v tleap >/dev/null; then
  echo "tleap not found: source \\$AMBERHOME/amber.sh first" >&2; exit 1
fi
tleap -f leap.in > leap.log
grep -q "Errors = 0" leap.log || { echo "tleap reported errors, see leap.log" >&2; exit 1; }
echo "topology: system.parm7 system.rst7"
if [ "${1:-}" != "--no-min" ]; then
  sander -O -i min.in -p system.parm7 -c system.rst7 -ref system.rst7 \\
         -r min.rst7 -o min.out
  ambpdb -p system.parm7 -c min.rst7 > min.pdb 2>/dev/null || true
  echo "minimized: min.rst7 min.pdb"
fi
"""


def readme_text(res: BuildResult, opts: AmberOptions) -> str:
    m = res.model
    return f"""dsXNAbuilder {__version__}
===============================================================

Fiber model : #{m.id} {m.family} - {m.description}
Strand 1    : 5'-{res.duplex.top.raw.strip()}-3'
Strand 2    : 3'-{res.duplex.bottom.raw.strip()}-5'

Files
  structure.pdb   built duplex (heavy atoms; tleap adds hydrogens)
  leap.in         tleap input ({DNA_FF[opts.dna_ff]}, {RNA_FF[opts.rna_ff]}{', modXNA' if res.modified() else ''})
  min.in          restrained minimization (sander)
  run.sh          ./run.sh   -> system.parm7/rst7 and min.rst7
  manifest.json   inputs, provenance, residue names, warnings, file hashes

Requirements: AmberTools (tleap, sander). Run: source $AMBERHOME/amber.sh; ./run.sh

Please cite
  3DNA fiber models: Lu XJ, Olson WK (2003) Nucleic Acids Res 31:5108-5121;
                     Lu XJ, Olson WK (2008) Nat Protoc 3:1213-1227.
  modXNA (if modified residues are used): Love O et al. (2024)
                     J Chem Theory Comput 20:9354-9363.

Fiber templates derived from 3DNA data (CC BY-NC 4.0); modXNA libraries and
frcmod.modxna are GPL-3.0. Your structures are yours.
"""


def manifest(res: BuildResult, opts: AmberOptions, files: Dict[str, bytes],
             fiber_provenance: dict, extra: Optional[dict] = None) -> dict:
    m = res.model
    chains = {}
    for r in res.residues:
        chains.setdefault(r.chain, []).append({
            "resseq": r.resseq, "resname": r.resname, "token": r.spec.token,
            "triplet": list(r.spec.triplet()), "cap": r.spec.cap,
            "column": r.spec.column + 1, "modxna": bool(r.flags.get("modxna")),
            "base_swapped": bool(r.flags.get("swapped")),
            "o2_changed": bool(r.flags.get("chem_changed"))})
    return {
        "schema": "dsxnabuilder/manifest/1",
        "dsxnabuilder_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {"strand1_5to3": res.duplex.top.raw, "strand2_3to5": res.duplex.bottom.raw,
                  "model": m.id, "amber": asdict(opts)},
        "reproduce": (f"dsxnabuilder build --top {json.dumps(res.duplex.top.raw)} "
                      f"--bottom {json.dumps(res.duplex.bottom.raw)} --model {m.id} "
                      f"-o package.zip"),
        "model": {"id": m.id, "family": m.family, "description": m.description,
                  "repeat_bp": m.k, "twist_per_bp": m.bp_twist, "rise_per_bp": m.bp_rise,
                  "handedness": m.handedness, "references": m.references,
                  "orientation_flipped": res.flipped, "repeat_phase": res.phase},
        "fiber_data": fiber_provenance,
        "residue_name_map": res.name_map,
        "chains": chains,
        "issues": [i.to_json() for i in res.issues],
        "qc": {"o3p_min": round(min(d for *_, d in res.o3p), 3) if res.o3p else None,
               "o3p_max": round(max(d for *_, d in res.o3p), 3) if res.o3p else None,
               "clashes": res.n_clashes},
        "files": {name: {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
                  for name, data in sorted(files.items())},
        **(extra or {}),
    }


def package(res: BuildResult, opts: AmberOptions, fiber_provenance: dict,
            libs: Optional[Dict[str, str]] = None,
            frcmod: Optional[str] = None) -> Tuple[Dict[str, bytes], bytes]:
    """All files of the package and the zip bytes."""
    opts.validate()
    libs = libs or {}
    files: Dict[str, bytes] = {
        "structure.pdb": structure_pdb(res).encode(),
        "leap.in": leap_script(res, opts, sorted(libs)).encode(),
        "min.in": min_input(res, opts).encode(),
        "run.sh": RUN_SH.encode(),
        "README.txt": readme_text(res, opts).encode(),
    }
    for name, text in libs.items():
        files[f"libs/{name}"] = text.encode()
    if libs:
        if frcmod is None:
            raise ValueError("modified residues need frcmod.modxna")
        files["frcmod.modxna"] = frcmod.encode()
    files["manifest.json"] = json.dumps(
        manifest(res, opts, files, fiber_provenance), indent=1).encode()

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in sorted(files):
            info = zipfile.ZipInfo(f"dsxnabuilder/{name}", date_time=(2026, 1, 1, 0, 0, 0))
            info.external_attr = (0o755 if name.endswith(".sh") else 0o644) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, files[name])
    return files, buf.getvalue()


def _ranges(nums: List[int]) -> str:
    out, nums = [], sorted(set(nums))
    start = prev = nums[0]
    for x in nums[1:]:
        if x == prev + 1:
            prev = x
            continue
        out.append(f"{start}-{prev}" if start != prev else str(start))
        start = prev = x
    out.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(out)
