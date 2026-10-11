"""
dsxnabuilder.amberio
====================

Advanced, optional output: an Amber/modXNA topology package.

The default output of dsXNAbuilder is a PDB file (``pdbout``), usable with
any program. For Amber users this package adds what tleap needs to turn it
into a topology, including the modXNA parameters of the modified residues:

    structure.pdb   the same PDB as the default output
    leap.in         tleap input (force field, modXNA libraries, optional solvent)
    libs/*.lib      modXNA residues used, renamed to their PDB residue names
    frcmod.modxna   modXNA parameters (only when modified residues are present)
    run.sh          runs tleap -> system.parm7 / system.rst7
    manifest.json   inputs, model provenance, residue map, issues, file hashes
    README.txt

dsXNAbuilder is a structure builder: minimization, equilibration and
production are left to the user. The tleap sources mirror what modxna.sh uses
to parametrize its residues (leaprc.protein.ff19SB, leaprc.DNA.OL15,
leaprc.RNA.OL3, leaprc.gaff, frcmod.modxna), so the modXNA units see the same
atom types they were built with; the DNA force field can be switched among
the leaprc files shipped with AmberTools.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from . import __version__
from .build import BuildResult
from .pdbout import structure_pdb

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
    water: str = "none"          # none (vacuum / implicit solvent) | opc | tip3p | spce
    box: str = "octahedral"      # octahedral | cubic
    buffer: float = 10.0
    neutralize: bool = True

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
    lines += ["", "mol = loadpdb structure.pdb", "check mol", "charge mol"]
    if water:
        solv = "solvateOct" if opts.box == "octahedral" else "solvateBox"
        if opts.neutralize:
            lines += ["addIons mol Na+ 0", "addIons mol Cl- 0"]
        lines.append(f"{solv} mol {water[1]} {opts.buffer:.1f}")
    lines += ["saveamberparm mol system.parm7 system.rst7", "savepdb mol system.pdb",
              "quit", ""]
    return "\n".join(lines)


RUN_SH = """#!/bin/bash
# dsXNAbuilder package: Amber topology with tleap (AmberTools).
# Minimization, equilibration and production are up to you.
set -euo pipefail
cd "$(dirname "$0")"
if ! command -v tleap >/dev/null; then
  echo "tleap not found: source \\$AMBERHOME/amber.sh first" >&2; exit 1
fi
tleap -f leap.in > leap.log
grep -q "Errors = 0" leap.log || { echo "tleap reported errors, see leap.log" >&2; exit 1; }
echo "topology: system.parm7 system.rst7 (system.pdb)"
"""


def readme_text(res: BuildResult, opts: AmberOptions) -> str:
    m = res.model
    mod = " + modXNA (frcmod.modxna, libs/)" if res.modified() else ""
    water = WATER[opts.water]
    solv = (f"{opts.water.upper()} water, {opts.box} box, {opts.buffer:.0f} A buffer"
            + (", neutralized with Na+/Cl-" if opts.neutralize else "")) if water else "no solvent"
    return f"""dsXNAbuilder {__version__} - Amber/modXNA topology package
===============================================================

Fiber model : #{m.id} {m.family} - {m.description}
Strand 1    : 5'-{res.duplex.top.raw.strip()}-3'
Strand 2    : 3'-{res.duplex.bottom.raw.strip()}-5'
Force field : {DNA_FF[opts.dna_ff]}, {RNA_FF[opts.rna_ff]}{mod}
Solvent     : {solv}

Run (AmberTools): source $AMBERHOME/amber.sh; ./run.sh   ->  system.parm7 / system.rst7

dsXNAbuilder builds the structure and this topology; minimization, equilibration
and production are up to you. Modified nucleotides keep the sugar pucker of the
fiber model; see manifest.json for any warning before simulating.

Please cite
  3DNA fiber models: Lu XJ, Olson WK (2003) Nucleic Acids Res 31:5108-5121;
                     Lu XJ, Olson WK (2008) Nat Protoc 3:1213-1227.
  modXNA (if modified residues are used): Love O et al. (2024)
                     J Chem Theory Comput 20:9354-9363.

Fiber templates derived from 3DNA data (CC BY-NC 4.0); modXNA libraries and
frcmod.modxna are GPL-3.0. Your structures are yours.
"""


def manifest(res: BuildResult, opts: AmberOptions, files: Dict[str, bytes],
             fiber_provenance: dict, residue_index: Optional[dict] = None) -> dict:
    m = res.model
    chains = {}
    for r in res.residues:
        chains.setdefault(r.chain, []).append({
            "resseq": r.resseq, "resname": r.resname, "token": r.spec.token,
            "triplet": list(r.spec.triplet()), "cap": r.spec.cap,
            "column": r.spec.column + 1, "modxna": bool(r.flags.get("modxna")),
            "base_swapped": bool(r.flags.get("swapped")),
            "o2_changed": bool(r.flags.get("chem_changed"))})
    out = {
        "schema": "dsxnabuilder/manifest/2",
        "dsxnabuilder_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {"strand1_5to3": res.duplex.top.raw, "strand2_3to5": res.duplex.bottom.raw,
                  "model": m.id, "amber": asdict(opts)},
        "reproduce": (f"dsxnabuilder build --top {json.dumps(res.duplex.top.raw)} "
                      f"--bottom {json.dumps(res.duplex.bottom.raw)} --model {m.id} "
                      f"--amber package.zip"),
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
    }
    if residue_index and res.modified():
        out["modxna"] = {k: residue_index.get(k) for k in
                         ("modxna_version", "modxna_commit", "cpptraj", "cpptraj_commit",
                          "frcmod_sha256")}
    return out


def package(res: BuildResult, opts: AmberOptions, fiber_provenance: dict,
            libs: Optional[Dict[str, str]] = None, frcmod: Optional[str] = None,
            residue_index: Optional[dict] = None, pdb: Optional[str] = None
            ) -> Tuple[Dict[str, bytes], bytes]:
    """All files of the Amber package and the zip bytes."""
    opts.validate()
    libs = libs or {}
    files: Dict[str, bytes] = {
        "structure.pdb": (pdb or structure_pdb(res, fiber_provenance, residue_index)).encode(),
        "leap.in": leap_script(res, opts, sorted(libs)).encode(),
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
        manifest(res, opts, files, fiber_provenance, residue_index), indent=1).encode()

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in sorted(files):
            info = zipfile.ZipInfo(f"dsxnabuilder/{name}", date_time=(2026, 1, 1, 0, 0, 0))
            info.external_attr = (0o755 if name.endswith(".sh") else 0o644) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, files[name])
    return files, buf.getvalue()
