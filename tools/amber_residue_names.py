#!/usr/bin/env python3
"""Extract the canonical Amber nucleic-acid units (names, charges, geometry).

    python tools/amber_residue_names.py [--amberhome /opt/amber26]

Writes data/amber_residues.json. Per unit: the heavy-atom names and net charge
(tests and the build preflight check that every atom written to the PDB
exists in the unit tleap will match, without needing AmberTools at test time)
and all atoms with elements, coordinates and bonds, from which
dsxnabuilder.hydrogens places the hydrogens of canonical nucleotides under the
names tleap expects. The source libraries are recorded with their sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from dsxnabuilder.offlib import parse_off  # noqa: E402

LIBS = {"DNA": "DNA.OL15.lib", "RNA": "RNA.lib"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--amberhome", default=os.environ.get("AMBERHOME", "/opt/amber26"))
    a = ap.parse_args(argv)
    out = {"source": {}, "units": {}}
    for chem, lib in LIBS.items():
        path = os.path.join(a.amberhome, "dat", "leap", "lib", lib)
        with open(path, "rb") as fh:
            raw = fh.read()
        units = parse_off(raw.decode("utf-8"))
        out["source"][chem] = {"file": lib, "sha256": hashlib.sha256(raw).hexdigest()}
        for name, u in sorted(units.items()):
            out["units"][name] = {
                "heavy": u.heavy_atoms(), "charge": round(u.net_charge, 4),
                "atoms": u.atom_names, "elements": u.elements,
                "xyz": [[round(float(x), 4) for x in p] for p in u.xyz],
                "bonds": [list(b) for b in u.bonds]}
    dest = os.path.join(ROOT, "data", "amber_residues.json")
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"{len(out['units'])} units -> {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
