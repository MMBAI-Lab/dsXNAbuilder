#!/usr/bin/env python3
"""Extract heavy-atom names of the canonical Amber nucleic-acid units.

    python tools/amber_residue_names.py [--amberhome /opt/amber26]

Writes data/amber_residues.json, used by tests (and the build preflight) to
check that every atom written to structure.pdb exists in the unit tleap will
match, without needing AmberTools at test time.
"""

from __future__ import annotations

import argparse
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
        with open(path, encoding="utf-8") as fh:
            units = parse_off(fh.read())
        out["source"][chem] = lib
        for name, u in sorted(units.items()):
            out["units"][name] = {"heavy": u.heavy_atoms(), "charge": round(u.net_charge, 4)}
    dest = os.path.join(ROOT, "data", "amber_residues.json")
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"{len(out['units'])} units -> {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
