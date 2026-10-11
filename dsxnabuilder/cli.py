"""Command line: the same build as the web app, for scripting and reproduction.

    dsxnabuilder models
    dsxnabuilder build --top "5'-CGCGAATTCGCG-3'" --bottom "3'-GCGCTTAAGCGC-5'" \\
                       --model 4 -o package.zip [--pdb out.pdb] [--water opc]
"""

from __future__ import annotations

import argparse
import json
import sys

from .service import Session


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dsxnabuilder", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("models", help="list the buildable fiber models")
    b = sub.add_parser("build", help="build a duplex and write the package")
    b.add_argument("--top", required=True, help="strand 1, 5'->3'")
    b.add_argument("--bottom", required=True, help="strand 2, 3'->5', aligned under strand 1")
    b.add_argument("--model", type=int, default=4, help="3DNA fiber model id (default 4, B-DNA)")
    b.add_argument("-o", "--out", default="dsxnabuilder.zip")
    b.add_argument("--pdb", help="also write the structure PDB here")
    b.add_argument("--dna-ff", default="OL15")
    b.add_argument("--water", default="none")
    b.add_argument("--json", action="store_true", help="print issues as JSON")
    a = ap.parse_args(argv)

    s = Session(residue_library=_residue_library())
    if a.cmd == "models":
        info = s.models()
        for m in info["models"]:
            print(f"{m['id']:>3}  {m['family']:<26} {m['twist']:>7.2f} deg  "
                  f"{m['rise']:.3f} A  {m['handedness'][0].upper()}  k={m['k']:<2} "
                  f"{m['policy']:<8} {m['description'][:60]}")
        for m in info["excluded"]:
            print(f"{m['id']:>3}  (excluded) {'; '.join(m.get('issues') or [])[:90]}")
        return 0

    req = {"top": a.top, "bottom": a.bottom, "model": a.model,
           "amber": {"dna_ff": a.dna_ff, "water": a.water}}
    out = s.build(req)
    if a.json:
        print(json.dumps(out["issues"], indent=1))
    else:
        for i in out["issues"]:
            print(f"[{i['level']}] {i['code']}: {i['message']}", file=sys.stderr)
    if not out["ok"]:
        return 1
    with open(a.out, "wb") as fh:
        fh.write(out["zip"])
    if a.pdb:
        with open(a.pdb, "w", encoding="utf-8") as fh:
            fh.write(out["pdb"])
    sm = out["summary"]
    print(f"{sm['n_bp']} bp on model #{sm['model']} ({sm['family']}): "
          f"{sm['n_atoms']} heavy atoms, O3'-P {sm['o3p']}, clashes {sm['clashes']} -> {a.out}")
    return 0


def _residue_library():
    try:
        from .reslib import ResidueLibrary
    except ImportError:
        return None
    return ResidueLibrary.load_default()


if __name__ == "__main__":
    sys.exit(main())
