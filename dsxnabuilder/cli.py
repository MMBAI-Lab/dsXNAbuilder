"""Command line: the same build as the web app, for scripting and reproduction.

    dsxnabuilder models
    dsxnabuilder build --top "5'-CGCGAATTCGCG-3'" --bottom "3'-GCGCTTAAGCGC-5'" \\
                       --model 4 -o duplex.pdb [--amber package.zip --water opc]
"""

from __future__ import annotations

import argparse
import json
import sys

from .amberio import BOXES, BUFFER_RANGE, DNA_FF, RNA_FF, WATER, AmberOptions
from .reslib import ResidueLibrary
from .service import Session

_DEFAULTS = AmberOptions()


def main(argv=None) -> int:
    """Run ``dsxnabuilder models`` or ``dsxnabuilder build``; returns the exit status."""
    ap = _parser()
    a = ap.parse_args(argv)
    s = Session(residue_library=ResidueLibrary.load_default())
    if a.cmd == "models":
        _print_models(s)
        return 0
    bad = [n for n in a.name if "=" not in n]
    if bad:
        ap.error(f"--name expects KEY=NAME, got {', '.join(bad)}")
    return _build(s, a)


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="dsxnabuilder", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("models", help="list the buildable fiber models")
    b = sub.add_parser("build", help="build a duplex and write its PDB")
    b.add_argument("--top", required=True, help="strand 1, 5'->3'")
    b.add_argument("--bottom", required=True, help="strand 2, 3'->5', aligned under strand 1")
    b.add_argument("--model", type=int, default=4, help="3DNA fiber model id (default 4, B-DNA)")
    b.add_argument("-o", "--out", default="dsxnabuilder.pdb", help="output PDB (default %(default)s)")
    b.add_argument("--name", action="append", default=[], metavar="KEY=NAME",
                   help="residue name for a modXNA residue, e.g. PS1/OME/RAA/int=OMA "
                        "(repeatable; the others are named X01, X02, ...)")
    b.add_argument("--json", action="store_true", help="print issues as JSON")
    p = b.add_argument_group("Amber/modXNA package (advanced)")
    p.add_argument("--amber", metavar="ZIP", help="also write the topology package here")
    p.add_argument("--dna-ff", choices=sorted(DNA_FF), default=_DEFAULTS.dna_ff,
                   help="DNA force field (default %(default)s)")
    p.add_argument("--rna-ff", choices=sorted(RNA_FF), default=_DEFAULTS.rna_ff,
                   help="RNA force field (default %(default)s)")
    p.add_argument("--water", choices=list(WATER), default=_DEFAULTS.water,
                   help="explicit water model (default %(default)s)")
    p.add_argument("--box", choices=list(BOXES), default=_DEFAULTS.box,
                   help="solvent box (default %(default)s)")
    p.add_argument("--buffer", type=float, default=_DEFAULTS.buffer,
                   help="solute to box edge, %g-%g A (default %%(default)s)" % BUFFER_RANGE)
    p.add_argument("--no-neutralize", dest="neutralize", action="store_false",
                   help="do not add Na+/Cl- counterions to the solvated system")
    return ap


def _print_models(s: Session) -> None:
    info = s.models()
    for m in info["models"]:
        print(f"{m['id']:>3}  {m['family']:<26} {m['twist']:>7.2f} deg  "
              f"{m['rise']:.3f} A  {m['handedness'][0].upper()}  k={m['k']:<2} "
              f"{m['policy']:<8} {m['description'][:60]}")
    for m in info["excluded"]:
        print(f"{m['id']:>3}  (excluded) {'; '.join(m.get('issues') or [])[:90]}")


def _build(s: Session, a: argparse.Namespace) -> int:
    req = {"top": a.top, "bottom": a.bottom, "model": a.model, "package": bool(a.amber),
           "names": dict(n.rsplit("=", 1) for n in a.name),
           "amber": {"dna_ff": a.dna_ff, "rna_ff": a.rna_ff, "water": a.water, "box": a.box,
                     "buffer": a.buffer, "neutralize": a.neutralize}}
    out = s.build(req)
    if a.json:
        print(json.dumps(out["issues"], indent=1))
    else:
        for i in out["issues"]:
            print(f"[{i['level']}] {i['code']}: {i['message']}", file=sys.stderr)
    if not out["ok"]:
        return 1
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(out["pdb"])
    if a.amber:
        if "zip" not in out:
            print(f"no Amber package written to {a.amber} (see the error above)", file=sys.stderr)
            return 1
        with open(a.amber, "wb") as fh:
            fh.write(out["zip"])
    sm = out["summary"]
    print(f"{sm['n_bp']} bp on model #{sm['model']} ({sm['family']}): "
          f"{sm['n_atoms']} atoms ({sm['n_heavy']} heavy), O3'-P {sm['o3p']}, "
          f"clashes {sm['clashes']} -> {a.out}"
          + (f" + {a.amber}" if a.amber else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
