"""Command line: the same build as the web app, for scripting and reproduction.

    dsxnabuilder models
    dsxnabuilder build --top "5'-CGCGAATTCGCG-3'" --bottom "3'-GCGCTTAAGCGC-5'" \\
                       --model 4 -o duplex.pdb [--amber package.zip --water opc]
    dsxnabuilder build --high-throughput duplexes.txt -o structures.zip
"""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .amberio import BOXES, BUFFER_RANGE, DNA_FF, RNA_FF, WATER, AmberOptions
from .batch import (CLI_LIMIT, DEFAULT_MODEL, Record, batch_files, build_records,
                    parse_batch, write_directory, zip_bytes)
from .reslib import ResidueLibrary
from .service import Session

_DEFAULTS = AmberOptions()


def main(argv=None) -> int:
    """Run ``dsxnabuilder models`` or ``dsxnabuilder build``; returns the exit status."""
    ap = _parser()
    a = ap.parse_args(argv)
    if a.cmd == "build":
        _check_build_args(ap, a)
    s = Session(residue_library=ResidueLibrary.load_default())
    if a.cmd == "models":
        _print_models(s)
        return 0
    if a.high_throughput:
        return _build_batch(s, a)
    return _build(s, a)


def _check_build_args(ap: argparse.ArgumentParser, a: argparse.Namespace) -> None:
    """Options that cannot go together end in a usage error (exit status 2)."""
    if a.high_throughput:
        given = [f for f, v in (("--top", a.top), ("--bottom", a.bottom), ("--amber", a.amber),
                                ("--name", a.name)) if v]
        if given:
            ap.error(f"--high-throughput builds PDB files only; drop {', '.join(given)}")
        if a.out and a.out.lower().endswith(".pdb"):
            ap.error("with --high-throughput, -o is a .zip file or a directory")
        return
    if not a.top or not a.bottom:
        ap.error("build needs --top and --bottom, or --high-throughput FILE")
    bad = [n for n in a.name if "=" not in n]
    if bad:
        ap.error(f"--name expects KEY=NAME, got {', '.join(bad)}")


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="dsxnabuilder", description=__doc__.split("\n\n")[0])
    ap.add_argument("--version", action="version", version=f"dsxnabuilder {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("models", help="list the buildable fiber models")
    b = sub.add_parser("build", help="build a duplex and write its PDB")
    b.add_argument("--top", help="strand 1, 5'->3'")
    b.add_argument("--bottom", help="strand 2, 3'->5', aligned under strand 1")
    b.add_argument("--high-throughput", metavar="FILE",
                   help=f"build every duplex of a FASTA-like file, up to {CLI_LIMIT} ('-' reads "
                        "standard input): one '>name=...; model=...' header line, then strand "
                        "1 and strand 2, per duplex")
    b.add_argument("--model", type=int, default=DEFAULT_MODEL,
                   help="3DNA fiber model id (default 4, B-DNA); with --high-throughput, the "
                        "model of the duplexes whose header has no model=")
    b.add_argument("-o", "--out",
                   help="output PDB (default dsxnabuilder.pdb); with --high-throughput a .zip "
                        "file or a directory (default dsxnabuilder_batch.zip)")
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
    a.out = a.out or "dsxnabuilder.pdb"
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


def _build_batch(s: Session, a: argparse.Namespace) -> int:
    """--high-throughput: every duplex of the file, PDB files into a zip or a directory."""
    try:
        if a.high_throughput == "-":
            text = sys.stdin.read()
        else:
            with open(a.high_throughput, encoding="utf-8") as fh:
                text = fh.read()
    except OSError as exc:
        print(f"cannot read {a.high_throughput}: {exc.strerror}", file=sys.stderr)
        return 1
    records, issues = parse_batch(text, CLI_LIMIT, a.model)
    if issues:
        for i in issues:
            print(f"[error] {i.code}: {i.message}", file=sys.stderr)
        return 1
    width = len(str(len(records)))

    def report(i: int, r: Record, out: dict) -> None:
        if out["ok"]:
            sm = out["summary"]
            line = f"built, {sm['n_bp']} bp on model #{sm['model']}, {sm['n_atoms']} atoms"
        else:
            err = next(x for x in out["issues"] if x["level"] == "error")
            line = f"failed: {err['code']}: {err['message']}"
        print(f"[{i + 1:>{width}}/{len(records)}] {r.name}: {line}", file=sys.stderr)

    outcomes = build_records(s, records, report)
    files = batch_files(records, outcomes, text, a.model)
    out_path = a.out or "dsxnabuilder_batch.zip"
    if out_path.lower().endswith(".zip"):
        with open(out_path, "wb") as fh:
            fh.write(zip_bytes(files))
    else:
        write_directory(files, out_path)
    if a.json:
        print(json.dumps([{"name": r.name, "model": r.model, "ok": o["ok"],
                           "issues": o["issues"]} for r, o in zip(records, outcomes)], indent=1))
    failed = sum(1 for o in outcomes if not o["ok"])
    print(f"{len(records) - failed} of {len(records)} duplexes built -> {out_path}"
          + (f"; {failed} failed (see report.tsv)" if failed else ""),
          file=sys.stderr if a.json else sys.stdout)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
