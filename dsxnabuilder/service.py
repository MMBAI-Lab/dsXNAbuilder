"""
dsxnabuilder.service
====================

Single entry point shared by the browser worker and the CLI. Everything is
JSON in, JSON (plus file bytes) out, and loading happens once per session.

Browser flow (two calls, so that Python stays synchronous):

1. ``validate(request)`` parses both strands, reports issues and the list of
   residue-library keys the build needs;
2. the page fetches those library files and calls ``build(request, libs)``.

A batch (``batch``) runs the same way: ``batch_validate`` checks the text and
every record, the page builds the records one by one with ``build`` (so it
can show progress and stop) and ``batch_zip`` packs the results.
"""

from __future__ import annotations

import os
from dataclasses import fields
from typing import Dict, List, Optional

from . import __version__
from .amberio import AmberOptions, package
from .batch import WEB_LIMIT, batch_files, parse_batch, zip_bytes
from .build import BuildError, BuildResult, build_duplex
from .catalog import Catalog, default_chem
from .fibers import FiberLibrary, load_library
from .fragments import pieces
from .issues import Issue
from .nucleotide import pair_kind
from .paths import DATA
from .pdbout import describe, structure_pdb
from .seqlang import AliasTable, DuplexSpec, SequenceParseError, parse_duplex
from .structure import element_of


class Session:
    """Data loaded once (catalogue, aliases, fiber models) and the requests.

    ``residue_library`` is a reslib.ResidueLibrary, or None to build only
    canonical duplexes.
    """

    def __init__(self, data_dir: Optional[str] = None, residue_library=None):
        d = data_dir or DATA
        self.data_dir = d
        self.catalog = Catalog.from_json(os.path.join(d, "modxna_catalog.json"))
        self.aliases = AliasTable.load(os.path.join(d, "aliases.json"))
        self.fibers: FiberLibrary = load_library(os.path.join(d, "fiber_models.json"))
        self.residue_library = residue_library

    # ------------------------------------------------------------------
    def models(self) -> dict:
        """Model list for the UI, grouped by family, buildable ones first."""
        out = []
        for m in sorted(self.fibers.models.values(), key=lambda m: m.id):
            out.append({
                "id": m.id, "family": m.family, "description": m.description,
                "k": m.k, "twist": m.bp_twist, "rise": m.bp_rise,
                "handedness": m.handedness, "policy": m.sequence_policy,
                "kind": m.kind,
                "strands": [m.strand_chem("I"), m.strand_chem("II")],
                "native": ["".join(s.I.native for s in m.slots),
                           "".join(s.II.native for s in m.slots)],
                "pucker": [(m.slots[0].I.pucker or {}).get("domain"),
                           (m.slots[0].II.pucker or {}).get("domain")],
                "syn": any(s.I.syn or s.II.syn for s in m.slots),
                "pairing": sorted({s.pairing for s in m.slots}),
                "references": m.references})
        excluded = [dict(v, id=k) for k, v in sorted(self.fibers.excluded.items())]
        return {"version": __version__, "models": out, "excluded": excluded,
                "references": self.fibers.references,
                "provenance": self.fibers.provenance}

    def pieces(self) -> dict:
        """Sugars, bases and linkages of the prebuilt library (modXNA library tab)."""
        return pieces(self.data_dir)

    def validate(self, request: dict) -> dict:
        """Parse and check without building coordinates."""
        try:
            model = self.fibers.get(int(request.get("model", 4)))
        except (KeyError, ValueError) as exc:
            return {"ok": False, "issues": [Issue("error", None, "UNKNOWN_MODEL",
                                                  str(exc), {}).to_json()]}
        defaults = (default_chem(model.strand_chem("I")), default_chem(model.strand_chem("II")))
        try:
            dup = parse_duplex(request.get("top", ""), request.get("bottom", ""),
                               self.catalog, self.aliases, defaults)
        except SequenceParseError as exc:
            return {"ok": False, "issues": [exc.issue.to_json()]}
        try:
            res = self._build(request, dry_run=True)
        except BuildError as exc:
            return {"ok": False, "issues": [i.to_json() for i in exc.issues],
                    "columns": _columns(dup)}
        return {"ok": True, "issues": [i.to_json() for i in res.issues],
                "columns": _columns(dup), "needed_libs": self._needed(res),
                "n_bp": dup.n_columns}

    def build(self, request: dict, libs: Optional[Dict[str, str]] = None) -> dict:
        """Build the duplex. Always returns the ``pdb`` text (the default
        output); with ``request["package"]`` true also the Amber/modXNA
        topology package as ``zip`` bytes (advanced option)."""
        opts = None
        if request.get("package"):
            try:
                opts = amber_options(request.get("amber", {}))
            except (TypeError, ValueError) as exc:
                return {"ok": False, "issues": [Issue("error", None, "AMBER_OPTIONS",
                                                      f"Amber package: {exc}", {}).to_json()]}
        if self.residue_library is not None:
            self.residue_library.provide(libs, request.get("frcmod"))
        try:
            res = self._build(request)
        except BuildError as exc:
            return {"ok": False, "issues": [i.to_json() for i in exc.issues]}
        index = self.residue_library.index if self.residue_library is not None else None
        pdb = structure_pdb(res, self.fibers.provenance, index)
        out = {"ok": True, "issues": [i.to_json() for i in res.issues], "pdb": pdb,
               "summary": summary(res)}
        if opts is not None and any(r.flags.get("abasic") for r in res.residues):
            out["issues"].append(Issue(
                "error", None, "ABASIC_NO_PARAMETERS",
                "the Amber/modXNA package cannot include abasic sites: modXNA has no abasic "
                "fragment yet. The PDB is still available.", {}).to_json())
        elif opts is not None:
            lib_texts, frcmod = None, None
            if self.residue_library is not None and res.modified():
                lib_texts, frcmod = self.residue_library.package_libs(res, libs)
            _, out["zip"] = package(res, opts, self.fibers.provenance, lib_texts, frcmod,
                                    index, pdb)
        return out

    def batch_validate(self, request: dict, limit: int = WEB_LIMIT) -> dict:
        """Check a batch text (``request["text"]``, see dsxnabuilder.batch) and
        validate each record. ``ok`` when the format is right and at least one
        record can be built; the others are reported and skipped."""
        records, issues = parse_batch(request.get("text", ""), limit)
        if any(i.code == "BATCH_LIMIT" for i in issues):
            records = []                    # not worth validating
        rows = []
        for r in records:
            v = self.validate(r.request())
            rows.append({"name": r.name, "model": r.model, "line": r.line, "top": r.top,
                         "bottom": r.bottom, "ok": v["ok"], "n_bp": v.get("n_bp"),
                         "issues": v["issues"], "needed_libs": v.get("needed_libs", [])})
        return {"ok": not issues and any(row["ok"] for row in rows),
                "issues": [i.to_json() for i in issues], "records": rows, "limit": limit}

    def batch_zip(self, request: dict, limit: int = WEB_LIMIT) -> dict:
        """The zip of a batch built record by record: ``request["text"]`` is the
        batch and ``request["results"]`` maps record names to build outputs."""
        records, issues = parse_batch(request.get("text", ""), limit)
        if issues:
            raise ValueError("; ".join(i.message for i in issues))
        results = request.get("results") or {}
        missing = {"ok": False, "issues": [Issue("error", None, "NOT_BUILT",
                                                 "not built (stopped)", {}).to_json()]}
        outcomes = [results.get(r.name) or missing for r in records]
        built = sum(1 for o in outcomes if o.get("ok"))
        return {"zip": zip_bytes(batch_files(records, outcomes, request.get("text", ""))),
                "built": built, "failed": len(records) - built}

    # ------------------------------------------------------------------
    def _build(self, request: dict, dry_run: bool = False) -> BuildResult:
        """build_duplex with the request's strands, model and residue names."""
        return build_duplex(request.get("top", ""), request.get("bottom", ""),
                            int(request.get("model", 4)), self.fibers, self.catalog,
                            self.aliases, self.residue_library, dry_run=dry_run,
                            names=request.get("names"))

    def _needed(self, res: BuildResult) -> List[str]:
        """Residue-library keys of the modified residues (files the page fetches)."""
        if self.residue_library is None:
            return []
        return sorted({r.flags.get("key") or self.residue_library.key_for(r.spec)
                       for r in res.modified()})


def _columns(dup: DuplexSpec) -> List[dict]:
    """The alignment, column by column, for the page's pair table."""
    out = []
    for i, (a, b) in enumerate(dup.pairs()):
        kind = None
        if a is not None and b is not None:
            kind = pair_kind(a.letter, b.letter)
        out.append({"i": i, "top": a.token if a else None, "bottom": b.token if b else None,
                    "top_mod": bool(a and not a.is_canonical),
                    "bottom_mod": bool(b and not b.is_canonical), "pair": kind})
    return out


def amber_options(d: dict) -> AmberOptions:
    """AmberOptions from the request (unknown keys ignored), validated."""
    names = {f.name for f in fields(AmberOptions)}
    opts = AmberOptions(**{k: v for k, v in (d or {}).items() if k in names})
    opts.validate()
    return opts


def summary(res: BuildResult) -> dict:
    """Figures of a build shown by the page and printed by the CLI."""
    return {"model": res.model.id, "family": res.model.family,
            "n_bp": res.duplex.n_columns,
            "n_residues": len(res.residues), "n_atoms": sum(len(r.names) for r in res.residues),
            "n_heavy": sum(element_of(n) != "H" for r in res.residues for n in r.names),
            "n_modified": len(res.modified()), "flipped": res.flipped, "phase": res.phase,
            "o3p": [round(min(d for *_, d in res.o3p), 3), round(max(d for *_, d in res.o3p), 3)]
            if res.o3p else None,
            "clashes": res.n_clashes, "name_map": res.name_map,
            # the residues that carry a generated name, for the renaming table
            "named": _named(res),
            # per residue, for the viewer's base and modification colouring:
            # [chain, resseq, pairing letter, canonical | modified | abasic]
            "residues": [[r.chain, r.resseq, "-" if r.flags.get("abasic") else r.spec.letter,
                          "abasic" if r.flags.get("abasic") else
                          "canonical" if r.spec.is_canonical else "modified"]
                         for r in res.residues]}


def _named(res: BuildResult) -> List[dict]:
    """Residues that carry a generated name, once per residue key."""
    out, seen = [], set()
    for r in res.residues:
        key = r.spec.key()
        if key not in res.name_map or key in seen:
            continue
        seen.add(key)
        abasic = bool(r.flags.get("abasic"))
        out.append({"key": key, "name": res.name_map[key], "triplet": list(r.spec.triplet()),
                    "cap": r.spec.cap, "custom": key in res.custom_names,
                    "description": "abasic site" if abasic
                    else describe(r.spec.triplet(), r.spec.cap)})
    return out
