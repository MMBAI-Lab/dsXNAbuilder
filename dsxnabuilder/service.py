"""
dsxnabuilder.service
====================

Single entry point shared by the browser worker and the CLI. Everything is
JSON in, JSON (plus file bytes) out, and loading happens once per session.

Browser flow (two calls, so that Python stays synchronous):

1. ``validate(request)`` parses both strands, reports issues and the list of
   residue-library keys the build needs;
2. the page fetches those library files and calls ``build(request, libs)``.
"""

from __future__ import annotations

import os
from dataclasses import fields
from typing import Dict, List, Optional

from . import __version__
from .amberio import AmberOptions, package
from .build import BuildError, BuildResult, build_duplex
from .catalog import Catalog, Issue
from .fibers import FiberLibrary, load_library
from .seqlang import AliasTable, SequenceParseError, parse_duplex

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


class Session:
    def __init__(self, data_dir: Optional[str] = None, residue_library=None):
        d = data_dir or DATA
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

    def validate(self, request: dict) -> dict:
        """Parse and check without building coordinates."""
        try:
            model = self.fibers.get(int(request.get("model", 4)))
        except (KeyError, ValueError) as exc:
            return {"ok": False, "issues": [Issue("error", None, "UNKNOWN_MODEL",
                                                  str(exc), {}).to_json()]}
        defaults = (_chem(model.strand_chem("I")), _chem(model.strand_chem("II")))
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
        """Build and package. Returns dict with ``pdb`` text and ``zip`` bytes."""
        if self.residue_library is not None:
            self.residue_library.provide(libs, request.get("frcmod"))
        try:
            res = self._build(request)
        except BuildError as exc:
            return {"ok": False, "issues": [i.to_json() for i in exc.issues]}
        opts = amber_options(request.get("amber", {}))
        lib_texts, frcmod = None, None
        if self.residue_library is not None and res.modified():
            lib_texts, frcmod = self.residue_library.package_libs(res, libs)
        files, zbytes = package(res, opts, self.fibers.provenance, lib_texts, frcmod)
        return {"ok": True, "issues": [i.to_json() for i in res.issues],
                "pdb": files["structure.pdb"].decode(), "zip": zbytes,
                "summary": summary(res)}

    # ------------------------------------------------------------------
    def _build(self, request: dict, dry_run: bool = False) -> BuildResult:
        return build_duplex(request.get("top", ""), request.get("bottom", ""),
                            int(request.get("model", 4)), self.fibers, self.catalog,
                            self.aliases, self.residue_library, dry_run=dry_run)

    def _needed(self, res: BuildResult) -> List[str]:
        if self.residue_library is None:
            return []
        return sorted({r.flags.get("key") or self.residue_library.key_for(r.spec)
                       for r in res.modified()})


def _chem(c: str) -> str:
    return c if c in ("DNA", "RNA") else "DNA"


def _columns(dup) -> List[dict]:
    out = []
    from .nucleotide import pair_kind
    for i, (a, b) in enumerate(dup.pairs()):
        kind = None
        if a is not None and b is not None:
            kind = pair_kind(a.letter, b.letter)
        out.append({"i": i, "top": a.token if a else None, "bottom": b.token if b else None,
                    "top_mod": bool(a and not a.is_canonical),
                    "bottom_mod": bool(b and not b.is_canonical), "pair": kind})
    return out


def amber_options(d: dict) -> AmberOptions:
    names = {f.name for f in fields(AmberOptions)}
    opts = AmberOptions(**{k: v for k, v in (d or {}).items() if k in names})
    opts.validate()
    return opts


def summary(res: BuildResult) -> dict:
    return {"model": res.model.id, "family": res.model.family,
            "n_bp": res.duplex.n_columns,
            "n_residues": len(res.residues), "n_atoms": sum(len(r.names) for r in res.residues),
            "n_modified": len(res.modified()), "flipped": res.flipped, "phase": res.phase,
            "o3p": [round(min(d for *_, d in res.o3p), 3), round(max(d for *_, d in res.o3p), 3)]
            if res.o3p else None,
            "clashes": res.n_clashes, "name_map": res.name_map}
