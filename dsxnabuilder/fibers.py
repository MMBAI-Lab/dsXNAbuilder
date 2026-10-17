# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.fibers
===================

Fiber models (3DNA, Lu & Olson 2003 NAR 31:5108) as screw-symmetric repeats.

A model is a repeating unit of ``k`` base pairs ("slots") and the screw
operator S(Omega, h) that maps one repeat onto the next one along strand I
5'->3'. All models are stored normalized so that strand I runs 5'->3' along +z
(``tools/derive_fiber_models.py`` flips the ones that 3DNA writes the other way
round). For every slot and strand we keep

* the base standard reference frame in helix coordinates, and
* one or more nucleotide templates expressed in that frame (one per base
  letter for 3DNA's per-base "generic" models, the model's own residue
  otherwise).

Base pair ``i`` uses slot ``i mod k`` transformed by S^(i div k). Frames are
placed directly by the screw, never through step parameters, so a built
duplex reproduces the fiber model exactly and the backbone (O3'-P) stays as
continuous as in the source repeat.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .geometry import Frame, rot_z
from .nucleotide import BACKBONE_ATOMS, base_swap, to_dna, to_rna
from .paths import DATA


@dataclass
class NucleotideTemplate:
    """Atom names and coordinates in the standard reference frame of the base."""

    names: List[str]
    xyz: np.ndarray

    @classmethod
    def from_json(cls, d: dict) -> "NucleotideTemplate":
        return cls(list(d["names"]), np.array(d["xyz"], float))


@dataclass
class SlotStrand:
    """One strand of one slot of the repeat."""

    frame: Frame
    native: str              # letter of the model's own base (A/C/G/T/U/I/X)
    chem: str                # DNA | RNA
    templates: Dict[str, NucleotideTemplate]
    syn: bool = False
    pucker: Optional[dict] = None


@dataclass
class Slot:
    I: SlotStrand
    II: SlotStrand
    pairing: str = "wc"      # wc | hoogsteen | other


@dataclass
class FiberModel:
    """A fiber model: k slots and the screw that maps one repeat onto the next."""

    id: int
    name: str
    family: str
    description: str
    k: int
    twist: float             # degrees per repeat (normalized)
    rise: float              # angstrom per repeat (normalized, > 0)
    slots: List[Slot]
    kind: str = "repeat"     # per-base | repeat
    sequence_policy: str = "any"   # any | swap | pattern
    handedness: str = "right"
    bp_twist: float = 0.0
    bp_rise: float = 0.0
    references: List[str] = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def strand_chem(self, strand: str) -> str:
        """DNA or RNA for model strand "I" or "II"; "mixed" if its slots differ."""
        chems = {getattr(s, strand).chem for s in self.slots}
        return chems.pop() if len(chems) == 1 else "mixed"

    def screw(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        """Rotation and translation of the screw applied n times."""
        return rot_z(n * self.twist), np.array([0.0, 0.0, n * self.rise])

    def slot_strand(self, i: int, strand: str) -> Tuple[SlotStrand, Frame]:
        """Slot data and base frame (helix coordinates) of base pair ``i``."""
        j, n = i % self.k, i // self.k
        ss: SlotStrand = getattr(self.slots[j], strand)
        R, t = self.screw(n)
        return ss, Frame(R @ ss.frame.origin + t, R @ ss.frame.triad)

    def nucleotide(self, i: int, strand: str, letter: str, chem: str,
                   standard: Dict[str, NucleotideTemplate]
                   ) -> Tuple[List[str], np.ndarray, dict]:
        """Atom names and global coordinates of one nucleotide.

        ``letter`` "-" gives an abasic site. Also returns flags: ``swapped``
        (base taken from the standard bases), ``chem_changed`` (O2' added or
        removed) and, for an abasic site, ``abasic``.
        """
        ss, frame = self.slot_strand(i, strand)
        flags = {"swapped": False, "chem_changed": False}
        if letter == "-":
            # abasic site: the slot's sugar-phosphate without its base (O2' is
            # set from the template while the base still defines the face)
            t = next(iter(ss.templates.values()))
            names, xyz = list(t.names), t.xyz.copy()
            names, xyz = to_rna(names, xyz) if chem == "RNA" else to_dna(names, xyz)
            keep = [k for k, n in enumerate(names) if n in BACKBONE_ATOMS]
            return [names[k] for k in keep], frame.to_global(xyz[keep]), dict(flags, abasic=True)
        if letter in ss.templates:
            t = ss.templates[letter]
            names, xyz = list(t.names), t.xyz.copy()
        else:
            t = next(iter(ss.templates.values()))
            std = standard[letter]
            names, xyz = base_swap(t.names, t.xyz, std.names, std.xyz)
            flags["swapped"] = True
        if chem == "RNA" and "O2'" not in names:
            names, xyz = to_rna(names, xyz)
            flags["chem_changed"] = True
        elif chem == "DNA" and "O2'" in names:
            names, xyz = to_dna(names, xyz)
            flags["chem_changed"] = True
        return names, frame.to_global(xyz), flags


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def _frame(d: dict) -> Frame:
    return Frame(np.array(d["o"], float), np.array(d["R"], float))


def _slot_strand(d: dict) -> SlotStrand:
    return SlotStrand(frame=_frame(d["frame"]), native=d["native"], chem=d["chem"],
                      templates={k: NucleotideTemplate.from_json(v)
                                 for k, v in d["templates"].items()},
                      syn=d.get("syn", False), pucker=d.get("pucker"))


def model_from_json(d: dict) -> FiberModel:
    """A FiberModel from its entry in data/fiber_models.json."""
    slots = [Slot(I=_slot_strand(s["I"]), II=_slot_strand(s["II"]),
                  pairing=s.get("pairing", "wc")) for s in d["slots"]]
    return FiberModel(
        id=d["id"], name=d["name"], family=d["family"],
        description=d["description"], k=d["k"], twist=d["twist"], rise=d["rise"],
        slots=slots, kind=d.get("kind", "repeat"),
        sequence_policy=d.get("sequence_policy", "any"),
        handedness=d.get("handedness", "right"), bp_twist=d.get("bp_twist", 0.0),
        bp_rise=d.get("bp_rise", 0.0), references=d.get("references", []),
        meta=d.get("meta", {}))


@dataclass
class FiberLibrary:
    """The buildable models, the standard bases, the models left out (with the
    reason) and the provenance of the data."""

    models: Dict[int, FiberModel]
    standard: Dict[str, NucleotideTemplate]
    references: Dict[str, str]
    excluded: Dict[int, dict]
    provenance: dict

    def get(self, model_id: int) -> FiberModel:
        """A buildable model; KeyError says why an excluded one is left out."""
        try:
            return self.models[int(model_id)]
        except KeyError:
            if int(model_id) in self.excluded:
                why = "; ".join(self.excluded[int(model_id)].get("issues", []))
                raise KeyError(f"fiber model {model_id} is excluded: {why}") from None
            raise KeyError(f"unknown fiber model {model_id}") from None


def load_library(path: Optional[str] = None) -> FiberLibrary:
    """Read data/fiber_models.json (or ``path``); the sha256 of the file goes
    into the provenance, which the PDB header and the manifest report."""
    with open(path or os.path.join(DATA, "fiber_models.json"), encoding="utf-8") as fh:
        text = fh.read()
    raw = json.loads(text)
    models, excluded = {}, {}
    for key, d in raw["models"].items():
        if d.get("status") == "ok":
            models[int(key)] = model_from_json(d)
        else:
            excluded[int(key)] = {k: d.get(k) for k in
                                  ("id", "name", "family", "description", "issues")}
    standard = {k: NucleotideTemplate.from_json(v) for k, v in raw["standard_bases"].items()}
    prov = dict(raw.get("provenance", {}))
    prov["fiber_models_sha256"] = hashlib.sha256(text.encode()).hexdigest()
    return FiberLibrary(models, standard, raw.get("references", {}), excluded, prov)
