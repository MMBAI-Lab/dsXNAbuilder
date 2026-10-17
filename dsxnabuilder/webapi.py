# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.webapi
===================

JSON dispatcher for the browser worker.

``call(name, payload_json) -> result_json`` with ``name`` in models | pieces
| validate | build | complement | batch_validate | batch_zip. Binary outputs
(the zips) travel base64-encoded: a few MB at most.
"""

from __future__ import annotations

import base64
import json
import traceback
from typing import Optional

from .catalog import default_chem
from .reslib import ResidueLibrary
from .seqlang import parse_line
from .service import Session

_session: Optional[Session] = None


def session() -> Session:
    """The worker's session, created on the first call."""
    global _session
    if _session is None:
        _session = Session(residue_library=ResidueLibrary.load_default())
    return _session


def complement(payload: dict) -> dict:
    """Watson-Crick complement of strand 1, written 3'->5' under it.

    Only the letters are complemented (canonical chemistry of the model's
    strand II); modifications are for the user to add on either strand. An
    abasic site has no partner: the complement keeps a "-" in its column for
    the user to replace with the base they want opposite it.
    """
    s = session()
    model = s.fibers.get(int(payload.get("model", 4)))
    pair = {"A": "U" if model.strand_chem("II") == "RNA" else "T",
            "T": "A", "U": "A", "G": "C", "C": "G"}
    strand, _ = parse_line(payload.get("top", ""), True, s.catalog, s.aliases,
                           default_chem(model.strand_chem("I")))
    out = "".join("-" if c is None or c.is_abasic else pair[c.letter] for c in strand.columns)
    return {"bottom": out}


def call(name: str, payload_json: str) -> str:
    """Run one request; errors come back as ``{"ok": false, "error": ...}``."""
    try:
        payload = json.loads(payload_json or "{}")
        s = session()
        if name == "models":
            result = s.models()
        elif name == "pieces":
            result = s.pieces()
        elif name == "validate":
            result = s.validate(payload)
        elif name == "complement":
            result = complement(payload)
        elif name in ("build", "batch_zip"):
            out = s.build(payload, payload.get("libs")) if name == "build" else s.batch_zip(payload)
            if out.get("zip") is not None:
                out["zip_b64"] = base64.b64encode(out.pop("zip")).decode()
            result = out
        elif name == "batch_validate":
            result = s.batch_validate(payload)
        else:
            raise ValueError(f"unknown call {name!r}")
        return json.dumps({"ok": True, "result": result})
    except Exception as exc:  # report, never crash the worker
        return json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}",
                           "trace": traceback.format_exc(limit=4)})
