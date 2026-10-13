"""JSON dispatcher for the browser worker.

``call(name, payload_json) -> result_json`` with ``name`` in models | pieces
| validate | build | complement. Binary outputs (the zip) travel base64-encoded; the
package is a few hundred kB at most.
"""

from __future__ import annotations

import base64
import json
import traceback

from .service import Session

_session = None


def session() -> Session:
    global _session
    if _session is None:
        lib = None
        try:
            from .reslib import ResidueLibrary
            lib = ResidueLibrary.load_default()
        except Exception:
            lib = None
        _session = Session(residue_library=lib)
    return _session


def complement(payload: dict) -> dict:
    """Watson-Crick complement of strand 1, written 3'->5' under it.

    Only the letters are complemented (canonical chemistry of the model's
    strand II); modifications are for the user to add on either strand.
    """
    from .seqlang import parse_line
    s = session()
    model = s.fibers.get(int(payload.get("model", 4)))
    chem2 = model.strand_chem("II")
    strand, _ = parse_line(payload.get("top", ""), True, s.catalog, s.aliases,
                           model.strand_chem("I") if model.strand_chem("I") in ("DNA", "RNA")
                           else "DNA")
    pair = {"A": "U" if chem2 == "RNA" else "T", "T": "A", "U": "A", "G": "C", "C": "G"}
    out = "".join("-" if c is None else pair[c.letter] for c in strand.columns)
    return {"bottom": out}


def call(name: str, payload_json: str) -> str:
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
        elif name == "build":
            out = s.build(payload, payload.get("libs"))
            if out.get("zip") is not None:
                out["zip_b64"] = base64.b64encode(out.pop("zip")).decode()
            result = out
        else:
            raise ValueError(f"unknown call {name!r}")
        return json.dumps({"ok": True, "result": result})
    except Exception as exc:  # report, never crash the worker
        return json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}",
                           "trace": traceback.format_exc(limit=4)})
