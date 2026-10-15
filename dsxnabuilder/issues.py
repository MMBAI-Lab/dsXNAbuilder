"""
dsxnabuilder.issues
===================

Findings reported to the user: errors stop a build, warnings and notes do not.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class Issue:
    """A user-facing finding.

    ``code`` is stable (scripts and tests match on it), ``message`` is the
    English text the web app and the CLI show, and ``params`` carry the values
    in the message for programmatic use. ``position`` is the alignment column
    (0-based) and ``strand`` "top" or "bottom" when the finding concerns one
    residue.
    """

    level: str                          # "error" | "warning" | "info"
    position: Optional[int]
    code: str
    message: str
    params: Dict[str, object] = field(default_factory=dict)
    strand: Optional[str] = None

    def to_json(self) -> dict:
        return {"level": self.level, "code": self.code, "strand": self.strand,
                "position": self.position, "message": self.message,
                "params": self.params}
