#!/usr/bin/env python3
"""Convert the human-edited YAML tables in data/ to the JSON read at runtime.

The browser runtime (Pyodide) ships without PyYAML, so the package reads only
JSON. The YAML files stay the editable source (they carry the comments that
justify each entry); tests/test_data_sync.py fails if the two drift apart.

    python tools/build_data.py
"""

from __future__ import annotations

import json
import os
import sys

import yaml

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TABLES = ["aliases", "base_pairing"]


def convert(name: str) -> dict:
    with open(os.path.join(ROOT, "data", f"{name}.yaml"), encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def main() -> int:
    for name in TABLES:
        out = os.path.join(ROOT, "data", f"{name}.json")
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(convert(name), fh, indent=1, sort_keys=True)
            fh.write("\n")
        print(f"data/{name}.yaml -> data/{name}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
