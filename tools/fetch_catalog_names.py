#!/usr/bin/env python3
"""Names and notes of the modXNA fragments, from the modXNA catalogue page.

    python tools/fetch_catalog_names.py [--html saved_page.html]

The fragment files of the modXNA repository (data/modxna_catalog.json) carry
codes but no chemical names; the catalogue page of the modXNA website lists,
per section (backbone, sugar/linker, A-, C-, G-, T-, U-based and other base
modifications), each fragment code with its name and notes. This tool reads
those tables and writes data/modxna_names.json with the page URL, the date it
was read and the sha256 of the HTML, so the names can be traced to the page
version they came from. ``--html`` reads a saved copy instead of the network.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import html
import json
import os
import re
import sys
import urllib.request
from html.parser import HTMLParser

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
URL = "https://modxna.chpc.utah.edu/catalog/"
# section heading on the page -> (moiety, parent base)
SECTIONS = {
    "Backbone": ("backbone", None),
    "Sugar / Linker": ("sugar", None),
    "A-based modification": ("base", "A"),
    "C-based modification": ("base", "C"),
    "G-based modification": ("base", "G"),
    "T-based modification": ("base", "T"),
    "U-based modification": ("base", "U"),
    "Other base modification": ("base", "other"),
}


class _Tables(HTMLParser):
    """Rows of the table that follows each <h3> heading, as lists of cell texts."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.sections = []          # [(heading, [[cell, ...], ...])]
        self._h3 = None
        self._heading = None
        self._row = None
        self._cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "h3":
            self._h3 = []
        elif tag == "table" and self._heading:
            self.sections.append((self._heading, []))
            self._heading = None
        elif tag == "tr" and self.sections:
            self._row = []
        elif tag == "td" and self._row is not None:
            self._cell = []
        elif tag == "br" and self._cell is not None:
            self._cell.append("\n")

    def handle_endtag(self, tag):
        if tag == "h3" and self._h3 is not None:
            self._heading = "".join(self._h3).strip()
            self._h3 = None
        elif tag == "td" and self._cell is not None:
            self._row.append(_clean("".join(self._cell)))
            self._cell = None
        elif tag == "tr" and self._row is not None:
            self.sections[-1][1].append(self._row)
            self._row = None

    def handle_data(self, data):
        if self._h3 is not None:
            self._h3.append(data)
        if self._cell is not None:
            self._cell.append(data)


def _clean(s: str) -> str:
    lines = [re.sub(r"\s+", " ", html.unescape(x)).strip() for x in s.split("\n")]
    return "\n".join(x for x in lines if x)


def parse(text: str) -> list:
    p = _Tables()
    p.feed(text)
    out = []
    for heading, rows in p.sections:
        if heading not in SECTIONS:
            continue
        moiety, parent = SECTIONS[heading]
        entries = []
        for row in rows:
            if len(row) < 4 or row[0].upper().startswith("FRAGMENT"):
                continue
            codes = [c.strip() for c in row[0].split("/") if c.strip()]
            entries.append({"codes": codes, "name": row[2].replace("\n", " "),
                            "notes": row[3].split("\n") if row[3] else []})
        out.append({"title": heading, "moiety": moiety, "parent": parent,
                    "entries": entries})
    missing = set(SECTIONS) - {s["title"] for s in out}
    if missing:
        raise SystemExit(f"sections not found on the page: {sorted(missing)}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", help="saved copy of the catalogue page")
    a = ap.parse_args(argv)
    if a.html:
        with open(a.html, "rb") as fh:
            raw = fh.read()
    else:
        req = urllib.request.Request(URL, headers={"User-Agent": "dsXNAbuilder"})
        with urllib.request.urlopen(req, timeout=60) as r:
            raw = r.read()
    sections = parse(raw.decode("utf-8"))
    out = {"source": {"url": URL, "read": datetime.date.today().isoformat(),
                      "sha256": hashlib.sha256(raw).hexdigest()},
           "sections": sections}
    dest = os.path.join(ROOT, "data", "modxna_names.json")
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    for s in sections:
        print(f"{s['title']:26s} {len(s['entries']):4d}")
    print(f"-> {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
