"""
dsxnabuilder.batch
==================

Many duplexes at once, from a FASTA-like text with one record per duplex:

    >name=sirna_1; model=20
    mG*mA*fC*mU*fA*mG*fC*mU*mG*mAdGdG
    C*U*GAUCGAC*U--

A record starts with a ``>`` header of ``key=value`` fields separated by
``;``: ``name`` (letters, digits, '.', '_' and '-'; it names the PDB file)
and ``model`` (3DNA fiber model; DEFAULT_MODEL when absent). A field without
'=' is the name, so ``>sirna_1; model=20`` and ``>sirna_1`` work too. The two
lines that follow are strand 1 (5'->3') and strand 2 (3'->5'), written as on
the Build tab. Blank lines and lines starting with '#' are ignored.

Every record is built as a single duplex, so its PDB is the file the Build
tab gives for it. The PDB files go into one zip (or directory) with a report
(report.tsv), the input and a README with the command that rebuilds them. A
record that cannot be built is reported and the others are still built; a
mistake in the format stops the whole batch. The web app takes up to
WEB_LIMIT records at a time and the command line, which runs on local
resources, up to CLI_LIMIT.
"""

from __future__ import annotations

import io
import os
import re
import zipfile
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import __version__
from .issues import Issue

WEB_LIMIT = 50
CLI_LIMIT = 500
DEFAULT_MODEL = 4                  # B-DNA, as on the Build tab
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
FIELDS = ("name", "model")
ZIP_ROOT = "dsxnabuilder_batch"
REPORT_COLUMNS = ("name", "model", "status", "base_pairs", "atoms", "heavy_atoms",
                  "modxna_residues", "o3p_min", "o3p_max", "clashes", "warnings", "info",
                  "errors")


@dataclass
class Record:
    """One duplex of a batch; ``line`` is the line of its header (1-based)."""

    name: str
    model: int
    top: str
    bottom: str
    line: int

    def request(self) -> dict:
        """The build request of the record (service.Session.build)."""
        return {"top": self.top, "bottom": self.bottom, "model": self.model}


def parse_batch(text: str, limit: int = WEB_LIMIT,
                default_model: int = DEFAULT_MODEL) -> Tuple[List[Record], List[Issue]]:
    """The records of a batch text and the mistakes in its format.

    Records with a mistake are left out; the batch should be built only when
    the issue list is empty.
    """
    blocks, issues = _blocks(text)
    records: List[Record] = []
    for index, (line, header, seqs) in enumerate(blocks, start=1):
        record, found = _record(index, line, header, seqs, default_model, records)
        if found:
            issues += [_format_issue(line, msg) for msg in found]
        else:
            records.append(record)
    if not blocks and not issues:
        issues.append(Issue("error", None, "BATCH_EMPTY",
                            "no duplexes: each one starts with a '>' header line", {}))
    if len(blocks) > limit:
        issues.insert(0, Issue("error", None, "BATCH_LIMIT",
                               f"{len(blocks)} duplexes: at most {limit} at a time",
                               {"count": len(blocks), "max": limit}))
    return records, issues


def _blocks(text: str) -> Tuple[List[Tuple[int, str, List[str]]], List[Issue]]:
    """(header line number, header text, sequence lines) of every record."""
    blocks: List[Tuple[int, str, List[str]]] = []
    issues: List[Issue] = []
    for n, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            blocks.append((n, line[1:], []))
        elif blocks:
            blocks[-1][2].append(line)
        elif not issues:
            issues.append(_format_issue(n, "sequence before the first '>' header line"))
    return blocks, issues


def _format_issue(line: int, message: str) -> Issue:
    return Issue("error", None, "BATCH_FORMAT", f"line {line}: {message}", {"line": line})


def _record(index: int, line: int, header: str, seqs: List[str], default_model: int,
            taken: Sequence[Record]) -> Tuple[Optional[Record], List[str]]:
    """The record of one header and its lines, or the mistakes found."""
    found: List[str] = []
    values: Dict[str, str] = {}
    for part in (p.strip() for p in header.split(";")):
        if not part:
            continue
        key, value = ([x.strip() for x in part.split("=", 1)] if "=" in part
                      else ["name", part])
        key = key.lower()
        if key not in FIELDS:
            found.append(f"unknown header field {key!r} (use name= and model=)")
        elif key in values:
            found.append(f"the header gives {key} twice")
        else:
            values[key] = value
    name = values.get("name") or f"duplex_{index:03d}"
    if not NAME_RE.match(name):
        found.append(f"name {name!r}: use letters, digits, '.', '_' or '-' "
                     "(at most 64, starting with a letter or digit)")
    elif name.lower() in {r.name.lower() for r in taken}:
        found.append(f"name {name!r} is used twice (names become file names)")
    model = default_model
    if "model" in values:
        try:
            model = int(values["model"])
        except ValueError:
            found.append(f"model {values['model']!r} is not a 3DNA fiber model number")
    if len(seqs) != 2:
        found.append(f"{name}: expected two lines after the header (strand 1 5'->3', "
                     f"strand 2 3'->5'), found {len(seqs)}")
    if found:
        return None, found
    return Record(name, model, seqs[0], seqs[1], line), []


def build_records(session, records: Sequence[Record],
                  on_built: Optional[Callable[[int, Record, dict], None]] = None) -> List[dict]:
    """Build every record with ``session`` (service.Session): one build output
    per record, in order. ``on_built(index, record, output)`` follows each one."""
    outcomes = []
    for i, r in enumerate(records):
        outcomes.append(session.build(r.request()))
        if on_built is not None:
            on_built(i, r, outcomes[-1])
    return outcomes


def batch_files(records: Sequence[Record], outcomes: Sequence[dict], text: str,
                default_model: int = DEFAULT_MODEL) -> Dict[str, bytes]:
    """The output files: pdb/<name>.pdb for every duplex built, report.tsv,
    input.txt and README.txt."""
    files: Dict[str, bytes] = {}
    for r, out in zip(records, outcomes):
        if out.get("ok") and out.get("pdb"):
            files[f"pdb/{r.name}.pdb"] = out["pdb"].encode()
    files["report.tsv"] = report_tsv(records, outcomes).encode()
    files["input.txt"] = text.encode()
    files["README.txt"] = readme_text(records, outcomes, default_model).encode()
    return files


def zip_bytes(files: Dict[str, bytes]) -> bytes:
    """The files in a zip under ZIP_ROOT/; fixed timestamps and order, so the
    same batch gives the same bytes."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name in sorted(files):
            info = zipfile.ZipInfo(f"{ZIP_ROOT}/{name}", date_time=(2026, 1, 1, 0, 0, 0))
            info.external_attr = 0o644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, files[name])
    return buf.getvalue()


def write_directory(files: Dict[str, bytes], path: str) -> None:
    """The files under directory ``path`` (created when missing)."""
    for name, data in files.items():
        target = os.path.join(path, *name.split("/"))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, "wb") as fh:
            fh.write(data)


def report_rows(records: Sequence[Record], outcomes: Sequence[dict]) -> List[dict]:
    """One row per record: what was built and what was found."""
    rows = []
    for r, out in zip(records, outcomes):
        issues = out.get("issues", [])
        sm = out.get("summary") or {}
        o3p = sm.get("o3p") or ["", ""]
        rows.append({
            "name": r.name, "model": r.model, "status": "built" if out.get("ok") else "failed",
            "base_pairs": sm.get("n_bp", ""), "atoms": sm.get("n_atoms", ""),
            "heavy_atoms": sm.get("n_heavy", ""), "modxna_residues": sm.get("n_modified", ""),
            "o3p_min": o3p[0], "o3p_max": o3p[1], "clashes": sm.get("clashes", ""),
            "warnings": _codes(issues, "warning"), "info": _codes(issues, "info"),
            "errors": "; ".join(f"{i['code']}: {i['message']}" for i in issues
                                if i["level"] == "error")})
    return rows


def _codes(issues: Sequence[dict], level: str) -> str:
    """Codes of the issues of one level, once each."""
    return ", ".join(dict.fromkeys(i["code"] for i in issues if i["level"] == level))


def report_tsv(records: Sequence[Record], outcomes: Sequence[dict]) -> str:
    """report.tsv: REPORT_COLUMNS, one line per record."""
    lines = ["\t".join(REPORT_COLUMNS)]
    for row in report_rows(records, outcomes):
        lines.append("\t".join(re.sub(r"[\t\r\n]+", " ", str(row[c]))
                               for c in REPORT_COLUMNS))
    return "\n".join(lines) + "\n"


def readme_text(records: Sequence[Record], outcomes: Sequence[dict],
                default_model: int = DEFAULT_MODEL) -> str:
    """README.txt of a batch."""
    built = sum(1 for o in outcomes if o.get("ok"))
    failed = len(records) - built
    command = "dsxnabuilder build --high-throughput input.txt -o structures.zip"
    if default_model != DEFAULT_MODEL:
        command += f" --model {default_model}"
    return f"""dsXNAbuilder {__version__} - batch of {len(records)} duplexes
===============================================================

{built} built, {failed} failed (see report.tsv).

pdb/<name>.pdb  one PDB per duplex built: the same file the Build tab gives for it,
                with the strands, the fiber model and the data versions in its header
report.tsv      per duplex: model, size, O3'-P range, heavy-atom contacts, warnings,
                and the errors of the duplexes that could not be built
input.txt       the input, as given

Rebuild with the dsXNAbuilder command-line tool:
    {command}

Please cite
  3DNA fiber models: Lu XJ, Olson WK (2003) Nucleic Acids Res 31:5108-5121;
                     Lu XJ, Olson WK (2008) Nat Protoc 3:1213-1227.
  modXNA (if modified residues are used): Love O et al. (2024)
                     J Chem Theory Comput 20:9354-9363.

Fiber templates derived from 3DNA data (CC BY-NC 4.0); modXNA residues are
GPL-3.0. Your structures are yours.
"""
