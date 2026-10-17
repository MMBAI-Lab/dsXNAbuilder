# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Many duplexes at once: the FASTA-like format, the zip and the command line."""

import base64
import io
import json
import os
import zipfile

import pytest

from dsxnabuilder.batch import (CLI_LIMIT, WEB_LIMIT, batch_files, build_records, parse_batch,
                                zip_bytes)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
needs_lib = pytest.mark.skipif(not os.path.isdir(os.path.join(ROOT, "web", "lib", "residues")),
                               reason="residue library not built (tools/build_residue_library.py)")

EXAMPLE = """# dsXNAbuilder high-throughput input
>name=sequence1; model=20
mG*mA*fC*mU*fA*mG*fC*mU*mG*mAdGdG
C*U*GAUCGAC*U--

>name=sequence2; model=4
mG*mA*fC*mU*fA*mG*fC*mU*mG*mA
C*U*GAUCGAC*U

>dickerson
5'-CGCGAATTCGCG-3'
3'-GCGCTTAAGCGC-5'
"""


def codes(issues):
    return [i.code for i in issues]


def test_records_and_header_forms():
    records, issues = parse_batch(EXAMPLE)
    assert issues == []
    assert [(r.name, r.model, r.line) for r in records] == [
        ("sequence1", 20, 2), ("sequence2", 4, 6), ("dickerson", 4, 10)]
    assert records[0].top == "mG*mA*fC*mU*fA*mG*fC*mU*mG*mAdGdG"
    assert records[0].bottom == "C*U*GAUCGAC*U--"
    # a bare field is the name; spacing, case of keys, CRLF and missing names
    text = ">sirna ; MODEL = 22\r\nACGU\r\nUGCA\r\n>model=20\nACGU\nUGCA\n"
    records, issues = parse_batch(text, default_model=7)
    assert issues == []
    assert [(r.name, r.model) for r in records] == [("sirna", 22), ("duplex_002", 20)]
    records, _ = parse_batch(">x\nACGT\nTGCA\n", default_model=7)
    assert records[0].model == 7


@pytest.mark.parametrize("text, fragment", [
    (">a; colour=red\nAC\nTG\n", "unknown header field 'colour'"),
    (">a; model=B\nAC\nTG\n", "model 'B' is not a 3DNA fiber model number"),
    (">name=a; name=b\nAC\nTG\n", "gives name twice"),
    (">my duplex\nAC\nTG\n", "name 'my duplex'"),
    (">a\nACGT\n", "expected two lines after the header"),
    (">a\nAC\nTG\nCA\n", "found 3"),
    (">a\nAC\nTG\n>A\nAC\nTG\n", "name 'A' is used twice"),
    ("ACGT\n>a\nAC\nTG\n", "line 1: sequence before the first '>' header line"),
])
def test_format_mistakes(text, fragment):
    records, issues = parse_batch(text)
    assert issues and all(i.level == "error" and i.code == "BATCH_FORMAT" for i in issues)
    assert any(fragment in i.message for i in issues), [i.message for i in issues]
    assert all("line" in i.params for i in issues)


def test_empty_and_limits():
    assert codes(parse_batch("# nothing here\n\n")[1]) == ["BATCH_EMPTY"]
    many = "".join(f">d{i}\nACGT\nACGT\n" for i in range(WEB_LIMIT + 1))
    records, issues = parse_batch(many)
    assert codes(issues) == ["BATCH_LIMIT"] and issues[0].params == {"count": WEB_LIMIT + 1,
                                                                      "max": WEB_LIMIT}
    assert parse_batch(many, CLI_LIMIT)[1] == []


def test_each_pdb_is_the_single_build(session):
    """Every duplex of a batch gives the PDB of the same single build; failures
    are reported and the rest are built."""
    text = EXAMPLE.replace("model=20", "model=4") + ">bad; model=4\n+C+G+C\nGCG\n"
    records, issues = parse_batch(text)
    assert issues == []
    outcomes = build_records(session, records)
    # without the residue library only the canonical duplex builds
    assert [o["ok"] for o in outcomes] == [False, False, True, False]
    files = batch_files(records, outcomes, text)
    assert sorted(f for f in files if f.startswith("pdb/")) == ["pdb/dickerson.pdb"]
    single = session.build({"top": "5'-CGCGAATTCGCG-3'", "bottom": "3'-GCGCTTAAGCGC-5'",
                            "model": 4})
    assert files["pdb/dickerson.pdb"].decode() == single["pdb"]
    report = files["report.tsv"].decode().splitlines()
    head = report[0].split("\t")
    rows = [dict(zip(head, line.split("\t"))) for line in report[1:]]
    assert [r["status"] for r in rows] == ["failed", "failed", "built", "failed"]
    assert rows[2]["atoms"] == "758" and rows[2]["heavy_atoms"] == "486"
    assert "PUCKER_LOCKED" in rows[3]["errors"]
    assert files["input.txt"].decode() == text
    readme = files["README.txt"].decode()
    assert "1 built, 3 failed" in readme
    assert "dsxnabuilder build --high-throughput input.txt" in readme
    # the same batch gives the same bytes
    z = zip_bytes(files)
    assert z == zip_bytes(batch_files(records, outcomes, text))
    assert "dsxnabuilder_batch/pdb/dickerson.pdb" in zipfile.ZipFile(io.BytesIO(z)).namelist()


@needs_lib
def test_web_calls():
    """batch_validate lists the records with the library files each one needs;
    batch_zip packs the builds the page made one by one."""
    from dsxnabuilder import webapi

    def call(name, payload):
        return json.loads(webapi.call(name, json.dumps(payload)))

    v = call("batch_validate", {"text": EXAMPLE})
    assert v["ok"] and v["result"]["ok"], v
    rows = v["result"]["records"]
    assert [r["name"] for r in rows] == ["sequence1", "sequence2", "dickerson"]
    assert rows[0]["needed_libs"] and rows[2]["needed_libs"] == []
    results = {}
    for r in rows[:2]:          # the page stopped before the third
        out = call("build", {"top": r["top"], "bottom": r["bottom"], "model": r["model"]})
        results[r["name"]] = out["result"]
    z = call("batch_zip", {"text": EXAMPLE, "results": results})["result"]
    assert (z["built"], z["failed"]) == (2, 1)
    names = zipfile.ZipFile(io.BytesIO(base64.b64decode(z["zip_b64"]))).namelist()
    assert "dsxnabuilder_batch/pdb/sequence2.pdb" in names
    assert "dsxnabuilder_batch/pdb/dickerson.pdb" not in names
    bad = call("batch_validate", {"text": ">a; colour=red\nAC\nTG\n"})["result"]
    assert not bad["ok"] and bad["issues"][0]["code"] == "BATCH_FORMAT"


@needs_lib
def test_cli_high_throughput(tmp_path, capsys):
    from dsxnabuilder import cli
    src = tmp_path / "duplexes.txt"
    src.write_text(EXAMPLE)
    out = tmp_path / "structures.zip"
    assert cli.main(["build", "--high-throughput", str(src), "-o", str(out)]) == 0
    z = zipfile.ZipFile(out)
    assert sorted(n for n in z.namelist() if n.endswith(".pdb")) == [
        "dsxnabuilder_batch/pdb/dickerson.pdb", "dsxnabuilder_batch/pdb/sequence1.pdb",
        "dsxnabuilder_batch/pdb/sequence2.pdb"]
    err = capsys.readouterr().err
    assert "[1/3] sequence1: built" in err
    # a directory instead of a zip; a failing record gives exit status 1
    src.write_text(EXAMPLE + ">lna_bdna; model=4\n+C+G+C\nGCG\n")
    assert cli.main(["build", "--high-throughput", str(src), "-o", str(tmp_path / "out")]) == 1
    assert (tmp_path / "out" / "pdb" / "dickerson.pdb").exists()
    assert not (tmp_path / "out" / "pdb" / "lna_bdna.pdb").exists()
    assert "lna_bdna\t4\tfailed" in (tmp_path / "out" / "report.tsv").read_text()
    # a mistake in the format: nothing is built
    src.write_text(">a; colour=red\nAC\nTG\n")
    assert cli.main(["build", "--high-throughput", str(src), "-o", str(tmp_path / "x.zip")]) == 1
    assert not (tmp_path / "x.zip").exists()
    assert "BATCH_FORMAT" in capsys.readouterr().err


def test_cli_high_throughput_options(tmp_path):
    from dsxnabuilder import cli
    src = tmp_path / "d.txt"
    src.write_text(">a\nACGT\nTGCA\n")
    for extra in (["--top", "AC"], ["--amber", "p.zip"], ["--name", "K=V"], ["-o", "x.pdb"]):
        with pytest.raises(SystemExit) as exc:
            cli.main(["build", "--high-throughput", str(src)] + extra)
        assert exc.value.code == 2, extra
    with pytest.raises(SystemExit) as exc:
        cli.main(["build", "--bottom", "AC"])          # neither --top nor a batch
    assert exc.value.code == 2
    many = "".join(f">d{i}\nACGT\nTGCA\n" for i in range(CLI_LIMIT + 1))
    src.write_text(many)
    assert cli.main(["build", "--high-throughput", str(src), "-o", str(tmp_path / "m.zip")]) == 1
    assert not (tmp_path / "m.zip").exists()
