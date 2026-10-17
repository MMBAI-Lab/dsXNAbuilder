#!/usr/bin/env python3
# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""End-to-end check of the static web app in headless Firefox.

    scripts/build_web.sh && .venv/bin/python scripts/e2e.py [--shots DIR]

Serves web/ locally (scripts/serve.py), loads the page, waits for Pyodide,
runs the one-click examples and checks the summary, the issues, the viewer
controls, the highlighted generic models, the high-throughput tab and the
downloads. Needs selenium and geckodriver.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
import zipfile

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PORT = 8767
# Stand-in for components/AppFrame.tsx of web-danslab: same-origin frame with a 1px
# border, content-box, no scrolling, sized from the app's "app:height" messages.
HOST = """
const url = arguments[0];
document.open();
document.write(`<!doctype html><html><body style="margin:0">
<iframe id="f" src="${url}" scrolling="no"
  style="width:1100px;height:560px;border:1px solid #888;box-sizing:content-box"></iframe>
<script>
addEventListener('message', (e) => {
  if (e.origin !== location.origin || !e.data || e.data.type !== 'app:height') return;
  document.getElementById('f').style.height =
    Math.min(Math.max(Math.round(e.data.height), 560), 20000) + 'px';
});
<\/script></body></html>`);
document.close();
"""
# inside the frame: [html box height, frame viewport height, scrollHeight, body overflow-y]
FRAME_METRICS = ("const r = document.documentElement; return [r.getBoundingClientRect().height, "
                 "innerHeight, r.scrollHeight, getComputedStyle(document.body).overflowY];")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shots", help="directory for screenshots")
    ap.add_argument("--url", help="test a deployed URL instead of a local server")
    a = ap.parse_args()
    server = None
    url = a.url
    if not url:
        server = subprocess.Popen([sys.executable, os.path.join(ROOT, "scripts", "serve.py"),
                                   str(PORT)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        url = f"http://localhost:{PORT}/"
        time.sleep(1.0)
    # snap Firefox cannot write to the host /tmp: download inside the repo (build/ is ignored)
    os.makedirs(os.path.join(ROOT, "build"), exist_ok=True)
    dl = tempfile.mkdtemp(prefix="e2e-", dir=os.path.join(ROOT, "build"))
    opts = Options()
    opts.add_argument("-headless")
    opts.set_preference("browser.download.folderList", 2)
    opts.set_preference("browser.download.dir", dl)
    opts.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/zip,chemical/x-pdb")
    opts.set_preference("intl.accept_languages", "en-US")
    # software WebGL so that the NGL viewer also runs headless
    opts.set_preference("webgl.force-enabled", True)
    opts.set_preference("webgl.disabled", False)
    opts.set_preference("gfx.webrender.software", True)
    gecko = "/snap/bin/geckodriver" if os.path.exists("/snap/bin/geckodriver") else "geckodriver"
    drv = webdriver.Firefox(options=opts, service=Service(gecko))
    drv.set_window_size(1400, 1100)
    fails = []

    def check(cond, msg):
        print(("ok   " if cond else "FAIL ") + msg)
        if not cond:
            fails.append(msg)

    def wait_file(path, tries=50):
        for _ in range(tries):
            if os.path.exists(path) and not os.path.exists(path + ".part"):
                return path
            time.sleep(0.2)
        return None

    def shot(name):
        if a.shots:
            os.makedirs(a.shots, exist_ok=True)
            drv.save_screenshot(os.path.join(a.shots, f"{name}.png"))

    try:
        t0 = time.time()
        drv.get(url)
        boot = drv.find_element(By.ID, "boot")
        check(boot.is_displayed() and boot.find_elements(By.CSS_SELECTOR, ".helix span"),
              f"loading animation while the engine starts ({boot.text!r})")
        WebDriverWait(drv, 180).until(
            lambda d: d.find_element(By.ID, "model").find_elements(By.TAG_NAME, "option"))
        check(True, f"engine ready in {time.time() - t0:.1f} s")
        WebDriverWait(drv, 10).until(lambda d: not d.find_element(By.ID, "boot").is_displayed())
        check(True, "loading animation gone once the engine is ready")
        n_opt = len(drv.find_element(By.ID, "model").find_elements(By.TAG_NAME, "option"))
        check(n_opt == 40, f"{n_opt} models in the selector")
        generic = drv.find_elements(By.CSS_SELECTOR, "#model option.generic")
        # the seven generic 3DNA models and C-DNA #8, which also has all four bases
        check(len(generic) == 8 and all(o.text.startswith("\u2605") for o in generic),
              f"{len(generic)} generic models highlighted in the selector")

        def run_example(ex, expect_bp):
            drv.find_element(By.CSS_SELECTOR, f'[data-ex="{ex}"]').click()
            WebDriverWait(drv, 60).until(
                lambda d: d.find_element(By.ID, "summary").is_displayed()
                and f">{expect_bp}<" in d.find_element(By.ID, "summary").get_attribute("innerHTML")
                and not d.find_element(By.ID, "dlzip").get_attribute("disabled"))
            time.sleep(0.5)
            return drv.find_element(By.ID, "summary").text, \
                drv.find_element(By.ID, "issues").text

        s, iss = run_example("dickerson", 12)
        check("758" in s and "486 heavy" in s and "1.60" in s,
              "Dickerson: 758 atoms (486 heavy), O3'-P 1.60")
        check(not drv.find_elements(By.CSS_SELECTOR, "#viewer .busy")
              and drv.find_element(By.ID, "build").text == "Build duplex",
              "build animation removed after the build")
        pressed = [b.get_attribute("data-v") for b in
                   drv.find_elements(By.CSS_SELECTOR, '#repr button[aria-pressed="true"], '
                                     '#color button[aria-pressed="true"]')]
        check(pressed == ["sticks", "base"], f"viewer defaults {pressed}")
        has_canvas = bool(drv.find_elements(By.CSS_SELECTOR, "#viewer canvas"))
        print(("ok   " if has_canvas else "note ") + "3D viewer " +
              ("rendered" if has_canvas else "not available (no WebGL in this browser)"))
        if has_canvas:      # the legend belongs to the viewer
            legend = drv.find_element(By.ID, "legend").text
            check(all(x in legend for x in ("A", "C", "G", "T", "sugar\u2013phosphate")),
                  "base colour legend")
        shot("dickerson")
        s, iss = run_example("zdna", 6)
        check("SEQ_PATTERN" not in iss, "Z-DNA (CG)3 fits the repeat")
        s, iss = run_example("hybrid", 10)
        check("ORIENTATION_FLIP" not in iss, "RNA:DNA hybrid on model 22")
        s, iss = run_example("mismatch", 14)
        check("MISMATCH" in iss and "OVERHANG" in iss, "mismatch + overhang reported")
        shot("mismatch")
        s, iss = run_example("abasic", 12)
        check("ABASIC" in iss, "abasic site built and reported")
        s, iss = run_example("lna", 16)
        check("modXNA residues" in s and "TRANSPLANT_STRAIN" not in iss,
              "LNA gapmer : RNA built with modXNA residues")
        drv.find_element(By.CSS_SELECTOR, '#color button[data-v="modified"]').click()
        if has_canvas:
            legend = drv.find_element(By.ID, "legend").text
            check("canonical" in legend and "modified" in legend, "modified colour legend")
        shot("lna")
        drv.find_element(By.CSS_SELECTOR, '#color button[data-v="base"]').click()
        for ex, n in (("sirna", 10), ("moe", 10), ("mc", 12)):
            s, iss = run_example(ex, n)
            check("modXNA residues" in s, f"example {ex} built")
        # explicit triplet example, then rename its residues
        s, iss = run_example("modxna", 14)
        rows = drv.find_elements(By.CSS_SELECTOR, "#nametable tbody tr")
        check("ABASIC" in iss and len(rows) == 3 and "PS2/OME/M5C" in rows[0].text,
              f"modXNA code example: {len(rows)} named residues")
        inp = drv.find_element(By.CSS_SELECTOR, "#nametable input.resname")
        inp.clear()
        inp.send_keys("omc")
        drv.find_element(By.ID, "applynames").click()
        WebDriverWait(drv, 60).until(lambda d: "names=" in d.current_url)
        time.sleep(0.5)
        first = drv.find_element(By.CSS_SELECTOR, "#nametable input.resname").get_attribute("value")
        check(first == "OMC" and "PS2%2FOME%2FM5C%2Fint%3DOMC" in drv.current_url,
              "residue renamed to OMC (kept in the reproducible link)")
        run_example("mc", 12)        # the package below needs a duplex without abasic sites
        drv.find_element(By.CSS_SELECTOR, "#advanced summary").click()
        drv.find_element(By.ID, "dlzip").click()
        path = wait_file(os.path.join(dl, "dsxnabuilder-amber.zip"))
        names = zipfile.ZipFile(path).namelist() if path else []
        check("dsxnabuilder/frcmod.modxna" in names
              and any(n.startswith("dsxnabuilder/libs/X") for n in names)
              and "dsxnabuilder/min.in" not in names,
              "advanced package: libs/ + frcmod.modxna, topology only")
        if path:
            os.remove(path)
        drv.find_element(By.CSS_SELECTOR, "#advanced summary").click()

        # live validation error
        top = drv.find_element(By.ID, "top")
        top.clear()
        top.send_keys("ACGT")
        bot = drv.find_element(By.ID, "bottom")
        bot.clear()
        bot.send_keys("TGC")
        WebDriverWait(drv, 20).until(lambda d: "LENGTH_MISMATCH" in d.find_element(By.ID, "issues").text)
        check(drv.find_element(By.ID, "build").get_attribute("disabled") is not None,
              "length mismatch disables Build")

        check(not drv.find_elements(By.CSS_SELECTOR, "#lang, #theme"),
              "no language or theme switch (English only, theme from danslab.xyz)")

        # default download: the PDB
        run_example("dickerson", 12)
        drv.find_element(By.ID, "dlpdb").click()
        path = wait_file(os.path.join(dl, "dsxnabuilder.pdb"))
        text = open(path).read() if path else ""
        atoms = [line for line in text.splitlines() if line.startswith("ATOM")]
        n_h = sum(1 for line in atoms if line[76:78].strip() == "H")
        check(len(atoms) == 758 and n_h == 272 and "REMARK 999 strand 1" in text,
              "PDB downloaded: 758 atoms with 272 hydrogens and self-describing header")

        # Reset (top right): back to the starting state
        drv.find_element(By.ID, "model").send_keys("#20")
        drv.find_element(By.ID, "reset").click()
        time.sleep(0.5)
        check(drv.find_element(By.ID, "top").get_attribute("value") == ""
              and drv.find_element(By.ID, "bottom").get_attribute("value") == ""
              and drv.find_element(By.ID, "model").get_attribute("value") == "4"
              and not drv.find_element(By.ID, "summary").is_displayed()
              and not drv.find_element(By.ID, "names").is_displayed()
              and drv.find_element(By.ID, "dlpdb").get_attribute("disabled") is not None
              and "?" not in drv.current_url
              and "appears here" in drv.find_element(By.ID, "viewer").text,
              "Reset: empty strands, model #4, no structure, clean URL")
        run_example("dickerson", 12)
        check("758" in drv.find_element(By.ID, "summary").text, "builds again after Reset")

        # high-throughput tab: the example, its zip, a record opened in the Build tab
        drv.find_element(By.ID, "tab-batch").click()
        drv.find_element(By.ID, "batchexample").click()
        WebDriverWait(drv, 30).until(
            lambda d: len(d.find_elements(By.CSS_SELECTOR, "#batchtable tbody tr")) == 4
            and not d.find_element(By.ID, "batchbuild").get_attribute("disabled"))
        check(True, "high-throughput: example checked, 4 duplexes listed")
        drv.find_element(By.ID, "batchbuild").click()
        path = wait_file(os.path.join(dl, "dsxnabuilder_batch.zip"), tries=300)
        names = zipfile.ZipFile(path).namelist() if path else []
        pdbs = sorted(n.split("/")[-1] for n in names if n.endswith(".pdb"))
        check(pdbs == ["dickerson.pdb", "lna_gapmer.pdb", "sequence1.pdb", "sequence2.pdb"]
              and "dsxnabuilder_batch/report.tsv" in names,
              f"high-throughput: zip with {len(pdbs)} PDB files and the report")
        WebDriverWait(drv, 10).until(
            lambda d: not d.find_element(By.ID, "batchprogress").is_displayed())
        sts = [e.text for e in drv.find_elements(By.CSS_SELECTOR, "#batchtable .st")]
        status = drv.find_element(By.ID, "batchstatus").text
        check(sts == ["built"] * 4 and status.startswith("4 of 4 built")
              and not drv.find_element(By.ID, "batchzip").get_attribute("disabled"),
              f"high-throughput: {status}")
        shot("batch")
        if path:
            os.remove(path)
        drv.find_element(By.CSS_SELECTOR, '#batchtable [data-view="2"]').click()     # dickerson
        WebDriverWait(drv, 60).until(lambda d: "758" in d.find_element(By.ID, "summary").text)
        check(drv.find_element(By.ID, "tab-build").get_attribute("aria-selected") == "true",
              "high-throughput: a record opens in the Build tab")
        drv.find_element(By.ID, "tab-batch").click()
        box = drv.find_element(By.ID, "batchtext")
        box.clear()
        box.send_keys(">a; colour=red\nACGT\nTGCA\n")
        WebDriverWait(drv, 20).until(lambda d: "BATCH_FORMAT" in d.find_element(By.ID, "batchissues").text)
        check(drv.find_element(By.ID, "batchbuild").get_attribute("disabled") is not None,
              "high-throughput: a format mistake disables the build")
        drv.find_element(By.ID, "batchclear").click()

        # models tab
        drv.find_element(By.ID, "tab-models").click()
        rows = drv.find_elements(By.CSS_SELECTOR, "#modeltable tbody tr")
        check(len(rows) == 40, f"models table has {len(rows)} rows")
        grows = drv.find_elements(By.CSS_SELECTOR, "#modeltable tbody tr.generic")
        check(len(grows) == 8 and all("generic" in r.text for r in grows),
              f"{len(grows)} generic models highlighted in the table")
        shot("models")
        # modXNA library tab: click a sugar, a base and a linkage, add it to strand 1
        drv.find_element(By.ID, "tab-library").click()
        WebDriverWait(drv, 30).until(lambda d: d.find_elements(By.CSS_SELECTOR, "#c-sugar button"))
        for sel in ('[data-sugar="5"]', '[data-base="5"]', '[data-link="1"]'):     # +, 5mC, *1
            drv.find_element(By.CSS_SELECTOR, f"#composer {sel}").click()
        top_code = drv.find_element(By.ID, "c-top-code").text
        bottom_code = drv.find_element(By.ID, "c-bottom-code").text
        n_pieces = len(drv.find_elements(By.CSS_SELECTOR, "#piecetable tbody tr"))
        check(top_code == "*1+C^5mC" and bottom_code == "+C^5mC*1" and n_pieces == 15,
              f"composer: {top_code} / {bottom_code}, {n_pieces} pieces listed")
        shot("library")
        drv.find_element(By.ID, "c-top").click()
        top = drv.find_element(By.ID, "top").get_attribute("value")
        check(drv.find_element(By.ID, "tab-build").get_attribute("aria-selected") == "true"
              and top.endswith("*1+C^5mC"), "composer: nucleotide added to strand 1")
        drv.find_element(By.ID, "tab-docs").click()
        shot("docs")

        # embedded in danslab.xyz: one scrollbar, the site's. The frame must take the
        # app's height (no overflow inside it), grow with the help and shrink back.
        drv.get(url.rstrip("/") + "/style.css")      # any same-origin page shown inline
        drv.execute_script(HOST, url)
        drv.switch_to.frame(drv.find_element(By.ID, "f"))
        WebDriverWait(drv, 180).until(
            lambda d: d.find_element(By.ID, "model").find_elements(By.TAG_NAME, "option"))

        def settled():
            time.sleep(1.5)     # ResizeObserver + the 1 s fallback report
            return drv.execute_script(FRAME_METRICS)

        heights = {}
        for tab in ("build", "batch", "docs", "models", "library", "build"):
            drv.find_element(By.ID, f"tab-{tab}").click()
            if tab == "library":
                WebDriverWait(drv, 30).until(
                    lambda d: d.find_elements(By.CSS_SELECTOR, "#piecetable tbody tr"))
            box, view, scroll, ovf = settled()
            heights.setdefault(tab, []).append(view)
            check(abs(view - box) <= 2 and scroll <= view and ovf == "visible",
                  f"embedded, {tab}: frame {view} px for a {box:.0f} px page, no inner scroll")
            if tab == "models":
                inner = drv.execute_script(
                    "return [...document.querySelectorAll('.tablewrap')]"
                    ".every(t => t.scrollHeight <= t.clientHeight + 1);")
                check(inner, "embedded: model tables flow with the page (no scroll box)")
        check(heights["docs"][0] > heights["build"][0] + 500
              and abs(heights["build"][1] - heights["build"][0]) <= 2,
              f"embedded: frame grows for the help and shrinks back ({heights})")
        drv.switch_to.default_content()
        logs = []
        try:
            logs = drv.get_log("browser")
        except Exception:
            pass
        for entry in logs:
            print("console:", entry)
    finally:
        drv.quit()
        if server:
            server.terminate()
    print(f"\n{len(fails)} failure(s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
