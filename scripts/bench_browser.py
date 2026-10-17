#!/usr/bin/env python3
# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""Wall time of a build in the browser (headless Firefox), click to result.

    scripts/build_web.sh && .venv/bin/python scripts/bench_browser.py [--out timings.json]

Canonical B-DNA (model #4) and fully modified siRNA-like duplexes (2'-OMe /
2'-F alternating, phosphorothioate ends; model #20) of increasing length.
Includes validation, fetching the modXNA library files and the build; the
first, cold build of each series is reported separately.
"""
import argparse, json, os, random, subprocess, sys, time
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PAIR_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}
PAIR_RNA = {"A": "U", "U": "A", "G": "C", "C": "G"}


def canonical(n, rng):
    s = "".join(rng.choice("ACGT") for _ in range(n))
    return s, "".join(PAIR_DNA[c] for c in s), 4


def modified(n, rng):
    s = "".join(rng.choice("ACGU") for _ in range(n))
    def mod(seq, first):
        toks = [(("m" if (i % 2 == 0) == (first == "m") else "f") + c) for i, c in enumerate(seq)]
        return toks[0] + "".join(("*" if i <= 2 or i >= len(toks) - 2 else "") + t
                                 for i, t in enumerate(toks) if i > 0)
    comp = "".join(PAIR_RNA[c] for c in s)          # 3'->5' under s
    return mod(s, "m"), mod(comp, "f"), 20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "paper", "data", "bench_browser.json"))
    ap.add_argument("--repeats", type=int, default=5)
    a = ap.parse_args()
    port = 8791
    srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "scripts", "serve.py"), str(port)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(1)
    o = Options(); o.add_argument("-headless")
    d = webdriver.Firefox(options=o, service=Service("/snap/bin/geckodriver"))
    rng = random.Random(1)
    out = {"canonical": {}, "modified": {}}
    try:
        t0 = time.time()
        d.get(f"http://localhost:{port}/")
        WebDriverWait(d, 180).until(lambda d: d.find_element(By.ID, "model").find_elements(By.TAG_NAME, "option"))
        out["load_s"] = round(time.time() - t0, 2)
        for kind, gen in (("canonical", canonical), ("modified", modified)):
            for n in (12, 24, 48, 96, 192):
                times = []
                for rep in range(a.repeats + 1):
                    top, bottom, model = gen(n, rng)
                    d.execute_script(
                        "document.body.dataset.buildMs='';"
                        "const s=(id,v)=>{const e=document.getElementById(id); e.value=v;"
                        "e.dispatchEvent(new Event('input'));e.dispatchEvent(new Event('change'));};"
                        "s('top',arguments[0]); s('bottom',arguments[1]); s('model',arguments[2]);",
                        top, bottom, str(model))
                    time.sleep(0.4)
                    d.find_element(By.ID, "build").click()
                    WebDriverWait(d, 120).until(lambda d: d.find_element(By.TAG_NAME, "body").get_attribute("data-build-ms"))
                    ms = int(d.find_element(By.TAG_NAME, "body").get_attribute("data-build-ms"))
                    times.append(ms)
                out[kind][n] = {"first_ms": times[0], "ms": times[1:]}
                print(kind, n, times)
    finally:
        d.quit(); srv.terminate()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print("load", out["load_s"], "s ->", a.out)


if __name__ == "__main__":
    main()
