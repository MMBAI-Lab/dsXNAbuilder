#!/usr/bin/env python3
"""Regenera ``data/modxna_catalog.json`` desde un arbol de modXNA.

    python tools/build_catalog.py --modxna /opt/modxna [-o data/modxna_catalog.json]

Lee los ``mol2`` de ``dat/lib_{backbone,sugar,base}``, extrae los metadatos de
conexion de la linea de titulo, suma las cargas parciales y calcula la
pseudorrotacion de cada azucar con el mismo codigo que usa el servidor
(``xnafiber.geometry.pseudorotation``), de modo que la anotacion de pucker del
catalogo y la del control de calidad no pueden divergir.

Ejecutar siempre que se actualice el catalogo de modXNA. El fichero resultante
lleva el hash del contenido del arbol de origen; la cache de residuos NO se
invalida sola al cambiar el catalogo (su clave depende del sha256 de
``modxna.sh``, no de los fragmentos), asi que si cambia ``content_sha256``
conviene decidir si procede vaciarla.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from typing import List, Optional

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from xnafiber.geometry import pseudorotation

MOIETIES = ("backbone", "sugar", "base")
GLYCOSIDIC = {"N9": "purine", "N1": "pyrimidine"}


def parse_title(line: str) -> Optional[dict]:
    """Interpreta la linea de titulo modXNA de un mol2."""
    line = line.strip()
    for prefix, cap in (("modXNA-5'-fragment:", "5prime"),
                        ("modXNA-3'-fragment:", "3prime"),
                        ("modXNA-fragment:", None)):
        if line.startswith(prefix):
            body = line[len(prefix):]
            break
    else:
        return None
    toks = body.split(":")
    meta = {"name": toks[0], "cap": cap}
    for i in range(1, len(toks) - 1, 2):
        meta[toks[i]] = toks[i + 1]
    return meta


def read_mol2(path: str):
    """Devuelve (atomos, linea de titulo). Atomo = {name, xyz, type, q}."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        i = next(k for k, l in enumerate(lines) if l.startswith("@<TRIPOS>ATOM"))
    except StopIteration:
        raise ValueError(f"{path}: sin seccion ATOM")
    try:
        j = next(k for k, l in enumerate(lines) if k > i and l.startswith("@<TRIPOS>"))
    except StopIteration:
        j = len(lines)
    atoms = []
    for l in lines[i + 1:j]:
        p = l.split()
        if len(p) < 6:
            continue
        atoms.append({"name": p[1],
                      "xyz": [float(p[2]), float(p[3]), float(p[4])],
                      "type": p[5],
                      "q": float(p[8]) if len(p) > 8 else 0.0})
    return atoms, (lines[1].strip() if len(lines) > 1 else "")


def build(modxna_root: str) -> dict:
    dat = os.path.join(modxna_root, "dat")
    if not os.path.isdir(dat):
        raise SystemExit(f"no existe {dat!r}: pase --modxna con la raiz del arbol modXNA")

    fragments: List[dict] = []
    digest = hashlib.sha256()
    for moiety in MOIETIES:
        d = os.path.join(dat, f"lib_{moiety}")
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".mol2"):
                continue
            path = os.path.join(d, fn)
            with open(path, "rb") as fh:
                digest.update(fh.read())
            atoms, title = read_mol2(path)
            meta = parse_title(title)
            if meta is None:
                print(f"aviso: {fn} sin cabecera modXNA reconocible, se omite",
                      file=sys.stderr)
                continue
            code = fn[:-5]
            row = {
                "path": f"dat/lib_{moiety}/{fn}",
                "source": "dat",
                "moiety": moiety,
                "code": code,
                "meta": meta,
                "natoms": len(atoms),
                "qsum": round(sum(a["q"] for a in atoms), 4),
                "atom_names": [a["name"] for a in atoms],
                "thumb_png": f"web/lib_{moiety}/{code}-2D.png",
                "sdf": f"web/lib_{moiety}/{code}.col6_fixed.sdf",
            }
            if moiety == "sugar":
                row["sugar_geom"] = pseudorotation({a["name"]: a["xyz"] for a in atoms})
            elif moiety == "base":
                head = (meta.get("HEAD01") or "").lstrip("@")
                row["base_class"] = GLYCOSIDIC.get(head, "other")
            elif moiety == "backbone":
                row["has_sulfur"] = any(a["type"].split(".")[0] == "S" for a in atoms)
            fragments.append(row)

    counts = {m: sum(1 for f in fragments if f["moiety"] == m) for m in MOIETIES}
    sugars = [f for f in fragments if f["moiety"] == "sugar" and f.get("sugar_geom")]
    north = sum(1 for f in sugars if f["sugar_geom"]["domain"] == "North")
    return {
        "generated_from": os.path.abspath(modxna_root),
        "content_sha256": digest.hexdigest()[:32],
        "pucker_method": ("Altona-Sundaralingam via xnafiber.geometry.pseudorotation "
                          "(acepta anillos con O4' o S4')"),
        "n_fragments": len(fragments),
        "counts": counts,
        "sugar_pucker_summary": {"with_furanose": len(sugars), "north": north,
                                 "south": len(sugars) - north,
                                 "no_furanose": counts["sugar"] - len(sugars)},
        "fragments": fragments,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--modxna", required=True, help="raiz del arbol modXNA")
    ap.add_argument("-o", "--output",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "data", "modxna_catalog.json"))
    args = ap.parse_args()
    blob = build(args.modxna)
    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump(blob, fh, indent=1)
    print(f"{blob['n_fragments']} fragmentos -> {args.output}")
    print(f"  {blob['counts']}")
    print(f"  azucares: {blob['sugar_pucker_summary']}")


if __name__ == "__main__":
    main()
