"""
xnafiber.catalog
================

Modelo de datos del catalogo de fragmentos modXNA y reglas de compatibilidad.

Cada fragmento modXNA es un fichero mol2 cuya linea de titulo codifica los
puntos de conexion, con la gramatica

    modXNA[-5'|-3']-fragment:<NOMBRE>:HEAD01:@<at>:HEAD01STRIP:@<at,...>
                             [:TAIL01:@<at>:TAIL01STRIP:@<at,...>]
                             [:HEAD03:@<at>:HEAD03STRIP:@<at,...>]

Semantica de los puntos de conexion (verificada sobre el arbol
github.com/modxna/modxna):

  backbone : HEAD01 = @P    (enlaza con O3' del residuo 5' precedente)
             TAIL01 = @O5'  (enlaza con C5' del azucar)
  sugar    : HEAD01 = @C5'  (enlaza con TAIL01 del backbone)
             TAIL01 = @O3'  (enlaza con P del residuo 3' siguiente)
             HEAD03 = @C1'  (enlaza con la base)
  base     : HEAD01 = @N9 (purinas) / @N1 (pirimidinas)

Los conjuntos ``*STRIP`` son los atomos de bloqueo (capping) que modXNA
elimina al fusionar los fragmentos; son exactamente los que ``cpptraj``
descarta en la operacion ``sequence`` tras ``dataset connect``.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

Moiety = str  # "backbone" | "sugar" | "base"

CANONICAL: Dict[Tuple[str, str], str] = {
    # (quimica, letra) -> codigo de base modXNA
    ("DNA", "A"): "DAA", ("DNA", "C"): "DCC", ("DNA", "G"): "DGG", ("DNA", "T"): "DTT",
    ("RNA", "A"): "RAA", ("RNA", "C"): "RCC", ("RNA", "G"): "RGG", ("RNA", "U"): "RUU",
}
CANONICAL_SUGAR = {"DNA": "DC2", "RNA": "RC3"}      # 2'-desoxirribosa S / ribosa N
CANONICAL_BACKBONE = {"DNA": "DPO", "RNA": "RPO"}   # fosfodiester
CAP_5_OH = "5PO"

WATSON_CRICK = {"A": "T", "T": "A", "G": "C", "C": "G", "U": "A"}
WATSON_CRICK_RNA = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"}


def canonical_base(chem: str, letter: str) -> str:
    """Codigo de base canonica para una quimica de azucar y una letra.

    La identidad de la base es independiente del azucar: el prefijo D/R de los
    codigos modXNA indica el contexto (desoxi o ribo) en el que se derivaron
    las cargas, no una base distinta. Timina solo existe como ``DTT`` y uracilo
    solo como ``RUU``, de modo que una T sobre un azucar de tipo ribo (p. ej.
    LNA-T, 2'-OMe-T) resuelve a ``DTT`` y una U sobre 2'-desoxirribosa resuelve
    a ``RUU``. Es el comportamiento correcto para los ASO reales, donde ambas
    combinaciones son habituales.
    """
    code = CANONICAL.get((chem, letter))
    if code is not None:
        return code
    other = "DNA" if chem == "RNA" else "RNA"
    code = CANONICAL.get((other, letter))
    if code is None:
        raise KeyError(f"not a canonical base letter: {letter!r}")
    return code


@dataclass
class Fragment:
    code: str
    moiety: Moiety
    path: str
    head01: Optional[str] = None
    head01_strip: Optional[str] = None
    tail01: Optional[str] = None
    tail01_strip: Optional[str] = None
    head03: Optional[str] = None
    head03_strip: Optional[str] = None
    cap: Optional[str] = None            # "5prime" | "3prime" | None
    natoms: int = 0
    qsum: float = 0.0
    atom_names: List[str] = field(default_factory=list)
    sugar_geom: Optional[dict] = None    # P, nu_max, dominio, envelope
    base_class: Optional[str] = None     # purine | pyrimidine | other
    has_sulfur: bool = False
    thumb_png: str = ""
    sdf: str = ""

    # -- utilidades ------------------------------------------------------
    @property
    def pucker_domain(self) -> Optional[str]:
        return (self.sugar_geom or {}).get("domain")

    @property
    def glycosidic_atom(self) -> Optional[str]:
        """Atomo de la base que forma el enlace glicosidico (N9 o N1)."""
        return (self.head01 or "").lstrip("@") or None

    def strip_atoms(self, which: str) -> List[str]:
        raw = getattr(self, f"{which}_strip") or ""
        raw = raw.lstrip("@")
        if not raw or raw.upper() == "NONE":
            return []
        return [a for a in raw.split(",") if a]


class Catalog:
    """Indice en memoria del catalogo modXNA."""

    def __init__(self, fragments: Sequence[Fragment], meta: Optional[dict] = None):
        self.fragments = list(fragments)
        self.meta = meta or {}
        self._by_code: Dict[str, Fragment] = {f.code: f for f in self.fragments}

    # -- constructores ---------------------------------------------------
    @classmethod
    def from_json(cls, path: Optional[str] = None) -> "Catalog":
        if path is None:
            path = os.path.join(os.path.dirname(__file__), "..", "data",
                                "modxna_catalog.json")
        with open(path, "r", encoding="utf-8") as fh:
            blob = json.load(fh)
        frags = []
        for row in blob["fragments"]:
            m = row.get("meta") or {}
            frags.append(Fragment(
                code=row["code"], moiety=row["moiety"], path=row["path"],
                head01=m.get("HEAD01"), head01_strip=m.get("HEAD01STRIP"),
                tail01=m.get("TAIL01"), tail01_strip=m.get("TAIL01STRIP"),
                head03=m.get("HEAD03"), head03_strip=m.get("HEAD03STRIP"),
                cap=m.get("cap"), natoms=row.get("natoms", 0),
                qsum=row.get("qsum", 0.0), atom_names=row.get("atom_names", []),
                sugar_geom=row.get("sugar_geom"), base_class=row.get("base_class"),
                has_sulfur=row.get("has_sulfur", False),
                thumb_png=row.get("thumb_png", ""), sdf=row.get("sdf", "")))
        return cls(frags, {k: v for k, v in blob.items() if k != "fragments"})

    # -- consultas -------------------------------------------------------
    def __getitem__(self, code: str) -> Fragment:
        try:
            return self._by_code[code]
        except KeyError:
            raise KeyError(f"unknown modXNA fragment: {code!r}") from None

    def __contains__(self, code: str) -> bool:
        return code in self._by_code

    def by_moiety(self, moiety: Moiety) -> List[Fragment]:
        return [f for f in self.fragments if f.moiety == moiety]

    def search(self, *, moiety: Optional[Moiety] = None,
               pucker: Optional[str] = None,
               base_class: Optional[str] = None,
               sulfur: Optional[bool] = None,
               text: Optional[str] = None) -> List[Fragment]:
        out = self.fragments
        if moiety:
            out = [f for f in out if f.moiety == moiety]
        if pucker:
            out = [f for f in out if f.pucker_domain == pucker]
        if base_class:
            out = [f for f in out if f.base_class == base_class]
        if sulfur is not None:
            out = [f for f in out if f.has_sulfur == sulfur]
        if text:
            t = text.lower()
            out = [f for f in out if t in f.code.lower()]
        return out

    def stereo_partner(self, code: str) -> Optional[str]:
        """Devuelve el diastereomero de P asociado (sufijos 1/2), si existe.

        En el catalogo modXNA los backbones modificados en el fosforo vienen
        en pares ``XX1``/``XX2`` que difieren en cual de los oxigenos
        proquirales (OP1/OP2) porta la modificacion. La asignacion absoluta
        Rp/Sp debe tomarse de la documentacion de modXNA; aqui solo se expone
        la relacion de pareja y el oxigeno retenido.
        """
        if len(code) == 3 and code[2] in "12":
            partner = code[:2] + ("2" if code[2] == "1" else "1")
            return partner if partner in self._by_code else None
        return None

    def retained_pro_oxygen(self, code: str) -> Optional[str]:
        f = self._by_code.get(code)
        if f is None or f.moiety != "backbone":
            return None
        for name in ("OP1", "OP2"):
            if name in f.atom_names:
                return name
        return None


# --------------------------------------------------------------------------
# Reglas de compatibilidad
# --------------------------------------------------------------------------
@dataclass
class Issue:
    """A user-facing finding. ``code`` is stable and drives the UI translation
    (web/i18n); ``message`` is the English fallback; ``params`` fill the
    placeholders of the translated message."""

    level: str   # "error" | "warning" | "info"
    position: Optional[int]
    code: str
    message: str
    params: Dict[str, object] = field(default_factory=dict)
    strand: Optional[str] = None     # "top" | "bottom" | None

    def to_json(self) -> dict:
        return {"level": self.level, "code": self.code, "strand": self.strand,
                "position": self.position, "message": self.message,
                "params": self.params}


def check_triplet(cat: Catalog, backbone: str, sugar: str, base: str,
                  form_pucker: str, position: Optional[int] = None,
                  pucker_policy: str = "fragment-native") -> List[Issue]:
    """Valida un triplete (backbone, sugar, base) contra una conformacion.

    ``pucker_policy``
        ``fragment-native``  acepta el pucker nativo del fragmento y lo
                             reporta (comportamiento por defecto: las
                             modificaciones 2' de los ASO *imponen* su propio
                             pucker y ese es el resultado deseado).
        ``match-form``       exige repuckerizar al dominio de la conformacion
                             mediante minimizacion restringida.
        ``strict``           rechaza cualquier fragmento cuyo dominio de
                             pucker no coincida con el de la conformacion.
    """
    issues: List[Issue] = []
    for code, moiety in ((backbone, "backbone"), (sugar, "sugar"), (base, "base")):
        if code not in cat:
            issues.append(Issue("error", position, "UNKNOWN_FRAGMENT",
                                f"{moiety} {code!r} is not in the modXNA catalog",
                                {"moiety": moiety, "fragment": code}))
            return issues
        if cat[code].moiety != moiety:
            issues.append(Issue("error", position, "WRONG_MOIETY",
                                f"{code!r} is a {cat[code].moiety} fragment, "
                                f"expected {moiety}",
                                {"fragment": code, "moiety": moiety}))
    if issues:
        return issues

    sug = cat[sugar]
    dom = sug.pucker_domain
    if dom is None:
        issues.append(Issue("warning", position, "NO_FURANOSE",
                            f"sugar {sugar!r} has no furanose ring (acyclic or "
                            "morpholino)", {"sugar": sugar}))
    elif dom != form_pucker:
        if pucker_policy == "strict":
            issues.append(Issue("error", position, "PUCKER_MISMATCH",
                                f"{sugar!r} is {dom} and the model requires "
                                f"{form_pucker}",
                                {"sugar": sugar, "native": dom, "model": form_pucker}))
        elif pucker_policy == "match-form":
            issues.append(Issue("info", position, "REPUCKER",
                                f"{sugar!r} ({dom}, P={sug.sugar_geom['P']} deg) "
                                f"takes the {form_pucker} pucker of the model",
                                {"sugar": sugar, "native": dom, "model": form_pucker}))
        else:
            issues.append(Issue("warning", position, "PUCKER_NATIVE",
                                f"{sugar!r} keeps its native {dom} pucker "
                                f"(P={sug.sugar_geom['P']} deg) in a "
                                f"{form_pucker} model",
                                {"sugar": sugar, "native": dom, "model": form_pucker}))

    bb = cat[backbone]
    if bb.cap is not None:
        issues.append(Issue("info", position, "CAP_FRAGMENT",
                            f"{backbone!r} is a {bb.cap} terminal fragment",
                            {"fragment": backbone, "cap": bb.cap}))
    partner = cat.stereo_partner(backbone)
    if partner:
        issues.append(Issue("info", position, "P_STEREOISOMER",
                            f"{backbone!r} and {partner!r} are P diastereomers "
                            f"(retained pro-chiral oxygen: "
                            f"{cat.retained_pro_oxygen(backbone)})",
                            {"fragment": backbone, "partner": partner}))
    return issues


def estimated_residue_charge(cat: Catalog, backbone: str, sugar: str, base: str) -> float:
    """Carga formal estimada del residuo ensamblado.

    Es solo una estimacion a partir de la suma de cargas de los fragmentos;
    la carga autoritativa es la del fichero ``.lib`` generado por modXNA tras
    la eliminacion de los grupos de bloqueo y las correcciones de union.
    """
    return round(sum(cat[c].qsum for c in (backbone, sugar, base)), 3)
