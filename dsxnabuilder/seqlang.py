"""
dsxnabuilder.seqlang
====================

Sequence language for both strands of a duplex, written the way oligo
chemists write them, resolved to modXNA triplets (backbone, sugar, base).

The user types strand 1 5'->3' and, aligned under it, strand 2 3'->5'.
Column ``i`` of both lines is base pair ``i``; mismatches are allowed. A
``-`` at either end of a line leaves that column without a nucleotide
(overhang); a ``-`` between nucleotides is an abasic site: the strand keeps
its sugar-phosphate backbone there but has no base.

Grammar (one line)
------------------
::

    line     := wrap5? item* wrap3?
    item     := residue | linkage | group | '-'
    residue  := explicit | shorthand
    explicit := '{' BB '/' SUGAR '/' BASE (':' LETTER)? '}'   e.g. {PS1/MOE/RUU}
    shorthand:= sugar_prefix? LETTER ('^' BASECODE)?          e.g. +A, mU, dC, C^5mC
    linkage  := '*' | '*1' | '*2' | '*R' | '*S' | '.' | '{' BB '}'
                ('.' is an explicit phosphodiester)
    group    := ('d' | 'r')? '(' item* ')'                    chemistry of the group
    '-'      := no nucleotide at the ends of a line (overhang); between
                nucleotides, an abasic site (sugar-phosphate without base)
    wrap5    := "5'" '-'?        (top line)    |  "3'" '-'?  (bottom line)
    wrap3    := '-'? "3'"        (top line)    |  '-'? "5'"  (bottom line)

Ends are recognized only with the prime, so any other ``-`` is a gap.
Prefixes, linkages and base suffixes come from ``data/aliases.json``, the
pairing letters of modified bases from ``data/base_pairing.json``.

Backbone convention
-------------------
In modXNA the backbone fragment of a residue is its 5' phosphate, so a
linkage written between two residues belongs to the residue on its 3' side.
The bottom line is written 3'->5': it is tokenized as typed, the token list is
reversed into 5'->3' order and only then are linkages assigned. In the
bottom line as typed, a linkage therefore belongs to the residue on its left.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .catalog import CANONICAL_BACKBONE, CANONICAL_SUGAR, CAP_5_OH, Catalog, canonical_base
from .issues import Issue
from .paths import DATA


class SequenceParseError(ValueError):
    """A line that cannot be read; ``issue`` holds the error for the user."""

    def __init__(self, message: str, code: str = "PARSE", strand: Optional[str] = None,
                 position: Optional[int] = None, params: Optional[dict] = None):
        super().__init__(message)
        self.issue = Issue("error", position, code, message, params or {}, strand)


# --------------------------------------------------------------------------
# Alias table
# --------------------------------------------------------------------------
@dataclass
class AliasTable:
    """The notation: sugar prefixes (and the strand chemistry each implies),
    base suffixes and linkage symbols, mapped to modXNA codes; aliases whose
    assignment is not verified; and the pairing letter of modified bases."""

    sugar: Dict[str, str] = field(default_factory=dict)
    sugar_chem: Dict[str, str] = field(default_factory=dict)
    base: Dict[str, str] = field(default_factory=dict)
    linkage: Dict[str, str] = field(default_factory=dict)
    unverified: List[str] = field(default_factory=list)
    pairing: Dict[str, str] = field(default_factory=dict)

    @classmethod
    def load(cls, path: Optional[str] = None) -> "AliasTable":
        """aliases.json and, from the same directory, base_pairing.json."""
        path = path or os.path.join(DATA, "aliases.json")
        with open(path, encoding="utf-8") as fh:
            raw = json.load(fh)
        with open(os.path.join(os.path.dirname(path), "base_pairing.json"),
                  encoding="utf-8") as fh:
            pairing = json.load(fh).get("pairing", {})
        return cls(sugar=raw.get("sugar", {}), sugar_chem=raw.get("sugar_chem", {}),
                   base=raw.get("base", {}), linkage=raw.get("linkage", {}),
                   unverified=raw.get("unverified", []), pairing=dict(pairing))


# --------------------------------------------------------------------------
# Residues
# --------------------------------------------------------------------------
@dataclass
class ResidueSpec:
    """One nucleotide as written: its modXNA triplet and where it sits."""

    backbone: str
    sugar: str
    base: str
    letter: str                 # Watson-Crick letter A/C/G/T/U, or "-" for an abasic site
    chem: str                   # DNA | RNA (sugar family)
    token: str = ""
    column: int = -1            # alignment column (0-based)
    cap: str = "int"            # int | 5cap | 3cap
    explicit_backbone: bool = False

    @property
    def is_abasic(self) -> bool:
        return self.letter == "-"

    def triplet(self) -> Tuple[str, str, str]:
        return (self.backbone, self.sugar, self.base)

    @property
    def is_canonical(self) -> bool:
        """Standard Amber residue (DA.., A..): canonical sugar, base, phosphate."""
        if self.is_abasic:
            return False
        if self.sugar != CANONICAL_SUGAR[self.chem]:
            return False
        if self.base != canonical_base(self.chem, self.letter):
            return False
        if self.chem == "DNA" and self.letter == "U":
            return False                      # deoxy-U: no Amber unit
        if self.chem == "RNA" and self.letter == "T":
            return False                      # ribo-T: no Amber unit
        return self.backbone in (CANONICAL_BACKBONE[self.chem], CAP_5_OH)

    def key(self) -> str:
        """BACKBONE/SUGAR/BASE/cap: the residue's identity as users write it
        (residue names, the manifest). The residue library files use the
        same fields joined by "-" (reslib.residue_key)."""
        return "/".join(self.triplet()) + f"/{self.cap}"


@dataclass
class Strand:
    """One line resolved to 5'->3' order (gaps kept as ``None``)."""

    columns: List[Optional[ResidueSpec]]   # in the order typed on screen
    residues: List[ResidueSpec]            # 5'->3', no gaps
    raw: str


# --------------------------------------------------------------------------
# Tokenizer
# --------------------------------------------------------------------------
_EXPLICIT = re.compile(
    r"\{(?P<bb>[A-Za-z0-9]{3})/(?P<sug>[A-Za-z0-9]{3})/(?P<base>[A-Za-z0-9]{3})"
    r"(?::(?P<letter>[ACGTU]))?\}")
_EXPLICIT_BB = re.compile(r"\{(?P<bb>[A-Za-z0-9]{3})\}")
_PRIMES = "'\u2032\u2019"     # ', prime, right single quotation mark


@dataclass
class _Tok:
    kind: str              # res | link | gap
    pos: int               # character position in the cleaned line
    text: str = ""
    linkage: Optional[str] = None
    res: Optional[ResidueSpec] = None
    chem: str = "DNA"      # chemistry in effect (for '-' that become abasic sites)


def _strip_ends(text: str, top: bool) -> str:
    s = re.sub(r"\s+", "", text)
    for p in _PRIMES[1:]:
        s = s.replace(p, "'")
    first, last = ("5", "3") if top else ("3", "5")
    other_first, other_last = last, first
    m = re.match(r"^(\d)'-?", s)
    if m:
        if m.group(1) == other_first:
            raise SequenceParseError(
                f"line starts with {m.group(1)}' but this line is written "
                f"{first}'->{last}'", "WRONG_END_LABEL", "top" if top else "bottom")
        if m.group(1) == first:
            s = s[m.end():]
    m = re.search(r"-?(\d)'$", s)
    if m:
        if m.group(1) == other_last:
            raise SequenceParseError(
                f"line ends with {m.group(1)}' but this line is written "
                f"{first}'->{last}'", "WRONG_END_LABEL", "top" if top else "bottom")
        if m.group(1) == last:
            s = s[:m.start()]
    return s


def _tokenize(body: str, catalog: Catalog, aliases: AliasTable, default_chem: str,
              strand: str, issues: List[Issue]) -> List[_Tok]:
    toks: List[_Tok] = []
    stack = [default_chem]
    i = 0
    while i < len(body):
        ch = body[i]
        chem = stack[-1]
        if ch in "dr" and body[i + 1:i + 2] == "(":
            stack.append({"d": "DNA", "r": "RNA"}[ch])
            i += 2
            continue
        if ch == "(":
            stack.append(chem)
            i += 1
            continue
        if ch == ")":
            if len(stack) == 1:
                raise SequenceParseError(f"unbalanced ')' at character {i + 1}",
                                         "UNBALANCED", strand, i, {"char": i + 1})
            stack.pop()
            i += 1
            continue
        if ch == "-":
            toks.append(_Tok("gap", i, "-", chem=chem))
            i += 1
            continue
        if ch == "*":
            suf = body[i + 1:i + 2]
            tok = f"*{suf}" if suf in ("R", "S", "1", "2") else "*"
            code = aliases.linkage.get(tok)
            if tok == "*":
                issues.append(Issue("warning", None, "PS_STEREO_UNSPECIFIED",
                                    f"phosphorothioate without stereo label; {code} is used "
                                    "(write *1 or *2 to choose the retained oxygen)",
                                    {"fragment": code}, strand))
            elif tok in aliases.unverified:
                issues.append(Issue("warning", None, "UNVERIFIED_ALIAS",
                                    f"{tok} -> {code} is not verified against modXNA "
                                    "documentation; *1/*2 are unambiguous",
                                    {"alias": tok, "fragment": code}, strand))
            toks.append(_Tok("link", i, tok, linkage=code))
            i += len(tok)
            continue
        if ch == ".":
            toks.append(_Tok("link", i, ".", linkage=None))
            i += 1
            continue
        if ch == "{":
            mm = _EXPLICIT.match(body, i)
            if mm:
                bb, sug, base = mm.group("bb"), mm.group("sug"), mm.group("base")
                letter = mm.group("letter") or _letter_of(base, aliases, strand, i)
                # the sugar decides the chemistry, whatever group it is written in
                if sug in catalog:
                    rchem = "RNA" if catalog[sug].is_ribo_type else "DNA"
                else:
                    rchem = chem
                toks.append(_Tok("res", i, mm.group(0), res=ResidueSpec(
                    backbone=bb, sugar=sug, base=base, letter=letter, chem=rchem,
                    token=mm.group(0), explicit_backbone=True)))
                i = mm.end()
                continue
            mb = _EXPLICIT_BB.match(body, i)
            if mb:
                toks.append(_Tok("link", i, mb.group(0), linkage=mb.group("bb")))
                i = mb.end()
                continue
            raise SequenceParseError(f"malformed {{...}} block at character {i + 1}",
                                     "BAD_BLOCK", strand, i, {"char": i + 1})
        # shorthand residue
        prefix = ""
        for plen in (3, 2, 1):
            cand = body[i:i + plen]
            if cand in aliases.sugar and body[i + plen:i + plen + 1].upper() in tuple("ACGTU"):
                prefix = cand
                break
        start = i
        i += len(prefix)
        letter = body[i:i + 1].upper()
        if letter not in ("A", "C", "G", "T", "U"):
            raise SequenceParseError(f"unexpected character {body[i:i + 1]!r} at "
                                     f"character {i + 1}", "BAD_CHAR", strand, i,
                                     {"char": i + 1, "text": body[i:i + 1]})
        i += 1
        suffix = ""
        if body[i:i + 1] == "^":
            mm = re.match(r"\^([A-Za-z0-9]{3})", body[i:])
            if not mm:
                raise SequenceParseError(f"malformed base modifier at character {i + 1}",
                                         "BAD_BASE_MOD", strand, i, {"char": i + 1})
            suffix = mm.group(1)
            i += mm.end()
        rchem = aliases.sugar_chem.get(prefix, chem) if prefix else chem
        sugar = aliases.sugar.get(prefix, CANONICAL_SUGAR[rchem])
        if suffix:
            base = aliases.base.get(suffix, suffix)
            if base not in catalog or catalog[base].moiety != "base":
                raise SequenceParseError(f"unknown base modifier ^{suffix}",
                                         "UNKNOWN_BASE_MOD", strand, start,
                                         {"suffix": suffix})
        else:
            base = canonical_base(rchem, letter)
        toks.append(_Tok("res", start, body[start:i], res=ResidueSpec(
            backbone=CANONICAL_BACKBONE[rchem], sugar=sugar, base=base, letter=letter,
            chem=rchem, token=body[start:i])))
    if len(stack) != 1:
        raise SequenceParseError("unbalanced '(': missing ')'", "UNBALANCED", strand)
    return toks


def _letter_of(base: str, aliases: AliasTable, strand: str, pos: int) -> str:
    """Watson-Crick letter of a modXNA base, from the verified pairing table.

    Guessing the pairing of a base analogue from its three-letter code would
    build silently wrong duplexes, so other bases must state the letter.
    """
    if base in aliases.pairing:
        return aliases.pairing[base]
    raise SequenceParseError(
        f"base {base!r} has no verified Watson-Crick letter; write "
        f"{{BB/SUGAR/{base}:X}} with X = A, C, G, T or U",
        "UNKNOWN_PAIRING", strand, pos, {"base": base})


# --------------------------------------------------------------------------
# Lines and duplex
# --------------------------------------------------------------------------
def parse_line(text: str, top: bool, catalog: Catalog, aliases: AliasTable,
               default_chem: str = "DNA") -> Tuple[Strand, List[Issue]]:
    """One line: strand 1 (``top``, typed 5'->3') or strand 2 (typed 3'->5').

    ``default_chem`` is the chemistry of plain letters, that of the model
    strand the line is built on. Raises SequenceParseError.
    """
    strand = "top" if top else "bottom"
    issues: List[Issue] = []
    body = _strip_ends(text, top)
    if not body:
        raise SequenceParseError("empty strand", "EMPTY", strand)
    toks = _tokenize(body, catalog, aliases, default_chem, strand, issues)

    # '-' between nucleotides is an abasic site; at the ends, no nucleotide
    res_idx = [k for k, t in enumerate(toks) if t.kind == "res"]
    if not res_idx:
        raise SequenceParseError("strand has no nucleotides", "EMPTY", strand)
    for k in range(res_idx[0] + 1, res_idx[-1]):
        t = toks[k]
        if t.kind == "gap":
            t.kind = "res"
            t.res = ResidueSpec(backbone=CANONICAL_BACKBONE[t.chem],
                                sugar=CANONICAL_SUGAR[t.chem], base="ABS", letter="-",
                                chem=t.chem, token="-")

    # columns in typed order
    columns: List[Optional[ResidueSpec]] = []
    for t in toks:
        if t.kind == "res":
            t.res.column = len(columns)
            columns.append(t.res)
        elif t.kind == "gap":
            columns.append(None)

    # 5'->3' token order; linkages belong to the following residue
    order = toks if top else list(reversed(toks))
    residues: List[ResidueSpec] = []
    pending: Optional[_Tok] = None
    for t in order:
        if t.kind == "link":
            if pending is not None:
                issues.append(Issue("warning", None, "DOUBLE_LINKAGE",
                                    f"two linkages in a row near {t.text!r}; the last "
                                    "one is used", {"text": t.text}, strand))
            pending = t
        elif t.kind == "res":
            r = t.res
            if pending is not None:
                if not residues:
                    issues.append(Issue("warning", r.column, "LEADING_LINKAGE",
                                        "linkage before the 5'-terminal residue is ignored "
                                        "(5' end is 5'-OH)", {}, strand))
                elif pending.linkage is None:
                    if not r.explicit_backbone:
                        r.backbone = CANONICAL_BACKBONE[r.chem]
                else:
                    r.backbone = pending.linkage
                pending = None
            residues.append(r)
    if pending is not None:
        issues.append(Issue("warning", None, "DANGLING_LINKAGE",
                            "linkage after the 3'-terminal residue is ignored", {}, strand))
    if len(residues) < 2:
        raise SequenceParseError("a strand needs at least two nucleotides",
                                 "TOO_SHORT", strand)
    residues[0].cap = "5cap"
    if residues[0].explicit_backbone and residues[0].backbone not in (
            CAP_5_OH, CANONICAL_BACKBONE[residues[0].chem]):
        issues.append(Issue("warning", residues[0].column, "FIVE_PRIME_BACKBONE",
                            f"backbone {residues[0].backbone} of the 5'-terminal residue "
                            "is replaced by the 5'-OH cap", {"fragment": residues[0].backbone},
                            strand))
    residues[0].backbone = CAP_5_OH
    residues[-1].cap = "3cap"
    return Strand(columns=columns, residues=residues, raw=text), issues


@dataclass
class DuplexSpec:
    """Both strands, aligned column by column."""

    top: Strand
    bottom: Strand
    issues: List[Issue]

    @property
    def n_columns(self) -> int:
        return len(self.top.columns)

    def pairs(self) -> List[Tuple[Optional[ResidueSpec], Optional[ResidueSpec]]]:
        return list(zip(self.top.columns, self.bottom.columns))


def parse_duplex(top: str, bottom: str, catalog: Catalog, aliases: AliasTable,
                 default_chem: Tuple[str, str] = ("DNA", "DNA")) -> DuplexSpec:
    """Both lines, checked to align: same number of columns, none empty on
    both strands. ``default_chem`` holds the chemistry of plain letters on
    each strand. Raises SequenceParseError."""
    t, it = parse_line(top, True, catalog, aliases, default_chem[0])
    b, ib = parse_line(bottom, False, catalog, aliases, default_chem[1])
    issues = it + ib
    if len(t.columns) != len(b.columns):
        raise SequenceParseError(
            f"the two strands have {len(t.columns)} and {len(b.columns)} columns; "
            "align them (use '-' at the ends for overhangs)", "LENGTH_MISMATCH",
            None, None, {"top": len(t.columns), "bottom": len(b.columns)})
    for i, (x, y) in enumerate(zip(t.columns, b.columns)):
        if x is None and y is None:
            raise SequenceParseError(f"column {i + 1} is empty on both strands",
                                     "EMPTY_COLUMN", None, i, {"column": i + 1})
    return DuplexSpec(t, b, issues)
