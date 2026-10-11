"""
dsxnabuilder.build
==================

From two typed strands and a fiber model to an atomic duplex.

Column ``i`` of the alignment is base pair ``i`` of the model. Strand 1 (top
line) sits on the model's strand I and strand 2 on strand II, unless the
chemistry fits better the other way round (DNA:RNA hybrid models, whose two
strands differ): then strand 1 goes on strand II, the column order is
reversed and the result is turned by Rx(180) so that chain A still runs
5'->3' along +z.

Every nucleotide is placed in its slot frame by the screw (``fibers``); base
letters missing from the model come from the standard bases (base swap) and
O2' is added or removed for the strand chemistry. None of this moves the
sugar-phosphate atoms of the template, so the backbone stays as continuous as
in the fiber repeat.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .catalog import Catalog, Issue
from .fibers import FiberLibrary, FiberModel
from .geometry import rot_x
from .naming import ResidueNamer
from .nucleotide import PHOSPHATE_ATOMS, PURINES, pair_kind
from .seqlang import AliasTable, DuplexSpec, ResidueSpec, SequenceParseError, parse_duplex
from .structure import Atom, element_of, o3p_distances

MAX_COLUMNS = 200
CLASH_DISTANCE = 2.2


@dataclass
class BuiltResidue:
    chain: str
    resseq: int
    spec: ResidueSpec
    resname: str
    names: List[str]
    xyz: np.ndarray
    flags: dict = field(default_factory=dict)


@dataclass
class BuildResult:
    model: FiberModel
    duplex: DuplexSpec
    residues: List[BuiltResidue]
    issues: List[Issue]
    flipped: bool
    name_map: Dict[str, str]
    o3p: List[Tuple[str, int, float]]
    n_clashes: int = 0
    phase: int = 0

    @property
    def ok(self) -> bool:
        return not any(i.level == "error" for i in self.issues)

    def atoms(self) -> List[Atom]:
        out = []
        for r in self.residues:
            for n, p in zip(r.names, r.xyz):
                out.append(Atom(n, r.resname, r.resseq, r.chain, p, element_of(n)))
        return out

    def modified(self) -> List[BuiltResidue]:
        return [r for r in self.residues if r.flags.get("modxna")]


class BuildError(Exception):
    def __init__(self, issues: List[Issue]):
        super().__init__("; ".join(i.message for i in issues if i.level == "error"))
        self.issues = issues


def choose_orientation(model: FiberModel, duplex: DuplexSpec) -> bool:
    """True when strand 1 fits the model's strand II better (flip)."""
    chem_i, chem_ii = model.strand_chem("I"), model.strand_chem("II")
    if chem_i == chem_ii:
        return False
    top = [r.chem for r in duplex.top.residues]
    bot = [r.chem for r in duplex.bottom.residues]
    keep = sum(c == chem_i for c in top) + sum(c == chem_ii for c in bot)
    flip = sum(c == chem_ii for c in top) + sum(c == chem_i for c in bot)
    return flip > keep


def build_duplex(top: str, bottom: str, model_id: int, lib: FiberLibrary,
                 catalog: Catalog, aliases: AliasTable,
                 residue_library=None, dry_run: bool = False) -> BuildResult:
    """Build the duplex. Raises :class:`BuildError` with the issues on failure.

    ``residue_library`` provides the modXNA residues; without it any
    non-canonical residue is reported as ``NOT_IN_LIBRARY``. ``dry_run``
    checks that they exist without transplanting them (the browser validates
    before it has fetched the library files).
    """
    issues: List[Issue] = []
    try:
        model = lib.get(model_id)
    except KeyError as exc:
        raise BuildError([Issue("error", None, "UNKNOWN_MODEL", str(exc.args[0]),
                                {"model": model_id})])
    chem_defaults = (_default_chem(model.strand_chem("I")),
                     _default_chem(model.strand_chem("II")))
    try:
        duplex = parse_duplex(top, bottom, catalog, aliases, chem_defaults)
    except SequenceParseError as exc:
        raise BuildError([exc.issue])
    issues += _dedupe(duplex.issues)
    n = duplex.n_columns
    if n > MAX_COLUMNS:
        raise BuildError([Issue("error", None, "TOO_LONG",
                                f"at most {MAX_COLUMNS} base pairs", {"max": MAX_COLUMNS})])

    flipped = choose_orientation(model, duplex)
    slot_top, slot_bot = ("II", "I") if flipped else ("I", "II")
    phase = choose_phase(model, duplex, flipped)

    def bp_of(column: int) -> int:
        return ((n - 1 - column) if flipped else column) + phase
    if flipped:
        issues.append(Issue("info", None, "ORIENTATION_FLIP",
                            "strand 1 was placed on the model's strand II to match "
                            "its chemistry", {"model": model.id}))

    issues += sequence_issues(model, duplex, flipped, phase)

    # -- place nucleotides -------------------------------------------------
    namer = ResidueNamer()
    res_issues: List[Issue] = []
    placed: Dict[Tuple[str, int], Tuple[List[str], np.ndarray, dict]] = {}
    for strand_name, strand, slot in (("top", duplex.top, slot_top),
                                      ("bottom", duplex.bottom, slot_bot)):
        for r in strand.residues:
            bp = bp_of(r.column)
            names, xyz, flags = model.nucleotide(bp, slot, r.letter, r.chem, lib.standard)
            ss, _ = model.slot_strand(bp, slot)
            res_issues += residue_issues(model, ss, r, strand_name, catalog)
            placed[(strand_name, r.column)] = (names, xyz, flags)

    issues += _group(res_issues)

    # -- modXNA residues (phase 2) -----------------------------------------
    modx = _modxna_positions(duplex)
    for strand_name, strand in (("top", duplex.top), ("bottom", duplex.bottom)):
        for r in strand.residues:
            if (strand_name, r.column) not in modx:
                continue
            if residue_library is None or not residue_library.has(r):
                issues.append(Issue("error", r.column, "NOT_IN_LIBRARY",
                                    f"{r.token} ({'/'.join(r.triplet())}, {r.cap}) is "
                                    "not available in the residue library",
                                    {"token": r.token, "triplet": "/".join(r.triplet()),
                                     "cap": r.cap}, strand_name))
                continue
            names, xyz, flags = placed[(strand_name, r.column)]
            if dry_run:
                flags = dict(flags, modxna=True, key=residue_library.key_for(r))
            else:
                names, xyz, extra = residue_library.transplant(r, names, xyz)
                issues += extra.pop("issues", [])
                flags = dict(flags, modxna=True, **extra)
            placed[(strand_name, r.column)] = (names, xyz, flags)

    if any(i.level == "error" for i in issues):
        raise BuildError(issues)

    # -- chains, terminal phosphate, names ---------------------------------
    residues: List[BuiltResidue] = []
    for chain, strand_name, strand in (("A", "top", duplex.top), ("B", "bottom", duplex.bottom)):
        for seq, r in enumerate(strand.residues, start=1):
            names, xyz, flags = placed[(strand_name, r.column)]
            if r.cap == "5cap":
                keep = [i for i, nm in enumerate(names) if nm not in PHOSPHATE_ATOMS]
                names, xyz = [names[i] for i in keep], xyz[keep]
            resname = namer.name_for(r, (strand_name, r.column) in modx)
            residues.append(BuiltResidue(chain, seq, r, resname, names, xyz, flags))

    # -- orientation and centring -----------------------------------------
    allxyz = np.vstack([r.xyz for r in residues])
    R = rot_x(180.0) if flipped else np.eye(3)
    centre = np.array([0.0, 0.0, float(np.mean((R @ allxyz.T).T[:, 2]))])
    for r in residues:
        r.xyz = (R @ r.xyz.T).T - centre

    result = BuildResult(model, duplex, residues, issues, flipped, namer.mapping, [])
    result.phase = phase
    result.o3p = o3p_distances(result.atoms())
    bad = [(c, i, d) for c, i, d in result.o3p if not (1.45 <= d <= 1.75)]
    if bad:
        result.issues.append(Issue("warning", None, "BACKBONE_GAP",
                                   f"{len(bad)} O3'-P distance(s) outside 1.45-1.75 A",
                                   {"count": len(bad)}))
    result.n_clashes = count_clashes(residues)
    if result.n_clashes:
        result.issues.append(Issue("warning", None, "CLASHES",
                                   f"{result.n_clashes} heavy-atom contact(s) shorter than "
                                   f"{CLASH_DISTANCE} A between residues; minimize before MD",
                                   {"count": result.n_clashes, "cutoff": CLASH_DISTANCE}))
    return result


def choose_phase(model: FiberModel, duplex: DuplexSpec, flipped: bool) -> int:
    """Offset of the repeat (0..k-1) that best fits strand 1's purine/pyrimidine
    pattern, so that e.g. CGCG on Z-DNA (repeat GC) starts on the C slot."""
    if model.k == 1:
        return 0
    n = duplex.n_columns
    slot = "II" if flipped else "I"
    best = (None, 0)
    for p in range(model.k):
        miss = 0
        for a in duplex.top.residues:
            if a.is_abasic:
                continue
            bp = ((n - 1 - a.column) if flipped else a.column) + p
            ss = getattr(model.slots[bp % model.k], slot)
            miss += (ss.native in PURINES) != (a.letter in PURINES)
        if best[0] is None or miss < best[0]:
            best = (miss, p)
    return best[1]


def _group(issues: List[Issue]) -> List[Issue]:
    """One issue per (code, strand) for per-residue findings, listing columns."""
    groups: Dict[Tuple[str, Optional[str]], List[Issue]] = {}
    for i in issues:
        groups.setdefault((i.code, i.strand), []).append(i)
    out = []
    for (code, strand), items in groups.items():
        if len(items) == 1:
            out.append(items[0])
            continue
        cols = [i.position + 1 for i in items if i.position is not None]
        first = items[0]
        params = dict(first.params, columns=_short(cols), count=len(items))
        level = "error" if any(i.level == "error" for i in items) else first.level
        out.append(Issue(level, None, code,
                         f"{len(items)} residues (columns {_short(cols)}): {first.message}",
                         params, strand))
    return out


def _default_chem(c: str) -> str:
    return c if c in ("DNA", "RNA") else "DNA"


def _dedupe(issues: List[Issue]) -> List[Issue]:
    seen, out = set(), []
    for i in issues:
        key = (i.code, i.strand, i.position if i.code not in ("PS_STEREO_UNSPECIFIED",) else None)
        if key in seen:
            continue
        seen.add(key)
        out.append(i)
    return out


def _modxna_positions(duplex: DuplexSpec) -> set:
    """Positions that must be modXNA units.

    Every non-canonical residue, plus both terminal residues of any strand
    whose 5' or 3' end is a modXNA residue: Amber and modXNA split the
    terminal charge differently (-0.3079/-0.6921 vs -0.320348/-0.679652), and
    mixing them would leave the strand with a non-integer charge.
    """
    out = set()
    for name, strand in (("top", duplex.top), ("bottom", duplex.bottom)):
        mod = [r for r in strand.residues if not r.is_canonical and not r.is_abasic]
        out |= {(name, r.column) for r in mod}
        ends = (strand.residues[0], strand.residues[-1])
        if any(not r.is_canonical for r in ends):
            out |= {(name, r.column) for r in ends}
    return out


# --------------------------------------------------------------------------
# Issues
# --------------------------------------------------------------------------
def sequence_issues(model: FiberModel, duplex: DuplexSpec, flipped: bool,
                    phase: int = 0) -> List[Issue]:
    out: List[Issue] = []
    n = duplex.n_columns
    for i, (a, b) in enumerate(duplex.pairs()):
        if a is None or b is None:
            continue
        if a.is_abasic or b.is_abasic:
            other = b if a.is_abasic else a
            what = "abasic sites on both strands" if (a.is_abasic and b.is_abasic) else \
                f"abasic site opposite {other.letter}"
            out.append(Issue("info", i, "ABASIC",
                             f"column {i + 1}: {what} (sugar-phosphate kept, no base)",
                             {"column": i + 1, "opposite": "-" if a.is_abasic and b.is_abasic
                              else other.letter}))
            continue
        kind = pair_kind(a.letter, b.letter)
        if kind == "mismatch":
            out.append(Issue("warning", i, "MISMATCH",
                             f"column {i + 1}: {a.letter}:{b.letter} mismatch placed in the "
                             "model's base-pair frame; expect clashes, minimize",
                             {"column": i + 1, "pair": f"{a.letter}:{b.letter}"}))
        elif kind == "wobble":
            out.append(Issue("info", i, "WOBBLE",
                             f"column {i + 1}: {a.letter}:{b.letter} wobble pair built in a "
                             "Watson-Crick frame", {"column": i + 1,
                                                    "pair": f"{a.letter}:{b.letter}"}))
    over = [i for i, (a, b) in enumerate(duplex.pairs()) if a is None or b is None]
    if over:
        out.append(Issue("info", None, "OVERHANG",
                         f"{len(over)} unpaired column(s) at the ends (overhangs); unpaired "
                         "nucleotides keep the stacked fiber geometry", {"count": len(over)}))
    if model.sequence_policy == "swap":
        out.append(Issue("info", None, "MODEL_SEQUENCE",
                         f"model #{model.id} was derived for a specific sequence; other "
                         "bases are placed on its backbone geometry", {"model": model.id}))
    if model.sequence_policy == "pattern":
        slot_top = "II" if flipped else "I"
        wrong = []
        for a in duplex.top.residues:
            if a.is_abasic:
                continue
            bp = ((n - 1 - a.column) if flipped else a.column) + phase
            ss, _ = model.slot_strand(bp, slot_top)
            if (ss.native in PURINES) != (a.letter in PURINES):
                wrong.append(a.column + 1)
        if wrong:
            out.append(Issue("warning", None, "SEQ_PATTERN",
                             f"model #{model.id} repeats {_native_top(model, slot_top)}; "
                             f"columns {_short(wrong)} break its purine/pyrimidine pattern",
                             {"model": model.id, "columns": _short(wrong),
                              "pattern": _native_top(model, slot_top)}))
    return out


def residue_issues(model: FiberModel, ss, r: ResidueSpec, strand: str,
                   catalog: Catalog) -> List[Issue]:
    out: List[Issue] = []
    if r.is_abasic:
        return out
    if ss.syn and r.letter not in PURINES:
        out.append(Issue("warning", r.column, "SYN_PYRIMIDINE",
                         f"{r.token}: pyrimidine placed in a syn slot of model #{model.id}",
                         {"token": r.token, "model": model.id}, strand))
    dom = (ss.pucker or {}).get("domain")
    sugar_dom = catalog[r.sugar].pucker_domain if r.sugar in catalog else None
    # 2'-deoxyribose populates both domains (B: South, A: North), so only
    # 2'-substituted sugars are flagged.
    if dom and sugar_dom and dom != sugar_dom and r.sugar != "DC2":
        hard = r.sugar in ("LNA", "CET")
        code = "PUCKER_LOCKED" if hard else "PUCKER_FORCED"
        out.append(Issue("error" if hard else "warning", r.column, code,
                         f"{r.token}: {r.sugar} prefers {sugar_dom} but model #{model.id} "
                         f"has {dom} sugars" + ("; its bicyclic ring cannot adopt it" if hard
                                                else "; it takes the model's pucker"),
                         {"token": r.token, "sugar": r.sugar, "native": sugar_dom,
                          "model_pucker": dom, "model": model.id}, strand))
    return out


def _native_top(model: FiberModel, slot: str) -> str:
    return "".join(getattr(s, slot).native for s in model.slots)


def _short(cols: Sequence[int], limit: int = 8) -> str:
    s = ", ".join(str(c) for c in cols[:limit])
    return s + (", ..." if len(cols) > limit else "")


def count_clashes(residues: Sequence[BuiltResidue], cutoff: float = CLASH_DISTANCE) -> int:
    """Inter-residue heavy-atom contacts below ``cutoff`` (spatial hashing).

    Pairs bonded across the phosphodiester linkage (O3'/C3' of i with
    P/OP1/OP2/O5' of i+1) are not contacts.
    """
    pts, owner, names = [], [], []
    for k, r in enumerate(residues):
        for nm, p in zip(r.names, r.xyz):
            pts.append(p)
            owner.append(k)
            names.append(nm)
    if not pts:
        return 0
    pts = np.array(pts)
    cell = {}
    keys = np.floor(pts / cutoff).astype(int)
    for idx, key in enumerate(map(tuple, keys)):
        cell.setdefault(key, []).append(idx)
    linked = set()
    for k in range(len(residues) - 1):
        a, b = residues[k], residues[k + 1]
        if a.chain == b.chain:
            linked.add((k, k + 1))
    count = 0
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for key, members in cell.items():
        for off in offs:
            other = cell.get((key[0] + off[0], key[1] + off[1], key[2] + off[2]))
            if not other:
                continue
            for i in members:
                for j in other:
                    if j <= i or owner[i] == owner[j]:
                        continue
                    if np.sum((pts[i] - pts[j]) ** 2) >= cutoff * cutoff:
                        continue
                    lo, hi = sorted((owner[i], owner[j]))
                    if (lo, hi) in linked:
                        ni, nj = (names[i], names[j]) if owner[i] == lo else (names[j], names[i])
                        if ni in ("O3'", "C3'") and nj in ("P", "OP1", "OP2", "O5'"):
                            continue
                    count += 1
    return count
