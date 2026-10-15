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

from .catalog import Catalog, default_chem
from .fibers import FiberLibrary, FiberModel
from .geometry import rot_x
from .graph import far_from
from .hydrogens import Template, add_hydrogens, amber_templates, amber_unit
from .issues import Issue
from .naming import ResidueNamer, check_names
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
    custom_names: Dict[str, str] = field(default_factory=dict)

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
                 residue_library=None, dry_run: bool = False,
                 names: Optional[Dict[str, str]] = None) -> BuildResult:
    """Build the duplex. Raises :class:`BuildError` with the issues on failure.

    ``residue_library`` provides the modXNA residues; without it any
    non-canonical residue is reported as ``NOT_IN_LIBRARY``. ``dry_run``
    checks that they exist without transplanting them (the browser validates
    before it has fetched the library files). ``names`` maps residue keys
    (``BACKBONE/SUGAR/BASE/cap``) to residue names chosen by the user.
    """
    issues: List[Issue] = []
    bad_names = check_names({k: v.strip().upper() for k, v in (names or {}).items()
                             if v and v.strip()})
    if bad_names:
        raise BuildError([Issue("error", None, "RESIDUE_NAME",
                                "residue names: " + "; ".join(bad_names),
                                {"detail": "; ".join(bad_names)})])
    try:
        model = lib.get(model_id)
    except KeyError as exc:
        raise BuildError([Issue("error", None, "UNKNOWN_MODEL", str(exc.args[0]),
                                {"model": model_id})])
    chem_defaults = (default_chem(model.strand_chem("I")),
                     default_chem(model.strand_chem("II")))
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
    namer = ResidueNamer(names)
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
                names, xyz, extra = residue_library.install(r, names, xyz, strand_name)
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

    if not dry_run:
        _add_hydrogens(residues, residue_library)
    if not dry_run and residue_library is not None:
        _orient_2o_chains(residues, residue_library)

    # -- orientation and centring -----------------------------------------
    allxyz = np.vstack([r.xyz for r in residues])
    R = rot_x(180.0) if flipped else np.eye(3)
    centre = np.array([0.0, 0.0, float(np.mean((R @ allxyz.T).T[:, 2]))])
    for r in residues:
        r.xyz = (R @ r.xyz.T).T - centre

    result = BuildResult(model, duplex, residues, issues, flipped, namer.mapping, [])
    result.phase = phase
    result.custom_names = namer.used_custom
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


CHAIN_CONTACT = 2.9        # heavy-atom distance a 2'-O substituent should keep
CHAIN_H_CONTACT = 2.0      # the same for pairs involving a hydrogen
CHAIN_STEP = 30.0          # torsion scan step (degrees)


def _orient_2o_chains(residues: Sequence[BuiltResidue], residue_library) -> None:
    """Choose the torsions of flexible 2'-O substituents (2'-OMe, MOE) in context.

    The transplant copies the substituent from the modXNA residue, minimized as
    an isolated nucleotide, so in a duplex the MOE chain can run into the next
    nucleotide (C6' 1.2 A from its C5'). The library conformation is kept when
    it makes no contact shorter than CHAIN_CONTACT; otherwise the chain
    torsions are scanned in CHAIN_STEP steps, two bonds at a time from O2'
    outwards, and the contact-free conformation closest to the library one is
    taken. Rotations are rigid about each bond, so bond lengths and angles stay
    those of the library. The LNA bridge closes a ring and is not touched.

    Hydrogens are already placed and turn with their chain atoms; contacts
    involving a hydrogen count below CHAIN_H_CONTACT, which also staggers the
    terminal methyl against the sugar. A bare 2'-OH is left as placed.
    """
    allxyz = np.vstack([r.xyz for r in residues])
    owner = np.concatenate([[k] * len(r.names) for k, r in enumerate(residues)])
    is_h = np.array([element_of(n) == "H" for r in residues for n in r.names])
    for k, r in enumerate(residues):
        if not r.flags.get("modxna") or "O2'" not in r.names:
            continue
        unit = residue_library.unit(r.flags["key"])
        present = set(r.names)
        unb = unit.neighbours()
        nb = {a: [b for b in unb.get(a, []) if b in present] for a in r.names}
        tree = _substituent_tree(nb)
        if not tree:
            continue
        parent, chain = tree
        if all(element_of(a) == "H" for a in chain):
            continue
        idx = {n: i for i, n in enumerate(r.names)}
        # rotatable bonds from O2' outwards and the atoms each one moves
        desc = {a: [b for b in chain if _ancestor(parent, b, a)] for a in ["O2'"] + chain}
        bonds = [(parent[a] if a != "O2'" else "C2'", a, desc[a]) for a in ["O2'"] + chain
                 if desc[a]]
        # environment: other residues near O2', and own atoms >= 5 bonds away from
        # each chain atom (1-5 contacts such as CM2...C4' at 2.9 A belong to the
        # minimized modXNA residue)
        o2 = r.xyz[idx["O2'"]]
        near = (owner != k) & (np.linalg.norm(allxyz - o2, axis=1) < 10.0)
        own = [n for n in r.names if n not in chain]
        env = np.vstack([allxyz[near], r.xyz[[idx[n] for n in own]]])
        env_h = np.concatenate([is_h[near], [element_of(n) == "H" for n in own]])
        far = {a: far_from(nb, [a], 5) for a in chain}
        counted = np.array([[True] * int(near.sum()) + [b in far[a] for b in own]
                            for a in chain])
        chain_h = np.array([element_of(a) == "H" for a in chain])
        cut = np.where(chain_h[:, None] | env_h[None, :], CHAIN_H_CONTACT, CHAIN_CONTACT)
        cut = np.where(counted, cut, 0.0)
        chain_far = [(i, j, CHAIN_H_CONTACT if chain_h[i] or chain_h[j] else CHAIN_CONTACT)
                     for i, a in enumerate(chain) for j, b in enumerate(chain)
                     if i < j and b in far[a]]
        start = {a: r.xyz[idx[a]].copy() for a in r.names}

        def build(angles):
            pos = dict(start)
            for (u, v, moving), ang in zip(bonds, angles):
                if ang:
                    pts = np.array([pos[m] for m in moving])
                    pts = _rotate(pts, pos[v], pos[v] - pos[u], ang)
                    for m, p in zip(moving, pts):
                        pos[m] = p
            return pos

        def penalty(pos):
            c = np.array([pos[a] for a in chain])
            d = np.linalg.norm(c[:, None, :] - env[None, :, :], axis=2)
            p = float(np.sum(np.clip(cut - d, 0.0, None) ** 2))
            for i, j, lim in chain_far:
                p += max(0.0, lim - float(np.linalg.norm(c[i] - c[j]))) ** 2
            return p

        best = [0.0] * len(bonds)
        best_p = penalty(build(best))
        if best_p == 0.0:
            continue
        grid = [CHAIN_STEP * i for i in range(int(round(360.0 / CHAIN_STEP)))]
        cost = lambda a: sum(min(x % 360.0, 360.0 - x % 360.0) for x in a)
        for _ in range(2):
            for s in range(0, len(bonds), 2):
                pair = range(s, min(s + 2, len(bonds)))
                cands = [[]]
                for _j in pair:
                    cands = [c + [g] for c in cands for g in grid]
                # closest to the library first (stable, so ties keep the grid order):
                # once a contact-free candidate is reached no later one can improve,
                # which ends most scans after a few trials (a staggered methyl)
                cands.sort(key=cost)
                for c in cands:
                    trial = list(best)
                    for j, g in zip(pair, c):
                        trial[j] = g
                    if best_p <= 1e-9 and cost(trial) >= cost(best):
                        break
                    p = penalty(build(trial))
                    if p < best_p - 1e-9 or (abs(p - best_p) <= 1e-9 and cost(trial) < cost(best)):
                        best, best_p = trial, p
        pos = build(best)
        for a in chain:
            r.xyz[idx[a]] = pos[a]
        r.flags["chain_torsions"] = [round(x, 1) for x in best]


def _add_hydrogens(residues: Sequence[BuiltResidue], residue_library) -> None:
    """Hydrogens from the template tleap will match each residue to (see
    dsxnabuilder.hydrogens). Runs before the 2'-O chains are oriented, so that
    the torsion scan sees them and they turn with the chain."""
    amber = amber_templates()
    modx: Dict[str, Template] = {}
    for r in residues:
        abasic = bool(r.flags.get("abasic"))
        if r.flags.get("modxna"):
            key = r.flags["key"]
            if key not in modx:
                modx[key] = Template.from_off(residue_library.unit(key))
            tpl = modx[key]
        else:
            # an abasic site keeps the sugar-phosphate of its strand chemistry;
            # any canonical unit of that chemistry provides its hydrogens
            unit = amber_unit("A" if abasic else r.spec.letter, r.spec.chem, r.spec.cap)
            if unit not in amber:
                raise ValueError(f"no Amber unit {unit} for the hydrogens of {r.spec.token}")
            tpl = amber[unit]
        r.names, r.xyz = add_hydrogens(r.names, r.xyz, tpl, abasic)


def _substituent_tree(nb: Dict[str, List[str]]):
    """Atoms hanging from O2' (not through C2') as a tree: (parent, atoms in BFS
    order), or None when there are none or they close a ring (LNA/cEt)."""
    parent = {"O2'": None}
    order: List[str] = []
    queue = ["O2'"]
    while queue:
        a = queue.pop(0)
        for b in nb.get(a, []):
            if a == "O2'" and b == "C2'":
                continue
            if b in parent:
                if parent[a] != b:
                    return None
                continue
            parent[b] = a
            order.append(b)
            queue.append(b)
    return (parent, order) if order else None


def _ancestor(parent: Dict[str, Optional[str]], atom: str, anc: str) -> bool:
    a = parent.get(atom)
    while a is not None:
        if a == anc:
            return True
        a = parent.get(a)
    return False


def _rotate(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, deg: float) -> np.ndarray:
    """Rotate points by ``deg`` about the axis through ``origin`` (Rodrigues)."""
    k = axis / np.linalg.norm(axis)
    t = math.radians(deg)
    p = points - origin
    return (origin + p * math.cos(t) + np.cross(k, p) * math.sin(t)
            + np.outer(p @ k, k) * (1.0 - math.cos(t)))


def count_clashes(residues: Sequence[BuiltResidue], cutoff: float = CLASH_DISTANCE) -> int:
    """Inter-residue heavy-atom contacts below ``cutoff`` (spatial hashing).

    Hydrogens are left out: they are placed from the residue templates after
    the heavy atoms and the 2.2 A criterion is for heavy atoms.

    Pairs bonded across the phosphodiester linkage (O3'/C3' of i with
    P/OP1/OP2/O5' of i+1) are not contacts.
    """
    pts, owner, names = [], [], []
    for k, r in enumerate(residues):
        for nm, p in zip(r.names, r.xyz):
            if element_of(nm) == "H":
                continue
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
