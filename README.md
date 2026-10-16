# dsXNAbuilder

[![tests](https://github.com/MMBAI-Lab/dsXNAbuilder/actions/workflows/test.yml/badge.svg)](https://github.com/MMBAI-Lab/dsXNAbuilder/actions/workflows/test.yml)

**Double-stranded DNA, RNA and DNA:RNA builder with modXNA chemistry.**

dsXNAbuilder builds duplexes in any of the regular helical forms of the
[3DNA](https://doi.org/10.1093/nar/gkg680) fiber library: A, B, C, D, Z, B′, S and B*-DNA,
A- and A′-RNA, and DNA:RNA hybrids. You can put [modXNA](https://modxna.chpc.utah.edu/)
sugar, backbone and base modifications on either strand. The output is a PDB file usable
with any modeling or simulation program; Amber users can also get a topology package with
the modXNA parameters.

**Web server: <https://danslab.xyz/dsxnabuilder/>** (no registration). It is a static web app:
Python runs in the browser through Pyodide, and nothing is sent to a server. The same Python
package works from the command line, for scripting and for batches of up to 500 duplexes.

## What it does

- **Both strands are typed explicitly.** Strand 1 is written 5′→3′ and strand 2 3′→5′ under
  it, aligned column by column. That makes it possible to modify either strand, build
  mismatches and leave overhangs (`-` at the ends).
- **40 fiber models from 3DNA** (Lu & Olson 2003, Table 4), derived directly from the 3DNA
  repeating units:
  - every nucleotide is placed by the screw symmetry of its model;
  - a built duplex reproduces the fiber model to within 1e-3 Å, with O3′–P at 1.60 Å and no
    relaxation;
  - generic models accept any sequence (marked ★ in the web app);
  - sequence-specific models (Z-DNA, D-DNA, poly(dA)·poly(dT) forms, …) keep their
    backbone, and bases are swapped in the standard reference frame. Warnings report when a
    sequence breaks the model's purine/pyrimidine or syn/anti pattern.
- **Modified nucleotides.** 244 modXNA residues are prebuilt with `modxna.sh`:
  - sugars: 2′-deoxyribose, ribose, 2′-OMe, 2′-F, 2′-MOE and LNA;
  - bases: A, C, G, T, U and 5mC on each sugar;
  - linkages: phosphodiester and the two phosphorothioates;
  - internal, 5′-terminal and 3′-terminal versions of each.

  Every combination of them is prebuilt, and the web app's **modXNA library** tab composes
  them by clicking. Modified residues are named X01, X02, … unless the user names them
  (`--name KEY=NAME` on the command line).

  They are *transplanted* onto the fiber scaffold. Shared atoms keep the fiber coordinates,
  so the helix and the backbone stay exact; only the modification atoms are placed from
  the modXNA geometry. LNA bridges are relaxed to their library geometry.
- **Output.** One self-describing PDB file (all atoms, hydrogens included; chain A = strand 1,
  chain B = strand 2). Modified nucleotides are named X01, X02, …; `REMARK 999` records list
  what each one is, the input strands, the model and the data versions. Hydrogens come from
  the residue template tleap will match (Amber DNA.OL15/RNA.lib units for canonical
  nucleotides, the modXNA residue for modified ones), placed on the built heavy atoms under
  the names tleap expects.
- **Many duplexes at once.** The **High-throughput** tab takes up to 50 duplexes in a
  FASTA-like text (a `>name=…; model=…` header, then strand 1 and strand 2) and downloads
  their PDB files in one zip, with a report. The command line takes up to 500
  (`--high-throughput FILE`).
- **Advanced: Amber/modXNA topology package.** It contains `structure.pdb`, `leap.in`
  (force field, modXNA libraries, optional water box), `libs/`, `frcmod.modxna`, `run.sh`
  (tleap only) and `manifest.json`. dsXNAbuilder is a builder: minimization and simulation
  are up to the user.

## Notation

| write | meaning |
|---|---|
| `ACGTU` | bases; plain letters take the chemistry of the model strand |
| `d(…)`, `r(…)` | DNA / RNA group |
| `m`, `f`, `e` | 2′-O-methyl, 2′-fluoro, 2′-O-methoxyethyl sugar prefixes (`mA`, `fU`, `eT`) |
| `+` (or `l`) | LNA, locked nucleic acid: 2′-O,4′-C-methylene-bridged ribose (`+G`, `+T`) |
| `*`, `*1`, `*2` | phosphorothioate (PS1 keeps OP2, PS2 keeps OP1); belongs to the residue on its 3′ side |
| `C^5mC` | base modification |
| `{PS1/MOE/RUU}` | explicit modXNA backbone/sugar/base triplet |
| `-` | at the ends of a line: no nucleotide (overhang); between nucleotides: abasic site (sugar-phosphate kept, no base; 2′-deoxy abasic = PDB component 3DR) |

Example, an LNA gapmer annealed to its RNA target on the hybrid model #22:

```
strand 1 (5'->3'):  r(GACUAGCUGA)
strand 2 (3'->5'):  +C*+T*+G*dA*dT*dC*dG*+A*+C*+T
```

## Command line

Requirements: Python 3.10 or later and NumPy. AmberTools is needed only to run the optional
topology package and the offline tools.

```bash
git clone https://github.com/MMBAI-Lab/dsXNAbuilder && cd dsXNAbuilder
python -m venv .venv && .venv/bin/pip install -e '.[dev]'
.venv/bin/dsxnabuilder models
.venv/bin/dsxnabuilder build --top "5'-CGCGAATTCGCG-3'" --bottom "3'-GCGCTTAAGCGC-5'" \
                             --model 4 -o dickerson.pdb
# advanced, for Amber users: topology package (tleap)
.venv/bin/dsxnabuilder build --top "5'-CGCGAATTCGCG-3'" --bottom "3'-GCGCTTAAGCGC-5'" \
                             --model 4 -o dickerson.pdb --amber dickerson.zip --water opc
unzip dickerson.zip && cd dsxnabuilder && ./run.sh     # tleap -> system.parm7/rst7
```

`dsxnabuilder build -h` lists the package options (force fields, water model, box, buffer,
counterions). The package's `manifest.json` holds the `reproduce` command that rebuilds it.

### Many duplexes at once

One record per duplex: a header line, then strand 1 (5′→3′) and strand 2 (3′→5′), written as
in the web app. `name` (letters, digits, `.`, `_`, `-`) names the PDB file; `model` is the
3DNA fiber model (B-DNA #4 when absent). A header field without `=` is the name. Blank lines
and lines starting with `#` are ignored.

```
>name=sirna_1; model=20
mG*mA*fC*mU*fA*mG*fC*mU*mG*mAdGdG
C*U*GAUCGAC*U--

>dickerson
CGCGAATTCGCG
GCGCTTAAGCGC
```

```bash
.venv/bin/dsxnabuilder build --high-throughput duplexes.txt -o structures.zip
```

Up to 500 duplexes per run (50 in the web app). The zip, or the directory given by `-o DIR`,
holds `pdb/<name>.pdb` (the same file a single build gives), `report.tsv`, `input.txt` and
`README.txt`. `--model N` sets the model of records without `model=`, and `-` reads the input
from standard input. A mistake in the format stops the run. A duplex that cannot be built is
reported in `report.tsv` while the others are still built, and the exit status is then 1.

## Tests

```bash
.venv/bin/python -m pytest -q     # tests that need external/3dna (tools/fetch_3dna.py) are skipped without it
```

The suite covers the geometry (checked against CPPTRAJ `nastruct`), the fiber models, the
sequence language, construction, the modXNA transplant, the hydrogens and the high-throughput
format. `tools/smoke_tleap.py` runs every library residue through tleap and sander, and
`scripts/e2e.py` drives the web app in headless Firefox.

## Repository

| path | content |
|---|---|
| `dsxnabuilder/` | the package (pure Python + numpy; the same code runs in Pyodide) |
| `data/` | `fiber_models.json` (derived 3DNA models), `residue_index.json`, `amber_residues.json` (Amber nucleotide units, for names and hydrogens), modXNA catalog, aliases |
| `web/` | static app (`index.html`, `app.js`, `worker.js`, `style.css`); `web/lib/` holds the prebuilt modXNA residues |
| `tools/` | offline tools (see below) |
| `tests/` | pytest suite |
| `scripts/` | web build, local server, end-to-end browser test, danslab deployment |
| `docs/paper/` | material of the paper: the figures and their scripts, the inputs, scripts and analysis of the use-case simulations, the browser benchmark and the built and simulated structures |

### Offline tools

| tool | role |
|---|---|
| `tools/fetch_3dna.py` | fetches the 3DNA fiber repeating units and standard bases from a pinned mirror of 3DNA (the forum downloads no longer work), with sha256 per file |
| `tools/derive_fiber_models.py` | derives `data/fiber_models.json`: screw direction from backbone closure, orientation normalization, OP1/OP2 chirality, pairing by H-bond contacts, per-model QC gate |
| `tools/update_cpptraj.sh` | builds cpptraj from the `Amber-MD/cpptraj` master branch against AmberTools |
| `tools/build_residue_library.py` | builds the modXNA residue library with `modxna.sh` (run under `LC_ALL=C`; see the note in the script) |
| `tools/smoke_tleap.py` | runs tleap (and sander) on every library residue and on random modified duplexes |
| `tools/build_data.py`, `tools/amber_residue_names.py`, `tools/build_catalog.py` | regenerate the JSON tables |
| `tools/fetch_catalog_names.py` | fragment names and notes from the modXNA catalogue page → `data/modxna_names.json` |

### Rebuilding everything

```bash
python tools/fetch_3dna.py
python tools/derive_fiber_models.py
tools/update_cpptraj.sh                       # AMBERHOME=/opt/amber26 by default
git clone https://github.com/modxna/modxna external/modxna
python tools/build_residue_library.py --jobs 8
python tools/smoke_tleap.py                   # 260 systems with AmberTools
python -m pytest -q
scripts/build_web.sh && python scripts/serve.py   # http://localhost:8766/
python scripts/e2e.py                          # headless Firefox
```

## Versions

Releases are tagged `vX.Y.Z` and listed in [CHANGELOG.md](CHANGELOG.md) and on the
[releases page](https://github.com/MMBAI-Lab/dsXNAbuilder/releases); between releases the
version carries a `-dev` suffix. `dsxnabuilder --version` prints it, every PDB file records it in
its header, and the web app shows it in its footer (with the commit when it is not a release).
Maintainers release with `scripts/release.sh X.Y.Z` (checks, tests, version, changelog, tag) and
then `scripts/release.sh X.Y.Z --publish` (push and GitHub release).

## Citation

A paper describing dsXNAbuilder is in preparation. Until it is published, please cite this
repository ([CITATION.cff](CITATION.cff); GitHub's "Cite this repository") and the works below.

## Credits and licenses

- Code: GPL-3.0-or-later (`LICENSE`).
- Fiber model data derive from 3DNA (CC BY-NC 4.0); modXNA residues and `frcmod.modxna` are
  GPL-3.0. Details are in `LICENSES/README.md`.
- Please cite:
  - Lu XJ, Olson WK (2003) *Nucleic Acids Res* 31:5108–5121;
  - Lu XJ, Olson WK (2008) *Nat Protoc* 3:1213–1227;
  - and, for modified nucleotides, Love O, Galindo-Murillo R, Roe DR, Dans PD, Cheatham TE III,
    Bergonzo C (2024) *J Chem Theory Comput* 20:9354–9363.

Dans Lab, Department of Biological Sciences, CENUR Litoral Norte, Universidad de la República.
