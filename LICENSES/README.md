# Licenses by resource

| Resource | License | Origin |
|---|---|---|
| Code (`dsxnabuilder/`, `tools/`, `tests/`, `scripts/`, `web/*.js`, `web/*.html`, `web/*.css`) | GPL-3.0-or-later (`LICENSE`) | this project |
| `data/fiber_models.json` (templates and frames of the fiber models, standard bases) | CC BY-NC 4.0 | derived from the 3DNA data files (`fiber/strNN`, `config/Atomic_*.pdb`) |
| `data/modxna_catalog.json`, residue libraries (`web/lib/residues/*.lib`), `frcmod.modxna` | GPL-3.0 | derived from modXNA (github.com/modxna/modxna) |
| `data/amber_residues.json` | public domain, as the `dat/leap` files of AmberTools under its `LICENSE` | atom names, charges, coordinates and bonds of the nucleic-acid units of AmberTools 26 (`DNA.OL15.lib`, `RNA.lib`, with sha256); they give the hydrogens of canonical nucleotides |
| `data/modxna_names.json` | factual data (fragment names and notes), with attribution | modXNA catalogue page (modxna.chpc.utah.edu/catalog), with URL, date read and sha256 of the HTML (`tools/fetch_catalog_names.py`) |
| `data/modxna_catalog_review.json` | GPL-3.0 (this project) | issues of the modXNA catalogue found by panmer (`web-panmer` `data/modxna_catalog_review.csv` @ddda1bd) |
| Structures and packages built by users | no restrictions: they belong to the user | — |

## 3DNA

The fiber models come from 3DNA (Lu & Olson 2003, *Nucleic Acids Res* 31:5108–5121; Lu &
Olson 2008, *Nat Protoc* 3:1213–1227), whose v2.4 distribution is licensed CC BY-NC 4.0
(https://creativecommons.org/licenses/by-nc/4.0/). The coordinates of the repeating units
come in turn from Arnott (1999), Chandrasekaran & Arnott (1989), Alexeev et al. (1987),
van Dam & Levitt (2000) and the papers of Premilat & Albiser (1983–2001); the standard bases
from Olson et al. (2001). `data/fiber_models.json` re-expresses them (standard reference
frames, normalized orientation, Amber atom names) and is redistributed with attribution, for
non-commercial use.

The 3DNA files were obtained from a third-party mirror pinned to commit
`1363f200eaa75a7931817029f5f65af15bb0ac0b` (see `tools/fetch_3dna.py`), which is the source
of this project because the download links of the 3DNA forum no longer work. The license of
that mirror cannot relicense 3DNA, so the 3DNA terms apply.

## modXNA

Love O, Galindo-Murillo R, Roe DR, Dans PD, Cheatham TE III, Bergonzo C (2024)
*J Chem Theory Comput* 20:9354–9363, doi:10.1021/acs.jctc.4c01164. Repository under GPL-3.0.
