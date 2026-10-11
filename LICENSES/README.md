# Licencias por recurso

| Recurso | Licencia | Origen |
|---|---|---|
| Código (`dsxnabuilder/`, `tools/`, `tests/`, `web/*.js`, `web/*.html`, `web/*.css`, `scripts/`) | GPL-3.0-or-later (`LICENSE`) | este proyecto |
| `data/fiber_models.json` (plantillas y marcos de los modelos de fibra, bases estándar) | CC BY-NC 4.0 | derivado de los ficheros de datos de 3DNA (`fiber/strNN`, `config/Atomic_*.pdb`) |
| `data/modxna_catalog.json`, librerías de residuos (`web/lib/residues/*.lib`), `frcmod.modxna` | GPL-3.0 | derivado de modXNA (github.com/modxna/modxna) |
| `data/amber_residues.json` | nombres de átomo y cargas extraídos de AmberTools (`DNA.OL15.lib`, `RNA.lib`) | AmberTools |
| Estructuras y paquetes generados por los usuarios | sin restricciones: son del usuario | — |

## 3DNA

Los modelos de fibra proceden de 3DNA (Lu & Olson 2003, *Nucleic Acids Res* 31:5108–5121;
Lu & Olson 2008, *Nat Protoc* 3:1213–1227), cuya distribución v2.4 se licencia
CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/). Las coordenadas de las
unidades repetitivas tienen a su vez su origen en Arnott (1999), Chandrasekaran & Arnott
(1990), Alexeev et al. (1987), van Dam & Levitt (2000) y Premilat & Albiser (2001); las
bases estándar, en Olson et al. (2001). `data/fiber_models.json` las reexpresa (marcos de
referencia estándar, normalización de orientación, nombres de Amber) y se redistribuye con
atribución y para uso no comercial. Pendiente: confirmar con X.-J. Lu la redistribución.

Los ficheros de 3DNA se obtuvieron de un espejo de terceros fijado al commit
`1363f200eaa75a7931817029f5f65af15bb0ac0b` (ver `tools/fetch_3dna.py`); la licencia de ese
espejo no puede relicenciar 3DNA, así que rige la de 3DNA.

## modXNA

Love O, Galindo-Murillo R, Roe DR, Dans PD, Cheatham TE III, Bergonzo C (2024)
*J Chem Theory Comput* 20:9354–9363, doi:10.1021/acs.jctc.4c01164. Repositorio bajo GPL-3.0.
