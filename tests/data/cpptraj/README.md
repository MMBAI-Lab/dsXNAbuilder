# Valores de referencia de cpptraj `nastruct` (convenio 3DNA)

Dúplex construidos por dsXNAbuilder sobre ocho modelos de fibra de 3DNA (#1 A-DNA, #20 A-RNA,
#4 B-DNA, #7 C-DNA, #14 D-DNA, #22 híbrido, #49 S-DNA, #15 Z-DNA). Se les aplicó un movimiento
rígido arbitrario, porque cpptraj no detecta un par cuyo marco coincide con la identidad.
Después se analizaron con:

    parm X.pdb
    trajin X.pdb
    nastruct NA naout X.out
    run

- cpptraj: CPPTRAJ: Version V7.12.0 (GitHub) OpenMP, commit 27d70d8b4c.
- Generados el 2026-10-11.
- Los usa `tests/test_cpptraj_parameters.py`.
