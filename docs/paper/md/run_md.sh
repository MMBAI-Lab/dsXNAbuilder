#!/bin/bash
# MD use case of the dsXNAbuilder paper. This is what a user does with the
# structures dsXNAbuilder builds; none of it is part of the tool.
#
#   docs/paper/md/run_md.sh <system dir made by systems.py>
#
# 1. tleap from the Amber/modXNA package (leap.in: OL15, OL3, modXNA, OPC,
#    truncated octahedron 12 A, neutralizing Na+/Cl-), plus 150 mM NaCl added
#    with addIonsRand (count from the water count, 0.15 M x N_wat / 55.5).
# 2. Hydrogen mass repartitioning (cpptraj hmassrepartition, H mass 3.024) for 4 fs.
# 3. Minimization on CPU (restrained, then free), heating to 300 K (NVT, 200 ps),
#    NPT equilibration with solute restraints 5 -> 1 -> 0.1 kcal/mol/A2
#    (3 x 500 ps), and 100 ns NPT production (Langevin 1 ps-1, MC barostat,
#    PME, 9 A cutoff, SHAKE).
set -euo pipefail
D=$(cd "$1" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
set +u; source "${AMBERHOME:-/opt/amber26}/amber.sh"; set -u
export LC_ALL=C
cd "$D/dsxnabuilder"

# --- topology with 150 mM NaCl -------------------------------------------
tleap -f leap.in > leap0.log
NWAT=$(awk '$4=="WAT" && $3=="O"' system.pdb | wc -l)
NSALT=$(python3 -c "print(round(0.15*$NWAT/55.5))")
NSOLUTE=$(awk '/^ATOM/ && $4!="WAT" && $4!="Na+" && $4!="Cl-" {print substr($0,23,4)+0}' system.pdb | sort -n | tail -1)
sed "s#^saveamberparm mol system.parm7 system.rst7#addIonsRand mol Na+ $NSALT Cl- $NSALT\nsaveamberparm mol system.parm7 system.rst7#" leap.in > leap_salt.in
tleap -f leap_salt.in > leap.log
grep -q "Errors = 0" leap.log
echo "waters $NWAT, salt pairs $NSALT, solute residues $NSOLUTE" | tee "$D/setup.txt"

# HMR with cpptraj (the parmed of AmberTools 26 fails to import here)
CPPTRAJ_BIN=${CPPTRAJ_BIN:-$HERE/../../../external/cpptraj/bin/cpptraj}
"$CPPTRAJ_BIN" > hmr.log <<CP
parm system.parm7
hmassrepartition
parmwrite out system.hmr.parm7
CP

for f in "$HERE"/inputs/*.in; do
  sed "s#SOLUTE#:1-$NSOLUTE#" "$f" > "$(basename "$f")"
done

P=system.hmr.parm7
run() { pmemd.cuda -O -i $1.in -p $P -c $2 -r $1.rst7 -o $1.out ${3:-} ; }
# Minimization on CPU: with OPC (4-point water) the GPU minimizer drove two
# waters to O-O = 1.86 A in one system (LINMIN failure, NaN on heating);
# pmemd.MPI converges normally on the same input.
cpu() { mpirun -np ${NP_MIN:-8} pmemd.MPI -O -i $1.in -p $P -c $2 -r $1.rst7 -o $1.out ${3:-} ; }
cpu min1 system.rst7 "-ref system.rst7"
cpu min2 min1.rst7
run heat min2.rst7 "-ref min2.rst7"
run eq_5.0 heat.rst7 "-ref heat.rst7"
run eq_1.0 eq_5.0.rst7 "-ref heat.rst7"
run eq_0.1 eq_1.0.rst7 "-ref heat.rst7"
pmemd.cuda -O -i prod.in -p $P -c eq_0.1.rst7 -r prod.rst7 -o prod.out -x prod.nc
echo done > "$D/DONE"
