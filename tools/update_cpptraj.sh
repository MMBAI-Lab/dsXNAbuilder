#!/bin/bash
# Clone or update cpptraj (github.com/Amber-MD/cpptraj, master) in
# external/cpptraj and build it against the AmberTools libraries.
#
#   tools/update_cpptraj.sh            # uses $AMBERHOME (default /opt/amber26)
#
# The commit built is written to external/cpptraj/COMMIT; the residue library
# records it, so rebuilt residues are traceable to the exact cpptraj.
# Periodic update (not installed automatically), e.g. weekly:
#   0 3 * * 1  cd /path/to/dsXNAbuilder && tools/update_cpptraj.sh >> build/cpptraj.log 2>&1
set -euo pipefail
cd "$(dirname "$0")/.."
AMBERHOME=${AMBERHOME:-/opt/amber26}
DEST=external/cpptraj

if [ -d "$DEST/.git" ]; then
  git -C "$DEST" fetch --quiet origin master
  if [ "$(git -C "$DEST" rev-parse HEAD)" = "$(git -C "$DEST" rev-parse origin/master)" ] \
     && [ -x "$DEST/bin/cpptraj" ]; then
    echo "cpptraj up to date at $(git -C "$DEST" rev-parse --short HEAD)"; exit 0
  fi
  git -C "$DEST" reset --quiet --hard origin/master
else
  mkdir -p external
  git clone --quiet https://github.com/Amber-MD/cpptraj.git "$DEST"
fi

cd "$DEST"
# shellcheck disable=SC1091
set +u; source "$AMBERHOME/amber.sh"; set -u
make clean >/dev/null 2>&1 || true
./configure -amberlib -nohdf5 -openmp gnu
make -j"$(nproc)" install
# -openmp installs bin/cpptraj.OMP; modxna.sh calls plain "cpptraj" from PATH
[ -x bin/cpptraj ] || ln -sf cpptraj.OMP bin/cpptraj
git rev-parse HEAD > COMMIT
./bin/cpptraj --version
