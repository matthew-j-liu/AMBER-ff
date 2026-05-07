#!/bin/bash

# AMBER setup
export AMBERHOME=/Users/jahnavigandhi/Desktop/chem_279/project/AMBER_benchmarking/ambertools25
export PATH=$AMBERHOME/bin:$PATH
export DYLD_LIBRARY_PATH=/usr/lib:$DYLD_LIBRARY_PATH
source $AMBERHOME/amber.sh

# 10 benchmark molecules — c, c2, c3, c1, hc, ha, ho, oh, cl only (no c6/ce)
BENCHMARK_MOLECULES=(
  "ethane"
  "propane"
  "n-butane"
  "n-pentane"
  "n-hexane"
)

run_molecule() {
  local MOLECULE="$1"
  # Sanitize name: replace spaces with underscores for all file/dir names
  local SAFE="${MOLECULE// /_}"
  local INPUT="input_molecules_copy/${MOLECULE}.xyz"
  local WORKDIR="amber_runs/${SAFE}"

  if [[ ! -f "$INPUT" ]]; then
    echo "Warning: '$INPUT' not found, skipping."
    return 1
  fi

  rm -rf "$WORKDIR"
  mkdir -p "$WORKDIR"

  echo ""
  echo "=== $(echo "$MOLECULE" | tr '[:lower:]' '[:upper:]') PIPELINE START ==="

  # Step 1: xyz -> mol2
  echo "Step 1: Open Babel"
  obabel "$INPUT" -O "$WORKDIR/${SAFE}.raw.mol2" || return 1

  # Step 2: antechamber (atom type assignment)
  echo "Step 2: antechamber"
  antechamber -i "$WORKDIR/${SAFE}.raw.mol2" \
              -fi mol2 \
              -o "$WORKDIR/${SAFE}.mol2" \
              -fo mol2 \
              -c bcc \
              -at gaff2 \
              -pf y \
              > "$WORKDIR/antechamber.log" 2>&1 || return 1

  # Step 3: parmchk2 (verify all required parameters are available)
  echo "Step 3: parmchk2"
  parmchk2 -i "$WORKDIR/${SAFE}.mol2" \
           -f mol2 \
           -o "$WORKDIR/${SAFE}.frcmod" || return 1

  # Step 4: tleap (generate topology and coordinate files)
  echo "Step 4: tleap"
  cat > "$WORKDIR/tleap.in" <<EOF
source leaprc.gaff2
mol = loadmol2 ${SAFE}.mol2
loadamberparams ${SAFE}.frcmod
saveamberparm mol ${SAFE}.prmtop ${SAFE}.inpcrd
quit
EOF
  (cd "$WORKDIR" && tleap -f tleap.in > tleap.log 2>&1) || return 1

  # Step 5: sander single point energy
  echo "Step 5: sander"
  cat > "$WORKDIR/sp.in" <<EOF
Single point
&cntrl
  imin=1,
  maxcyc=0,
  ntb=0,
  cut=999.0,
/
EOF
  (cd "$WORKDIR" && sander -O -i sp.in -p "${SAFE}.prmtop" -c "${SAFE}.inpcrd" -o sp.out -r sp.rst) || return 1

  # Step 6: extract energy
  echo "Step 6: energy"
  grep -A 9 "FINAL RESULTS" "$WORKDIR/sp.out" || echo "Energy not found (check $WORKDIR/sp.out)"

  echo "=== DONE: $MOLECULE ==="
}

# Usage: ./run_AMBER_benchmarking.sh [molecule]
# With no argument, loops all benchmark molecules.
if [[ $# -gt 0 ]]; then
  run_molecule "$1"
else
  FAILED=()
  for mol in "${BENCHMARK_MOLECULES[@]}"; do
    run_molecule "$mol" || FAILED+=("$mol")
  done

  echo ""
  echo "=== BENCHMARK COMPLETE ==="
  if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Failed molecules: ${FAILED[*]}"
  else
    echo "All molecules completed successfully."
  fi
fi

rm -f sqm.pdb sqm.out sqm.in
