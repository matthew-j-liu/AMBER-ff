#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
LATEST_FILE="$RESULTS_DIR/latest_benchmark.txt"

# --load / -l: print the last saved benchmark results and exit
if [[ "${1:-}" == "--load" || "${1:-}" == "-l" ]]; then
  if [[ -f "$LATEST_FILE" ]]; then
    echo "=== Last saved benchmark results ==="
    cat "$LATEST_FILE"
  else
    echo "No saved results found. Run without --load to generate them."
  fi
  exit 0
fi

# AMBER setup
export AMBERHOME=/Users/jahnavigandhi/Desktop/chem_279/project/AMBER_benchmarking/ambertools25
export PATH=$AMBERHOME/bin:$PATH
export DYLD_LIBRARY_PATH=/usr/lib:$DYLD_LIBRARY_PATH
source $AMBERHOME/amber.sh

FUNCTIONAL_GROUPS=()
for _d in "$SCRIPT_DIR/input_molecules_copy"/*/; do
  [[ -d "$_d" ]] && FUNCTIONAL_GROUPS+=("$(basename "$_d")")
done
IFS=$'\n' FUNCTIONAL_GROUPS=($(sort <<<"${FUNCTIONAL_GROUPS[*]}")); unset IFS

run_molecule() {
  local INPUT="$1"
  local REL="${INPUT#input_molecules_copy/}"
  local SAFE="${REL%.xyz}"; SAFE="${SAFE// /_}"
  local MOLECULE; MOLECULE=$(basename "$SAFE")
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
  obabel "$INPUT" -O "$WORKDIR/${MOLECULE}.raw.mol2" || return 1

  # Step 2: antechamber (atom type assignment)
  echo "Step 2: antechamber"
  antechamber -i "$WORKDIR/${MOLECULE}.raw.mol2" \
              -fi mol2 \
              -o "$WORKDIR/${MOLECULE}.mol2" \
              -fo mol2 \
              -c bcc \
              -at gaff2 \
              -pf y \
              > "$WORKDIR/antechamber.log" 2>&1 || return 1

  # Step 3: parmchk2 (verify all required parameters are available)
  echo "Step 3: parmchk2"
  parmchk2 -i "$WORKDIR/${MOLECULE}.mol2" \
           -f mol2 \
           -o "$WORKDIR/${MOLECULE}.frcmod" || return 1

  # Step 4: tleap (generate topology and coordinate files)
  echo "Step 4: tleap"
  cat > "$WORKDIR/tleap.in" <<EOF
source leaprc.gaff2
mol = loadmol2 ${MOLECULE}.mol2
loadamberparams ${MOLECULE}.frcmod
saveamberparm mol ${MOLECULE}.prmtop ${MOLECULE}.inpcrd
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
  (cd "$WORKDIR" && sander -O -i sp.in -p "${MOLECULE}.prmtop" -c "${MOLECULE}.inpcrd" -o sp.out -r sp.rst) || return 1

  # Step 6: extract energy and save to results file
  echo "Step 6: energy"
  local energy_block
  energy_block=$(grep -A 9 "FINAL RESULTS" "$WORKDIR/sp.out" 2>/dev/null || echo "Energy not found (check $WORKDIR/sp.out)")
  echo "$energy_block"
  {
    echo ""
    echo "--- $MOLECULE ---"
    echo "$energy_block"
  } >> "$RESULTS_FILE"

  echo "=== DONE: $MOLECULE ==="
}

# Usage: ./run_AMBER_benchmarking.sh [xyz_path]
# With no argument, loops all molecules by functional group.
mkdir -p "$RESULTS_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULTS_FILE="$RESULTS_DIR/benchmark_${TIMESTAMP}.txt"
{
  echo "AMBER Benchmark Results — $(date)"
  echo "========================================"
} > "$RESULTS_FILE"

if [[ $# -gt 0 ]]; then
  run_molecule "$1"
else
  FAILED=()
  for group in "${FUNCTIONAL_GROUPS[@]}"; do
    echo ""
    echo "====== GROUP: $(echo "$group" | tr '[:lower:]' '[:upper:]') ======"
    echo "" >> "$RESULTS_FILE"
    echo "====== GROUP: $(echo "$group" | tr '[:lower:]' '[:upper:]') ======" >> "$RESULTS_FILE"
    for xyz in input_molecules_copy/${group}/*.xyz; do
      [[ -f "$xyz" ]] || continue
      run_molecule "$xyz" || FAILED+=("$xyz")
    done
  done

  echo ""
  echo "=== BENCHMARK COMPLETE ==="
  if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Failed molecules: ${FAILED[*]}"
    echo "" >> "$RESULTS_FILE"
    echo "FAILED: ${FAILED[*]}" >> "$RESULTS_FILE"
  else
    echo "All molecules completed successfully."
    echo "" >> "$RESULTS_FILE"
    echo "All molecules completed successfully." >> "$RESULTS_FILE"
  fi
fi

cp "$RESULTS_FILE" "$LATEST_FILE"
echo ""
echo "Results saved to: $RESULTS_FILE"
echo "Latest symlink:   $LATEST_FILE"

rm -f sqm.pdb sqm.out sqm.in
