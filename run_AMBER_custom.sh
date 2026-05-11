#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

FUNCTIONAL_GROUPS=()
for _d in "$SCRIPT_DIR/input_molecules_copy"/*/; do
  [[ -d "$_d" ]] && FUNCTIONAL_GROUPS+=("$(basename "$_d")")
done
IFS=$'\n' FUNCTIONAL_GROUPS=($(sort <<<"${FUNCTIONAL_GROUPS[*]}")); unset IFS

usage() {
  echo "Usage: $0 [molecule] [params]"
  echo "  molecule  molecule name or full path to .xyz file (default: loops all benchmark molecules)"
  echo "  params    params filename or full path to .dat file (default: selected_atoms.dat)"
  exit 1
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage

run_molecule() {
  local MOLECULE_ARG="$1"
  local PARAMS_ARG="${2:-selected_atoms.dat}"

  local MOLECULE
  if [[ "$MOLECULE_ARG" == /* || "$MOLECULE_ARG" == *.xyz ]]; then
    MOLECULE="$MOLECULE_ARG"
  else
    local NAME="${MOLECULE_ARG%.xyz}"
    # Search subdirectories first, fall back to flat layout
    local found
    found=$(find "$SCRIPT_DIR/input_molecules_copy" -name "${NAME}.xyz" | head -1)
    MOLECULE="${found:-$SCRIPT_DIR/input_molecules_copy/${NAME}.xyz}"
  fi

  local PARAMS
  if [[ "$PARAMS_ARG" == /* ]]; then
    PARAMS="$PARAMS_ARG"
  else
    PARAMS="$SCRIPT_DIR/params/${PARAMS_ARG}"
  fi

  if [[ ! -f "$MOLECULE" ]]; then
    echo "Warning: molecule file not found: $MOLECULE, skipping."
    return 1
  fi

  if [[ ! -f "$PARAMS" ]]; then
    echo "Error: params file not found: $PARAMS"
    return 1
  fi

  echo ""
  echo "==> Running: amber_ff $(basename "$MOLECULE") with $(basename "$PARAMS")"
  "$BUILD_DIR/amber_ff" "$MOLECULE" "$PARAMS"
}

# Usage: ./run_AMBER_custom.sh [molecule] [params]
# With no arguments, loops all molecules by functional group using selected_atoms.dat.
if [[ $# -ge 1 ]]; then
  run_molecule "${1}" "${2:-selected_atoms.dat}"
else
  FAILED=()
  for group in "${FUNCTIONAL_GROUPS[@]}"; do
    echo ""
    echo "====== GROUP: $(echo "$group" | tr '[:lower:]' '[:upper:]') ======"
    for xyz in "$SCRIPT_DIR/input_molecules_copy/${group}"/*.xyz; do
      [[ -f "$xyz" ]] || continue
      run_molecule "$xyz" "selected_atoms.dat" || FAILED+=("$xyz")
    done
  done

  echo ""
  echo "=== LOOP COMPLETE ==="
  if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "Failed molecules: ${FAILED[*]}"
  else
    echo "All molecules completed successfully."
  fi
fi
