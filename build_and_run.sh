#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

usage() {
  echo "Usage: $0 [molecule] [params]"
  echo "  molecule  molecule name (e.g. ethane) or full path to .xyz file (default: ethane)"
  echo "  params    params filename (e.g. gaff2.dat) or full path to .dat file (default: gaff2.dat)"
  exit 1
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage

# Resolve molecule: accept name, name+extension, or full path
MOLECULE_ARG="${1:-ethane}"
if [[ "$MOLECULE_ARG" == /* ]]; then
  MOLECULE="$MOLECULE_ARG"
else
  NAME="${MOLECULE_ARG%.xyz}"
  MOLECULE="$SCRIPT_DIR/input_molecules/${NAME}.xyz"
fi

# Resolve params: accept filename or full path
PARAMS_ARG="${2:-hydrocarbons.dat}"
if [[ "$PARAMS_ARG" == /* ]]; then
  PARAMS="$PARAMS_ARG"
else
  PARAMS="$SCRIPT_DIR/params/${PARAMS_ARG}"
fi

if [[ ! -f "$MOLECULE" ]]; then
  echo "Error: molecule file not found: $MOLECULE"
  exit 1
fi

if [[ ! -f "$PARAMS" ]]; then
  echo "Error: params file not found: $PARAMS"
  exit 1
fi

echo "==> Configuring..."
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release

echo "==> Building..."
cmake --build "$BUILD_DIR" --parallel

echo "==> Running: amber_ff $MOLECULE $PARAMS"
echo ""
"$BUILD_DIR/amber_ff" "$MOLECULE" "$PARAMS"
