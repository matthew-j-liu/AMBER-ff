#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

MOLECULE="${1:-$SCRIPT_DIR/input_molecules/ethane.xyz}"
PARAMS="${2:-$SCRIPT_DIR/params/gaff2.dat}"

echo "==> Configuring..."
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release

echo "==> Building..."
cmake --build "$BUILD_DIR" --parallel

echo "==> Running: amber_ff $MOLECULE $PARAMS"
echo ""
"$BUILD_DIR/amber_ff" "$MOLECULE" "$PARAMS"
