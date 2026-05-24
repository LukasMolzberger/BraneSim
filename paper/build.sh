#!/usr/bin/env bash
# Build the LaTeX paper into the build/ subdirectory.
#
# Usage:
#   ./build.sh         compile paper.tex -> build/paper.pdf
#   ./build.sh clean   remove build/

set -euo pipefail

cd "$(dirname "$0")"

BUILD_DIR="build"
MAIN="paper"

if [[ "${1:-}" == "clean" ]]; then
  rm -rf "$BUILD_DIR"
  echo "Removed $BUILD_DIR/"
  exit 0
fi

if ! command -v latexmk >/dev/null 2>&1; then
  echo "error: latexmk not found in PATH (install TeX Live / MacTeX)" >&2
  exit 1
fi

mkdir -p "$BUILD_DIR"

latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -output-directory="$BUILD_DIR" \
  "$MAIN.tex"

echo ""
echo "Built: $(pwd)/$BUILD_DIR/$MAIN.pdf"
