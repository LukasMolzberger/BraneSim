#!/usr/bin/env bash
# Build the manuscript into a dedicated build/ directory.
# Usage:
#   ./build.sh          build paper.pdf into build/
#   ./build.sh clean    remove the build/ directory
# Output: paper/build/paper.pdf
#
# Adapted from the multi-paper build.sh of the BraneSim project for this
# single-paper layout, where Definitions/ and figures/ live inside paper/.

set -euo pipefail

PAPER_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PAPER_DIR"

if [ "${1:-}" = "clean" ]; then
    latexmk -C -output-directory=build paper.tex >/dev/null 2>&1 || true
    rm -rf build
    echo "==> cleaned build/"
    exit 0
fi

mkdir -p build

# Prepend parent dir so Definitions/{mdpi.cls,mdpi.bst,journalnames.tex} are found.
# The trailing colon preserves the default TeX/BibTeX search paths.
TEXINPUTS="..:${TEXINPUTS:-}" \
BSTINPUTS="..:${BSTINPUTS:-}" \
latexmk \
    -pdf \
    -bibtex \
    -output-directory=build \
    -interaction=nonstopmode \
    -halt-on-error \
    paper.tex

echo "==> $PAPER_DIR/build/paper.pdf"
