#!/usr/bin/env bash
# Build one or more papers to their build/ subdirectory.
# Usage: ./build.sh [paper1 paper2 ...]   (default: all papers)
# Output: papers/<name>/build/paper.pdf

set -euo pipefail

PAPERS_DIR="$(cd "$(dirname "$0")" && pwd)"
ALL_PAPERS=(core lorentz_gravity gauge_color bell matter_mass)

if [ $# -gt 0 ]; then
    PAPERS=("$@")
else
    PAPERS=("${ALL_PAPERS[@]}")
fi

for paper in "${PAPERS[@]}"; do
    PAPER_DIR="$PAPERS_DIR/$paper"
    if [ ! -f "$PAPER_DIR/paper.tex" ]; then
        echo "WARNING: $paper/paper.tex not found, skipping" >&2
        continue
    fi

    echo "=== Building $paper ==="
    mkdir -p "$PAPER_DIR/build"

    (
        cd "$PAPER_DIR"
        # Prepend parent dir so Definitions/{mdpi.cls,mdpi.bst} are found.
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
    )

    echo "==> $PAPER_DIR/build/paper.pdf"
done