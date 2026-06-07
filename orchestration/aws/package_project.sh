#!/usr/bin/env bash
set -euo pipefail

OUT_PATH="${1:-/tmp/branesim-project.tar.gz}"
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"

cd "${ROOT_DIR}"

tar \
  --exclude=".git" \
  --exclude=".idea" \
  --exclude="__pycache__" \
  --exclude="*.pyc" \
  --exclude="test-runs" \
  --exclude="runs" \
  --exclude="vortex-out" \
  --exclude=".pytest_cache" \
  --exclude=".cache" \
  -czf "${OUT_PATH}" \
  .

echo "Wrote ${OUT_PATH}"
