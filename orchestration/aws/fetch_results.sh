#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 s3://bucket/prefix/results local_dir"
  exit 1
fi

S3_URI="$1"
LOCAL_DIR="$2"

mkdir -p "${LOCAL_DIR}"
aws s3 sync "${S3_URI}" "${LOCAL_DIR}"

echo "Results downloaded to ${LOCAL_DIR}"
