#!/usr/bin/env bash
# AWS remote command: run the full dispersion sweep on 96^3 x 512 configs.
#
# This script is base64-encoded and passed as --remote-command to
# orchestration/aws/launch_branesim_job.py.  It runs inside the EC2 instance
# after pip install -e . completes.
#
# Memory budget:
#   One worldvolume (96^3 x 512): (512+1) * 96^3 * 4 * 8 bytes = 14.5 GB
#   bvp_chiral is a direct march — no JFNK Krylov stack.
#   r7i.4xlarge: 128 GB RAM. We run configs SEQUENTIALLY (one at a time),
#   so peak memory ~ 14.5 GB (well within budget).
#   Total AWS configs: 38 (see generate_aws_configs.py).
#   Estimated wall time per config: ~5 min (96^3 * 512 Verlet steps).
#   Total estimate: ~3 hours.

set -euo pipefail

RESULTS_DIR="${BRANESIM_RESULTS_DIR:-/opt/brane-job/results}"
CONFIG_DIR="orchestration/configs/dispersion_sweep"
OUT_DIR="${RESULTS_DIR}/dispersion_4d_bvp"
mkdir -p "${OUT_DIR}"

# Generate the AWS configs (in case they weren't in the archive)
python orchestration/configs/dispersion_sweep/generate_aws_configs.py

# Run each AWS config in sequence
CONFIGS=$(ls "${CONFIG_DIR}"/aws_*.json 2>/dev/null | sort)
TOTAL=$(echo "${CONFIGS}" | wc -l | tr -d ' ')
COUNT=0

for CFG in ${CONFIGS}; do
    COUNT=$((COUNT + 1))
    LABEL=$(basename "${CFG}" .json)
    CFG_OUT="${OUT_DIR}/${LABEL}"
    mkdir -p "${CFG_OUT}"
    echo "[${COUNT}/${TOTAL}] Running ${LABEL}..."
    python -m branesim.run_experiment \
        --config "${CFG}" \
        --output-dir "${CFG_OUT}" \
        2>&1 | tee "${CFG_OUT}/run.log"
    echo "  Done: ${LABEL}"
done

# Run the diagnostic to extract omega(k) and compute observables
echo ""
echo "Running dispersion diagnostic..."
python test-runs/dispersion_4d_bvp/run_dispersion_4d_bvp.py \
    --config-dir "${CONFIG_DIR}" \
    --pattern "aws_*" \
    --output-dir "${OUT_DIR}/diagnostics" \
    --linearity-check

echo ""
echo "All done. Results in ${OUT_DIR}"
