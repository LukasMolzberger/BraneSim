#!/usr/bin/env bash
# =============================================================================
# Sprint 4b -- Corrected Skyrme-twisted hedgehog sweep
# VSH channel decomposition §2.5a: alpha/u0/w triples fixed by derivation.
#
# Uses branesim.run_experiment (NOT the legacy orchestration/run_pipeline.py).
# Diagnostics via branesim.diagnostics.confinement (strain weight mode).
#
# Usage:
#   cd /path/to/BraneSim
#   bash test-runs/sprint4b_skyrme_corrected/run_sweep.sh [--dry-run]
#
# Output:
#   test-runs/sprint4b_skyrme_corrected/runs/<stem>/worldvolume.zip
#   test-runs/sprint4b_skyrme_corrected/runs/<stem>/summary.json
#   test-runs/sprint4b_skyrme_corrected/runs/<stem>/stdout.log
#   test-runs/sprint4b_skyrme_corrected/runs/<stem>/stderr.log
#   test-runs/sprint4b_skyrme_corrected/sweep.csv
#
# Sizing (all 5 runs identical grid):
#   grid_shape=[64,64,64], n_slices=500, m_ambient=4
#   worldvolume per run: 4.203 GB (within 6 GB budget)
#   t_phys = 25.0, box_crossings = 0.39 (< "several")
#
# IMPORTANT: The IVP march holds the entire worldvolume (n_slices+1, n_nodes, m)
# float64 array in memory.  The current march() / WorldVolumeWriter have no
# stride / checkpoint-subsampling support.  At 64^3, the maximum n_slices within
# 6 GB is 714, giving t_phys=35.7, which is only 0.56 box-crossings at c_L=1
# (need several >= 3 crossings = t_phys >= 192 for 64^3).  The 48^3 alternative
# gives t_phys=84.7 (1.76 crossings) but box/w_max=4.8 < 6, masking dispersion
# for the w=10 run.  There is NO (grid, n_slices) within 6 GB that simultaneously
# satisfies box/w >= 6 AND t_phys >= 3*box WITHOUT checkpoint subsampling.
#
# ACTION REQUIRED before production: add stride/subsample parameter to
# branesim/solver/ivp.py march() and branesim/io/contracts.py WorldVolumeWriter
# so only every k-th slice is stored (e.g. stride=10 gives t_phys=357 at 64^3
# within the same 6 GB budget).
#
# These runs are valid for:
#   (a) Confirming the seed builds without error.
#   (b) Checking early-time confinement metrics (strain weight mode).
#   (c) Extracting short-time radius as a proxy for the R_h ratio test.
# They are NOT sufficient to confirm multi-crossing soliton stability.
# =============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SWEEP_DIR="${REPO_ROOT}/test-runs/sprint4b_skyrme_corrected"
CONFIGS_DIR="${REPO_ROOT}/orchestration/configs/baryon_skyrme_corrected"
RUNS_DIR="${SWEEP_DIR}/runs"
CSV="${SWEEP_DIR}/sweep.csv"
ANALYZE_PY="${SWEEP_DIR}/analyze_run.py"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    echo "[dry-run] Will print commands but not execute."
fi

mkdir -p "${RUNS_DIR}"

# -----------------------------------------------------------------------
# Write CSV header via the analyze helper
# -----------------------------------------------------------------------
python3 "${ANALYZE_PY}" --header > "${CSV}"

# -----------------------------------------------------------------------
# Config list (fixed order: R_h(alpha) trace, then R_h(u0) trace)
# -----------------------------------------------------------------------
CONFIGS=(
    "${CONFIGS_DIR}/a0p5_u10_w5.json"
    "${CONFIGS_DIR}/a0p7_u10_w8.json"
    "${CONFIGS_DIR}/a0p8_u10_w10.json"
    "${CONFIGS_DIR}/a0p7_u6_w5.json"
    "${CONFIGS_DIR}/a0p7_u3_w2p5.json"
)

echo "Found ${#CONFIGS[@]} config files."
echo "Output directory: ${RUNS_DIR}"
echo "CSV: ${CSV}"
echo ""

TOTAL=0
PASSED=0
FAILED=0

for CFG in "${CONFIGS[@]}"; do
    STEM=$(basename "${CFG}" .json)
    RUN_DIR="${RUNS_DIR}/${STEM}"
    LOG_OUT="${RUN_DIR}/stdout.log"
    LOG_ERR="${RUN_DIR}/stderr.log"

    echo "----------------------------------------------------------------------"
    echo "Config: ${STEM}"
    echo "  config: ${CFG}"
    echo "  run dir: ${RUN_DIR}"

    mkdir -p "${RUN_DIR}"

    CMD=(
        python3 -m branesim.run_experiment
        --config "${CFG}"
        --output-dir "${RUN_DIR}"
    )

    echo "  cmd: ${CMD[*]}"

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "  [dry-run] skipping branesim.run_experiment execution"
        echo "  [dry-run] skipping analyze_run.py execution"
        # Emit a dry-run placeholder row to CSV
        python3 - "${CFG}" "${STEM}" "${RUN_DIR}" <<'PYEOF'
import json, sys
cfg_path, stem, run_dir = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(cfg_path))
seed = d.get("seed", {})
action = d.get("action", {})
lattice = d.get("lattice", {})
alpha = action.get("alpha", "")
u0 = seed.get("u0", "")
w = seed.get("w", "")
prof = seed.get("profile_shape", "")
steep = seed.get("tanh_steepness", "")
ns = action.get("n_slices", "")
dt = action.get("dt", 0.05)
t_phys = float(ns)*float(dt) if ns != "" else ""
gs = lattice.get("grid_shape", [])
gs_str = "x".join(str(v) for v in gs)
print(
    f"{stem},{alpha},{u0},{w},{prof},{steep},"
    f"{ns},,,{gs_str},{t_phys},"
    f"dry-run,dry-run,"
    f"dry-run,dry-run,dry-run,"
    f"dry-run,dry-run,dry-run,"
    f"dry-run,dry-run,dry-run"
)
PYEOF
        TOTAL=$(( TOTAL + 1 ))
        continue
    fi

    T_START=$(date +%s)
    EXIT_CODE=0
    set +e
    python3 -m branesim.run_experiment \
        --config "${CFG}" \
        --output-dir "${RUN_DIR}" \
        >"${LOG_OUT}" 2>"${LOG_ERR}"
    EXIT_CODE=$?
    set -e
    T_END=$(date +%s)
    WALL_SECS=$(( T_END - T_START ))

    if [[ ${EXIT_CODE} -eq 0 ]]; then
        STATUS="ok"
        PASSED=$(( PASSED + 1 ))
        echo "  PASSED  (wall=${WALL_SECS}s)"
    else
        STATUS="failed"
        FAILED=$(( FAILED + 1 ))
        echo "  FAILED  (exit=${EXIT_CODE}, wall=${WALL_SECS}s)"
        echo "  --- stderr tail ---"
        tail -5 "${LOG_ERR}" 2>/dev/null || true
        echo "  ---"
    fi

    # Run confinement diagnostics (always — partial runs still emit a row).
    python3 "${ANALYZE_PY}" "${RUN_DIR}" --weight-mode strain >> "${CSV}" || \
        echo "${STEM},,,,,,,,,,,,,,,,,,,,,diag-script-error" >> "${CSV}"

    TOTAL=$(( TOTAL + 1 ))
done

echo ""
echo "======================================================================"
echo "Sweep complete: ${TOTAL} configs, ${PASSED} passed, ${FAILED} failed."
echo "CSV: ${CSV}"
echo "======================================================================"