#!/usr/bin/env bash
# =============================================================================
# Sprint 4 -- Skyrme-twisted hedgehog sweep
# Phase 2 Candidate 2 existence test at alpha=0.20
#
# Iterates over all 20 configs in orchestration/configs/baryon_candidates/,
# runs the full pipeline for each, captures stdout/stderr, and aggregates
# diagnostics into summary.csv.
#
# Usage:
#   cd /path/to/BraneSim
#   bash test-runs/sprint4_skyrme_hedgehog_focused/run_sweep.sh [--dry-run]
#
# Output:
#   test-runs/sprint4_skyrme_hedgehog_focused/runs/<config_stem>/  (per-run artifacts)
#   test-runs/sprint4_skyrme_hedgehog_focused/sweep.csv            (aggregated)
# =============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SWEEP_DIR="${REPO_ROOT}/test-runs/sprint4_skyrme_hedgehog_focused"
CONFIGS_DIR="${REPO_ROOT}/orchestration/configs/baryon_candidates"
RUNS_DIR="${SWEEP_DIR}/runs"
CSV="${SWEEP_DIR}/sweep.csv"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    echo "[dry-run] Will print commands but not execute."
fi

mkdir -p "${RUNS_DIR}"

# Write CSV header
cat > "${CSV}" <<'CSVHDR'
config,u0,w,profile,curve_tag,grid_size,num_steps,run_dir,exit_code,wall_seconds,final_energy_note,status
CSVHDR

# -----------------------------------------------------------------------
# Helper: extract a JSON field from a config file (requires python3)
# -----------------------------------------------------------------------
json_field() {
    python3 -c "import json,sys; d=json.load(open('$1')); print($2)" 2>/dev/null || echo "N/A"
}

# -----------------------------------------------------------------------
# Main sweep loop
# -----------------------------------------------------------------------
CONFIGS=("${CONFIGS_DIR}"/sthh_*.json)
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

    # Parse metadata from config name: sthh_<u0>_<w>_<profile>[_<tag>]
    # Pattern: sthh_0p006_10_tanh  or  sthh_0p006_5_power2_off_low
    U0_RAW=$(echo "${STEM}" | sed -E 's/sthh_([^_]+)_.*/\1/' | sed 's/p/./')
    W_VAL=$(echo  "${STEM}" | sed -E 's/sthh_[^_]+_([0-9]+)_.*/\1/')
    PROFILE=$(echo "${STEM}" | sed -E 's/sthh_[^_]+_[0-9]+_([a-z0-9]+).*/\1/')
    CURVE_TAG="on-curve"
    if echo "${STEM}" | grep -q "off_"; then
        CURVE_TAG=$(echo "${STEM}" | grep -o 'off_[a-z0-9]*$' || echo "off_curve")
    fi
    GRID_SIZE=$(python3 -c "import json; d=json.load(open('${CFG}')); print(d['initialization'].get('grid_size', 'N/A'))" 2>/dev/null || echo "N/A")
    NUM_STEPS=$(python3 -c "import json; d=json.load(open('${CFG}')); print(d['initialization'].get('num_steps', 'N/A'))" 2>/dev/null || echo "N/A")

    echo "----------------------------------------------------------------------"
    echo "Config: ${STEM}"
    echo "  u0=${U0_RAW}  w=${W_VAL}  profile=${PROFILE}  curve=${CURVE_TAG}"
    echo "  grid=${GRID_SIZE}  steps=${NUM_STEPS}"
    echo "  run dir: ${RUN_DIR}"

    mkdir -p "${RUN_DIR}"

    CMD=(
        python3 -m orchestration.run_pipeline
        --config "${CFG}"
        --output-dir "${RUN_DIR}"
    )

    echo "  cmd: ${CMD[*]}"

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "  [dry-run] skipping execution"
        echo "${STEM},${U0_RAW},${W_VAL},${PROFILE},${CURVE_TAG},${GRID_SIZE},${NUM_STEPS},${RUN_DIR},dry-run,0,dry-run,dry-run" >> "${CSV}"
        continue
    fi

    T_START=$(date +%s)
    EXIT_CODE=0
    set +e
    python3 -m orchestration.run_pipeline \
        --config "${CFG}" \
        --output-dir "${RUN_DIR}" \
        >"${LOG_OUT}" 2>"${LOG_ERR}"
    EXIT_CODE=$?
    set -e
    T_END=$(date +%s)
    WALL_SECS=$(( T_END - T_START ))

    # Extract energy info from diagnostics if available
    DIAG_JSON="${RUN_DIR}/diagnostics/diagnostics_summary.json"
    ENERGY_NOTE="no-diag"
    if [[ -f "${DIAG_JSON}" ]]; then
        ENERGY_NOTE=$(python3 /dev/stdin "${DIAG_JSON}" 2>/dev/null <<'PYEOF'
import json, sys
path = sys.argv[1]
d = json.load(open(path))
frames = d.get("frames", [])
if len(frames) >= 2:
    e0 = frames[0].get("total_energy", None)
    ef = frames[-1].get("total_energy", None)
    if e0 is not None and ef is not None and abs(e0) > 1e-30:
        drift = abs((ef - e0) / e0) * 100
        print(f"drift={drift:.2f}pct,E0={e0:.4e},Ef={ef:.4e}")
    else:
        print("energy-unavailable")
else:
    print("too-few-frames")
PYEOF
        ) || ENERGY_NOTE="parse-error"
    fi

    if [[ ${EXIT_CODE} -eq 0 ]]; then
        STATUS="ok"
        PASSED=$(( PASSED + 1 ))
        echo "  PASSED  (wall=${WALL_SECS}s, ${ENERGY_NOTE})"
    else
        STATUS="failed"
        FAILED=$(( FAILED + 1 ))
        echo "  FAILED  (exit=${EXIT_CODE}, wall=${WALL_SECS}s)"
        echo "  --- stderr tail ---"
        tail -5 "${LOG_ERR}" 2>/dev/null || true
        echo "  ---"
    fi

    echo "${STEM},${U0_RAW},${W_VAL},${PROFILE},${CURVE_TAG},${GRID_SIZE},${NUM_STEPS},${RUN_DIR},${EXIT_CODE},${WALL_SECS},${ENERGY_NOTE},${STATUS}" >> "${CSV}"

    TOTAL=$(( TOTAL + 1 ))
done

echo ""
echo "======================================================================"
echo "Sweep complete: ${TOTAL} configs, ${PASSED} passed, ${FAILED} failed."
echo "CSV: ${CSV}"
echo "======================================================================"