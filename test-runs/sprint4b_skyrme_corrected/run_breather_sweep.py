"""Baryon soliton search via time-periodic eigen-BVP (breather solver).

Physics rationale
-----------------
The forward-Verlet IVP reintroduces a causal time-direction (violating the
time-symmetry stance; OPEN_PROBLEMS A) and cannot find a bound eigenstate
without artificial damping — a non-eigenstate radiates without dissipation,
so IVP "dispersion" was partly a solver artifact.

The correct vehicle: the cyclic time-collocation JFNK root-find in
``branesim/solver/breather.py::solve_breather(mode="topological")``.
One period (P slices), R^P == R^0, saddle discipline (objective = residual_norm).

Decisive observables (per NOTES.md §breather-eigen-BVP):
  - converged AND Floquet stable (radius <= 1+tol) AND radiationless.
  - R_h(alpha)/R_h(0.5) = sqrt(alpha/(1-alpha)) ratio test.
  - R_h proportional to u0 amplitude test.

Usage
-----
Run a single smoke config (16^3, default):
    python run_breather_sweep.py --smoke

Run the full bracket (64^3, flags 64³ walltime risk):
    python run_breather_sweep.py --full

Individual config from the BRACKET list (0-indexed):
    python run_breather_sweep.py --idx 0

STEP 1 — validate the analytic band-top against the numeric at 16^3:
    python run_breather_sweep.py --validate-band-top

STEP 2 — stability-trend sweep (32^3 / 40^3, decides AWS go/no-go):
    python run_breather_sweep.py --trend

STEP 3 — convergence-robustness spot-check (alpha=0.7, u0=4, w=4, 32^3):
    python run_breather_sweep.py --spot-check
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

# Ensure the branesim package is on the path when running from this directory
_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.solver.breather import (
    BreatherOpts,
    floquet_multipliers,
    harmonic_resonance_check,
    omega_band_top_analytic,
    phonon_band_top,
    solve_breather,
)


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

OUT_DIR = _HERE / "breather_runs"
OUT_DIR.mkdir(exist_ok=True)

CSV_PATH = _HERE / "breather_sweep.csv"

# ---------------------------------------------------------------------------
# Physics bracket
# ---------------------------------------------------------------------------
# Common: k_s=rho=a=1, m_ambient=4, P=16 (even), open boundaries.
# mass = rho * a^dim = 1.0 * 1.0^3 = 1.0
#
# Full bracket grid: 64^3 (box >= 6*max(w) = 60; fits soliton).
# WARNING: 64^3 with P=16 => system size = 16*262144*4+1 = 16_777_217 unknowns.
# Each Newton-Krylov outer step calls O(inner_maxiter) matvecs, each O(n_nodes).
# Extrapolated from 16^3 smoke: ~wall_factor * (64/16)^3 = 64x per outer Newton step.
# Expect hours per bracket point on a laptop; recommend AWS/HPC for full bracket.
#
# Smoke config: 16^3, alpha=0.7, u0=3.0, skyrme_w=2.5, P=16.

FULL_BRACKET = [
    # (label, alpha, u0, w, grid_size)
    # R_h(alpha) trace at fixed u0=10
    ("a0p5_u10_w5",  0.5, 10.0,  5.0, 64),
    ("a0p7_u10_w8",  0.7, 10.0,  8.0, 64),
    ("a0p8_u10_w10", 0.8, 10.0, 10.0, 64),
    # R_h(u0) trace at fixed alpha=0.7
    ("a0p7_u6_w5",   0.7,  6.0,  5.0, 64),
    ("a0p7_u3_w2p5", 0.7,  3.0,  2.5, 64),
]

SMOKE_CONFIG = [
    # 16^3 quick smoke run: grid small enough for dense Floquet + quick JFNK
    ("smoke_a0p7_u3_w2p5", 0.7, 3.0, 2.5, 16),
]

# Fixed physics constants
K_S = 1.0
RHO = 1.0
A = 1.0        # lattice spacing
M_AMBIENT = 4
P = 16         # temporal slices per period; must be even

# ---------------------------------------------------------------------------
# BreatherOpts: smoke vs full
# ---------------------------------------------------------------------------

SMOKE_OPTS = BreatherOpts(
    tol=1e-6,
    max_iter=200,
    inner_maxiter=500,
    method="lgmres",
    verbose=True,
)

FULL_OPTS = BreatherOpts(
    tol=1e-8,
    max_iter=1000,
    inner_maxiter=2000,
    method="lgmres",
    verbose=True,
)

# ---------------------------------------------------------------------------
# STEP 1 — Analytic band-top validation grid  (16^3, quick numeric eigensystem)
# ---------------------------------------------------------------------------
# Use the SAME 16^3 config as the smoke run so the comparison is apples-to-apples.
BAND_TOP_VALIDATE_GRID = 16

# ---------------------------------------------------------------------------
# STEP 2 — Stability-trend sweep
# ---------------------------------------------------------------------------
# Common physics: k_s=rho=a=1, m_ambient=4, P=16, BreatherOpts(tol=1e-6,
# max_iter=150, inner_maxiter=1000).  Soliton-fit floor: box/w >= 6.
#
# Series A (alpha-trend, fixed moderate size): grid 32^3, u0=4.0, w=4.0
#   box=32, w=4 → ratio 8 ≥ 6  ✓
# Series B (size-trend, fixed alpha=0.7): varying (u0, w, grid)
#   (2.5,2.5,32^3) → 32/2.5=12.8 ≥ 6  ✓
#   (4.0,4.0,32^3) → shared with Series A  ✓
#   (6.0,6.0,40^3) → 40/6≈6.7 ≥ 6  ✓

TREND_OPTS = BreatherOpts(
    tol=1e-6,
    max_iter=100,
    inner_maxiter=1000,
    method="lgmres",
    verbose=True,
)

# (label, alpha, u0, w, grid_n)
TREND_SERIES_A = [
    # alpha-trend, fixed u0=4.0, w=4.0, 32^3
    ("trend_a0p5_u4_w4_32", 0.5, 4.0, 4.0, 32),
    ("trend_a0p7_u4_w4_32", 0.7, 4.0, 4.0, 32),
    ("trend_a0p8_u4_w4_32", 0.8, 4.0, 4.0, 32),
]

TREND_SERIES_B_EXTRA = [
    # size-trend, fixed alpha=0.7 — only the non-shared points
    ("trend_a0p7_u2p5_w2p5_32", 0.7, 2.5, 2.5, 32),
    ("trend_a0p7_u6_w6_40",     0.7, 6.0, 6.0, 40),
]

# Full trend bracket: 5 unique runs (Series A includes the shared alpha=0.7 point)
TREND_BRACKET = TREND_SERIES_A + TREND_SERIES_B_EXTRA

TREND_CSV_PATH = _HERE / "trend_sweep.csv"

# ---------------------------------------------------------------------------
# STEP 3 — Convergence-robustness spot-check
# ---------------------------------------------------------------------------
# Alpha=0.7, u0=4, w=4, 32^3 at inner_maxiter 1000 vs 2000.
# We re-use the "trend_a0p7_u4_w4_32" result from STEP 2 as the 1000-matvec run.
# STEP 3 re-runs at inner_maxiter=2000 and compares Floquet radii.

SPOT_CHECK_OPTS_2000 = BreatherOpts(
    tol=1e-6,
    max_iter=100,
    inner_maxiter=2000,
    method="lgmres",
    verbose=True,
)

SPOT_CHECK_CONFIG = ("spotcheck_a0p7_u4_w4_32", 0.7, 4.0, 4.0, 32)


# ---------------------------------------------------------------------------
# STEP 1: Analytic band-top validation
# ---------------------------------------------------------------------------


def validate_band_top() -> dict[str, float]:
    """Compare omega_band_top_analytic vs numeric phonon_band_top at 16^3.

    Uses: alpha=0.7, k_s=rho=a=1, m_ambient=4 (the smoke-run config).
    Reports both values and their ratio.  The analytic formula is
    sqrt(4*dim*k_s/m) which is an over-estimate (upper bound) of the true
    finite-lattice band top (see breather.py docstring for derivation).

    This is a one-time validation certificate, not new physics.
    """
    grid_n = BAND_TOP_VALIDATE_GRID
    alpha_val = 0.7
    grid_shape = (grid_n, grid_n, grid_n)
    lp = LatticeParams(grid_shape=grid_shape, spacing=A)
    lattice = SpacelikeLattice(lp)
    params = ActionParams(
        k_s=K_S,
        alpha=alpha_val,
        rho=RHO,
        dt=0.1,
        n_slices=1,
        m_ambient=M_AMBIENT,
    )
    mass = RHO * (A ** lattice.dim)

    analytic = omega_band_top_analytic(K_S, mass, lattice.dim, A)

    print(f"\n{'='*70}")
    print("  STEP 1 — Analytic vs Numeric Band-Top Validation")
    print(f"  Grid: {grid_n}^3, alpha={alpha_val}, k_s={K_S}, m={mass}, dim={lattice.dim}")
    print(f"  Analytic  omega_band_top = sqrt(4*{lattice.dim}*{K_S}/{mass}) = {analytic:.6f}")
    print(f"  Computing numeric phonon_band_top (dense eigensystem at {grid_n}^3)...")
    print(f"{'='*70}")

    t0 = time.perf_counter()
    numeric = phonon_band_top(lattice, params, mass)
    elapsed = time.perf_counter() - t0

    ratio = analytic / numeric
    pct_diff = 100.0 * (analytic - numeric) / numeric
    print(f"  Numeric   phonon_band_top = {numeric:.6f}  ({elapsed:.1f}s)")
    print(f"  Ratio analytic/numeric    = {ratio:.4f}  ({pct_diff:+.1f}%)")
    print(f"  Verdict: analytic is a {'OVER' if analytic > numeric else 'UNDER'}-estimate"
          f" ({abs(pct_diff):.1f}% off)")

    result = {
        "analytic": analytic,
        "numeric": numeric,
        "ratio": ratio,
        "pct_diff": pct_diff,
        "elapsed_s": elapsed,
    }
    json_path = OUT_DIR / "band_top_validation.json"
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {json_path}")
    return result


# ---------------------------------------------------------------------------
# R_h extraction
# ---------------------------------------------------------------------------

def fit_rh(slices: np.ndarray, lattice: SpacelikeLattice, m_ambient: int) -> float:
    """Estimate the soliton half-radius from the l=0 X⁴ profile.

    Takes the l=0 slice, computes |X⁴(r) - X⁴_far| along the l=0 Skyrme angle
    profile and finds the first radius where the angle F crosses pi/2,
    i.e. where the X⁴ component crosses u0*cos(pi/2) = 0.

    Operationally: scan the radially-binned mean of the X⁴ displacement,
    interpolate to find r where it equals 0 (the F=pi/2 crossing).
    Falls back to the RMS radius if no zero-crossing is found.

    Returns
    -------
    float  — R_h in lattice units.
    """
    dim = lattice.dim
    x4_comp = dim  # ambient component index for X⁴

    R0 = slices[0]  # (n_nodes, m_ambient)
    ref = lattice.reference_positions(m_ambient)

    # Spatial coords relative to box centre
    mi = lattice.multi_indices  # (n_nodes, dim)
    center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
    coords_c = (mi - center_mi) * lattice.params.spacing  # (n_nodes, dim)
    r = np.linalg.norm(coords_c, axis=1)  # (n_nodes,)

    # X⁴ displacement from reference (reference X⁴ = 0)
    x4_disp = R0[:, x4_comp] - ref[:, x4_comp]  # (n_nodes,)

    # Bin by radial distance: 20 bins up to max(r)
    n_bins = 20
    r_max = float(r.max())
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    r_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    x4_mean = np.full(n_bins, float("nan"))
    for b in range(n_bins):
        mask = (r >= bin_edges[b]) & (r < bin_edges[b + 1])
        if mask.sum() > 0:
            x4_mean[b] = float(x4_disp[mask].mean())

    # Find first zero-crossing of x4_mean (F = pi/2 crossing)
    valid = np.isfinite(x4_mean)
    if valid.sum() < 2:
        # Not enough data: fall back to RMS
        w2 = x4_disp ** 2
        if w2.sum() > 0:
            return float(math.sqrt((r ** 2 * w2).sum() / w2.sum()))
        return float("nan")

    r_valid = r_mid[valid]
    x4_valid = x4_mean[valid]

    # Linear interpolation for first sign change
    for i in range(len(x4_valid) - 1):
        if x4_valid[i] * x4_valid[i + 1] <= 0:
            # Linear crossing
            x0, x1 = x4_valid[i], x4_valid[i + 1]
            r0, r1 = r_valid[i], r_valid[i + 1]
            if abs(x1 - x0) > 1e-30:
                return float(r0 + (r1 - r0) * (-x0) / (x1 - x0))

    # No zero-crossing found: return RMS radius
    w2 = x4_disp ** 2
    if w2.sum() > 0:
        return float(math.sqrt((r ** 2 * w2).sum() / w2.sum()))
    return float("nan")


# ---------------------------------------------------------------------------
# Single-run executor
# ---------------------------------------------------------------------------

def run_one(
    label: str,
    alpha: float,
    u0: float,
    w: float,
    grid_n: int,
    opts: BreatherOpts,
    band_top_override: float | None = None,
) -> dict[str, Any]:
    """Execute one breather eigen-BVP solve and return a flat result dict.

    Parameters
    ----------
    band_top_override : float, optional
        If given, pass this pre-computed band_top to harmonic_resonance_check
        instead of computing the dense phonon eigensystem.  Use
        ``omega_band_top_analytic(k_s, mass, dim)`` for large grids where the
        full eigensystem is infeasible (O(N³) cost).  The analytic value is an
        over-estimate (conservative); see NOTES.md §analytic-band-top-validation.
    """
    print(f"\n{'='*70}")
    print(f"  Run: {label}")
    print(f"  alpha={alpha}, u0={u0}, w={w}, grid={grid_n}^3, P={P}")
    print(f"  tol={opts.tol}, max_iter={opts.max_iter}, inner_maxiter={opts.inner_maxiter}")
    if band_top_override is not None:
        print(f"  band_top (analytic override): {band_top_override:.6f}")
    print(f"{'='*70}")

    grid_shape = (grid_n, grid_n, grid_n)
    lp = LatticeParams(grid_shape=grid_shape, spacing=A)
    lattice = SpacelikeLattice(lp)

    params = ActionParams(
        k_s=K_S,
        alpha=alpha,
        rho=RHO,
        dt=0.1,       # ignored by breather solver; kept for ActionParams validity
        n_slices=1,   # ignored by breather solver
        m_ambient=M_AMBIENT,
    )

    mass = RHO * (A ** lattice.dim)  # = 1.0 for k_s=rho=a=1

    t_wall_start = time.perf_counter()

    result = solve_breather(
        lattice,
        params,
        mass,
        P=P,
        amplitude=u0,
        mode="topological",
        skyrme_w=w,
        skyrme_profile="power2",
        opts=opts,
    )

    t_wall_solve = time.perf_counter() - t_wall_start

    converged = result["converged"]
    res_init = result["residual_initial"]
    res_final = result["residual_norm"]
    T_sol = result["T"]
    omega_sol = result["omega"]

    print(f"\n  Solve complete: converged={converged}")
    print(f"  residual: {res_init:.3e} -> {res_final:.3e}")
    print(f"  T={T_sol:.6f}, omega={omega_sol:.6f}")
    print(f"  Solve walltime: {t_wall_solve:.1f}s")

    # ---- Floquet stability (on converged or best-available orbit) ----
    t_floquet_start = time.perf_counter()
    floquet = floquet_multipliers(
        result["slices"],
        result["T"],
        lattice,
        params,
        mass,
    )
    t_floquet = time.perf_counter() - t_floquet_start

    spectral_radius = floquet["spectral_radius"]
    floquet_stable = floquet["stable"]
    floquet_method = floquet["method"]
    n_state = floquet["n_state"]

    print(f"  Floquet: spectral_radius={spectral_radius:.4f}, stable={floquet_stable}"
          f"  (method={floquet_method}, n_state={n_state}, {t_floquet:.1f}s)")

    # ---- Harmonic resonance check ----
    # If band_top_override is provided (analytic), use it directly to avoid the
    # O(N³) dense eigensystem cost, which is infeasible at 32^3+ grids.
    t_res_start = time.perf_counter()
    resonance = harmonic_resonance_check(
        omega_sol,
        lattice,
        params,
        mass,
        band_top=band_top_override,  # None → compute dense; float → use analytic
    )
    t_res = time.perf_counter() - t_res_start

    radiationless = resonance["radiationless"]
    lowest_n_ge_2 = resonance["lowest_in_band_n_ge_2"]
    band_top = resonance["band_top"]
    transverse_top = resonance["transverse_top"]

    print(f"  Resonance: radiationless={radiationless}, "
          f"lowest_in_band_n_ge_2={lowest_n_ge_2}")
    print(f"  band_top={band_top:.4f}, transverse_top={transverse_top:.4f} ({t_res:.1f}s)")

    # ---- R_h from l=0 X⁴ profile ----
    R_h = fit_rh(result["slices"], lattice, M_AMBIENT)
    print(f"  R_h (fitted)={R_h:.4f} lattice units")

    total_walltime = time.perf_counter() - t_wall_start

    # ---- Save converged worldtube ----
    npz_path = OUT_DIR / f"{label}_worldtube.npz"
    np.savez_compressed(
        npz_path,
        slices=result["slices"].astype(np.float32),  # float32 to save space
        T=np.float64(T_sol),
        omega=np.float64(omega_sol),
        alpha=np.float64(alpha),
        u0=np.float64(u0),
        w=np.float64(w),
        grid_n=np.int32(grid_n),
        P=np.int32(P),
        residual_norm=np.float64(res_final),
        converged=np.bool_(converged),
    )
    print(f"  Worldtube saved: {npz_path}")

    # ---- Assemble flat result dict ----
    row = {
        "label": label,
        "alpha": alpha,
        "u0": u0,
        "w": w,
        "grid_n": grid_n,
        "P": P,
        "m_ambient": M_AMBIENT,
        "converged": converged,
        "residual_initial": res_init,
        "residual_final": res_final,
        "T": T_sol,
        "omega": omega_sol,
        "floquet_spectral_radius": spectral_radius,
        "floquet_stable": floquet_stable,
        "floquet_method": floquet_method,
        "radiationless": radiationless,
        "lowest_in_band_n_ge_2": lowest_n_ge_2,
        "band_top": band_top,
        "transverse_top": transverse_top,
        "R_h": R_h,
        "walltime_s": total_walltime,
    }

    # ---- Save per-run JSON ----
    json_path = OUT_DIR / f"{label}_result.json"
    # JSON-serialise: convert numpy scalars to Python types
    json_row = {
        k: (v.item() if hasattr(v, "item") else v)
        for k, v in row.items()
    }
    # Add full harmonic detail
    json_row["harmonics"] = resonance["harmonics"]
    json_row["floquet_dense"] = floquet["dense"]
    json_row["floquet_n_state"] = n_state
    with open(json_path, "w") as f:
        json.dump(json_row, f, indent=2)
    print(f"  Per-run JSON saved: {json_path}")

    print(f"\n  Total walltime: {total_walltime:.1f}s")
    return row


# ---------------------------------------------------------------------------
# CSV writer helpers
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "label", "alpha", "u0", "w", "grid_n", "P", "m_ambient",
    "converged", "residual_initial", "residual_final", "T", "omega",
    "floquet_spectral_radius", "floquet_stable", "floquet_method",
    "radiationless", "lowest_in_band_n_ge_2",
    "band_top", "transverse_top", "R_h", "walltime_s",
]


def append_csv(row: dict[str, Any], csv_path: Path = CSV_PATH) -> None:
    """Append one result row to the CSV, creating it with a header if needed."""
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _error_row(label: str, alpha: float, u0: float, w: float, grid_n: int) -> dict[str, Any]:
    """Return a sentinel error row for the CSV when a run throws an exception."""
    return {
        "label": label, "alpha": alpha, "u0": u0, "w": w,
        "grid_n": grid_n, "P": P, "m_ambient": M_AMBIENT,
        "converged": False, "residual_initial": float("nan"),
        "residual_final": float("nan"), "T": float("nan"),
        "omega": float("nan"),
        "floquet_spectral_radius": float("nan"),
        "floquet_stable": False, "floquet_method": "ERROR",
        "radiationless": False, "lowest_in_band_n_ge_2": None,
        "band_top": float("nan"), "transverse_top": float("nan"),
        "R_h": float("nan"), "walltime_s": float("nan"),
    }


def _print_summary(rows: list[dict[str, Any]], csv_path: Path) -> None:
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'label':<30} {'conv':>5} {'res_final':>12} {'Floq_rad':>9} {'radless':>7} {'R_h':>7}")
    for r in rows:
        print(
            f"{r['label']:<30} {str(r['converged']):>5} "
            f"{r['residual_final']:>12.3e} "
            f"{r['floquet_spectral_radius']:>9.4f} "
            f"{str(r['radiationless']):>7} "
            f"{r['R_h']:>7.3f}"
        )
    print(f"\nCSV: {csv_path}")
    print(f"Run JSON + worldtube NPZ: {OUT_DIR}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Baryon soliton eigen-BVP breather sweep (time-periodic JFNK)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--smoke",
        action="store_true",
        help="Run the 16^3 smoke config (alpha=0.7, u0=3, w=2.5) to confirm convergence",
    )
    group.add_argument(
        "--full",
        action="store_true",
        help="Run the full 64^3 bracket (WARNING: hours per point on a laptop)",
    )
    group.add_argument(
        "--idx",
        type=int,
        metavar="I",
        help="Run single FULL_BRACKET entry by 0-based index",
    )
    group.add_argument(
        "--validate-band-top",
        action="store_true",
        help=(
            "STEP 1: Compute numeric phonon_band_top at 16^3 and compare to the "
            "analytic omega_band_top_analytic.  Validates the closed-form formula "
            "before it is used on large grids where the dense eigensystem is infeasible."
        ),
    )
    group.add_argument(
        "--trend",
        action="store_true",
        help=(
            "STEP 2: Run 5-point stability-trend sweep at 32^3 / 40^3 "
            "(Series A: alpha-trend; Series B: size-trend).  Uses analytic band_top "
            "to skip the dense phonon eigensystem.  Writes trend_sweep.csv."
        ),
    )
    group.add_argument(
        "--spot-check",
        action="store_true",
        help=(
            "STEP 3: Convergence-robustness check for alpha=0.7, u0=4, w=4, 32^3 "
            "at inner_maxiter=2000 (compare Floquet radius vs the 1000-matvec run "
            "from --trend)."
        ),
    )
    group.add_argument(
        "--trend-idx",
        type=int,
        metavar="I",
        help="Run single TREND_BRACKET entry by 0-based index (0=alpha=0.5, 1=alpha=0.7, ...).",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # STEP 1: Analytic band-top validation
    # ------------------------------------------------------------------
    if args.validate_band_top:
        result = validate_band_top()
        print(f"\n[BAND-TOP VALIDATION COMPLETE]")
        print(f"  analytic = {result['analytic']:.6f}")
        print(f"  numeric  = {result['numeric']:.6f}")
        print(f"  ratio    = {result['ratio']:.4f}  ({result['pct_diff']:+.1f}%)")
        return

    # ------------------------------------------------------------------
    # STEP 2: Stability-trend sweep
    # ------------------------------------------------------------------
    if args.trend:
        print("\n[TREND SWEEP] 5-point stability-trend sweep")
        print("[TREND] Series A: alpha-trend (32^3, u0=4, w=4, alpha in {0.5,0.7,0.8})")
        print("[TREND] Series B: size-trend  (alpha=0.7, sizes: 32^3 x2 + 40^3)")
        print("[TREND] Using analytic band_top (over-estimate) — no dense eigensystem.")
        print("[TREND] BreatherOpts: tol=1e-6, max_iter=100 (32^3) / max_iter=50 (40^3), inner_maxiter=1000")

        # Walltime cap: 40^3 is ~(40/32)^3 ≈ 2× heavier per matvec than 32^3.
        # At TREND_OPTS (max_iter=150) the 40^3 run is estimated at ~76 min — above
        # the ~20 min local-run budget.  Cap max_iter=50 for 40^3 to stay bounded.
        # The Floquet radius after 50 outer Newton steps is still informative for
        # the trend verdict; the note below records the capped budget explicitly.
        _GRID40_OPTS = BreatherOpts(
            tol=1e-6,
            max_iter=50,
            inner_maxiter=1000,
            method="lgmres",
            verbose=True,
        )
        print("[TREND] NOTE: 40^3 run capped at max_iter=50 to stay within ~20 min local budget.")

        rows = []
        for entry in TREND_BRACKET:
            label, alpha, u0, w, grid_n = entry
            # Compute analytic band_top for THIS grid/alpha (k_s=1, mass=1, dim=3)
            mass = RHO * (A ** 3)  # = 1.0 for rho=a=1
            band_top_an = omega_band_top_analytic(K_S, mass, 3, A)
            run_opts = _GRID40_OPTS if grid_n >= 40 else TREND_OPTS
            if grid_n >= 40:
                print(f"[TREND] {label}: using max_iter=50 (40^3 walltime cap)")
            try:
                row = run_one(
                    label, alpha, u0, w, grid_n, run_opts,
                    band_top_override=band_top_an,
                )
            except Exception as exc:
                import traceback
                print(f"\n  ERROR in run {label}: {exc}")
                traceback.print_exc()
                row = _error_row(label, alpha, u0, w, grid_n)
            rows.append(row)
            append_csv(row, csv_path=TREND_CSV_PATH)
        _print_summary(rows, TREND_CSV_PATH)
        return

    # ------------------------------------------------------------------
    # STEP 3: Convergence-robustness spot-check
    # ------------------------------------------------------------------
    if args.spot_check:
        label, alpha, u0, w, grid_n = SPOT_CHECK_CONFIG
        mass = RHO * (A ** 3)
        band_top_an = omega_band_top_analytic(K_S, mass, 3, A)

        print("\n[SPOT-CHECK] Convergence-robustness: alpha=0.7, u0=4, w=4, 32^3")
        print("[SPOT-CHECK] inner_maxiter=2000 (vs 1000 in --trend)")
        print("[SPOT-CHECK] Comparing Floquet radii at two convergence levels.")

        # Check if 1000-matvec result already exists from --trend
        trend_json = OUT_DIR / "trend_a0p7_u4_w4_32_result.json"
        radius_1000: float | None = None
        if trend_json.exists():
            with open(trend_json) as f:
                trend_data = json.load(f)
            radius_1000 = trend_data.get("floquet_spectral_radius")
            res_1000 = trend_data.get("residual_final")
            print(f"[SPOT-CHECK] Loaded 1000-matvec result: Floquet radius={radius_1000:.4f}, "
                  f"residual={res_1000:.3e}")
        else:
            print("[SPOT-CHECK] No 1000-matvec result found; run --trend first for comparison.")

        try:
            row = run_one(
                label, alpha, u0, w, grid_n, SPOT_CHECK_OPTS_2000,
                band_top_override=band_top_an,
            )
        except Exception as exc:
            import traceback
            print(f"\n  ERROR in spot-check: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)

        radius_2000 = row["floquet_spectral_radius"]
        res_2000 = row["residual_final"]
        print(f"\n[SPOT-CHECK RESULT]")
        if radius_1000 is not None:
            delta = abs(radius_2000 - radius_1000)
            print(f"  Floquet radius @ inner_maxiter=1000: {radius_1000:.4f}  (res={res_1000:.3e})")
            print(f"  Floquet radius @ inner_maxiter=2000: {radius_2000:.4f}  (res={res_2000:.3e})")
            print(f"  Delta: {delta:.4f}  ({'stable' if delta < 0.1 else 'SIGNIFICANT shift'})")
        else:
            print(f"  Floquet radius @ inner_maxiter=2000: {radius_2000:.4f}  (res={res_2000:.3e})")
            print(f"  (No 1000-matvec baseline; run --trend first)")

        append_csv(row, csv_path=TREND_CSV_PATH)
        return

    # ------------------------------------------------------------------
    # Single TREND_BRACKET entry by index (resume after partial run)
    # ------------------------------------------------------------------
    if args.trend_idx is not None:
        idx = args.trend_idx
        if idx < 0 or idx >= len(TREND_BRACKET):
            print(f"ERROR: --trend-idx must be 0..{len(TREND_BRACKET)-1}; got {idx}")
            sys.exit(1)
        entry = TREND_BRACKET[idx]
        label, alpha, u0, w, grid_n = entry
        mass = RHO * (A ** 3)
        band_top_an = omega_band_top_analytic(K_S, mass, 3, A)
        run_opts = BreatherOpts(tol=1e-6, max_iter=50, inner_maxiter=1000,
                                method="lgmres", verbose=True) if grid_n >= 40 else TREND_OPTS
        print(f"\n[TREND-IDX={idx}] {label}  (grid={grid_n}^3, alpha={alpha}, u0={u0}, w={w})")
        try:
            row = run_one(label, alpha, u0, w, grid_n, run_opts,
                          band_top_override=band_top_an)
        except Exception as exc:
            import traceback
            print(f"\n  ERROR: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)
        append_csv(row, csv_path=TREND_CSV_PATH)
        _print_summary([row], TREND_CSV_PATH)
        return

    # ------------------------------------------------------------------
    # Legacy modes: --smoke, --full, --idx
    # ------------------------------------------------------------------
    if args.smoke:
        bracket = SMOKE_CONFIG
        opts = SMOKE_OPTS
        print("\n[SMOKE] 16^3 grid, alpha=0.7, u0=3.0, w=2.5, P=16")
        print("[SMOKE] tol=1e-6, max_iter=200, inner_maxiter=500")
    elif args.full:
        bracket = FULL_BRACKET
        opts = FULL_OPTS
        print(f"\n[FULL] {len(bracket)}-point bracket at 64^3.")
        print("[FULL] WARNING: JFNK at 64^3 (system size ~16.8M unknowns) is ~64x")
        print("[FULL] heavier than 16^3 per outer Newton step.  Expect hours per")
        print("[FULL] bracket point on a laptop.  Recommend AWS/HPC for this sweep.")
    else:
        idx = args.idx
        if idx < 0 or idx >= len(FULL_BRACKET):
            print(f"ERROR: --idx must be 0..{len(FULL_BRACKET)-1}; got {idx}")
            sys.exit(1)
        bracket = [FULL_BRACKET[idx]]
        opts = FULL_OPTS

    rows = []
    for entry in bracket:
        label, alpha, u0, w, grid_n = entry
        try:
            row = run_one(label, alpha, u0, w, grid_n, opts)
        except Exception as exc:
            print(f"\n  ERROR in run {label}: {exc}")
            row = _error_row(label, alpha, u0, w, grid_n)
        rows.append(row)
        append_csv(row)

    _print_summary(rows, CSV_PATH)


if __name__ == "__main__":
    main()
