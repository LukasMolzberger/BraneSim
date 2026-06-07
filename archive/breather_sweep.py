"""Baryon soliton search via time-periodic eigen-BVP (breather solver).

Packaged entrypoint for local and AWS remote runs.  Imports only branesim.*
APIs; no test-runs/ dependencies.

Physics rationale
-----------------
The cyclic time-collocation JFNK root-find in
``branesim.solver.breather.solve_breather(mode="topological")`` finds one
period (P slices), R^P == R^0.  Decisive observables:
  - converged AND Floquet stable (radius <= 1 + tol) AND radiationless.
  - R_h(alpha)/R_h(0.5) = sqrt(alpha/(1-alpha)) ratio test.
  - R_h proportional to u0 amplitude test.

CLI usage
---------
Validate analytic band-top (16^3, quick):
    branesim-breather --validate-band-top

Smoke run (16^3):
    branesim-breather --smoke

Stability-trend sweep (32^3 / 40^3):
    branesim-breather --trend

Single trend-bracket point by index:
    branesim-breather --trend-idx 0

AWS fat-soliton bracket (64^3 + 48^3, alpha in {0.8, 0.85, 0.9}):
    branesim-breather --aws-bracket --output-dir "$BRANESIM_RESULTS_DIR"

Single AWS bracket point by index:
    branesim-breather --aws-idx 0 --output-dir "$BRANESIM_RESULTS_DIR"

Full 64^3 bracket (warning: hours per point on a laptop):
    branesim-breather --full --output-dir /tmp/breather_out

Single full-bracket point:
    branesim-breather --idx 0 --output-dir /tmp/breather_out
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
# Fixed physics constants (k_s=rho=a=1, dimensionless; mass = rho*a^dim = 1)
# ---------------------------------------------------------------------------

K_S: float = 1.0
RHO: float = 1.0
A: float = 1.0          # lattice spacing
M_AMBIENT: int = 4
P: int = 16             # temporal slices per period; must be even

# ---------------------------------------------------------------------------
# Bracket definitions
# ---------------------------------------------------------------------------

# Full 64^3 bracket (original sprint4b sweep)
FULL_BRACKET: list[tuple[str, float, float, float, int]] = [
    # (label, alpha, u0, w, grid_n)
    ("a0p5_u10_w5",  0.5, 10.0,  5.0, 64),
    ("a0p7_u10_w8",  0.7, 10.0,  8.0, 64),
    ("a0p8_u10_w10", 0.8, 10.0, 10.0, 64),
    ("a0p7_u6_w5",   0.7,  6.0,  5.0, 64),
    ("a0p7_u3_w2p5", 0.7,  3.0,  2.5, 64),
]

SMOKE_CONFIG: list[tuple[str, float, float, float, int]] = [
    ("smoke_a0p7_u3_w2p5", 0.7, 3.0, 2.5, 16),
]

# ---------------------------------------------------------------------------
# AWS fat-soliton bracket
# ---------------------------------------------------------------------------
# Target: alpha in {0.8, 0.85, 0.9}, u0=10, fat soliton w~8.
# Constraint: box/w >= 6  =>  box >= 6*w.
# 64^3: box=64, w=8 -> ratio=8 >=6 (fits).  w=10 -> ratio=6.4 >=6 (fits).
# 48^3: box=48, w=8 -> ratio=6.0 (marginal; use w<=7 for clean fit, or accept w=8).
# AWS opts: tol=1e-6, max_iter=200, inner_maxiter=2000  (from problem statement).
#
# Points are INDEPENDENT — run in parallel on a multi-core instance via
# subprocess or just sequentially (6 points, each hours on 64^3).

# RE-SCOPED 2026-06-04 (alpha-scan -> SIZE-scan). The local 32^3 study found
# converged breather orbits at alpha=0.80/0.85/0.90 with Floquet rho = 3.59/3.47/3.70
# -- robustly unstable, FLAT in alpha (alpha ruled out as stabilizer). The one
# remaining lever is soliton SIZE: does fattening the soliton drive rho below 1?
# This bracket fixes alpha=0.85 and scans the physical width w (u0=w, near-stationary
# seed), all box/w >= 6. idx2 vs idx3 (w=7 at 48^3 vs 64^3) is the grid-control that
# separates a physical-size trend from a grid-resolution/containment artifact.
# Decisive read: does rho fall monotonically with w toward 1?
AWS_BRACKET: list[tuple[str, float, float, float, int]] = [
    # (label, alpha, u0, w, grid_n)   -- fixed alpha=0.85, ascending width
    ("aws_a85_w4_32",  0.85,  4.0,  4.0, 32),  # local reference (rho=3.47)
    ("aws_a85_w6_48",  0.85,  6.0,  6.0, 48),  # box/w = 8.0
    ("aws_a85_w7_48",  0.85,  7.0,  7.0, 48),  # box/w = 6.9
    ("aws_a85_w7_64",  0.85,  7.0,  7.0, 64),  # grid-control vs w7_48 (box/w = 9.1)
    ("aws_a85_w9_64",  0.85,  9.0,  9.0, 64),  # box/w = 7.1
    ("aws_a85_w10_64", 0.85, 10.0, 10.0, 64),  # box/w = 6.4 (fattest)
]

# BreatherOpts for the AWS bracket (per problem statement)
AWS_OPTS = BreatherOpts(
    tol=1e-6,
    max_iter=200,
    inner_maxiter=2000,
    method="lgmres",
    verbose=True,
)

# ---------------------------------------------------------------------------
# Stability-trend sweep (STEP 2: 32^3 / 40^3, pre-AWS go/no-go)
# ---------------------------------------------------------------------------

TREND_OPTS = BreatherOpts(
    tol=1e-6,
    max_iter=100,
    inner_maxiter=1000,
    method="lgmres",
    verbose=True,
)

TREND_BRACKET: list[tuple[str, float, float, float, int]] = [
    # Series A: alpha-trend, fixed u0=4.0, w=4.0, 32^3
    ("trend_a0p5_u4_w4_32",   0.5, 4.0, 4.0, 32),
    ("trend_a0p7_u4_w4_32",   0.7, 4.0, 4.0, 32),
    ("trend_a0p8_u4_w4_32",   0.8, 4.0, 4.0, 32),
    # Series B: size-trend, alpha=0.7
    ("trend_a0p7_u2p5_w2p5_32", 0.7, 2.5, 2.5, 32),
    ("trend_a0p7_u6_w6_40",     0.7, 6.0, 6.0, 40),
]

# Spot-check config (STEP 3: convergence robustness)
SPOT_CHECK_CONFIG: tuple[str, float, float, float, int] = (
    "spotcheck_a0p7_u4_w4_32", 0.7, 4.0, 4.0, 32
)

SPOT_CHECK_OPTS_2000 = BreatherOpts(
    tol=1e-6,
    max_iter=100,
    inner_maxiter=2000,
    method="lgmres",
    verbose=True,
)

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
# CSV fields
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "label", "alpha", "u0", "w", "grid_n", "P", "m_ambient",
    "converged", "residual_initial", "residual_final", "T", "omega",
    "floquet_spectral_radius", "floquet_stable", "floquet_method",
    "radiationless", "lowest_in_band_n_ge_2",
    "band_top", "transverse_top", "R_h", "walltime_s",
]

# ---------------------------------------------------------------------------
# R_h extraction
# ---------------------------------------------------------------------------


def fit_rh(slices: np.ndarray, lattice: SpacelikeLattice, m_ambient: int) -> float:
    """Estimate the soliton half-radius from the l=0 X^4 profile.

    Scans the radially-binned mean of the X^4 displacement and interpolates
    to find r where it crosses zero (the F=pi/2 crossing).  Falls back to
    the RMS radius if no zero-crossing is found.

    Returns
    -------
    float  — R_h in lattice units.
    """
    dim = lattice.dim
    x4_comp = dim  # ambient component index for X^4

    R0 = slices[0]  # (n_nodes, m_ambient)
    ref = lattice.reference_positions(m_ambient)

    mi = lattice.multi_indices  # (n_nodes, dim)
    center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
    coords_c = (mi - center_mi) * lattice.params.spacing  # (n_nodes, dim)
    r = np.linalg.norm(coords_c, axis=1)  # (n_nodes,)

    x4_disp = R0[:, x4_comp] - ref[:, x4_comp]  # (n_nodes,)

    n_bins = 20
    r_max = float(r.max())
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    r_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    x4_mean = np.full(n_bins, float("nan"))
    for b in range(n_bins):
        mask = (r >= bin_edges[b]) & (r < bin_edges[b + 1])
        if mask.sum() > 0:
            x4_mean[b] = float(x4_disp[mask].mean())

    valid = np.isfinite(x4_mean)
    if valid.sum() < 2:
        w2 = x4_disp ** 2
        if w2.sum() > 0:
            return float(math.sqrt((r ** 2 * w2).sum() / w2.sum()))
        return float("nan")

    r_valid = r_mid[valid]
    x4_valid = x4_mean[valid]

    for i in range(len(x4_valid) - 1):
        if x4_valid[i] * x4_valid[i + 1] <= 0:
            x0, x1 = x4_valid[i], x4_valid[i + 1]
            r0, r1 = r_valid[i], r_valid[i + 1]
            if abs(x1 - x0) > 1e-30:
                return float(r0 + (r1 - r0) * (-x0) / (x1 - x0))

    w2 = x4_disp ** 2
    if w2.sum() > 0:
        return float(math.sqrt((r ** 2 * w2).sum() / w2.sum()))
    return float("nan")


# ---------------------------------------------------------------------------
# STEP 1: Analytic band-top validation
# ---------------------------------------------------------------------------


def validate_band_top(out_dir: Path) -> dict[str, float]:
    """Compare omega_band_top_analytic vs numeric phonon_band_top at 16^3.

    The analytic formula sqrt(4*dim*k_s/m) is an over-estimate (upper bound)
    of the true finite-lattice band top.  This is a one-time validation
    certificate, not new physics.
    """
    grid_n = 16
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
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "band_top_validation.json"
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {json_path}")
    return result


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
    out_dir: Path,
    csv_path: Path,
    band_top_override: float | None = None,
) -> dict[str, Any]:
    """Execute one breather eigen-BVP solve and return a flat result dict.

    Parameters
    ----------
    band_top_override
        If given, pass this pre-computed band_top to harmonic_resonance_check
        instead of computing the dense phonon eigensystem.  Use
        omega_band_top_analytic(k_s, mass, dim) for large grids where the
        full eigensystem is infeasible.  The analytic value is an over-estimate
        (conservative).
    """
    print(f"\n{'='*70}")
    print(f"  Run: {label}")
    print(f"  alpha={alpha}, u0={u0}, w={w}, grid={grid_n}^3, P={P}")
    print(f"  tol={opts.tol}, max_iter={opts.max_iter}, inner_maxiter={opts.inner_maxiter}")
    if band_top_override is not None:
        print(f"  band_top (analytic override): {band_top_override:.6f}")
    print(f"{'='*70}")

    grid_shape = tuple(grid_n for _ in range(3))
    lp = LatticeParams(grid_shape=grid_shape, spacing=A)
    lattice = SpacelikeLattice(lp)

    params = ActionParams(
        k_s=K_S,
        alpha=alpha,
        rho=RHO,
        dt=0.1,      # unused by breather solver; required for ActionParams validity
        n_slices=1,  # unused by breather solver
        m_ambient=M_AMBIENT,
    )

    mass = RHO * (A ** lattice.dim)

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

    # Floquet stability
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

    # Harmonic resonance check
    t_res_start = time.perf_counter()
    resonance = harmonic_resonance_check(
        omega_sol,
        lattice,
        params,
        mass,
        band_top=band_top_override,
    )
    t_res = time.perf_counter() - t_res_start

    radiationless = resonance["radiationless"]
    lowest_n_ge_2 = resonance["lowest_in_band_n_ge_2"]
    band_top = resonance["band_top"]
    transverse_top = resonance["transverse_top"]

    print(f"  Resonance: radiationless={radiationless}, "
          f"lowest_in_band_n_ge_2={lowest_n_ge_2}")
    print(f"  band_top={band_top:.4f}, transverse_top={transverse_top:.4f} ({t_res:.1f}s)")

    # R_h from l=0 X^4 profile
    R_h = fit_rh(result["slices"], lattice, M_AMBIENT)
    print(f"  R_h (fitted)={R_h:.4f} lattice units")

    total_walltime = time.perf_counter() - t_wall_start

    # Save converged worldtube
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir / f"{label}_worldtube.npz"
    np.savez_compressed(
        npz_path,
        slices=result["slices"].astype(np.float32),
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

    row: dict[str, Any] = {
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

    # Per-run JSON
    json_path = out_dir / f"{label}_result.json"
    json_row = {
        k: (v.item() if hasattr(v, "item") else v)
        for k, v in row.items()
    }
    json_row["harmonics"] = resonance["harmonics"]
    json_row["floquet_dense"] = floquet["dense"]
    json_row["floquet_n_state"] = n_state
    with open(json_path, "w") as f:
        json.dump(json_row, f, indent=2)
    print(f"  Per-run JSON saved: {json_path}")

    # Append CSV row
    _append_csv(row, csv_path)

    print(f"\n  Total walltime: {total_walltime:.1f}s")
    return row


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------


def _append_csv(row: dict[str, Any], csv_path: Path) -> None:
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _error_row(label: str, alpha: float, u0: float, w: float, grid_n: int) -> dict[str, Any]:
    return {
        "label": label, "alpha": alpha, "u0": u0, "w": w,
        "grid_n": grid_n, "P": P, "m_ambient": M_AMBIENT,
        "converged": False,
        "residual_initial": float("nan"),
        "residual_final": float("nan"),
        "T": float("nan"),
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
    print(f"{'label':<35} {'conv':>5} {'res_final':>12} {'Floq_rad':>9} {'radless':>7} {'R_h':>7}")
    for r in rows:
        print(
            f"{r['label']:<35} {str(r['converged']):>5} "
            f"{r['residual_final']:>12.3e} "
            f"{r['floquet_spectral_radius']:>9.4f} "
            f"{str(r['radiationless']):>7} "
            f"{r['R_h']:>7.3f}"
        )
    print(f"\nCSV: {csv_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Baryon soliton eigen-BVP breather sweep (time-periodic JFNK). "
            "Entrypoint: branesim-breather (or python -m branesim.experiments.breather_sweep)."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Directory for CSV, JSON, and NPZ outputs.  Defaults to "
            "$BRANESIM_RESULTS_DIR if set, else ./breather_results."
        ),
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--validate-band-top",
        action="store_true",
        help=(
            "STEP 1: Compare numeric phonon_band_top (16^3) to analytic "
            "omega_band_top_analytic.  Validates the closed-form formula."
        ),
    )
    group.add_argument(
        "--smoke",
        action="store_true",
        help="Quick smoke run: 16^3, alpha=0.7, u0=3, w=2.5, P=16.",
    )
    group.add_argument(
        "--trend",
        action="store_true",
        help=(
            "STEP 2: 5-point stability-trend sweep at 32^3/40^3 "
            "(alpha-trend + size-trend).  Uses analytic band_top."
        ),
    )
    group.add_argument(
        "--spot-check",
        action="store_true",
        help=(
            "STEP 3: Convergence-robustness check for alpha=0.7, u0=4, w=4, "
            "32^3 at inner_maxiter=2000."
        ),
    )
    group.add_argument(
        "--trend-idx",
        type=int,
        metavar="I",
        help="Run single TREND_BRACKET entry by 0-based index.",
    )
    group.add_argument(
        "--aws-bracket",
        action="store_true",
        help=(
            "AWS fat-soliton bracket: alpha in {0.8, 0.85, 0.9}, u0=10, "
            "64^3 + 48^3 variants.  Uses analytic band_top.  "
            "BreatherOpts: tol=1e-6, max_iter=200, inner_maxiter=2000."
        ),
    )
    group.add_argument(
        "--aws-idx",
        type=int,
        metavar="I",
        help=(
            "Run single AWS_BRACKET entry by 0-based index "
            "(0=64^3/a=0.80, 1=64^3/a=0.85, 2=64^3/a=0.90, "
            "3=48^3/a=0.80, 4=48^3/a=0.85, 5=48^3/a=0.90)."
        ),
    )
    group.add_argument(
        "--full",
        action="store_true",
        help="Full 64^3 bracket (warning: hours per point on a laptop).",
    )
    group.add_argument(
        "--idx",
        type=int,
        metavar="I",
        help="Run single FULL_BRACKET entry by 0-based index.",
    )

    args = parser.parse_args()

    # Resolve output directory
    if args.output_dir is not None:
        base_out = Path(args.output_dir)
    elif os.environ.get("BRANESIM_RESULTS_DIR"):
        base_out = Path(os.environ["BRANESIM_RESULTS_DIR"])
    else:
        base_out = Path("./breather_results")

    base_out.mkdir(parents=True, exist_ok=True)
    csv_path = base_out / "breather_sweep.csv"

    # ------------------------------------------------------------------
    # STEP 1: Analytic band-top validation
    # ------------------------------------------------------------------
    if args.validate_band_top:
        result = validate_band_top(base_out)
        print(f"\n[BAND-TOP VALIDATION COMPLETE]")
        print(f"  analytic = {result['analytic']:.6f}")
        print(f"  numeric  = {result['numeric']:.6f}")
        print(f"  ratio    = {result['ratio']:.4f}  ({result['pct_diff']:+.1f}%)")
        return

    # ------------------------------------------------------------------
    # Smoke
    # ------------------------------------------------------------------
    if args.smoke:
        label, alpha, u0, w, grid_n = SMOKE_CONFIG[0]
        print(f"\n[SMOKE] 16^3 grid, alpha={alpha}, u0={u0}, w={w}, P={P}")
        try:
            row = run_one(label, alpha, u0, w, grid_n, SMOKE_OPTS, base_out, csv_path)
        except Exception as exc:
            import traceback
            print(f"\n  ERROR: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)
            _append_csv(row, csv_path)
        _print_summary([row], csv_path)
        return

    # ------------------------------------------------------------------
    # STEP 2: Stability-trend sweep
    # ------------------------------------------------------------------
    if args.trend:
        trend_csv = base_out / "trend_sweep.csv"
        _GRID40_OPTS = BreatherOpts(
            tol=1e-6, max_iter=50, inner_maxiter=1000, method="lgmres", verbose=True
        )
        rows = []
        for label, alpha, u0, w, grid_n in TREND_BRACKET:
            mass = RHO * (A ** 3)
            band_top_an = omega_band_top_analytic(K_S, mass, 3, A)
            run_opts = _GRID40_OPTS if grid_n >= 40 else TREND_OPTS
            try:
                row = run_one(label, alpha, u0, w, grid_n, run_opts, base_out,
                              trend_csv, band_top_override=band_top_an)
            except Exception as exc:
                import traceback
                print(f"\n  ERROR in {label}: {exc}")
                traceback.print_exc()
                row = _error_row(label, alpha, u0, w, grid_n)
                _append_csv(row, trend_csv)
            rows.append(row)
        _print_summary(rows, trend_csv)
        return

    # ------------------------------------------------------------------
    # STEP 3: Spot-check
    # ------------------------------------------------------------------
    if args.spot_check:
        trend_csv = base_out / "trend_sweep.csv"
        label, alpha, u0, w, grid_n = SPOT_CHECK_CONFIG
        mass = RHO * (A ** 3)
        band_top_an = omega_band_top_analytic(K_S, mass, 3, A)

        trend_json = base_out / "trend_a0p7_u4_w4_32_result.json"
        radius_1000: float | None = None
        if trend_json.exists():
            with open(trend_json) as f:
                trend_data = json.load(f)
            radius_1000 = trend_data.get("floquet_spectral_radius")
            res_1000 = trend_data.get("residual_final")
            print(f"[SPOT-CHECK] Loaded 1000-matvec result: "
                  f"Floquet radius={radius_1000:.4f}, residual={res_1000:.3e}")

        try:
            row = run_one(label, alpha, u0, w, grid_n, SPOT_CHECK_OPTS_2000,
                          base_out, trend_csv, band_top_override=band_top_an)
        except Exception as exc:
            import traceback
            print(f"\n  ERROR: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)
            _append_csv(row, trend_csv)

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
        return

    # ------------------------------------------------------------------
    # Single TREND_BRACKET entry
    # ------------------------------------------------------------------
    if args.trend_idx is not None:
        idx = args.trend_idx
        if idx < 0 or idx >= len(TREND_BRACKET):
            print(f"ERROR: --trend-idx must be 0..{len(TREND_BRACKET)-1}; got {idx}")
            sys.exit(1)
        trend_csv = base_out / "trend_sweep.csv"
        label, alpha, u0, w, grid_n = TREND_BRACKET[idx]
        mass = RHO * (A ** 3)
        band_top_an = omega_band_top_analytic(K_S, mass, 3, A)
        run_opts = (
            BreatherOpts(tol=1e-6, max_iter=50, inner_maxiter=1000,
                         method="lgmres", verbose=True)
            if grid_n >= 40 else TREND_OPTS
        )
        try:
            row = run_one(label, alpha, u0, w, grid_n, run_opts, base_out,
                          trend_csv, band_top_override=band_top_an)
        except Exception as exc:
            import traceback
            print(f"\n  ERROR: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)
            _append_csv(row, trend_csv)
        _print_summary([row], trend_csv)
        return

    # ------------------------------------------------------------------
    # AWS fat-soliton bracket (full or single point)
    # ------------------------------------------------------------------
    if args.aws_bracket or args.aws_idx is not None:
        if args.aws_idx is not None:
            idx = args.aws_idx
            if idx < 0 or idx >= len(AWS_BRACKET):
                print(f"ERROR: --aws-idx must be 0..{len(AWS_BRACKET)-1}; got {idx}")
                sys.exit(1)
            bracket = [AWS_BRACKET[idx]]
        else:
            bracket = AWS_BRACKET

        mass = RHO * (A ** 3)
        band_top_an = omega_band_top_analytic(K_S, mass, 3, A)

        print(f"\n[AWS BRACKET] {len(bracket)} point(s)")
        print(f"  opts: tol={AWS_OPTS.tol}, max_iter={AWS_OPTS.max_iter}, "
              f"inner_maxiter={AWS_OPTS.inner_maxiter}")
        print(f"  analytic band_top = {band_top_an:.6f}")
        print(f"  output_dir: {base_out}")

        rows = []
        for label, alpha, u0, w, grid_n in bracket:
            try:
                row = run_one(label, alpha, u0, w, grid_n, AWS_OPTS, base_out,
                              csv_path, band_top_override=band_top_an)
            except Exception as exc:
                import traceback
                print(f"\n  ERROR in {label}: {exc}")
                traceback.print_exc()
                row = _error_row(label, alpha, u0, w, grid_n)
                _append_csv(row, csv_path)
            rows.append(row)
        _print_summary(rows, csv_path)
        return

    # ------------------------------------------------------------------
    # Full 64^3 bracket or single point
    # ------------------------------------------------------------------
    if args.full:
        bracket = FULL_BRACKET
        opts = FULL_OPTS
        print(f"\n[FULL] {len(bracket)}-point bracket at 64^3.")
        print("[FULL] WARNING: expect hours per bracket point on a laptop.")
    else:
        idx = args.idx
        if idx < 0 or idx >= len(FULL_BRACKET):
            print(f"ERROR: --idx must be 0..{len(FULL_BRACKET)-1}; got {idx}")
            sys.exit(1)
        bracket = [FULL_BRACKET[idx]]
        opts = FULL_OPTS

    rows = []
    for label, alpha, u0, w, grid_n in bracket:
        try:
            row = run_one(label, alpha, u0, w, grid_n, opts, base_out, csv_path)
        except Exception as exc:
            import traceback
            print(f"\n  ERROR in {label}: {exc}")
            traceback.print_exc()
            row = _error_row(label, alpha, u0, w, grid_n)
            _append_csv(row, csv_path)
        rows.append(row)
    _print_summary(rows, csv_path)


if __name__ == "__main__":
    main()
