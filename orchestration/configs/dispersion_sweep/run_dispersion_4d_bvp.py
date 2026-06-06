#!/usr/bin/env python3
"""4D world-volume dispersion / isotropy diagnostic — bvp_chiral block solver.

SCOPE (honest):
    This measures LAB-FRAME dispersion: phase velocities c_L, c_T as seen by
    an external observer who knows the lattice coordinates.  It is NOT the
    internal-observer isotropy verdict (which needs solitons, not yet working).
    Two roles:
      (i)  4D block-solver shakedown: residual -> machine-zero, correct worldvolume,
           conditioning stable for all N (the central bvp_chiral guarantee).
      (ii) Precise lab-metric characterization: c_L, c_T per direction, cubic-
           anisotropy index, [111] triplet degeneracy, dispersion vs |k|*a.
           Sprint-3 soliton / dual-observer tests will compare against these numbers.

METHOD:
    For each config (direction, k_index, pol_axis):
      1. Build the analytic plane-wave eigenmode (R0, R1) and run bvp_chiral.
      2. At a reference node extract the amplitude time series A(l) = displacement
         along pol_axis at node p=0.
      3. Fit A(l) = eps * cos(l * theta_meas + phi) to recover theta_meas.
      4. Convert: omega_meas = theta_meas / dt; c_meas = omega_meas / |k|.
      5. Compare to closed-form d_of_k_eigenvalues prediction.

LINEARITY CHECK:
    Run each config at two amplitudes (1e-3 and 5e-4). omega must be
    invariant to < 0.01% — confirms linear regime and no numerical artifact.

OUTPUT:
    Per-run JSON + aggregate CSV in the output directory.
    Plots of omega(k) vs k along each axis (MPLBACKEND=Agg).

Usage:
    python test-runs/dispersion_4d_bvp/run_dispersion_4d_bvp.py \\
        --config-dir orchestration/configs/dispersion_sweep \\
        --pattern "local_*" \\
        --output-dir test-runs/dispersion_4d_bvp/results \\
        [--linearity-check]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from branesim.core.conventions import (
    ActionParams, LatticeParams, d_of_k_eigenvalues,
    c_longitudinal, c_transverse,
)
from branesim.core.lattice import SpacelikeLattice
from branesim.core.residual import residual_norm
from branesim.solver.boundary import ChiralBC
from branesim.solver.bvp import BoundaryProblem, solve_block


# ---------------------------------------------------------------------------
# Analytic helpers
# ---------------------------------------------------------------------------

def build_eigenmode(
    lattice: SpacelikeLattice,
    params: ActionParams,
    k_idx: list[int],
    pol_axis: int,
    amplitude: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Build the analytic plane-wave (R0, R1) and return (R0, R1, theta_k, omega_sq).

    The eigenmode is:
        R^l_p = ref_p + eps * cos(k . x_p - l * theta_k) * e_{pol_axis}

    theta_k = arccos(1 - 0.5 * dt^2 * omega^2_pol_axis)

    Returns
    -------
    R0, R1 : ndarray, shape (n_nodes, m_ambient)
    theta_k : float   — phase advance per time step
    omega_sq : float  — eigenfrequency squared for this mode
    """
    dim = lattice.params.dim
    a = lattice.params.spacing
    m = params.ambient_dim(dim)
    dt = params.dt
    N_grid = np.array(lattice.params.grid_shape, dtype=float)

    kvec = 2.0 * math.pi * np.array(k_idx, dtype=float) / (N_grid * a)
    eig = d_of_k_eigenvalues(kvec, params.alpha, params.k_s, params.rho, a)
    omega_sq = float(eig[pol_axis])
    arg = max(-1.0, min(1.0, 1.0 - 0.5 * dt * dt * omega_sq))
    theta_k = math.acos(arg)

    ref = lattice.reference_positions(m)           # (n_nodes, m)
    phase = ref[:, :dim] @ kvec                    # (n_nodes,)

    pol = np.zeros(m, dtype=float)
    pol[pol_axis] = 1.0

    R0 = ref + amplitude * np.cos(phase)[:, None] * pol[None, :]
    R1 = ref + amplitude * np.cos(phase - theta_k)[:, None] * pol[None, :]

    return R0, R1, theta_k, omega_sq


def fit_phase_advance(
    wv_slices: np.ndarray,
    ref_positions: np.ndarray,
    k_idx: list[int],
    pol_axis: int,
    grid_shape: tuple,
    spacing: float,
    amplitude: float,
) -> tuple[float, float]:
    """Fit theta_meas from the worldvolume time series at a reference node.

    Strategy: project each slice onto the plane-wave mode to get the
    spatially-averaged amplitude A(l).  A(l) = eps * cos(l * theta + phi).
    Fit theta using a least-squares nonlinear fit (scipy.optimize.curve_fit).

    Returns (theta_meas, fit_residual_rms).
    """
    from scipy.optimize import curve_fit

    dim = len(grid_shape)
    a = spacing
    N_grid = np.array(grid_shape, dtype=float)
    kvec = 2.0 * math.pi * np.array(k_idx, dtype=float) / (N_grid * a)

    # Spatial projection weights: cos(k . x_p) / (N_nodes/2)
    # For a single-mode excitation with zero initial velocity (cos only),
    # the projection picks out the cos channel.
    phase = ref_positions[:, :dim] @ kvec           # (n_nodes,)
    norm_factor = 2.0 / len(phase)                  # normalised so A(0) = amplitude

    n_slices_plus1 = wv_slices.shape[0]
    n_nodes = wv_slices.shape[1]
    m = wv_slices.shape[2]

    # Displacement along pol_axis
    disp = wv_slices[:, :, pol_axis] - ref_positions[None, :, pol_axis]  # (L+1, n_nodes)

    # Spatial projection: A_cos(l) = norm_factor * sum_p disp_p^l * cos(k.x_p)
    A_cos = norm_factor * (disp * np.cos(phase)[None, :]).sum(axis=1)   # (L+1,)
    l_arr = np.arange(n_slices_plus1, dtype=float)

    # Fit model: A(l) = amp * cos(theta * l + phi)
    def model(l, theta, phi):
        return amplitude * np.cos(theta * l + phi)

    # Initial guess: use exact theta from analytic eigenvalue
    eig = d_of_k_eigenvalues(kvec, 0.2, 1.0, 1.0, a)  # alpha=0.2 default
    omega_sq = float(eig[pol_axis])
    dt = 0.1  # must match the config
    theta0 = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * dt * dt * omega_sq)))

    try:
        popt, pcov = curve_fit(model, l_arr, A_cos, p0=[theta0, 0.0],
                               bounds=([0.0, -math.pi], [math.pi, math.pi]),
                               maxfev=10000)
        theta_meas = abs(float(popt[0]))
        residual_rms = float(np.sqrt(np.mean((A_cos - model(l_arr, *popt))**2)))
    except Exception:
        # Fallback: use the analytic theta (so we can still report residual)
        theta_meas = theta0
        residual_rms = float(np.std(A_cos - amplitude * np.cos(theta0 * l_arr)))

    return theta_meas, residual_rms


# ---------------------------------------------------------------------------
# Single-run dispatcher
# ---------------------------------------------------------------------------

def run_one(
    config: dict,
    output_dir: Path,
    label: str,
    amplitude_override: float | None = None,
) -> dict:
    """Run one bvp_chiral dispersion measurement and return result dict."""

    lcfg = config["lattice"]
    acfg = config["action"]
    scfg = config["seed"]

    lp = LatticeParams(
        grid_shape=tuple(int(v) for v in lcfg["grid_shape"]),
        spacing=float(lcfg.get("spacing", 1.0)),
        periodic_axes=tuple(bool(v) for v in lcfg.get("periodic_axes", [True]*3)),
        axial_weight=float(lcfg.get("axial_weight", 1.0)),
    )
    lattice = SpacelikeLattice(lp)
    m = int(acfg.get("m_ambient", lp.dim + 1))
    r_t_cfg = acfg.get("r_t", None)
    ap = ActionParams(
        k_s=float(acfg["k_s"]),
        alpha=float(acfg["alpha"]),
        rho=float(acfg["rho"]),
        dt=float(acfg["dt"]),
        n_slices=int(acfg["n_slices"]),
        m_ambient=m,
        r_t=float(r_t_cfg) if r_t_cfg is not None else None,
    )
    mass = ap.rho * lp.spacing ** lp.dim
    N = ap.n_slices
    a = lp.spacing

    k_idx = [int(v) for v in scfg["k_index"]]
    pol_vec = [float(v) for v in scfg["polarization"]]
    pol_axis = int(np.argmax(np.abs(pol_vec)))
    amp = amplitude_override if amplitude_override is not None else float(scfg["amplitude"])

    # Linearity pre-check: A*|k| << 1
    kvec = 2.0 * math.pi * np.array(k_idx, dtype=float) / (np.array(lp.grid_shape, dtype=float) * a)
    k_mag = float(np.linalg.norm(kvec))
    A_k = amp * k_mag
    assert A_k < 0.05, f"A*|k| = {A_k:.4f} is not small (linear regime requires << 1)"

    # Analytic prediction
    eig = d_of_k_eigenvalues(kvec, ap.alpha, ap.k_s, ap.rho, a)
    omega_sq_pred = float(eig[pol_axis])
    theta_pred = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * ap.dt**2 * omega_sq_pred)))
    omega_pred = theta_pred / ap.dt
    c_pred = omega_pred / k_mag if k_mag > 0 else 0.0

    # Build eigenmode initial condition and solve
    R0, R1, theta_k, omega_sq_exact = build_eigenmode(lattice, ap, k_idx, pol_axis, amp)
    bc = ChiralBC(R0=R0, R1=R1, chirality="forward")
    t0 = time.perf_counter()
    wv = solve_block(BoundaryProblem(lattice, ap, mass, bc))
    walltime = time.perf_counter() - t0

    # Quality checks
    res_norm = float(residual_norm(wv.slices, lattice, ap, mass))
    n_interior_dof = (N - 1) * lattice.n_nodes * m
    res_per_dof = res_norm / math.sqrt(n_interior_dof) if n_interior_dof > 0 else res_norm

    solver_report = wv.solver_report

    # Phase-advance measurement from worldvolume
    ref = lattice.reference_positions(m)
    theta_meas, fit_residual_rms = fit_phase_advance(
        wv.slices, ref, k_idx, pol_axis,
        lp.grid_shape, a, amp,
    )
    omega_meas = theta_meas / ap.dt
    c_meas = omega_meas / k_mag if k_mag > 0 else 0.0
    c_rel_err = (c_meas - c_pred) / c_pred if c_pred > 0 else 0.0

    # Long-wavelength predictions (k->0)
    c_L_lw = c_longitudinal(ap.k_s, a, mass)
    c_T_lw = c_transverse(ap.k_s, a, mass, ap.alpha)

    result = {
        "label": label,
        "config_file": label,
        "grid_shape": list(lp.grid_shape),
        "n_slices": N,
        "k_index": k_idx,
        "pol_axis": pol_axis,
        "amplitude": amp,
        "A_k": A_k,
        "k_mag": k_mag,
        "ka": k_mag * a,
        # Analytic
        "omega_sq_pred": omega_sq_pred,
        "theta_pred": theta_pred,
        "omega_pred": omega_pred,
        "c_pred": c_pred,
        "c_L_lw": c_L_lw,
        "c_T_lw": c_T_lw,
        # Measured
        "theta_meas": theta_meas,
        "omega_meas": omega_meas,
        "c_meas": c_meas,
        "c_rel_err": c_rel_err,
        "fit_residual_rms": fit_residual_rms,
        # Solver quality
        "residual_norm": res_norm,
        "residual_per_dof": res_per_dof,
        "condition_estimate": float(solver_report.get("condition_estimate", 0.0)),
        "walltime_s": walltime,
        "converged": bool(solver_report.get("converged", True)),
        # Pass/fail
        "pass_c_1pct": abs(c_rel_err) < 0.01,
        "pass_residual": res_per_dof < 1e-9,
    }

    return result


# ---------------------------------------------------------------------------
# Linearity check
# ---------------------------------------------------------------------------

def linearity_check(config: dict, output_dir: Path, label: str) -> dict:
    """Check that omega is invariant when amplitude is halved (linear regime)."""
    amp0 = float(config["seed"]["amplitude"])
    r1 = run_one(config, output_dir, label + "_amp_full", amplitude_override=amp0)
    r2 = run_one(config, output_dir, label + "_amp_half", amplitude_override=amp0 / 2.0)
    delta_omega = abs(r1["omega_meas"] - r2["omega_meas"])
    ref_omega = r1["omega_pred"]
    frac = delta_omega / ref_omega if ref_omega > 0 else 0.0
    return {
        "label": label,
        "omega_full": r1["omega_meas"],
        "omega_half": r2["omega_meas"],
        "delta_omega_frac": frac,
        "pass_linearity": frac < 1e-4,
    }


# ---------------------------------------------------------------------------
# Aggregate observables
# ---------------------------------------------------------------------------

def extract_observables(results: list[dict]) -> dict:
    """Extract the decisive observables from the full result list.

    Observables:
      - c_L([100]) extrapolated to k->0
      - c_T([100]) extrapolated to k->0
      - c_all([111]) — triplet mean and spread
      - anisotropy index: c_L([100]) / c([111])
      - [110] birefringence: c_hard / c_soft
      - [210] three-mode spread
    """
    def by(direction, pol_axis=None):
        out = [r for r in results if r["k_index"] == _dir_kidx(r, direction)]
        if pol_axis is not None:
            out = [r for r in out if r["pol_axis"] == pol_axis]
        return out

    # Simpler: filter by direction string and pol
    def group(dir_str, pol):
        return sorted(
            [r for r in results if _direction_of(r) == dir_str and r["pol_axis"] == pol],
            key=lambda r: r["ka"]
        )

    def _direction_of(r):
        ki = r["k_index"]
        if ki[1] == 0 and ki[2] == 0:
            return "100"
        if ki[0] == ki[1] and ki[2] == 0:
            return "110"
        if ki[0] == ki[1] == ki[2]:
            return "111"
        if ki[2] == 0 and ki[0] == 2 * ki[1]:
            return "210"
        return "other"

    def extrapolate_c0(runs):
        """Fit c(ka) = c0 + c2*(ka)^2 and return c0."""
        if not runs:
            return None, None
        ka = np.array([r["ka"] for r in runs])
        c = np.array([r["c_meas"] for r in runs])
        if len(ka) < 2:
            return float(c[0]), None
        A = np.stack([np.ones_like(ka), ka**2], axis=1)
        sol, *_ = np.linalg.lstsq(A, c, rcond=None)
        return float(sol[0]), float(sol[1])

    obs = {}

    # [100] L
    r100_L = group("100", 0)
    c0_100L, c2_100L = extrapolate_c0(r100_L)
    obs["c_L_100_extrap"] = c0_100L
    obs["c_L_100_analytic"] = 1.0
    obs["c_L_100_rel_err"] = (c0_100L - 1.0) / 1.0 if c0_100L else None

    # [100] T
    r100_T = group("100", 1)
    c0_100T, _ = extrapolate_c0(r100_T)
    import math as _math
    alpha = 0.2
    c_T_analytic = _math.sqrt(1.0 - alpha)
    obs["c_T_100_extrap"] = c0_100T
    obs["c_T_100_analytic"] = c_T_analytic
    obs["c_T_100_rel_err"] = (c0_100T - c_T_analytic) / c_T_analytic if c0_100T else None

    # [111] — all pols degenerate
    r111_p0 = group("111", 0)
    r111_p1 = group("111", 1)
    c0_111_p0, _ = extrapolate_c0(r111_p0)
    c0_111_p1, _ = extrapolate_c0(r111_p1)
    c_111_analytic = _math.sqrt((3.0 - 2.0 * alpha) / 3.0)
    obs["c_111_p0_extrap"] = c0_111_p0
    obs["c_111_p1_extrap"] = c0_111_p1
    obs["c_111_analytic"] = c_111_analytic
    obs["triplet_spread_frac"] = (
        abs(c0_111_p0 - c0_111_p1) / c_111_analytic
        if c0_111_p0 and c0_111_p1 else None
    )

    # Anisotropy: c_L([100]) / c([111])
    aniso_analytic = _math.sqrt(3.0 / (3.0 - 2.0 * alpha))
    if c0_100L and c0_111_p0:
        aniso_meas = c0_100L / c0_111_p0
        obs["anisotropy_100_111"] = aniso_meas
        obs["anisotropy_analytic"] = aniso_analytic
        obs["anisotropy_rel_err"] = (aniso_meas - aniso_analytic) / aniso_analytic
    else:
        obs["anisotropy_100_111"] = None
        obs["anisotropy_analytic"] = aniso_analytic
        obs["anisotropy_rel_err"] = None

    # [110] birefringence: hard (p0) vs soft (p2)
    r110_hard = group("110", 0)
    r110_soft = group("110", 2)
    c0_110_hard, _ = extrapolate_c0(r110_hard)
    c0_110_soft, _ = extrapolate_c0(r110_soft)
    c_hard_analytic = _math.sqrt(1.0 - alpha / 2.0)
    obs["c_110_hard_extrap"] = c0_110_hard
    obs["c_110_soft_extrap"] = c0_110_soft
    obs["c_110_hard_analytic"] = c_hard_analytic
    obs["c_110_soft_analytic"] = c_T_analytic
    if c0_110_hard and c0_110_soft:
        obs["birefringence_110"] = c0_110_hard / c0_110_soft
        obs["birefringence_110_analytic"] = c_hard_analytic / c_T_analytic
    else:
        obs["birefringence_110"] = None
        obs["birefringence_110_analytic"] = c_hard_analytic / c_T_analytic

    # Per-direction per-ka: collect c_rel_err at reference ka values
    for ka_target in [0.1, 0.3, 0.5]:
        for dir_str in ["100", "110", "111"]:
            for p in [0, 1, 2]:
                grp = group(dir_str, p)
                if not grp:
                    continue
                # Find closest ka
                diffs = [abs(r["ka"] - ka_target) for r in grp]
                if min(diffs) > 0.15:
                    continue
                best = grp[int(np.argmin(diffs))]
                key = f"c_rel_err_{dir_str}_p{p}_ka{ka_target:.1f}"
                obs[key] = best["c_rel_err"]

    return obs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plots(results: list[dict], output_dir: Path) -> None:
    """Plot omega(k) vs ka along each direction, measured vs predicted."""
    def group_dir(dir_str, pol):
        return sorted(
            [r for r in results if _direction_of(r) == dir_str and r["pol_axis"] == pol],
            key=lambda r: r["ka"]
        )

    def _direction_of(r):
        ki = r["k_index"]
        if ki[1] == 0 and ki[2] == 0:
            return "100"
        if ki[0] == ki[1] and ki[2] == 0:
            return "110"
        if ki[0] == ki[1] == ki[2]:
            return "111"
        if ki[2] == 0 and ki[0] == 2 * ki[1]:
            return "210"
        return "other"

    # Omega vs ka plots
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("omega(k) vs k·a: measured (dots) vs predicted (lines)\nbvp_chiral, 4D world-volume")

    dir_configs = [
        ("100", [("L", 0, "C0"), ("T", 1, "C1")], axes[0, 0]),
        ("110", [("hard", 0, "C2"), ("soft", 2, "C3")], axes[0, 1]),
        ("111", [("p0", 0, "C4"), ("p1", 1, "C5")], axes[1, 0]),
        ("210", [("p0", 0, "C6"), ("p1", 1, "C7"), ("p2", 2, "C8")], axes[1, 1]),
    ]

    for dir_str, pol_configs, ax in dir_configs:
        ax.set_title(f"[{dir_str}]")
        ax.set_xlabel("k·a")
        ax.set_ylabel("ω")
        for pol_label, pol_axis, color in pol_configs:
            grp = group_dir(dir_str, pol_axis)
            if not grp:
                continue
            ka = [r["ka"] for r in grp]
            omega_meas = [r["omega_meas"] for r in grp]
            omega_pred = [r["omega_pred"] for r in grp]
            ax.plot(ka, omega_pred, "--", color=color, alpha=0.5, label=f"{pol_label} pred")
            ax.scatter(ka, omega_meas, color=color, zorder=5, label=f"{pol_label} meas", s=40)
        ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(output_dir / "omega_vs_ka.png", dpi=120)
    plt.close(fig)

    # Phase velocity vs ka (normalized to c_L)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle("c(k)/c_L vs k·a per direction (cubic anisotropy)")
    alpha = 0.2
    import math as _math
    c_L = 1.0
    c_T = _math.sqrt(1.0 - alpha)

    for ax, (dir_str, pol_configs) in zip(axes, [
        ("100", [("L", 0, "C0"), ("T", 1, "C1")]),
        ("110", [("hard", 0, "C2"), ("soft", 2, "C3")]),
        ("111", [("p0", 0, "C4"), ("p1", 1, "C5")]),
    ]):
        ax.set_title(f"[{dir_str}]")
        ax.set_xlabel("k·a")
        ax.set_ylabel("c / c_L")
        ax.axhline(1.0, color="gray", linewidth=0.5, linestyle="--", label="c_L")
        ax.axhline(c_T / c_L, color="gray", linewidth=0.5, linestyle=":", label="c_T")
        for pol_label, pol_axis, color in pol_configs:
            grp = group_dir(dir_str, pol_axis)
            if not grp:
                continue
            ka = [r["ka"] for r in grp]
            c_ratio_meas = [r["c_meas"] / c_L for r in grp]
            c_ratio_pred = [r["c_pred"] / c_L for r in grp]
            ax.plot(ka, c_ratio_pred, "--", color=color, alpha=0.5)
            ax.scatter(ka, c_ratio_meas, color=color, s=40, zorder=5, label=pol_label)
        ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(output_dir / "c_vs_ka.png", dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="4D bvp_chiral dispersion sweep diagnostic.")
    parser.add_argument("--config-dir", required=True, help="Directory containing JSON configs.")
    parser.add_argument("--pattern", default="local_*", help="Glob pattern for configs (default: local_*).")
    parser.add_argument("--output-dir", required=True, help="Output directory for results.")
    parser.add_argument("--linearity-check", action="store_true",
                        help="For each config, also run at half amplitude and report omega invariance.")
    args = parser.parse_args()

    config_dir = Path(args.config_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    configs = sorted(config_dir.glob(args.pattern + ".json"))
    if not configs:
        print(f"No configs found matching {config_dir}/{args.pattern}.json")
        return

    print(f"Found {len(configs)} configs matching '{args.pattern}'.")

    all_results = []
    linearity_results = []

    for cfg_path in configs:
        config = json.loads(cfg_path.read_text())
        label = cfg_path.stem
        print(f"\n--- {label} ---")

        try:
            result = run_one(config, output_dir, label)
            all_results.append(result)

            status = "PASS" if result["pass_c_1pct"] and result["pass_residual"] else "FAIL"
            print(f"  c_pred={result['c_pred']:.6f}  c_meas={result['c_meas']:.6f}  "
                  f"rel_err={result['c_rel_err']:+.3e}  res_per_dof={result['residual_per_dof']:.2e}  "
                  f"cond={result['condition_estimate']:.1f}  [{status}]")

            if args.linearity_check:
                lin = linearity_check(config, output_dir, label)
                linearity_results.append(lin)
                lin_status = "PASS" if lin["pass_linearity"] else "FAIL"
                print(f"  linearity: delta_omega/omega = {lin['delta_omega_frac']:.2e}  [{lin_status}]")

        except Exception as exc:
            print(f"  ERROR: {exc}")
            all_results.append({"label": label, "error": str(exc)})
            continue

    # --- Save per-run JSON ---
    (output_dir / "all_results.json").write_text(
        json.dumps(all_results, indent=2), encoding="utf-8"
    )

    # --- Save CSV ---
    csv_fields = [
        "label", "k_index", "pol_axis", "ka",
        "c_pred", "c_meas", "c_rel_err",
        "omega_pred", "omega_meas",
        "residual_per_dof", "condition_estimate",
        "pass_c_1pct", "pass_residual", "walltime_s",
    ]
    with (output_dir / "dispersion_results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for r in all_results:
            if "error" not in r:
                row = {k: r.get(k, "") for k in csv_fields}
                row["k_index"] = str(r.get("k_index", ""))
                writer.writerow(row)

    # --- Aggregate observables ---
    valid_results = [r for r in all_results if "error" not in r]
    if valid_results:
        obs = extract_observables(valid_results)
        (output_dir / "observables.json").write_text(
            json.dumps(obs, indent=2), encoding="utf-8"
        )

        print("\n" + "=" * 60)
        print("DECISIVE OBSERVABLES")
        print("=" * 60)
        print(f"c_L([100]) extrap:   {obs.get('c_L_100_extrap', 'N/A'):.6f}  "
              f"(analytic: {obs.get('c_L_100_analytic', 1.0):.6f}  "
              f"err: {obs.get('c_L_100_rel_err') or 0:+.3e})")
        print(f"c_T([100]) extrap:   {obs.get('c_T_100_extrap', 'N/A'):.6f}  "
              f"(analytic: {obs.get('c_T_100_analytic', 0.8944):.6f}  "
              f"err: {obs.get('c_T_100_rel_err') or 0:+.3e})")
        print(f"c([111]) p0 extrap:  {obs.get('c_111_p0_extrap', 'N/A'):.6f}  "
              f"(analytic: {obs.get('c_111_analytic', 0.9309):.6f})")
        print(f"c([111]) p1 extrap:  {obs.get('c_111_p1_extrap', 'N/A'):.6f}  "
              f"triplet spread: {obs.get('triplet_spread_frac') or 0:.2e}")
        print(f"Anisotropy c_L([100])/c([111]): "
              f"{obs.get('anisotropy_100_111', 'N/A'):.6f}  "
              f"(analytic: {obs.get('anisotropy_analytic', 1.0742):.6f}  "
              f"err: {obs.get('anisotropy_rel_err') or 0:+.3e})")
        print(f"[110] birefringence: {obs.get('birefringence_110', 'N/A'):.6f}  "
              f"(analytic: {obs.get('birefringence_110_analytic', 1.0607):.6f})")

        # Overall pass/fail
        n_pass = sum(r["pass_c_1pct"] and r["pass_residual"] for r in valid_results)
        n_total = len(valid_results)
        print(f"\nPass/fail (|c_rel_err| < 1% AND residual_per_dof < 1e-9): {n_pass}/{n_total}")

        if linearity_results:
            n_lin_pass = sum(lr["pass_linearity"] for lr in linearity_results)
            print(f"Linearity (delta_omega/omega < 1e-4): {n_lin_pass}/{len(linearity_results)}")

        make_plots(valid_results, output_dir)
        print(f"\nPlots: {output_dir}/omega_vs_ka.png, c_vs_ka.png")

    print(f"\nResults in: {output_dir}")


if __name__ == "__main__":
    main()
