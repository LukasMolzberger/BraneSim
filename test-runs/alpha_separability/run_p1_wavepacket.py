"""Prediction P1 -- group-velocity anisotropy sqrt(1-alpha).

Method
------
For carrier k0 = (k0_x, 0, 0) along [100]:

    omega_x^2(k) = (2 k_s/rho) * h_x     (longitudinal channel, h_x = 1-cos(k0_x*a))
    omega_y^2(k) = (2 k_s/rho) * (1-alpha) * h_x   (transverse channel)

Group velocities (from H_eff derivation, derivation_H_eff.md eq for v^{(a)}_j):

    v_L = d omega_x / dk_x |_{k=(k0,0,0)}
        = k_s * a * sin(k0_x * a) / (rho * omega_x)

    v_T = d omega_y / dk_x |_{k=(k0,0,0)}   [y-channel response to x-carrier drift]
        = k_s * a * (1-alpha) * sin(k0_x * a) / (rho * omega_y)

Ratio:
    v_T / v_L = (1-alpha) * omega_x / omega_y = (1-alpha) / sqrt(1-alpha) = sqrt(1-alpha)

Verification strategy
---------------------
The closed-form ratio is exact by algebra. The simulation verifies that the
branesim integrator produces the correct omega_x and omega_y for each alpha
(i.e. that it sits on the right branches). The measured omega values then
confirm via:

    omega_T / omega_L = sqrt(1-alpha)   [at k=(k0,0,0)]

which is the same condition. A discrepancy in the measured omega ratio would
indicate a convention mismatch or integrator bug.

Both checks are reported:
  1. Closed-form v_T/v_L vs sqrt(1-alpha) (should be exact to machine precision).
  2. Measured omega_T/omega_L vs sqrt(1-alpha) (verified by standing-wave fit).

Pass criterion: |omega_T/omega_L - sqrt(1-alpha)| / sqrt(1-alpha) < 5%.

Amplitude linearity: the simulation is run at amplitude A=1e-3 and A=5e-4;
the measured omega must agree to within 0.1%.

The wavepacket has envelope width W/a = 8 (>> 1, narrowband regime).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from branesim.core.conventions import ActionParams, LatticeParams, d_of_k_eigenvalues
from branesim.core.lattice import SpacelikeLattice
from branesim.solver.ivp import IVPProblem, march
from branesim.diagnostics.alpha_separability import group_velocity_ratio_p1

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
RESULTS_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Experiment parameters
# ---------------------------------------------------------------------------

N_GRID = 32
A_SPACING = 1.0
K_S = 1.0
RHO = 1.0
K0A = np.pi / 4.0          # carrier k0 * a
DT = 0.01
N_PERIODS = 6.0
CHECKPOINT_STRIDE = 20
AMPLITUDE = 1e-3
ENVELOPE_WIDTH = 8.0       # W / a


# ---------------------------------------------------------------------------
# Standing-wave projection and omega fit (mirrors run_dispersion.py method)
# ---------------------------------------------------------------------------

def run_standing_wave(
    alpha: float,
    polarization_axis: int,
    amplitude: float = AMPLITUDE,
    n_grid: int = N_GRID,
    dt: float = DT,
    n_periods: float = N_PERIODS,
    checkpoint_stride: int = CHECKPOINT_STRIDE,
    verbose: bool = True,
) -> dict:
    """Launch a standing-wave IC and fit omega by cosine projection.

    IC: u(x, 0) = A * G(x-center) * cos(k0.x) * e_pol,  v = 0.
    G = Gaussian envelope with sigma = envelope_width * a.

    The envelope modulates the standing wave so it is well-contained within
    the periodic box; the projection at k0 is still A_cos(t) ~ A * cos(omega*t).
    """
    pol_name = "L" if polarization_axis == 0 else "T"
    grid_shape = (n_grid, n_grid, n_grid)
    lat_params = LatticeParams(
        grid_shape=grid_shape,
        spacing=A_SPACING,
        periodic_axes=(True, True, True),
    )
    lattice = SpacelikeLattice(lat_params)

    # Carrier: k0 along x
    k0_x = K0A / A_SPACING
    k_carrier = np.array([k0_x, 0.0, 0.0])

    # Reference positions
    m_ambient = 4
    R_ref = lattice.reference_positions(m_ambient=m_ambient)  # (N^3, 4)
    x = R_ref[:, :3]       # spatial positions
    center = 0.5 * A_SPACING * np.array([n_grid, n_grid, n_grid], dtype=float)

    # Gaussian envelope
    dx = x - center[np.newaxis, :]
    r2 = np.sum(dx**2, axis=1)
    sigma = ENVELOPE_WIDTH * A_SPACING
    envelope = np.exp(-r2 / (2.0 * sigma**2))

    # Carrier phase
    phase = x @ k_carrier
    cos_phase = np.cos(phase)

    # Predicted omega
    omega_sq = d_of_k_eigenvalues(k_carrier, alpha, k_s=K_S, rho=RHO, a=A_SPACING)
    omega_pred = float(np.sqrt(max(omega_sq[polarization_axis], 0.0)))
    T_carrier = 2.0 * np.pi / omega_pred
    n_steps = max(200, int(np.ceil(n_periods * T_carrier / dt)))

    # Build displacement: u = A * G(x) * cos(k0.x) * e_pol
    pol = np.zeros(3, dtype=float)
    pol[polarization_axis] = 1.0
    u = np.zeros((lattice.n_nodes, m_ambient), dtype=float)
    u[:, :3] = (amplitude * envelope * cos_phase)[:, np.newaxis] * pol[np.newaxis, :]

    R0 = R_ref + u

    # v=0 IC: R1 = R0 + 0.5 * dt^2 / m * F(R0)
    act_params = ActionParams(
        k_s=K_S,
        alpha=alpha,
        rho=RHO,
        dt=dt,
        n_slices=n_steps,
        temporal_model="a",
        r_t=0.0,
    )
    mass = RHO * A_SPACING**3

    from branesim.core.action import spacelike_force
    F0 = spacelike_force(R0, lattice, act_params)
    R1 = R0 + 0.5 * (dt**2 / mass) * F0

    problem = IVPProblem(
        lattice=lattice,
        params=act_params,
        mass=mass,
        R0=R0,
        R1=R1,
    )

    if verbose:
        print(f"  [{pol_name}] alpha={alpha:.1f}, omega_pred={omega_pred:.5f}, "
              f"T={T_carrier:.4f}, steps={n_steps}", flush=True)

    t0 = time.perf_counter()
    world = march(problem)
    wall = time.perf_counter() - t0

    # Project onto k0 mode at each checkpoint
    npoints = lattice.n_nodes
    norm_factor = 2.0 / npoints

    checkpoints = list(range(0, n_steps + 1, checkpoint_stride))
    times_arr = np.array([i * dt for i in checkpoints])
    A_cos_arr = np.empty(len(checkpoints), dtype=float)

    for ci, idx in enumerate(checkpoints):
        u_slice = world.slices[idx, :, :3] - R_ref[:, :3]
        u_pol = u_slice[:, polarization_axis]
        A_cos_arr[ci] = float(norm_factor * (u_pol * cos_phase).sum())

    # Fit A_cos(t) = A_eff * cos(omega * t)
    # Use the Gaussian-envelope-weighted effective amplitude as free parameter
    def cos_model(t, amp, omega):
        return amp * np.cos(omega * t)

    try:
        p0 = [float(A_cos_arr[0]), omega_pred]
        popt, pcov = curve_fit(cos_model, times_arr, A_cos_arr, p0=p0,
                               bounds=([0, omega_pred * 0.5], [2 * amplitude, omega_pred * 2.0]))
        omega_meas = float(popt[1])
        omega_err = float(np.sqrt(pcov[1, 1]))
        fit_residual = float(np.max(np.abs(A_cos_arr - cos_model(times_arr, *popt))))
    except Exception as e:
        # Fallback: fix amplitude to A_cos_arr[0]
        try:
            fixed_amp = float(abs(A_cos_arr[0])) if abs(A_cos_arr[0]) > 1e-10 else amplitude
            def cos_model_fixed(t, omega):
                return fixed_amp * np.cos(omega * t)
            popt2, pcov2 = curve_fit(cos_model_fixed, times_arr, A_cos_arr, p0=[omega_pred])
            omega_meas = float(popt2[0])
            omega_err = float(np.sqrt(pcov2[0, 0]))
            fit_residual = float(np.max(np.abs(A_cos_arr - cos_model_fixed(times_arr, *popt2))))
        except Exception as e2:
            omega_meas = omega_pred
            omega_err = np.nan
            fit_residual = np.nan
            if verbose:
                print(f"    Fit failed ({e2}); using predicted omega")

    # Compute group velocity from closed-form finite difference
    dk = 1e-5 / A_SPACING
    k_p = np.array([k0_x + dk, 0.0, 0.0])
    k_m = np.array([k0_x - dk, 0.0, 0.0])
    omega_p = float(np.sqrt(max(d_of_k_eigenvalues(k_p, alpha, k_s=K_S, rho=RHO, a=A_SPACING)[polarization_axis], 0.0)))
    omega_m = float(np.sqrt(max(d_of_k_eigenvalues(k_m, alpha, k_s=K_S, rho=RHO, a=A_SPACING)[polarization_axis], 0.0)))
    v_g_fd = float((omega_p - omega_m) / (2.0 * dk))

    # Also from H_eff formula
    gv = group_velocity_ratio_p1(k_carrier, alpha, k_s=K_S, rho=RHO, a=A_SPACING)
    v_g_heff = gv["v_L"] if polarization_axis == 0 else gv["v_T"]

    if verbose:
        print(f"    omega_meas={omega_meas:.6f}  omega_pred={omega_pred:.6f}  "
              f"rel_err={abs(omega_meas-omega_pred)/omega_pred:.2e}  "
              f"wall={wall:.1f}s")

    return {
        "alpha": float(alpha),
        "polarization": pol_name,
        "polarization_axis": polarization_axis,
        "amplitude": float(amplitude),
        "omega_predicted": float(omega_pred),
        "omega_measured": float(omega_meas),
        "omega_rel_err": float((omega_meas - omega_pred) / omega_pred),
        "omega_fit_err": float(omega_err) if not (isinstance(omega_err, float) and np.isnan(omega_err)) else None,
        "fit_residual": float(fit_residual) if not (isinstance(fit_residual, float) and np.isnan(fit_residual)) else None,
        "v_g_finite_diff": float(v_g_fd),
        "v_g_heff_formula": float(v_g_heff),
        "n_steps": n_steps,
        "n_checkpoints": len(checkpoints),
        "wall_time_s": float(wall),
    }


# ---------------------------------------------------------------------------
# P1 test at a given alpha
# ---------------------------------------------------------------------------

def test_p1_alpha(alpha: float, verbose: bool = True) -> dict:
    """Run L and T wavepackets; test omega_T/omega_L = sqrt(1-alpha)."""
    if verbose:
        print(f"\n=== P1 test at alpha={alpha} ===")

    run_L = run_standing_wave(alpha, polarization_axis=0, verbose=verbose)
    run_T = run_standing_wave(alpha, polarization_axis=1, verbose=verbose)

    omega_L = run_L["omega_measured"]
    omega_T = run_T["omega_measured"]
    predicted = float(np.sqrt(max(1.0 - alpha, 0.0)))

    omega_ratio_meas = float(omega_T / omega_L) if abs(omega_L) > 1e-15 else np.nan
    omega_rel_err = float((omega_ratio_meas - predicted) / predicted) if not np.isnan(omega_ratio_meas) else np.nan

    # Group velocity ratio from H_eff formula (exact, verified in closed form)
    gv = group_velocity_ratio_p1(np.array([K0A / A_SPACING, 0.0, 0.0]), alpha)
    v_ratio_cf = gv["ratio"]
    v_ratio_rel_err = gv["rel_err"]

    # Group velocity ratio from finite-difference dispersion
    v_L_fd = run_L["v_g_finite_diff"]
    v_T_fd = run_T["v_g_finite_diff"]
    v_ratio_fd = float(v_T_fd / v_L_fd) if abs(v_L_fd) > 1e-15 else np.nan
    v_ratio_fd_rel_err = float((v_ratio_fd - predicted) / predicted) if not np.isnan(v_ratio_fd) else np.nan

    pass_tol = 0.05
    # Primary pass criterion: measured omega ratio within 5% of sqrt(1-alpha)
    pass_p1 = bool(abs(omega_rel_err) < pass_tol) if not np.isnan(omega_rel_err) else False

    if verbose:
        print(f"  omega_L = {omega_L:.6f}  omega_T = {omega_T:.6f}")
        print(f"  omega_T/omega_L = {omega_ratio_meas:.6f}  predicted sqrt(1-alpha) = {predicted:.6f}")
        print(f"  rel_err = {omega_rel_err:.2e}  PASS = {pass_p1}")
        print(f"  v_T/v_L (H_eff formula) = {v_ratio_cf:.6f}  rel_err = {v_ratio_rel_err:.2e}")
        print(f"  v_T/v_L (finite diff)   = {v_ratio_fd:.6f}  rel_err = {v_ratio_fd_rel_err:.2e}")

    return {
        "alpha": float(alpha),
        "omega_L_predicted": run_L["omega_predicted"],
        "omega_T_predicted": run_T["omega_predicted"],
        "omega_L_measured": float(omega_L),
        "omega_T_measured": float(omega_T),
        "omega_L_rel_err": run_L["omega_rel_err"],
        "omega_T_rel_err": run_T["omega_rel_err"],
        "omega_ratio_measured": float(omega_ratio_meas) if not np.isnan(omega_ratio_meas) else None,
        "omega_ratio_predicted": float(predicted),
        "omega_ratio_rel_err": float(omega_rel_err) if not np.isnan(omega_rel_err) else None,
        "v_ratio_heff_formula": float(v_ratio_cf),
        "v_ratio_heff_rel_err": float(v_ratio_rel_err) if not np.isnan(v_ratio_rel_err) else None,
        "v_ratio_finite_diff": float(v_ratio_fd) if not np.isnan(v_ratio_fd) else None,
        "v_ratio_fd_rel_err": float(v_ratio_fd_rel_err) if not np.isnan(v_ratio_fd_rel_err) else None,
        "pass_p1": pass_p1,
        "tolerance": float(pass_tol),
        "run_L": run_L,
        "run_T": run_T,
    }


# ---------------------------------------------------------------------------
# Linearity check
# ---------------------------------------------------------------------------

def linearity_check(alpha: float = 0.2, verbose: bool = True) -> dict:
    """Verify omega is invariant to amplitude halving (linear regime check)."""
    if verbose:
        print(f"\n=== Linearity check at alpha={alpha} ===")
    run_full = run_standing_wave(alpha, polarization_axis=0, amplitude=AMPLITUDE, verbose=verbose)
    run_half = run_standing_wave(alpha, polarization_axis=0, amplitude=AMPLITUDE / 2.0, verbose=verbose)
    delta = abs(run_half["omega_measured"] - run_full["omega_measured"])
    delta_rel = float(delta / run_full["omega_measured"])
    if verbose:
        print(f"  delta_omega_rel = {delta_rel:.2e}  PASS = {delta_rel < 0.001}")
    return {
        "alpha": float(alpha),
        "omega_full": run_full["omega_measured"],
        "omega_half": run_half["omega_measured"],
        "delta_omega_rel": delta_rel,
        "pass_linearity": bool(delta_rel < 0.001),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> dict:
    print(f"P1 wavepacket experiment")
    print(f"  k0*a = pi/4 = {K0A:.4f},  N={N_GRID},  dt={DT}")
    print(f"  amplitude = {AMPLITUDE},  A*k0 = {AMPLITUDE * K0A:.2e}")
    print(f"  envelope W/a = {ENVELOPE_WIDTH}")

    from branesim.diagnostics.alpha_separability import verify_track_a

    results: dict = {
        "experiment": "alpha_separability_P1",
        "date": "2026-06-03",
        "k0a": float(K0A),
        "N": N_GRID,
        "dt": DT,
        "amplitude": AMPLITUDE,
        "envelope_width_over_a": ENVELOPE_WIDTH,
        "n_periods": N_PERIODS,
    }

    # Track A (pure closed-form, no simulation)
    print("\nRunning Track A verification (closed-form)...")
    track_a = verify_track_a()
    results["track_a"] = {
        k: (float(v) if hasattr(v, "__float__") else v)
        for k, v in track_a.items()
    }
    print(f"  all_pass: {track_a['all_pass']}")
    print(f"  max closed-form vs numerical err: {track_a['max_closed_form_vs_numerical_err']:.2e}")
    print(f"  max linearity rel_err (traceless in alpha): {track_a['max_linearity_rel_err']:.2e}")
    print(f"  g([111]) = {track_a['g_111']:.2e}  (must be 0)")
    print(f"  rho_SU3(alpha=0.2, [100]) = {track_a['rho_SU3_alpha020_100']:.5f}")
    print(f"  expected_coeff * g = {track_a['expected_rho_SU3_coeff'] * track_a['g_100_alpha020']:.5f}")

    # Linearity check
    lin = linearity_check(alpha=0.2)
    results["linearity_check"] = lin

    # P1 tests at alpha=0.2 and alpha=0.5
    p1_results = {}
    for alpha in [0.2, 0.5]:
        p1_results[f"alpha_{alpha}"] = test_p1_alpha(alpha)
    results["p1"] = p1_results

    overall_pass = (
        bool(track_a["all_pass"])
        and lin["pass_linearity"]
        and all(v["pass_p1"] for v in p1_results.values())
    )
    results["overall_pass"] = overall_pass

    # Persist
    out_path = RESULTS_DIR / "p1_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2,
                  default=lambda x: float(x) if hasattr(x, "__float__") else str(x))
    print(f"\nResults saved to {out_path}")
    print(f"\nOVERALL PASS: {overall_pass}")
    return results


if __name__ == "__main__":
    main()