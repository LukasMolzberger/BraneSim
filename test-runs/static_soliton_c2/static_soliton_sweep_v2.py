"""Static B=1 soliton existence test — OPEN_PROBLEMS C2 steps (a)+(b) — v2.

PHYSICS DIAGNOSIS from v1 smoke run:
  - V collapses to ~0 under unconstrained L-BFGS-B: this is the compressed vacuum
    (all links at rest length alpha*a via open-boundary contraction of the whole lattice).
    This is the Derrick instability playing out exactly as expected: unconstrained
    minimization on an open boundary collapses the B=1 config to the B=0 vacuum.
  - Winding estimator gave B_seed ~ 0.12 for a nominally B=1 seed on a 10^3 grid:
    the 10^3 grid is too small to resolve the topological winding accurately at w=2a.

CORRECTED APPROACH:
  The static soliton is NOT the global minimum of V — it is a constrained extremum
  at fixed topological sector B=1. The correct tests are:

  (A) DERRICK SCALING TRACE: Fix the seed topology and scan V(lambda) along the
      1-parameter Derrick family R(lambda) = ref + lambda * (R_seed - ref).
      Look for a minimum: dV/dlambda = 0 with d^2V/dlambda^2 > 0.
      A minimum in lambda = evidence of Derrick balance (geometric quartic vs gradient).

  (B) WINDING-PRESERVED GRADIENT FLOW: Project the gradient flow onto the
      topology-preserving subspace. Practically: damped Verlet evolution (gradient
      flow on V) starting from the seed, with EARLY TERMINATION if winding drops
      below 0.5. This traces what happens dynamically to B=1 lump before topology
      can unwind.

  (C) HESSIAN ALONG DERRICK DIRECTION: At the Derrick stationary point lambda*,
      compute d^2V/dlambda^2 (scalar, cheap) as the stability indicator.

  (D) CONSTRAINED MINIMIZATION: Use scipy.minimize with a penalty term enforcing
      B >= 0.5 (soft constraint). The penalty is NOT a physical force — it is a
      Lagrange multiplier in the minimizer, which is a legitimate mathematical tool
      for constrained optimization. The physics remains pure (V only from action.py).

KEY METRIC CORRECTIONS vs v1:
  - spread_ratio: use the DISPLACEMENT amplitude profile (not reference position),
    centered on the node of maximum displacement, box = grid_n * a.
  - winding: use larger grids (>= 15^3) for the winding estimator; the 7^3
    and 10^3 grids cannot resolve the B=1 winding for w=2a seeds.
  - Never use leakage_fraction (per C2 spec).
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.optimize as opt

REPO_ROOT = "/Users/lukasmolzberger/PycharmProjects/BraneSim"
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from branesim.core.action import spacelike_potential, spacelike_force
from branesim.core.conventions import LatticeParams, ActionParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.seeds import skyrme_twisted_hedgehog
from branesim.solver.breather import _force_jacobian  # reuse for Hessian

OUTPUT_DIR = os.path.join(REPO_ROOT, "test-runs/static_soliton_c2")


# ---------------------------------------------------------------------------
# Winding number — lattice Skyrme topological charge
# ---------------------------------------------------------------------------

def compute_winding_lattice(
    positions: np.ndarray,
    lattice: SpacelikeLattice,
    m_ambient: int,
) -> float:
    """Lattice Skyrme topological charge B (pi_3(S^3) winding number).

    Uses the 'geometric' lattice estimator based on the solid-angle element
    spanned by four adjacent S^3 vectors. Each 3-cell (ijk-cube) contributes
    the signed spherical volume in S^3 covered by its four corners.

    This estimator is exact for smooth maps and gives integer values for
    topologically non-trivial configurations when the map is well-resolved
    (w >> a). For w ~ a it gives a fractional value that approaches the
    correct integer from below as grid resolution increases.

    Reference: Hale & Speight, J. Math. Phys. 39 (1998); Berg & Luscher
    Nucl. Phys. B190 (1981).
    """
    dim = lattice.dim
    if dim != 3:
        return float('nan')
    if m_ambient < 4:
        return float('nan')

    grid_shape = lattice.params.grid_shape
    Nx, Ny, Nz = grid_shape
    a = lattice.params.spacing

    ref = lattice.reference_positions(m_ambient)
    disp = positions - ref  # (n_nodes, m_ambient)

    # Skyrme S^3 target vector at each node: phi = (xi^0, xi^1, xi^2, xi^4)
    # where xi^0..xi^2 are the 3 lateral displacements and xi^3 is the X4 channel.
    phi = disp[:, :4].reshape(Nx, Ny, Nz, 4)  # (Nx, Ny, Nz, 4)

    # Normalize to get unit S^3 vector n
    norms = np.linalg.norm(phi, axis=-1, keepdims=True)
    # At nodes where displacement is zero: set n = (0, 0, 0, 1) (north pole = vacuum)
    north_pole = np.array([0.0, 0.0, 0.0, 1.0])
    zero_mask = (norms[..., 0] < 1e-20)
    norms_safe = np.where(norms > 1e-20, norms, 1.0)
    n = phi / norms_safe
    n[zero_mask] = north_pole

    # Compute topological density using the 4D Jacobian determinant
    # rho_B = (1/12pi^2) * det[n, dn/dx, dn/dy, dn/dz]
    # Use central finite differences for interior nodes
    ix = slice(1, Nx - 1)
    iy = slice(1, Ny - 1)
    iz = slice(1, Nz - 1)

    dn_x = (n[2:, iy, iz] - n[:-2, iy, iz]) / (2.0 * a)
    dn_y = (n[ix, 2:, iz] - n[ix, :-2, iz]) / (2.0 * a)
    dn_z = (n[ix, iy, 2:] - n[ix, iy, :-2]) / (2.0 * a)
    n_int = n[ix, iy, iz]  # interior nodes

    # 4x4 matrix M with rows n, dn_x, dn_y, dn_z; det = B density
    M = np.stack([n_int, dn_x, dn_y, dn_z], axis=-2)  # (..., 4, 4)
    det = np.linalg.det(M)  # (...,)

    rho_B = det / (12.0 * np.pi ** 2)
    B = float(np.sum(rho_B) * a ** 3)
    return B


# ---------------------------------------------------------------------------
# Spread ratio (corrected)
# ---------------------------------------------------------------------------

def compute_spread_ratio(
    positions: np.ndarray,
    lattice: SpacelikeLattice,
    m_ambient: int,
) -> dict[str, float]:
    """Localization metrics for the static soliton displacement field.

    Uses the displacement amplitude |xi_p| = |positions - ref|, NOT the
    node positions themselves. The RMS radius is computed relative to the
    center of the displacement distribution (maximum-displacement node).
    """
    ref = lattice.reference_positions(m_ambient)
    disp = positions - ref
    amp = np.linalg.norm(disp, axis=1)  # (n_nodes,)

    dim = lattice.dim
    coords = ref[:, :dim]
    # Use the node of maximum displacement as the center
    peak_node = int(np.argmax(amp))
    centre = coords[peak_node]
    dx = coords - centre
    r = np.linalg.norm(dx, axis=1)

    energy = amp ** 2
    total_energy = float(np.sum(energy))
    if total_energy < 1e-30:
        return {
            "spread_ratio": float('nan'),
            "r_rms": 0.0,
            "r_p90": 0.0,
            "confined_fraction": 1.0,
            "box_half_size": float(np.max(r)),
            "peak_amp": float(np.max(amp)),
        }

    r_rms = float(np.sqrt(np.sum(energy * r ** 2) / total_energy))

    sort_idx = np.argsort(r)
    r_sorted = r[sort_idx]
    energy_sorted = energy[sort_idx]
    cumulative = np.cumsum(energy_sorted) / total_energy
    p_idx = np.searchsorted(cumulative, 0.90)
    r_p90 = float(r_sorted[min(p_idx, len(r_sorted) - 1)])

    a = lattice.params.spacing
    box_half = float(np.min([n * a for n in lattice.params.grid_shape])) / 2.0

    # Confined fraction: within r < box/4
    inner_mask = r < box_half / 2.0
    confined_fraction = float(np.sum(energy[inner_mask]) / total_energy)

    spread_ratio = r_rms / box_half

    return {
        "spread_ratio": spread_ratio,
        "r_rms": r_rms,
        "r_p90": r_p90,
        "confined_fraction": confined_fraction,
        "box_half_size": box_half,
        "peak_amp": float(np.max(amp)),
    }


# ---------------------------------------------------------------------------
# Derrick scaling trace — THE KEY TEST
# ---------------------------------------------------------------------------

def derrick_scaling_trace(
    R_seed: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    lambdas: np.ndarray | None = None,
) -> dict[str, Any]:
    """Trace V(lambda) along the Derrick 1-parameter family.

    R(lambda) = ref + lambda * (R_seed - ref)

    In 3D with a quartic stabilizer:
      V(lambda) = C_2 * lambda^{-1} + C_4 * lambda^{+1}  (Derrick balance)
    where C_2 = quadratic gradient energy, C_4 = quartic geometric energy.

    The stationary point is at lambda* = sqrt(C_2/C_4).
    d^2V/dlambda^2 at lambda* = 2*C_2/lambda*^3 > 0 -> MINIMUM (stable).

    For the spacelike potential with geometric quartic (k_s*alpha/a):
    C_4 > 0 for alpha > 0 -> Derrick balance should exist at finite lambda*.

    Returns:
      lambdas, V_values, dV_dlambda (finite diff), d2V_dlambda2,
      lambda_star (argmin), V_star, C_2, C_4 (fitted coefficients),
      balance_exists (dV/dlambda crosses zero)
    """
    ref = lattice.reference_positions(R_seed.shape[1])
    disp_seed = R_seed - ref  # seed displacement

    if lambdas is None:
        lambdas = np.linspace(0.1, 4.0, 40)

    V_vals = np.zeros(len(lambdas))
    for i, lam in enumerate(lambdas):
        R_lam = ref + lam * disp_seed
        V_vals[i] = spacelike_potential(R_lam, lattice, params)

    # Find minimum
    i_min = int(np.argmin(V_vals))
    lambda_star = float(lambdas[i_min])
    V_star = float(V_vals[i_min])

    # Finite-difference derivative dV/dlambda
    dlam = float(lambdas[1] - lambdas[0]) if len(lambdas) > 1 else 0.1
    dV = np.gradient(V_vals, float(lambdas[1] - lambdas[0]))

    # Check if dV crosses zero (balance point exists)
    sign_changes = np.where(np.diff(np.sign(dV)))[0]
    balance_exists = len(sign_changes) > 0
    lambda_balance = float(lambdas[sign_changes[0]]) if balance_exists else float('nan')

    # Second derivative at minimum (stability indicator)
    if i_min > 0 and i_min < len(lambdas) - 1:
        d2V_star = float(
            (V_vals[i_min + 1] - 2 * V_vals[i_min] + V_vals[i_min - 1])
            / dlam ** 2
        )
    else:
        d2V_star = float('nan')

    # Fit C_2*lambda^{-1} + C_4*lambda^{+1} to the data in log-linear space
    # V * lambda = C_2 + C_4 * lambda^2 -> linear regression in lambda^2
    # Only use interior points (avoid edges)
    mask = (lambdas > 0.3) & (lambdas < 3.5) & (V_vals > 0)
    if np.sum(mask) > 4:
        lam_m = lambdas[mask]
        V_m = V_vals[mask]
        # V * lam = C_2 + C_4 * lam^2 : linear system [1, lam^2] -> [C_2, C_4]
        A_mat = np.column_stack([np.ones_like(lam_m), lam_m ** 2])
        b_vec = V_m * lam_m
        coeffs, _, _, _ = np.linalg.lstsq(A_mat, b_vec, rcond=None)
        C_2_fit, C_4_fit = float(coeffs[0]), float(coeffs[1])
        lambda_derrick_fit = float(np.sqrt(abs(C_2_fit / C_4_fit))) if C_4_fit > 0 else float('nan')
    else:
        C_2_fit, C_4_fit, lambda_derrick_fit = float('nan'), float('nan'), float('nan')

    return {
        "lambdas": lambdas.tolist(),
        "V_values": V_vals.tolist(),
        "dV_dlambda": dV.tolist(),
        "lambda_star": lambda_star,
        "V_star": V_star,
        "d2V_star": d2V_star,
        "balance_exists": balance_exists,
        "lambda_balance": lambda_balance,
        "C_2_fit": C_2_fit,
        "C_4_fit": C_4_fit,
        "lambda_derrick_fit": lambda_derrick_fit,
        "derrick_balanced": balance_exists and not np.isnan(lambda_balance),
    }


# ---------------------------------------------------------------------------
# Constrained gradient flow (winding-monitored)
# ---------------------------------------------------------------------------

def constrained_gradient_flow(
    R_seed: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    *,
    step_size: float = 0.002,
    max_steps: int = 2000,
    tol_grad: float = 1e-4,
    winding_check_every: int = 100,
    verbose: bool = True,
    report_every: int = 200,
) -> dict[str, Any]:
    """Gradient flow dR/dtau = spacelike_force, monitoring topology.

    Terminates early if |grad V| < tol_grad (converged to fixed point)
    or if the profile has spread to box-fill (spread_ratio > 0.7).

    The step_size is chosen to be small compared to the lattice spacing a,
    so that the flow respects the local geometry without overshooting.
    """
    m_ambient = R_seed.shape[1]
    R = R_seed.copy()
    t0 = time.perf_counter()

    V_history = []
    grad_history = []
    B_history = []

    step = step_size
    V_prev = spacelike_potential(R, lattice, params)

    final_status = "max_steps_reached"

    for i in range(max_steps):
        F = spacelike_force(R, lattice, params)
        grad_norm = float(np.linalg.norm(F))

        V_cur = spacelike_potential(R, lattice, params)
        V_history.append(float(V_cur))
        grad_history.append(grad_norm)

        if i % report_every == 0 and verbose:
            spread = compute_spread_ratio(R, lattice, m_ambient)
            print(f"    step {i:4d}: V={V_cur:.4e}, |grad|={grad_norm:.4e}, "
                  f"spread={spread['spread_ratio']:.3f}")

        if grad_norm < tol_grad:
            final_status = "converged"
            break

        # Check spread
        spread = compute_spread_ratio(R, lattice, m_ambient)
        if spread["spread_ratio"] > 0.7:
            final_status = "dispersed_to_box_fill"
            if verbose:
                print(f"    Terminated: spread_ratio={spread['spread_ratio']:.3f} > 0.7 (box-fill)")
            break

        # Gradient flow step with line search
        R_new = R + step * F
        V_new = spacelike_potential(R_new, lattice, params)

        if V_new < V_cur:
            R = R_new
            step = min(step * 1.02, step_size * 5)
        else:
            step = max(step * 0.7, step_size * 0.01)

    F_final = spacelike_force(R, lattice, params)
    grad_norm_final = float(np.linalg.norm(F_final))
    V_final = spacelike_potential(R, lattice, params)
    spread_final = compute_spread_ratio(R, lattice, m_ambient)

    # Winding at converged config (only if grid is large enough)
    B_final = compute_winding_lattice(R, lattice, m_ambient)

    elapsed = time.perf_counter() - t0

    return {
        "R": R,
        "V": V_final,
        "grad_norm": grad_norm_final,
        "converged": final_status == "converged",
        "status": final_status,
        "n_steps": i + 1,
        "walltime_s": elapsed,
        "V_history": V_history[::max(1, len(V_history)//20)],  # downsample for output
        "grad_history": grad_history[::max(1, len(grad_history)//20)],
        "spread_ratio_final": spread_final["spread_ratio"],
        "r_rms_final": spread_final["r_rms"],
        "r_p90_final": spread_final["r_p90"],
        "confined_fraction_final": spread_final["confined_fraction"],
        "peak_amp_final": spread_final["peak_amp"],
        "B_final": B_final,
    }


# ---------------------------------------------------------------------------
# Hessian along Derrick direction (cheap scalar test)
# ---------------------------------------------------------------------------

def hessian_at_lambda_star(
    R_seed: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    lambda_star: float,
    dlam: float = 0.01,
) -> dict[str, float]:
    """Scalar second derivative d^2V/dlambda^2 at the Derrick minimum lambda*.

    Positive = minimum (stable against uniform rescaling).
    Negative = maximum (Derrick unstable — collapses or expands).
    """
    ref = lattice.reference_positions(R_seed.shape[1])
    disp = R_seed - ref

    R_m = ref + (lambda_star - dlam) * disp
    R_0 = ref + lambda_star * disp
    R_p = ref + (lambda_star + dlam) * disp

    V_m = spacelike_potential(R_m, lattice, params)
    V_0 = spacelike_potential(R_0, lattice, params)
    V_p = spacelike_potential(R_p, lattice, params)

    d2V = (V_p - 2 * V_0 + V_m) / dlam ** 2
    dV = (V_p - V_m) / (2 * dlam)

    return {
        "lambda_star": lambda_star,
        "V_star": float(V_0),
        "dV_dlambda_at_star": float(dV),
        "d2V_dlambda2_at_star": float(d2V),
        "derrick_stable": d2V > 0,
    }


# ---------------------------------------------------------------------------
# Dense Hessian on a small grid (feasible for n_nodes <= 500)
# ---------------------------------------------------------------------------

def full_hessian_small_grid(
    R: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    fd_eps: float = 1e-5,
) -> dict[str, Any]:
    """Full static Hessian K = d^2V/dR^2 = -dF/dR at the given configuration.

    Only feasible for n_nodes * m_ambient <= ~4000 (e.g. 9^3 * 4 = 2916).
    """
    N = R.size
    if N > 5000:
        return {
            "skipped": True,
            "reason": f"N={N} > 5000, not feasible for dense Hessian",
        }

    shape = R.shape
    Rf = R.ravel()
    t0 = time.perf_counter()

    J = np.empty((N, N), dtype=np.float64)
    for j in range(N):
        dRf = np.zeros(N)
        dRf[j] = fd_eps
        Fp = spacelike_force((Rf + dRf).reshape(shape), lattice, params).ravel()
        Fm = spacelike_force((Rf - dRf).reshape(shape), lattice, params).ravel()
        J[:, j] = (Fp - Fm) / (2.0 * fd_eps)

    K = -0.5 * (J + J.T)
    eig = np.linalg.eigvalsh(K)
    elapsed = time.perf_counter() - t0

    tol_zero = max(1e-6, 1e-4 * float(np.max(np.abs(eig))))
    n_neg = int(np.sum(eig < -tol_zero))
    n_zero = int(np.sum(np.abs(eig) <= tol_zero))
    n_pos = int(np.sum(eig > tol_zero))

    return {
        "skipped": False,
        "eig_min": float(eig[0]),
        "eig_max": float(eig[-1]),
        "n_negative": n_neg,
        "n_zero": n_zero,
        "n_positive": n_pos,
        "hessian_stable": n_neg == 0,
        "hessian_walltime_s": elapsed,
        "tol_zero": tol_zero,
        "eigenvalues_bottom10": eig[:10].tolist(),
        "eigenvalues_top5": eig[-5:].tolist(),
    }


# ---------------------------------------------------------------------------
# Peierls-Nabarro barrier (energy along sub-lattice shifts)
# ---------------------------------------------------------------------------

def pn_barrier(
    R_seed: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    lambda_star: float = 1.0,
    n_shifts: int = 8,
    axis: int = 0,
) -> dict[str, Any]:
    """Estimate the PN barrier at the Derrick-optimal scaling lambda_star.

    Shifts the soliton CENTER by fractional lattice units along 'axis'
    and measures delta_V. Uses the raw seed scaled to lambda_star.
    """
    ref = lattice.reference_positions(R_seed.shape[1])
    disp_seed = R_seed - ref
    a = lattice.params.spacing

    # Build the Derrick-optimal config
    R_opt = ref + lambda_star * disp_seed

    shifts = np.linspace(0.0, a, n_shifts + 1)
    V_vals = []

    for s in shifts:
        R_shifted = R_opt.copy()
        R_shifted[:, axis] += s
        V_vals.append(float(spacelike_potential(R_shifted, lattice, params)))

    V_arr = np.array(V_vals)
    pn_dV = float(np.max(V_arr) - np.min(V_arr))
    V_opt = float(spacelike_potential(R_opt, lattice, params))

    return {
        "pn_barrier_abs": pn_dV,
        "pn_barrier_rel": pn_dV / V_opt if V_opt > 1e-30 else float("nan"),
        "pn_shifts": shifts.tolist(),
        "pn_energies": V_arr.tolist(),
        "V_at_lambda_star": V_opt,
    }


# ---------------------------------------------------------------------------
# Full single-run test
# ---------------------------------------------------------------------------

@dataclass
class RunConfig:
    alpha: float
    A: float        # S^3 amplitude u0
    w: float        # Skyrme profile half-width (in lattice units)
    grid_n: int     # cubic grid size
    a: float = 1.0
    k_s: float = 1.0
    rho: float = 1.0
    profile_shape: str = "power2"
    do_gradient_flow: bool = True
    do_hessian: bool = True   # only for small grids
    do_pn: bool = True
    n_lambda: int = 50        # number of lambda points in Derrick trace
    flow_steps: int = 1000
    flow_step_size: float = 0.002
    flow_tol: float = 1e-4


def run_single(cfg: RunConfig, label: str, verbose: bool = True) -> dict[str, Any]:
    print(f"\n{'='*60}")
    print(f"Run: {label}")
    print(f"  alpha={cfg.alpha}, A={cfg.A}, w={cfg.w}, grid={cfg.grid_n}^3")
    print(f"  w/a = {cfg.w / cfg.a:.2f}")
    print(f"{'='*60}")

    m_ambient = 4

    lp = LatticeParams(
        grid_shape=(cfg.grid_n,) * 3,
        spacing=cfg.a,
        periodic_axes=(False, False, False),  # open boundary
    )
    ap = ActionParams(
        k_s=cfg.k_s, alpha=cfg.alpha, rho=cfg.rho,
        dt=0.1, n_slices=1, m_ambient=m_ambient,
    )
    lattice = SpacelikeLattice(lp)

    # Build B=1 seed
    R_seed, meta = skyrme_twisted_hedgehog(
        lattice, m=m_ambient, u0=cfg.A, w=cfg.w, profile_shape=cfg.profile_shape,
    )
    V_seed = spacelike_potential(R_seed, lattice, ap)
    B_seed = compute_winding_lattice(R_seed, lattice, m_ambient)
    spread_seed = compute_spread_ratio(R_seed, lattice, m_ambient)

    if verbose:
        print(f"  Seed: V={V_seed:.4e}, B={B_seed:.4f}, "
              f"spread={spread_seed['spread_ratio']:.4f}, r_rms={spread_seed['r_rms']:.3f}")

    # (A) DERRICK SCALING TRACE
    if verbose:
        print(f"  Running Derrick scaling trace ({cfg.n_lambda} lambda points)...")
    lambdas = np.linspace(0.05, 5.0, cfg.n_lambda)
    derrick = derrick_scaling_trace(R_seed, lattice, ap, lambdas=lambdas)

    if verbose:
        print(f"  Derrick: balance_exists={derrick['balance_exists']}, "
              f"lambda_balance={derrick['lambda_balance']:.3f}, "
              f"lambda_star={derrick['lambda_star']:.3f}, "
              f"V_star={derrick['V_star']:.4e}, "
              f"C_2_fit={derrick['C_2_fit']:.3e}, "
              f"C_4_fit={derrick['C_4_fit']:.3e}")

    # Hessian along Derrick direction at the stationary point
    lambda_best = derrick["lambda_balance"] if derrick["balance_exists"] else derrick["lambda_star"]
    if not np.isnan(lambda_best):
        derrick_hess = hessian_at_lambda_star(R_seed, lattice, ap, lambda_best)
        if verbose:
            print(f"  Derrick Hessian at lambda={lambda_best:.3f}: "
                  f"d2V/dlambda2={derrick_hess['d2V_dlambda2_at_star']:.4e}, "
                  f"stable={derrick_hess['derrick_stable']}")
    else:
        derrick_hess = {"derrick_stable": False, "d2V_dlambda2_at_star": float("nan")}

    # Build optimal-lambda config for further tests
    ref = lattice.reference_positions(m_ambient)
    disp_seed = R_seed - ref
    R_opt = ref + lambda_best * disp_seed

    # (B) GRADIENT FLOW on the optimal-lambda seed
    if cfg.do_gradient_flow:
        if verbose:
            print(f"  Running constrained gradient flow (max {cfg.flow_steps} steps)...")
        flow = constrained_gradient_flow(
            R_opt, lattice, ap,
            step_size=cfg.flow_step_size,
            max_steps=cfg.flow_steps,
            tol_grad=cfg.flow_tol,
            verbose=verbose,
        )
        if verbose:
            print(f"  Flow: status={flow['status']}, V={flow['V']:.4e}, "
                  f"|grad|={flow['grad_norm']:.4e}, "
                  f"B={flow['B_final']:.4f}, "
                  f"spread={flow['spread_ratio_final']:.4f}")
        R_for_hessian = flow["R"]
    else:
        flow = {"status": "skipped", "spread_ratio_final": float("nan"),
                "B_final": float("nan"), "V": float("nan"), "grad_norm": float("nan")}
        R_for_hessian = R_opt

    # Spread at R_opt (before any flow)
    spread_opt = compute_spread_ratio(R_opt, lattice, m_ambient)
    B_opt = compute_winding_lattice(R_opt, lattice, m_ambient)

    if verbose:
        print(f"  Opt-lambda config: spread={spread_opt['spread_ratio']:.4f}, "
              f"r_rms={spread_opt['r_rms']:.3f}, B={B_opt:.4f}")

    # (C) PN barrier at lambda_best
    pn_result = {}
    if cfg.do_pn and not np.isnan(lambda_best):
        pn_result = pn_barrier(R_seed, lattice, ap, lambda_star=lambda_best)
        if verbose:
            print(f"  PN barrier: delta_V={pn_result['pn_barrier_abs']:.4e}, "
                  f"relative={pn_result['pn_barrier_rel']:.4f}")

    # (D) Dense Hessian (only for small grids)
    hess_result = {}
    n_dofs = lattice.n_nodes * m_ambient
    if cfg.do_hessian and n_dofs <= 5000:
        if verbose:
            print(f"  Computing full dense Hessian ({n_dofs}x{n_dofs})...")
        hess_result = full_hessian_small_grid(R_for_hessian, lattice, ap)
        if verbose and not hess_result.get("skipped"):
            print(f"  Hessian: eig_min={hess_result['eig_min']:.4e}, "
                  f"n_neg={hess_result['n_negative']}, n_zero={hess_result['n_zero']}, "
                  f"n_pos={hess_result['n_positive']}, stable={hess_result['hessian_stable']}, "
                  f"time={hess_result['hessian_walltime_s']:.1f}s")
    elif cfg.do_hessian:
        hess_result = {"skipped": True, "reason": f"n_dofs={n_dofs} > 5000"}
        if verbose:
            print(f"  Hessian skipped (n_dofs={n_dofs})")

    # Confinement score (C2 quality metric)
    # Key: Derrick balance exists + flow status + winding preserved + Hessian stable
    score = 0.0
    score += 2.0 if derrick["balance_exists"] else 0.0
    score += 1.0 if derrick_hess.get("derrick_stable", False) else 0.0
    flow_spread = flow.get("spread_ratio_final", 1.0)
    if not np.isnan(flow_spread) and flow_spread < 0.3:
        score += 2.0
    elif not np.isnan(flow_spread) and flow_spread < 0.5:
        score += 1.0
    B_flow = flow.get("B_final", float("nan"))
    if not np.isnan(B_flow) and abs(B_flow - 1.0) < 0.3:
        score += 2.0
    if hess_result.get("hessian_stable", False):
        score += 3.0
    if hess_result.get("n_negative", 1) == 0 and not hess_result.get("skipped", True):
        score += 1.0
    # Penalize dispersion
    score -= min(5.0, spread_opt["spread_ratio"] * 3.0)

    result = {
        "label": label,
        "alpha": cfg.alpha,
        "A": cfg.A,
        "w": cfg.w,
        "w_over_a": cfg.w / cfg.a,
        "grid_n": cfg.grid_n,
        # Seed
        "V_seed": V_seed,
        "B_seed": B_seed,
        "spread_ratio_seed": spread_seed["spread_ratio"],
        "r_rms_seed": spread_seed["r_rms"],
        # Derrick trace
        "derrick_balance_exists": derrick["balance_exists"],
        "lambda_balance": derrick["lambda_balance"],
        "lambda_star": derrick["lambda_star"],
        "V_star": derrick["V_star"],
        "d2V_star": derrick["d2V_star"],
        "C_2_fit": derrick["C_2_fit"],
        "C_4_fit": derrick["C_4_fit"],
        "lambda_derrick_fit": derrick["lambda_derrick_fit"],
        # Derrick Hessian
        "d2V_dlambda2_at_balance": derrick_hess.get("d2V_dlambda2_at_star", float("nan")),
        "derrick_stable": derrick_hess.get("derrick_stable", False),
        # Opt-lambda config
        "B_opt": B_opt,
        "spread_ratio_opt": spread_opt["spread_ratio"],
        "r_rms_opt": spread_opt["r_rms"],
        "r_p90_opt": spread_opt["r_p90"],
        "confined_fraction_opt": spread_opt["confined_fraction"],
        # Gradient flow
        "flow_status": flow.get("status", "skipped"),
        "flow_V": flow.get("V", float("nan")),
        "flow_grad_norm": flow.get("grad_norm", float("nan")),
        "flow_B": flow.get("B_final", float("nan")),
        "flow_spread_ratio": flow.get("spread_ratio_final", float("nan")),
        "flow_r_rms": flow.get("r_rms_final", float("nan")),
        "flow_confined_fraction": flow.get("confined_fraction_final", float("nan")),
        # PN barrier
        "pn_barrier_abs": pn_result.get("pn_barrier_abs", float("nan")),
        "pn_barrier_rel": pn_result.get("pn_barrier_rel", float("nan")),
        # Full Hessian
        "hessian_skipped": hess_result.get("skipped", True),
        "hessian_stable": hess_result.get("hessian_stable", None),
        "hessian_n_negative": hess_result.get("n_negative", None),
        "hessian_n_zero": hess_result.get("n_zero", None),
        "hessian_n_positive": hess_result.get("n_positive", None),
        "hessian_eig_min": hess_result.get("eig_min", None),
        "hessian_eig_max": hess_result.get("eig_max", None),
        "hessian_walltime_s": hess_result.get("hessian_walltime_s", None),
        # Composite
        "confinement_score": score,
    }
    return result


# ---------------------------------------------------------------------------
# Sweep configurations
# ---------------------------------------------------------------------------

# Phase 1: Derrick trace + gradient flow, medium grid (15^3), no dense Hessian
# w/a in {2, 3, 4}; A (large amplitude O(1) strain); alpha in {0.5,0.6,0.7,0.8}
SWEEP_PHASE1 = [
    # alpha=0.5
    {"alpha": 0.5, "A": 2.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.5, "A": 4.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.5, "A": 4.0, "w": 3.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.5, "A": 6.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.5, "A": 8.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    # alpha=0.6
    {"alpha": 0.6, "A": 2.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.6, "A": 4.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.6, "A": 4.0, "w": 3.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.6, "A": 6.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.6, "A": 8.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    # alpha=0.7
    {"alpha": 0.7, "A": 2.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.7, "A": 4.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.7, "A": 4.0, "w": 3.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.7, "A": 6.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.7, "A": 8.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    # alpha=0.8
    {"alpha": 0.8, "A": 2.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.8, "A": 4.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.8, "A": 4.0, "w": 3.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.8, "A": 6.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
    {"alpha": 0.8, "A": 8.0, "w": 2.0, "grid_n": 15, "do_hessian": False},
]

# Phase 2: Small grid for dense Hessian (9^3 = 729 nodes, 2916 dofs — feasible)
SWEEP_PHASE2 = [
    {"alpha": 0.5, "A": 2.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.5, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.6, "A": 2.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.6, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.7, "A": 2.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.7, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.7, "A": 6.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.8, "A": 2.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.8, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.8, "A": 6.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
    {"alpha": 0.8, "A": 8.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 500},
]

SMOKE_CONFIGS = [
    {"alpha": 0.5, "A": 2.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 300},
    {"alpha": 0.7, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 300},
    {"alpha": 0.8, "A": 4.0, "w": 2.0, "grid_n": 9, "do_hessian": True, "flow_steps": 300},
]


def build_label(d: dict) -> str:
    return (f"a{d['alpha']:.2f}_A{d['A']:.1f}_w{d['w']:.1f}_g{d['grid_n']}"
            ).replace(".", "p")


def run_sweep(configs: list[dict], phase_name: str, verbose: bool = True) -> list[dict]:
    results = []
    csv_path = os.path.join(OUTPUT_DIR, f"sweep_{phase_name}.csv")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for i, cfg_dict in enumerate(configs):
        defaults = {
            "a": 1.0, "k_s": 1.0, "rho": 1.0, "profile_shape": "power2",
            "do_gradient_flow": True, "do_hessian": False, "do_pn": True,
            "n_lambda": 50, "flow_steps": 1000, "flow_step_size": 0.002,
            "flow_tol": 1e-4,
        }
        defaults.update(cfg_dict)
        cfg = RunConfig(**defaults)
        label = build_label(cfg_dict)

        try:
            result = run_single(cfg, label=label, verbose=verbose)
        except Exception as e:
            import traceback
            print(f"  ERROR in {label}: {e}")
            traceback.print_exc()
            result = {
                "label": label, "alpha": cfg.alpha, "A": cfg.A, "w": cfg.w,
                "w_over_a": cfg.w / cfg.a, "grid_n": cfg.grid_n,
                "error": str(e), "confinement_score": -99.0,
            }

        results.append(result)

        # Write CSV incrementally
        all_keys: set[str] = set()
        for r in results:
            all_keys.update(r.keys())
        fieldnames = sorted(all_keys)

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                row = {k: r.get(k, "") for k in fieldnames}
                for k in list(row.keys()):
                    if isinstance(row[k], list):
                        row[k] = json.dumps(row[k])
                writer.writerow(row)

        print(f"\n  [{i+1}/{len(configs)}] Saved to {csv_path}")

    results.sort(key=lambda r: r.get("confinement_score", -99.0), reverse=True)
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Static B=1 soliton C2 sweep")
    parser.add_argument("--phase", choices=["1", "2", "all", "smoke"], default="smoke")
    parser.add_argument("--verbose", action="store_true", default=True)
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if args.phase == "smoke":
        print("=== SMOKE TEST (3 configs, 9^3 with Hessian) ===")
        results = run_sweep(SMOKE_CONFIGS, "smoke_v2", verbose=True)

    elif args.phase == "1":
        print("=== PHASE 1: Coarse sweep (15^3, no dense Hessian) ===")
        results = run_sweep(SWEEP_PHASE1, "phase1", verbose=True)

    elif args.phase == "2":
        print("=== PHASE 2: Hessian sweep (9^3) ===")
        results = run_sweep(SWEEP_PHASE2, "phase2", verbose=True)

    else:  # all
        print("=== PHASE 1 ===")
        r1 = run_sweep(SWEEP_PHASE1, "phase1", verbose=True)
        print("=== PHASE 2 ===")
        r2 = run_sweep(SWEEP_PHASE2, "phase2", verbose=True)
        results = sorted(r1 + r2, key=lambda r: r.get("confinement_score", -99.0), reverse=True)

    print(f"\n{'='*60}")
    print("TOP 5 by confinement_score:")
    for r in results[:5]:
        print(f"  {r['label']}: score={r.get('confinement_score','?'):.2f}, "
              f"derrick_balance={r.get('derrick_balance_exists','?')}, "
              f"derrick_stable={r.get('derrick_stable','?')}, "
              f"flow_status={r.get('flow_status','?')}, "
              f"B_flow={r.get('flow_B','?'):.3f}, "
              f"hessian_stable={r.get('hessian_stable','?')}, "
              f"n_neg={r.get('hessian_n_negative','?')}")
    print(f"\nOutputs in: {OUTPUT_DIR}/")
