#!/usr/bin/env python3
"""Intermediate-w box-doubling study: decisive self-localization test for B=1 soliton.

CONTEXT
-------
Prior experiments established two negative results:
  C1: w >> a  (w=10a+) → no interior energy minimum; profile flattens.
  C2: w = 2a  → boundary-confined: V*/N_nodes constant to 4%, spread_ratio flat
      ~0.72, radius_rms ∝ grid_n across grid_n 13→20→26.

The ONLY untested window is a << w << box.  This experiment picks w = 5a
(justification: w/a=5 gives 5x lattice spacing separation, well into the
continuum corner; lattice effects fall as (a/w)^2 = 4%), holds alpha=0.7 fixed,
and varies the box: grid_n in {64, 96, 128}.

ISOLATION CHECK (power-2 profile F(r) = pi/(1+(r/w)^2)):
  Need F(r_edge) < 0.05 for the analytic B estimator to be reliable.
  r_edge = grid_n/2 (face center, conservative):
    grid_n=64: r_edge=32, F_edge = pi/(1+(32/5)^2) = pi/(1+40.96) = 0.0743  (marginal)
    grid_n=96: r_edge=48, F_edge = pi/(1+(48/5)^2) = pi/(1+92.16) = 0.0337  (clean)
    grid_n=128:r_edge=64, F_edge = pi/(1+(64/5)^2) = pi/(1+163.84) = 0.0190 (very clean)
  grid_n=64 is included as the fast reference point but is flagged marginal.
  grid_n=96 is the minimum for a clean verdict.  grid_n=128 confirms doubling.

MINIMIZER
---------
Plain gradient descent with backtracking (the C2 harness) did NOT converge even
at grid_n=13 in 8000 steps (grad_norm 0.05-0.09).  At grid_n 64-96 this leaves
the verdict muddy.  This harness implements FIRE (Fast Inertial Relaxation Engine,
Bitzek et al. 2006) for robust accelerated minimization:
  - When power (F · v) > 0: mix velocity toward force direction, grow step.
  - When power < 0: zero velocity, reduce step.
  - Converges 5-30x faster than plain gradient descent on rugged energy surfaces.
  - FIRE is a pure minimizer of V; it does not introduce any physics or damping
    beyond what gradient flow does.  It is equivalent to the overdamped Langevin
    limit with adaptive step — a legal minimization tool under principles.md §1.2.

VERDICT METRIC
--------------
  SELF-LOCALIZED: V* converges to a box-INDEPENDENT constant as box grows;
                  radius_rms flat; spread_ratio decreases as box grows.
  BOUNDARY-CONFINED: V*/N_nodes constant across grid_n; radius_rms ∝ grid_n;
                     spread_ratio flat ~0.7.
  ROUTE DEAD: if boundary-confined persists at w=5a, the substrate's geometric
              quartic cannot bind a self-localized B=1 soliton at any continuum w.

HARD RULES (principles.md enforced)
------------------------------------
- No confinement forces, no nonlinear saturation, no clamps.
- Physics from spacelike_potential / spacelike_force only.
- Dirichlet vacuum BC = legitimate constraint (boundary nodes frozen).
- FIRE = minimization tool only; no dynamical claim.

OUTPUTS (all small, CSV/JSON, no large arrays)
----------------------------------------------
  $BRANESIM_RESULTS_DIR/intermediate_w_sweep.csv
  $BRANESIM_RESULTS_DIR/intermediate_w_g{N}.json  (one per grid_n)

Usage:
  # local fast validation (w=3, grid_n=32, smoke test):
  python orchestration/configs/baryon_candidates/run_intermediate_w_box_doubling.py --validate

  # full production run (w=5, grid_n 64/96/128):
  python orchestration/configs/baryon_candidates/run_intermediate_w_box_doubling.py --run

  # single config (e.g. grid_n=64 only):
  python orchestration/configs/baryon_candidates/run_intermediate_w_box_doubling.py --run --grid-sizes 64

  # AWS: $BRANESIM_RESULTS_DIR is pre-set by ec2_user_data.sh.tmpl
  python orchestration/configs/baryon_candidates/run_intermediate_w_box_doubling.py --run
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Path setup (works both as a script in orchestration/ and when called from
# the project root as a module).
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from branesim.core.action import spacelike_potential, spacelike_force
from branesim.core.conventions import LatticeParams, ActionParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.confinement import confinement_metrics_per_slice
from branesim.initialization.seeds import skyrme_twisted_hedgehog


# ---------------------------------------------------------------------------
# Output directory: honour BRANESIM_RESULTS_DIR if set (AWS), else local.
# ---------------------------------------------------------------------------
def _output_dir() -> str:
    env = os.environ.get("BRANESIM_RESULTS_DIR", "").strip()
    if env:
        os.makedirs(env, exist_ok=True)
        return env
    # Local fallback: test-runs directory alongside the repo.
    d = os.path.join(_REPO_ROOT, "test-runs", "intermediate_w_study")
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Analytic B estimator (AUTHORITATIVE — not the Jacobian one)
# B = -(2/pi) * [antideriv(F_edge) - antideriv(F_inner)]
# F_inner = pi (continuum BC at r=0), F_edge from outermost node.
# ---------------------------------------------------------------------------

def _antideriv(Fv: float) -> float:
    return Fv / 2.0 - math.sin(2.0 * Fv) / 4.0


def compute_B_analytic(positions: np.ndarray, ref: np.ndarray) -> float:
    """Baryon number via the radial integral formula (exact for spherical hedgehog).

    Uses F(r=0)=pi as the inner boundary (true continuum BC) and reads F at
    the outermost lattice node as the outer boundary.  The integrand sin^2(F)*F'
    vanishes at F=pi, so the missing [0, r_center] strip contributes zero.
    """
    phi_4 = positions[:, :4] - ref[:, :4]
    norms = np.linalg.norm(phi_4, axis=1)
    eps = 1e-10 * (float(np.max(norms)) + 1e-30)
    valid = norms > eps
    cos_F = np.where(valid, phi_4[:, 3] / np.where(valid, norms, 1.0), 1.0)
    cos_F = np.clip(cos_F, -1.0, 1.0)
    F_nodes = np.arccos(cos_F)

    coords = ref[:, :3]
    r = np.linalg.norm(coords - coords.mean(axis=0), axis=1)
    i_edge = int(np.argmax(r))
    F_edge = float(F_nodes[i_edge])

    F_inner = math.pi  # continuum BC
    return -(2.0 / math.pi) * (_antideriv(F_edge) - _antideriv(F_inner))


# ---------------------------------------------------------------------------
# Boundary identification and Dirichlet vacuum BC
# ---------------------------------------------------------------------------

def identify_boundary(lattice: SpacelikeLattice, n_shells: int = 1) -> np.ndarray:
    nb = lattice.neighbors
    mask = np.any(nb < 0, axis=1)
    for _ in range(n_shells - 1):
        valid_nb = nb[mask][nb[mask] >= 0]
        mask[valid_nb] = True
    return mask


def apply_dirichlet(R: np.ndarray, ref: np.ndarray,
                    is_boundary: np.ndarray, u0: float) -> np.ndarray:
    """Freeze boundary nodes to north-pole vacuum: (0,0,0,u0) displacement."""
    R_out = R.copy()
    R_out[is_boundary, :4] = ref[is_boundary, :4]
    R_out[is_boundary, 3] = ref[is_boundary, 3] + u0
    return R_out


# ---------------------------------------------------------------------------
# FIRE minimizer (Bitzek et al. 2006)
# ---------------------------------------------------------------------------
# FIRE parameters (standard Bitzek 2006 values):
_FIRE_N_MIN = 5       # min steps after last power < 0 before accelerating
_FIRE_F_INC = 1.1     # step-size growth factor
_FIRE_F_DEC = 0.5     # step-size reduction factor
_FIRE_ALPHA_START = 0.1   # initial mixing fraction
_FIRE_F_ALPHA = 0.99  # mixing fraction decay


def fire_minimize(
    R_init: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    is_boundary: np.ndarray,
    u0: float,
    *,
    max_steps: int = 20000,
    dt_init: float = 2e-3,
    dt_max: float = 2e-2,
    tol_grad: float = 1e-4,
    report_every: int = 2000,
    verbose: bool = True,
) -> dict[str, Any]:
    """FIRE (Fast Inertial Relaxation Engine) constrained minimizer.

    Minimizes spacelike_potential(R) with boundary nodes frozen.
    No physics claim; pure V-minimization tool.

    The FIRE velocity vector is purely a minimization artifact; it has no
    dynamical interpretation and is discarded after convergence.

    Returns
    -------
    dict with keys: R, V, grad_norm, status, n_steps, walltime_s, min_dist
    """
    R = apply_dirichlet(R_init.copy(), ref, is_boundary, u0)
    n_nodes, m_ambient = R.shape

    vel = np.zeros_like(R)  # FIRE velocity (minimization artifact)
    alpha = _FIRE_ALPHA_START
    dt = dt_init
    n_pos_steps = 0  # steps since last power-negative event

    status = "max_steps_reached"
    t0 = time.perf_counter()

    for step in range(max_steps):
        # Force = -dV/dR
        F = spacelike_force(R, lattice, params)
        F[is_boundary] = 0.0  # boundary frozen

        grad_norm = float(np.linalg.norm(F))

        if step % report_every == 0 and verbose:
            V_cur = float(spacelike_potential(R, lattice, params))
            B_cur = compute_B_analytic(R, ref)
            print(f"    FIRE step {step:6d}: V={V_cur:.5e}, |grad|={grad_norm:.4e}, "
                  f"B={B_cur:.4f}, dt={dt:.3e}, alpha={alpha:.4f}")

        if grad_norm < tol_grad:
            status = "converged"
            if verbose:
                print(f"    FIRE converged at step {step}: |grad|={grad_norm:.4e}")
            break

        # FIRE update
        power = float(np.sum(F * vel))

        if power > 0.0:
            n_pos_steps += 1
            # Mix velocity toward force direction
            F_norm = float(np.linalg.norm(F))
            if F_norm > 1e-30:
                vel = (1.0 - alpha) * vel + alpha * (np.linalg.norm(vel) / F_norm) * F
            if n_pos_steps > _FIRE_N_MIN:
                dt = min(dt * _FIRE_F_INC, dt_max)
                alpha *= _FIRE_F_ALPHA
        else:
            # Power negative: reset velocity, reduce step
            vel[:] = 0.0
            alpha = _FIRE_ALPHA_START
            dt *= _FIRE_F_DEC
            n_pos_steps = 0

        # Velocity-Verlet half-step with FIRE velocity
        vel = vel + 0.5 * dt * F

        # Position update
        R_new = R + dt * vel

        # Boundary enforcement
        R_new = apply_dirichlet(R_new, ref, is_boundary, u0)

        # Update force at new position for second half-step
        F_new = spacelike_force(R_new, lattice, params)
        F_new[is_boundary] = 0.0

        vel = vel + 0.5 * dt * F_new
        vel[is_boundary] = 0.0  # boundary velocities irrelevant

        R = R_new

    # Final metrics
    V_final = float(spacelike_potential(R, lattice, params))
    F_final = spacelike_force(R, lattice, params)
    F_final[is_boundary] = 0.0
    grad_norm_final = float(np.linalg.norm(F_final))

    # Minimum bond distance (collapse watch)
    nb = lattice.neighbors
    min_dist = float("inf")
    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        valid = q_arr >= 0
        p_idx = np.where(valid)[0]
        if len(p_idx) == 0:
            continue
        q_idx = q_arr[p_idx]
        dists = np.linalg.norm(R[p_idx] - R[q_idx], axis=1)
        min_dist = min(min_dist, float(np.min(dists)))

    elapsed = time.perf_counter() - t0

    return {
        "R": R,
        "V": V_final,
        "grad_norm": grad_norm_final,
        "status": status,
        "n_steps": step + 1,
        "walltime_s": elapsed,
        "min_dist": min_dist,
    }


# ---------------------------------------------------------------------------
# Strain-weighted spread metrics (interior bonds only)
# ---------------------------------------------------------------------------

def compute_spread(positions: np.ndarray, ref: np.ndarray,
                   lattice: SpacelikeLattice,
                   is_boundary: np.ndarray) -> dict[str, float]:
    """Interior-only strain-weighted confinement metrics.

    Excludes bonds incident to boundary nodes so the frozen vacuum shell does
    not inflate spread estimates.
    """
    nb = lattice.neighbors.copy()
    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        bnd_nb = (q_arr >= 0) & is_boundary[q_arr]
        nb[bnd_nb, nb_idx] = -1
    nb[is_boundary] = -1

    return confinement_metrics_per_slice(
        positions, ref,
        dim=3,
        confinement_radius_factor=0.5,
        weight_mode="strain",
        _neighbor_table=nb,
    )


# ---------------------------------------------------------------------------
# F-edge calculation (used for isolation check and B estimator warning)
# ---------------------------------------------------------------------------

def f_edge_power2(grid_n: int, w: float, a: float = 1.0) -> float:
    r_edge = grid_n * a / 2.0
    return math.pi / (1.0 + (r_edge / w) ** 2)


# ---------------------------------------------------------------------------
# Single-box run
# ---------------------------------------------------------------------------

def run_box(
    grid_n: int,
    w: float = 5.0,
    alpha: float = 0.7,
    u0: float = 1.0,
    a: float = 1.0,
    k_s: float = 1.0,
    rho: float = 1.0,
    n_boundary_shells: int = 1,
    max_steps: int = 20000,
    dt_init: float = 2e-3,
    dt_max: float = 2e-2,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one box at fixed w, alpha, varying grid_n.

    Key design choices
    ------------------
    - w=5a places the soliton width firmly in the continuum corner (a/w=0.2).
    - n_boundary_shells=1: freeze only the surface layer; minimises boundary
      contamination of spread metrics while maintaining topology.
    - FIRE minimizer: converges reliably even when plain gradient descent stalls.
    - Hessian: SKIPPED here (grid_n >= 64 → n_dofs >> 10^5; dense Hessian not
      feasible; the Derrick verdict from V*/N_nodes + radius_rms is the primary
      localization metric and does not need the Hessian).
    """
    box_size = grid_n * a
    box_half = box_size / 2.0
    F_edge = f_edge_power2(grid_n, w, a)
    isolation_clean = F_edge < 0.05
    isolation_marginal = F_edge < 0.10

    print(f"\n{'='*70}")
    print(f"INTERMEDIATE-w BOX STUDY: grid_n={grid_n}, box={box_size:.0f}a, "
          f"w={w}a, w/a={w/a:.1f}")
    print(f"  alpha={alpha}, u0={u0}, a={a}")
    print(f"  box/w={box_size/w:.1f}, F(edge)={F_edge:.4f} "
          f"({'clean' if isolation_clean else 'marginal' if isolation_marginal else 'POOR ISOLATION'})")
    if not isolation_marginal:
        print(f"  WARNING: F(edge)={F_edge:.4f} > 0.10 -- box too small for w={w}; "
              f"B estimator unreliable, spread_ratio inflated by BCs.")
    print(f"{'='*70}")

    m_ambient = 4
    lp = LatticeParams(
        grid_shape=(grid_n,) * 3,
        spacing=a,
        periodic_axes=(False, False, False),
    )
    ap = ActionParams(
        k_s=k_s, alpha=alpha, rho=rho,
        dt=0.1, n_slices=1, m_ambient=m_ambient,
    )
    lattice = SpacelikeLattice(lp)
    ref = lattice.reference_positions(m_ambient)

    t_build = time.perf_counter()

    # Seed
    R_seed, _ = skyrme_twisted_hedgehog(
        lattice, m=m_ambient, u0=u0, w=w, profile_shape="power2"
    )
    B_seed = compute_B_analytic(R_seed, ref)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    print(f"  Seed: B_seed(analytic)={B_seed:.4f}, V_seed={V_seed:.4e}")
    if abs(B_seed - 1.0) > 0.10:
        print(f"  WARNING: B_seed={B_seed:.4f} deviates >0.1 from 1.0 -- "
              f"seed isolation poor (box too small for this w).")

    # Boundary
    is_boundary = identify_boundary(lattice, n_shells=n_boundary_shells)
    n_boundary = int(np.sum(is_boundary))
    n_interior = int(np.sum(~is_boundary))
    n_dofs = n_interior * m_ambient
    print(f"  Nodes: {lattice.n_nodes:,}  boundary={n_boundary} "
          f"({100*n_boundary/lattice.n_nodes:.1f}%)  interior={n_interior:,}  "
          f"n_dofs={n_dofs:,}")

    R_init = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)
    B_init = compute_B_analytic(R_init, ref)
    t_build_elapsed = time.perf_counter() - t_build
    print(f"  After Dirichlet BC: B_init={B_init:.4f}  (build time: {t_build_elapsed:.1f}s)")

    # FIRE minimization
    print(f"\n  FIRE minimizer (max {max_steps} steps, tol={tol_grad:.0e}, "
          f"dt_init={dt_init:.2e}, dt_max={dt_max:.2e})...")
    result = fire_minimize(
        R_init, ref, lattice, ap, is_boundary, u0,
        max_steps=max_steps,
        dt_init=dt_init,
        dt_max=dt_max,
        tol_grad=tol_grad,
        report_every=max(max_steps // 8, 500),
        verbose=verbose,
    )

    R_final = result["R"]
    V_star = result["V"]
    B_final = compute_B_analytic(R_final, ref)
    grad_norm_final = result["grad_norm"]
    not_converged = result["status"] != "converged"

    print(f"\n  Minimizer result:")
    print(f"    status={result['status']}, steps={result['n_steps']}, "
          f"walltime={result['walltime_s']:.1f}s")
    print(f"    V*={V_star:.5e}, grad_norm={grad_norm_final:.4e}, "
          f"B_final={B_final:.4f}, min_dist={result['min_dist']:.4f}a")
    if not_converged:
        print(f"    *** NOT CONVERGED (grad_norm={grad_norm_final:.4e} > tol={tol_grad:.0e}) ***")
        print(f"    V* and spread_ratio from a non-converged run.  "
              f"Interpret with caution; V* lower bound only.")

    # Spread metrics
    spread = compute_spread(R_final, ref, lattice, is_boundary)
    spread_ratio = spread["spread_ratio"]
    radius_rms = spread["radius_rms"]
    box_fill_radius = spread["box_fill_radius"]
    confined_fraction = spread["confined_fraction"]

    print(f"\n  Localization (interior-only strain weight):")
    print(f"    spread_ratio = {spread_ratio:.4f}  (<<0.5 = localized; ~0.7 = box-fill)")
    print(f"    radius_rms   = {radius_rms:.4f}a  (KEY: box-independent if self-localized)")
    print(f"    box_fill_rad = {box_fill_radius:.4f}a  (scales with grid_n)")
    print(f"    confined_frac= {confined_fraction:.4f}  (at 0.5*box_fill_radius)")

    # V* per node (normalization for Derrick scaling check)
    V_per_node = V_star / lattice.n_nodes
    V_per_interior = V_star / n_interior

    print(f"\n  Energy normalization:")
    print(f"    V*/N_nodes    = {V_per_node:.5e}  (constant across boxes = boundary-confined)")
    print(f"    V*/N_interior = {V_per_interior:.5e}")
    print(f"    V*/w^3        = {V_star / (w**3):.5e}  (box-independent if self-localized)")

    # Peak memory estimate (approximate)
    mem_estimate_mb = lattice.n_nodes * m_ambient * 8 * 12 / 1e6  # 12 arrays

    print(f"\n  SUMMARY grid_n={grid_n}:")
    print(f"    V*={V_star:.4e}  V*/N_nodes={V_per_node:.4e}  V*/w^3={V_star/(w**3):.4e}")
    print(f"    spread_ratio={spread_ratio:.4f}  radius_rms={radius_rms:.4f}a")
    print(f"    B_final={B_final:.4f}  grad_norm={grad_norm_final:.4e}")
    print(f"    F(edge)={F_edge:.4f}  isolation={'clean' if isolation_clean else 'marginal'}")
    print(f"    walltime={result['walltime_s']:.1f}s  peak_mem_est~{mem_estimate_mb:.0f}MB")
    if not_converged:
        print(f"    *** CONVERGENCE CAVEAT: results from non-converged run ***")

    return {
        # Config
        "grid_n": grid_n,
        "w": w,
        "w_over_a": w / a,
        "alpha": alpha,
        "u0": u0,
        "a": a,
        "k_s": k_s,
        "box_size": float(box_size),
        "box_half": float(box_half),
        "box_over_w": box_size / w,
        "F_edge_seed": F_edge,
        "isolation_clean": isolation_clean,
        "n_nodes": lattice.n_nodes,
        "n_boundary": n_boundary,
        "n_interior": n_interior,
        "n_dofs": n_dofs,
        # Seed
        "B_seed": B_seed,
        "V_seed": V_seed,
        "B_init": B_init,
        # Minimizer
        "minimizer": "FIRE",
        "descent_status": result["status"],
        "not_converged": not_converged,
        "n_steps": result["n_steps"],
        "walltime_s": result["walltime_s"],
        "peak_mem_est_mb": mem_estimate_mb,
        "min_dist_final": result["min_dist"],
        # Energy
        "V_star": V_star,
        "V_per_node": V_per_node,
        "V_per_interior": V_per_interior,
        "V_over_w3": V_star / (w ** 3),
        "V_positive": V_star > 0,
        # Winding
        "B_final": B_final,
        "B_preserved": abs(B_final - 1.0) < 0.15,
        "grad_norm_final": grad_norm_final,
        # Localization
        "spread_ratio": spread_ratio,
        "radius_rms": radius_rms,
        "box_fill_radius": box_fill_radius,
        "confined_fraction": confined_fraction,
        "localized": spread_ratio < 0.5,
    }


# ---------------------------------------------------------------------------
# Validation run (smoke test — fast, confirms harness and FIRE speedup)
# ---------------------------------------------------------------------------

def run_validate(verbose: bool = True) -> None:
    """Fast smoke test: grid_n=32, w=3a, alpha=0.7.

    w=3a on grid_n=32: r_edge=16, F_edge=pi/(1+(16/3)^2)=pi/29.7=0.106 (marginal
    but tractable for a validation run).  Demonstrates FIRE convergence speed.
    Also runs 200 steps of plain gradient descent for comparison.
    """
    print("\n" + "="*70)
    print("VALIDATION RUN: grid_n=32, w=3a, alpha=0.7")
    print("Goal: confirm harness runs, B=1 preserved, FIRE converges faster than GD")
    print("="*70)

    m_ambient = 4
    a = 1.0
    w = 3.0
    alpha = 0.7
    u0 = 1.0
    grid_n = 32

    F_edge = f_edge_power2(grid_n, w, a)
    print(f"\nIsolation check: F(edge)={F_edge:.4f} (marginal for validation; "
          f"not used for verdict)")

    lp = LatticeParams(
        grid_shape=(grid_n,) * 3,
        spacing=a,
        periodic_axes=(False, False, False),
    )
    ap = ActionParams(k_s=1.0, alpha=alpha, rho=1.0, dt=0.1, n_slices=1, m_ambient=m_ambient)
    lattice = SpacelikeLattice(lp)
    ref = lattice.reference_positions(m_ambient)

    R_seed, _ = skyrme_twisted_hedgehog(lattice, m=m_ambient, u0=u0, w=w, profile_shape="power2")
    B_seed = compute_B_analytic(R_seed, ref)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    print(f"Seed: B_seed={B_seed:.4f}, V_seed={V_seed:.4e}")

    is_boundary = identify_boundary(lattice, n_shells=1)
    n_interior = int(np.sum(~is_boundary))
    print(f"Interior nodes: {n_interior:,} (n_dofs={n_interior*m_ambient:,})")
    R_init = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)

    # --- FIRE run ---
    print("\n[A] FIRE minimizer (max 5000 steps, tol=1e-4):")
    t0 = time.perf_counter()
    res_fire = fire_minimize(
        R_init.copy(), ref, lattice, ap, is_boundary, u0,
        max_steps=5000,
        dt_init=2e-3,
        dt_max=2e-2,
        tol_grad=1e-4,
        report_every=500,
        verbose=verbose,
    )
    t_fire = time.perf_counter() - t0
    B_fire = compute_B_analytic(res_fire["R"], ref)
    print(f"  FIRE: status={res_fire['status']}, steps={res_fire['n_steps']}, "
          f"V*={res_fire['V']:.4e}, |grad|={res_fire['grad_norm']:.4e}, "
          f"B={B_fire:.4f}, time={t_fire:.1f}s")

    # --- Plain gradient descent comparison (200 steps only for timing) ---
    print("\n[B] Plain gradient descent (200 steps for timing comparison):")
    R_gd = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)
    V_gd = float(spacelike_potential(R_gd, lattice, ap))
    step_size = 5e-4
    t0 = time.perf_counter()
    for i in range(200):
        F_gd = spacelike_force(R_gd, lattice, ap)
        F_gd[is_boundary] = 0.0
        R_try = R_gd + step_size * F_gd
        V_try = float(spacelike_potential(R_try, lattice, ap))
        if V_try < V_gd:
            R_gd = R_try
            V_gd = V_try
            step_size = min(step_size * 1.05, 5e-3)
        else:
            step_size = max(step_size * 0.5, 5e-7)
        R_gd = apply_dirichlet(R_gd, ref, is_boundary, u0)
    F_gd_final = spacelike_force(R_gd, lattice, ap)
    F_gd_final[is_boundary] = 0.0
    t_gd = time.perf_counter() - t0
    print(f"  GD (200 steps): V={V_gd:.4e}, |grad|={np.linalg.norm(F_gd_final):.4e}, "
          f"time={t_gd:.1f}s")

    # Speedup comparison: note FIRE ran to convergence or 5000 steps
    if res_fire["status"] == "converged":
        print(f"\n  FIRE converged in {res_fire['n_steps']} steps vs GD still at "
              f"|grad|={np.linalg.norm(F_gd_final):.4e} after 200 steps.")
    else:
        print(f"\n  FIRE ran {res_fire['n_steps']} steps without full convergence "
              f"(tol=1e-4); |grad|={res_fire['grad_norm']:.4e}.")

    # Spread metrics
    spread = compute_spread(res_fire["R"], ref, lattice, is_boundary)
    print(f"\n  Spread (FIRE final):")
    print(f"    spread_ratio={spread['spread_ratio']:.4f}  radius_rms={spread['radius_rms']:.4f}a")
    print(f"    B_final={B_fire:.4f}  B_preserved={abs(B_fire-1.0)<0.15}")

    # Save validation result
    out_dir = _output_dir()
    val_path = os.path.join(out_dir, "validation_g32_w3.json")
    val_result = {
        "mode": "validation",
        "grid_n": grid_n, "w": w, "alpha": alpha,
        "B_seed": B_seed, "V_seed": V_seed,
        "F_edge": F_edge,
        "FIRE_status": res_fire["status"],
        "FIRE_steps": res_fire["n_steps"],
        "FIRE_V_star": res_fire["V"],
        "FIRE_grad_norm": res_fire["grad_norm"],
        "FIRE_B_final": float(B_fire),
        "FIRE_walltime_s": res_fire["walltime_s"],
        "GD_200steps_V": float(V_gd),
        "GD_200steps_grad_norm": float(np.linalg.norm(F_gd_final)),
        "GD_200steps_walltime_s": float(t_gd),
        "spread_ratio": spread["spread_ratio"],
        "radius_rms": spread["radius_rms"],
    }
    with open(val_path, "w") as f:
        json.dump(val_result, f, indent=2)
    print(f"\nValidation result saved: {val_path}")

    # PASS/FAIL
    b_ok = abs(B_fire - 1.0) < 0.15
    v_ok = res_fire["V"] > 0
    fire_faster = (res_fire["grad_norm"] < np.linalg.norm(F_gd_final)
                   or res_fire["status"] == "converged")
    print(f"\n  VALIDATION {'PASS' if (b_ok and v_ok and fire_faster) else 'PARTIAL'}:")
    print(f"    B=1 preserved: {'YES' if b_ok else 'NO'}")
    print(f"    V* > 0:        {'YES' if v_ok else 'NO'}")
    print(f"    FIRE faster:   {'YES' if fire_faster else 'check'} "
          f"(FIRE |grad|={res_fire['grad_norm']:.4e} vs GD |grad|={np.linalg.norm(F_gd_final):.4e})")


# ---------------------------------------------------------------------------
# Full production run
# ---------------------------------------------------------------------------

def run_production(
    w: float = 5.0,
    alpha: float = 0.7,
    grid_sizes: list[int] | None = None,
    max_steps: int = 20000,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> None:
    """Full box-doubling study at fixed w=5a, alpha=0.7.

    grid_n in {64, 96, 128} by default.
    """
    if grid_sizes is None:
        grid_sizes = [64, 96, 128]

    out_dir = _output_dir()
    results = []

    print(f"\n{'='*70}")
    print(f"INTERMEDIATE-w BOX-DOUBLING STUDY")
    print(f"  w={w}a, alpha={alpha}, grid_n in {grid_sizes}")
    print(f"  Verdict: V*/N_nodes + radius_rms vs grid_n")
    print(f"{'='*70}\n")

    for grid_n in grid_sizes:
        try:
            r = run_box(
                grid_n=grid_n,
                w=w,
                alpha=alpha,
                u0=1.0,
                a=1.0,
                k_s=1.0,
                rho=1.0,
                n_boundary_shells=1,
                max_steps=max_steps,
                dt_init=2e-3,
                dt_max=2e-2,
                tol_grad=tol_grad,
                verbose=verbose,
            )
        except Exception as exc:
            import traceback
            print(f"\n  ERROR for grid_n={grid_n}: {exc}")
            traceback.print_exc()
            r = {
                "grid_n": grid_n, "w": w, "alpha": alpha,
                "error": str(exc), "descent_status": "error",
            }
        results.append(r)

        # Save per-box JSON immediately (so partial runs are recoverable)
        json_path = os.path.join(out_dir, f"intermediate_w_g{grid_n}.json")
        r_json = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
                  for k, v in r.items() if k != "R"}
        with open(json_path, "w") as f:
            json.dump(r_json, f, indent=2, default=str)
        print(f"\n  JSON saved: {json_path}")

    # Write CSV
    csv_path = os.path.join(out_dir, "intermediate_w_sweep.csv")
    all_keys: set[str] = set()
    for r in results:
        all_keys.update(r.keys())
    fieldnames = sorted(k for k in all_keys if k != "R")

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in fieldnames}
            writer.writerow(row)
    print(f"\nCSV written: {csv_path}")

    # ---------- Verdict table ----------
    print(f"\n{'='*70}")
    print("BOX-DOUBLING VERDICT TABLE (intermediate w)")
    print(f"{'='*70}")
    print(f"{'grid_n':>8} {'box':>6} {'box/w':>6} {'F_edge':>8} "
          f"{'V*':>12} {'V*/N_nodes':>12} {'V*/w^3':>12} "
          f"{'radius_rms':>11} {'spread':>8} {'B_final':>8} "
          f"{'|grad|':>9} {'steps':>7} {'status':>20}")
    print("-" * 130)
    for r in results:
        if "error" in r:
            print(f"  grid_n={r['grid_n']}: ERROR - {r['error']}")
            continue
        status_str = r.get("descent_status", "?")
        if r.get("not_converged", False):
            status_str += "(!)"
        print(
            f"{r['grid_n']:>8} {r['box_size']:>6.0f}a {r['box_over_w']:>6.1f} "
            f"{r['F_edge_seed']:>8.4f} "
            f"{r['V_star']:>12.4e} {r['V_per_node']:>12.4e} {r['V_over_w3']:>12.4e} "
            f"{r['radius_rms']:>11.4f} {r['spread_ratio']:>8.4f} {r['B_final']:>8.4f} "
            f"{r['grad_norm_final']:>9.4e} {r['n_steps']:>7d} {status_str:>20}"
        )

    # ---------- Trend analysis ----------
    good = [r for r in results if "error" not in r and not r.get("not_converged", False)]
    print(f"\n  Converged runs: {len(good)}/{len(results)}")

    if len(good) >= 2:
        r0, r1 = good[0], good[-1]
        box_ratio = r1["box_size"] / r0["box_size"]
        radius_ratio = (r1["radius_rms"] / r0["radius_rms"]
                        if r0["radius_rms"] > 0 else float("nan"))
        V_ratio = r1["V_star"] / r0["V_star"] if r0["V_star"] > 0 else float("nan")
        V_per_node_ratio = (r1["V_per_node"] / r0["V_per_node"]
                            if r0["V_per_node"] > 0 else float("nan"))
        V_w3_ratio = (r1["V_over_w3"] / r0["V_over_w3"]
                      if r0["V_over_w3"] > 0 else float("nan"))
        spread_delta = r1["spread_ratio"] - r0["spread_ratio"]

        print(f"\n  Trend (largest/smallest converged box):")
        print(f"    box_ratio       = {box_ratio:.2f}x")
        print(f"    radius_rms ratio= {radius_ratio:.3f}  "
              f"(~1.0 = self-localized; ~{box_ratio:.1f} = boundary-confined)")
        print(f"    V* ratio        = {V_ratio:.3f}")
        print(f"    V*/N_nodes ratio= {V_per_node_ratio:.3f}  "
              f"(~1.0 = boundary-confined; << 1.0 = self-localized)")
        print(f"    V*/w^3 ratio    = {V_w3_ratio:.3f}  "
              f"(~1.0 = self-localized geometric scaling)")
        print(f"    spread_ratio delta = {spread_delta:+.4f}  "
              f"(negative = shrinks relative to box = self-localized)")
    else:
        print(f"\n  WARNING: fewer than 2 converged runs. Cannot compute trend ratios.")
        print(f"  Consider increasing --max-steps or reducing --tol.")
        non_conv = [r for r in results if "error" not in r and r.get("not_converged", False)]
        for r in non_conv:
            print(f"    grid_n={r['grid_n']}: grad_norm={r.get('grad_norm_final','?'):.4e}, "
                  f"V*={r.get('V_star','?'):.4e} (unconverged)")

    # ---------- Derrick verdict ----------
    print(f"\n{'='*70}")
    print("DERRICK VERDICT:")
    if len(good) < 2:
        if len(results) > 0:
            # Report what we have even if unconverged
            all_r = [r for r in results if "error" not in r]
            if len(all_r) >= 2:
                r0, r1 = all_r[0], all_r[-1]
                V_per_node_ratio = (r1["V_per_node"] / r0["V_per_node"]
                                    if r0.get("V_per_node", 0) > 0 else float("nan"))
                radius_ratio = (r1["radius_rms"] / r0["radius_rms"]
                                if r0.get("radius_rms", 0) > 0 else float("nan"))
                box_ratio = r1["box_size"] / r0["box_size"]
                print(f"  INCONCLUSIVE (convergence not achieved).")
                print(f"  Partial evidence: V*/N_nodes ratio={V_per_node_ratio:.3f}, "
                      f"radius_rms ratio={radius_ratio:.3f}, box ratio={box_ratio:.2f}x")
                if V_per_node_ratio > 0.7 or radius_ratio > 0.7 * box_ratio:
                    print(f"  Trend SUGGESTS boundary-confined, but run is NOT converged.")
                elif radius_ratio < 1.3 and not math.isnan(radius_ratio):
                    print(f"  Trend SUGGESTS self-localized, but run is NOT converged.")
            else:
                print(f"  INCONCLUSIVE: insufficient runs.")
    else:
        r_small = good[0]
        r_large = good[-1]
        radius_ratio = r_large["radius_rms"] / r_small["radius_rms"]
        V_per_node_ratio = r_large["V_per_node"] / r_small["V_per_node"]
        box_ratio = r_large["box_size"] / r_small["box_size"]

        if radius_ratio < 1.3 and V_per_node_ratio < 0.5:
            verdict = "SELF-LOCALIZED"
            detail = (
                f"radius_rms changed by {radius_ratio:.3f}x while box grew {box_ratio:.2f}x; "
                f"V*/N_nodes dropped by factor {V_per_node_ratio:.3f} -- "
                f"soliton width and energy are box-independent."
            )
        elif radius_ratio > 0.65 * box_ratio or abs(V_per_node_ratio - 1.0) < 0.15:
            verdict = "BOUNDARY-CONFINED (route dead at intermediate w)"
            detail = (
                f"radius_rms grew {radius_ratio:.3f}x tracking box growth of {box_ratio:.2f}x; "
                f"V*/N_nodes ratio={V_per_node_ratio:.3f} (~1.0). "
                f"The substrate's geometric quartic does NOT bind a self-localized B=1 soliton "
                f"at w={w}a, alpha={alpha}."
            )
        else:
            verdict = "AMBIGUOUS"
            detail = (
                f"radius_rms ratio={radius_ratio:.3f}, V*/N_nodes ratio={V_per_node_ratio:.3f}, "
                f"box ratio={box_ratio:.2f}x. "
                f"More boxes or tighter convergence needed."
            )

        print(f"  {verdict}")
        print(f"  {detail}")

    print(f"{'='*70}\n")

    print(f"All outputs in: {out_dir}")
    print(f"CSV: {csv_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Intermediate-w box-doubling: decisive B=1 soliton localization test"
    )
    parser.add_argument(
        "--validate", action="store_true",
        help="Smoke test: grid_n=32, w=3a; confirm harness + FIRE speedup (~2min)",
    )
    parser.add_argument(
        "--run", action="store_true",
        help="Full production: w=5a, grid_n in {64, 96, 128}",
    )
    parser.add_argument(
        "--w", type=float, default=5.0,
        help="Soliton profile half-width in lattice units (default 5.0)",
    )
    parser.add_argument(
        "--alpha", type=float, default=0.7,
        help="Prestress parameter (default 0.7)",
    )
    parser.add_argument(
        "--grid-sizes", type=int, nargs="+", default=None,
        help="Override grid_n list, e.g. --grid-sizes 64 96 (default: 64 96 128)",
    )
    parser.add_argument(
        "--max-steps", type=int, default=20000,
        help="Max FIRE steps per box (default 20000)",
    )
    parser.add_argument(
        "--tol", type=float, default=1e-4,
        help="FIRE gradient convergence tolerance (default 1e-4)",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Reduce per-step verbosity",
    )
    args = parser.parse_args()

    verbose = not args.quiet

    if args.validate:
        run_validate(verbose=verbose)
    elif args.run or not (args.validate):
        # Default when called with no flags (e.g. from AWS): run production.
        run_production(
            w=args.w,
            alpha=args.alpha,
            grid_sizes=args.grid_sizes,
            max_steps=args.max_steps,
            tol_grad=args.tol,
            verbose=verbose,
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
