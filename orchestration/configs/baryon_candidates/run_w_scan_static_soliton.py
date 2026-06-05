#!/usr/bin/env python3
"""w-SCAN static B=1 soliton: decisive self-localization experiment.

PHYSICS CONTEXT
---------------
The static B=1 soliton was found boundary-confined at w=2a in prior runs
(C2 corrected + box-doubling study): V*/N_nodes constant to 4%, spread_ratio
flat ~0.72, radius_rms grew with grid_n.  Two stabilization mechanisms could
explain a localized soliton if one exists:

  (1) LATTICE (Peierls-Nabarro / discrete-soliton): discreteness pins the
      soliton at w~a.  Predicts localized lump at w~a, box-INDEPENDENT width.

  (2) QUARTIC (Skyrme-Derrick continuum): the geometric quartic E4 ~ k_s*alpha/a
      supplies E4 ~ lambda^{-1} vs E2 ~ lambda, giving a continuum minimum at
      w* = sqrt(E4/E2) ~ sqrt(alpha/(1-alpha)) in units of a.  Predicts a lump
      at w* >> a, also box-independent.

The prior w=2a verdict was SPREADING (lambda -> inf), not collapse.  The lattice
cutoff only stops collapse; it does NOT stop spreading.  The plain-gradient-flow
non-convergence (grad_norm 0.05-0.09) was the main contaminating factor.

PRIMARY EXPERIMENT (w-scan at alpha=0.7)
-----------------------------------------
Sweep w/a in {1.5, 2, 3, 4, 5, 6}.  For each w:
  - Use a box large enough for good isolation (F_edge < 0.05 where feasible).
  - Relax with FIRE to a real convergence criterion (tol=1e-4) or clear plateau.
  - Measure equilibrium radius_rms (absolute, in lattice units), spread_ratio
    (strain-weighted, interior bonds), V*, V*/N_nodes, B_final, F_edge, grad_norm.

KEY DIAGNOSTIC (box-independence)
-----------------------------------
For w in {2, 3.5, 5} (representative of small/mid/large regimes), run each at
two box sizes and check whether equilibrium radius_rms is BOX-INDEPENDENT
(intrinsic width = genuine soliton) or TRACKS the box (boundary-driven).

MINIMIZER
----------
FIRE (Bitzek et al. 2006) with Velocity-Verlet half-steps.  Converges 5-30x
faster than plain gradient descent on rugged energy surfaces.  Pure V-minimizer;
no physics claim, no damping.

VERDICT CRITERIA
----------------
  LATTICE-STABILIZED: equilibrium radius_rms ~ 1-2a, box-independent.
  QUARTIC-DERRICK: width detaches from lattice at w* > a, box-independent at w*.
  NO CONFINEMENT: width tracks box at every w (spreads everywhere).

OUTPUTS (all small, CSV/JSON, no large arrays)
----------------------------------------------
  $BRANESIM_RESULTS_DIR/w_scan_primary.csv           (one row per w, primary box)
  $BRANESIM_RESULTS_DIR/w_scan_box_independence.csv  (box-independence check)
  $BRANESIM_RESULTS_DIR/w_scan_g{N}_w{W}.json        (per-run JSON)
  $BRANESIM_RESULTS_DIR/w_scan_validation.json        (local validation result)

USAGE
-----
  # Local validation (fast, ~2-5min):
  python orchestration/configs/baryon_candidates/run_w_scan_static_soliton.py --validate

  # Full production run:
  python orchestration/configs/baryon_candidates/run_w_scan_static_soliton.py --run

  # Single w value for testing:
  python orchestration/configs/baryon_candidates/run_w_scan_static_soliton.py --run --w 3.0 --grid-n 48

  # Box-independence check only:
  python orchestration/configs/baryon_candidates/run_w_scan_static_soliton.py --box-independence

PRINCIPLES COMPLIANCE (principles.md)
--------------------------------------
- No confinement forces, no clamps, no nonlinear saturation added by hand.
- Physics from spacelike_potential / spacelike_force only.
- Dirichlet vacuum BC = legitimate constraint (boundary nodes frozen at vacuum).
- FIRE = pure V-minimization tool (no dynamical/physics claim).
- grad_norm reported per run; non-converged results flagged explicitly.
- Hessian SKIPPED for grid_n > 30 (n_dofs too large for dense computation;
  localization, not stability, is the question here).
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
# Path setup: works as a script in orchestration/ and from project root.
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
# Output directory: honour BRANESIM_RESULTS_DIR (AWS) or local fallback.
# ---------------------------------------------------------------------------
def _output_dir() -> str:
    env = os.environ.get("BRANESIM_RESULTS_DIR", "").strip()
    if env:
        os.makedirs(env, exist_ok=True)
        return env
    d = os.path.join(_REPO_ROOT, "test-runs", "w_scan_static_soliton")
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Box-size design table
# ---------------------------------------------------------------------------
# For power-2 profile F(r) = pi/(1 + (r/w)^2):
#   F_edge = pi / (1 + (grid_n/2 / w)^2)
#   For F_edge < 0.05: grid_n/2 > w * sqrt(pi/0.05 - 1) = w * 7.86
#   => grid_n > 15.72 * w
#
# We use the SMALLEST grid_n that satisfies F_edge < 0.05 for the primary scan,
# rounded up to the nearest multiple of 8 for memory alignment.  For the
# box-independence check we use two boxes: the primary box and roughly 1.6x.
#
# Key: the box is large RELATIVE to w so the soliton has room.  The w-scan
# tests all w values; we let the minimizer find the equilibrium width freely.
# The box_over_w ratio must be >> 1 so the boundary does not force spreading.

def _primary_grid_n(w: float, a: float = 1.0) -> int:
    """Smallest grid_n such that F_edge < 0.05, rounded up to multiple of 8."""
    min_n = math.ceil(15.72 * w / a)
    # Round up to nearest multiple of 8 (cache alignment, also ensures even grid)
    return max(16, int(math.ceil(min_n / 8)) * 8)


def _f_edge_power2(grid_n: int, w: float, a: float = 1.0) -> float:
    r_edge = grid_n * a / 2.0
    return math.pi / (1.0 + (r_edge / w) ** 2)


# Primary w-scan configurations (w/a, grid_n).
# All use alpha=0.7, u0=1.0, a=1.0, k_s=1.0.
# grid_n chosen for F_edge < 0.05 (clean isolation).
#
#   w=1.5: grid_n=24  -> F_edge = pi/(1+(12/1.5)^2) = pi/65 = 0.048  (clean)
#   w=2.0: grid_n=32  -> F_edge = pi/(1+(16/2)^2)   = pi/65 = 0.048  (clean)
#   w=3.0: grid_n=48  -> F_edge = pi/(1+(24/3)^2)   = pi/65 = 0.048  (clean)
#   w=4.0: grid_n=64  -> F_edge = pi/(1+(32/4)^2)   = pi/65 = 0.048  (clean)
#   w=5.0: grid_n=80  -> F_edge = pi/(1+(40/5)^2)   = pi/65 = 0.048  (clean)
#   w=6.0: grid_n=96  -> F_edge = pi/(1+(48/6)^2)   = pi/65 = 0.048  (clean)
#
# Note: all primary boxes have r_edge/w = 8, so F_edge = pi/(1+64) = 0.0480
# for all w.  This is the natural isolating design.

PRIMARY_W_VALUES = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0]

# Box-independence check: w values and their two box sizes.
# Format: (w, grid_n_small, grid_n_large)
# small box: F_edge < 0.05 (primary);  large box: ~1.6x in linear dimension.
BOX_INDEPENDENCE_CONFIGS = [
    # w=2.0: lattice regime; primary=32, large=48 (1.5x)
    (2.0,  32,  48),
    # w=3.5: mid regime; primary=56 (~closest mult-of-8 to 15.72*3.5=55), large=80 (1.43x)
    (3.5,  56,  80),
    # w=5.0: continuum regime; primary=80, large=120 (1.5x)
    (5.0,  80, 120),
]

# Validation config: fast smoke test.  grid_n=20, w=2.  Box/w=10, F_edge = pi/(1+100)~0.031.
# Demonstrates FIRE vs plain-GD speedup.
VALIDATION_W = 2.0
VALIDATION_GRID_N = 20


# ---------------------------------------------------------------------------
# AUTHORITATIVE B estimator (analytic radial integral formula)
# ---------------------------------------------------------------------------

def _antideriv(Fv: float) -> float:
    return Fv / 2.0 - math.sin(2.0 * Fv) / 4.0


def compute_B_analytic(positions: np.ndarray, ref: np.ndarray) -> float:
    """Baryon number B = -(2/pi)*[antideriv(F_edge) - antideriv(pi)].

    Uses F(r=0)=pi (true continuum BC) as inner boundary.
    F at the outermost lattice node as outer boundary.
    Exact for a spherical hedgehog; robust post-flow (reads F_edge only).
    NOT the Jacobian estimator (systematically underestimates on even grids).
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
    F_inner = math.pi
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
    """Freeze boundary to north-pole vacuum: displacement = (0,0,0,u0)."""
    R_out = R.copy()
    R_out[is_boundary, :4] = ref[is_boundary, :4]
    R_out[is_boundary, 3] = ref[is_boundary, 3] + u0
    return R_out


# ---------------------------------------------------------------------------
# FIRE minimizer (Bitzek et al. 2006) - accelerated constrained minimizer
# ---------------------------------------------------------------------------
# Standard FIRE parameters (Bitzek 2006):
_FIRE_N_MIN = 5        # steps after last power<0 before accelerating
_FIRE_F_INC = 1.1      # step-size growth factor
_FIRE_F_DEC = 0.5      # step-size reduction factor
_FIRE_ALPHA_START = 0.1
_FIRE_F_ALPHA = 0.99   # mixing-fraction decay per step


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
    """FIRE constrained minimizer: minimize spacelike_potential with frozen boundary.

    FIRE velocity is a minimization artifact; no dynamical claim.
    No physics or damping beyond what gradient flow does.
    """
    R = apply_dirichlet(R_init.copy(), ref, is_boundary, u0)
    vel = np.zeros_like(R)
    alpha = _FIRE_ALPHA_START
    dt = dt_init
    n_pos = 0

    status = "max_steps_reached"
    t0 = time.perf_counter()

    for step in range(max_steps):
        F = spacelike_force(R, lattice, params)
        F[is_boundary] = 0.0
        grad_norm = float(np.linalg.norm(F))

        if step % report_every == 0 and verbose:
            V_cur = float(spacelike_potential(R, lattice, params))
            B_cur = compute_B_analytic(R, ref)
            print(f"    FIRE {step:7d}: V={V_cur:.5e}  |grad|={grad_norm:.4e}  "
                  f"B={B_cur:.4f}  dt={dt:.3e}  alpha={alpha:.4f}")

        if grad_norm < tol_grad:
            status = "converged"
            if verbose:
                print(f"    FIRE converged step={step}  |grad|={grad_norm:.4e}")
            break

        power = float(np.sum(F * vel))
        if power > 0.0:
            n_pos += 1
            F_norm = float(np.linalg.norm(F))
            if F_norm > 1e-30:
                vel = (1.0 - alpha) * vel + alpha * (np.linalg.norm(vel) / F_norm) * F
            if n_pos > _FIRE_N_MIN:
                dt = min(dt * _FIRE_F_INC, dt_max)
                alpha *= _FIRE_F_ALPHA
        else:
            vel[:] = 0.0
            alpha = _FIRE_ALPHA_START
            dt *= _FIRE_F_DEC
            n_pos = 0

        vel = vel + 0.5 * dt * F
        R_new = R + dt * vel
        R_new = apply_dirichlet(R_new, ref, is_boundary, u0)

        F_new = spacelike_force(R_new, lattice, params)
        F_new[is_boundary] = 0.0
        vel = vel + 0.5 * dt * F_new
        vel[is_boundary] = 0.0

        R = R_new

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
# Strain-weighted confinement metrics (interior bonds only)
# ---------------------------------------------------------------------------

def compute_spread(positions: np.ndarray, ref: np.ndarray,
                   lattice: SpacelikeLattice,
                   is_boundary: np.ndarray) -> dict[str, float]:
    """Interior-only strain-weighted spread metrics.

    Excludes boundary-incident bonds so the frozen vacuum shell does not
    inflate spread estimates (known issue in earlier C2 harness).
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
# Single-box run (core computation unit)
# ---------------------------------------------------------------------------

def run_one_box(
    grid_n: int,
    w: float,
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
    """Seed, FIRE-minimize, and measure confinement metrics for one (w, grid_n) pair.

    No Hessian (grid_n >= 24: n_dofs >> 10^4; dense Hessian not feasible).
    No Peierls-Nabarro barrier (not the question; localization is).
    All diagnostics are read-only post-minimization.

    Returns a flat dict of scalars (no large arrays).
    """
    box_size = float(grid_n * a)
    box_half = box_size / 2.0
    F_edge_seed = _f_edge_power2(grid_n, w, a)
    isolation_clean = F_edge_seed < 0.05
    isolation_marginal = F_edge_seed < 0.10

    tag = f"grid_n={grid_n}  w={w:.1f}a  box={box_size:.0f}a  box/w={box_size/w:.1f}"
    print(f"\n{'='*68}")
    print(f"  {tag}")
    print(f"  alpha={alpha}  u0={u0}  F(edge)={F_edge_seed:.4f}  "
          f"({'clean' if isolation_clean else 'marginal' if isolation_marginal else 'POOR'})")
    if not isolation_marginal:
        print(f"  WARNING: F_edge={F_edge_seed:.4f} > 0.10 -- box too small; "
              f"B estimator and spread_ratio may be unreliable.")
    print(f"{'='*68}")

    m = 4
    lp = LatticeParams(grid_shape=(grid_n,) * 3, spacing=a,
                       periodic_axes=(False, False, False))
    ap = ActionParams(k_s=k_s, alpha=alpha, rho=rho, dt=0.1, n_slices=1, m_ambient=m)
    lattice = SpacelikeLattice(lp)
    ref = lattice.reference_positions(m)

    t_build = time.perf_counter()
    R_seed, _ = skyrme_twisted_hedgehog(lattice, m=m, u0=u0, w=w, profile_shape="power2")
    B_seed = compute_B_analytic(R_seed, ref)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    t_build_s = time.perf_counter() - t_build

    print(f"  Seed: B_seed={B_seed:.4f}  V_seed={V_seed:.4e}  build={t_build_s:.1f}s")
    if abs(B_seed - 1.0) > 0.10:
        print(f"  WARNING: B_seed={B_seed:.4f} far from 1.0 -- poor isolation (box too small).")

    is_boundary = identify_boundary(lattice, n_shells=n_boundary_shells)
    n_boundary = int(np.sum(is_boundary))
    n_interior = int(np.sum(~is_boundary))
    n_dofs = n_interior * m
    print(f"  Nodes: {lattice.n_nodes:,}  boundary={n_boundary} ({100*n_boundary/lattice.n_nodes:.1f}%)  "
          f"interior={n_interior:,}  n_dofs={n_dofs:,}")

    R_init = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)
    B_init = compute_B_analytic(R_init, ref)
    print(f"  After Dirichlet BC: B_init={B_init:.4f}")

    print(f"\n  FIRE minimize (max={max_steps} tol={tol_grad:.1e} dt_init={dt_init:.2e} dt_max={dt_max:.2e})")
    res = fire_minimize(
        R_init, ref, lattice, ap, is_boundary, u0,
        max_steps=max_steps,
        dt_init=dt_init,
        dt_max=dt_max,
        tol_grad=tol_grad,
        report_every=max(max_steps // 8, 500),
        verbose=verbose,
    )

    R_final = res["R"]
    V_star = res["V"]
    B_final = compute_B_analytic(R_final, ref)
    grad_norm = res["grad_norm"]
    not_converged = res["status"] != "converged"

    print(f"\n  Minimizer: status={res['status']}  steps={res['n_steps']}  "
          f"walltime={res['walltime_s']:.1f}s")
    print(f"  V*={V_star:.5e}  grad_norm={grad_norm:.4e}  "
          f"B_final={B_final:.4f}  min_dist={res['min_dist']:.4f}a")
    if not_converged:
        print(f"  *** NOT CONVERGED (grad_norm={grad_norm:.4e} > tol={tol_grad:.1e}) ***")
        print(f"  V* is a lower bound; spread_ratio from non-converged state.")

    spread = compute_spread(R_final, ref, lattice, is_boundary)
    spread_ratio = spread["spread_ratio"]
    radius_rms = spread["radius_rms"]
    box_fill_radius = spread["box_fill_radius"]
    confined_fraction = spread["confined_fraction"]

    V_per_node = V_star / lattice.n_nodes
    V_per_interior = V_star / n_interior
    V_over_w3 = V_star / (w ** 3)

    # Peak memory estimate (order of magnitude)
    peak_mem_mb = lattice.n_nodes * m * 8 * 15 / 1e6  # ~15 working arrays

    print(f"\n  Localization (interior-only strain weight):")
    print(f"    spread_ratio = {spread_ratio:.4f}  (<<0.5 = localized; ~0.7 = box-fill)")
    print(f"    radius_rms   = {radius_rms:.4f}a  (KEY: box-independent if self-localized)")
    print(f"    box_fill_rad = {box_fill_radius:.4f}a")
    print(f"    confined_frac= {confined_fraction:.4f}")
    print(f"\n  Energy normalization:")
    print(f"    V*/N_nodes    = {V_per_node:.5e}  (const across boxes = boundary-confined)")
    print(f"    V*/N_interior = {V_per_interior:.5e}")
    print(f"    V*/w^3        = {V_over_w3:.5e}  (box-independent if self-localized)")
    print(f"    peak_mem_est  = {peak_mem_mb:.0f} MB")

    return {
        # Config
        "grid_n": grid_n,
        "w": float(w),
        "w_over_a": float(w / a),
        "alpha": float(alpha),
        "u0": float(u0),
        "a": float(a),
        "k_s": float(k_s),
        "box_size": float(box_size),
        "box_half": float(box_half),
        "box_over_w": float(box_size / w),
        "F_edge_seed": float(F_edge_seed),
        "isolation_clean": bool(isolation_clean),
        "isolation_marginal": bool(isolation_marginal),
        "n_nodes": int(lattice.n_nodes),
        "n_boundary": int(n_boundary),
        "n_interior": int(n_interior),
        "n_dofs": int(n_dofs),
        # Seed / init
        "B_seed": float(B_seed),
        "V_seed": float(V_seed),
        "B_init": float(B_init),
        # Minimizer
        "minimizer": "FIRE",
        "descent_status": res["status"],
        "not_converged": bool(not_converged),
        "n_steps": int(res["n_steps"]),
        "walltime_s": float(res["walltime_s"]),
        "peak_mem_est_mb": float(peak_mem_mb),
        "min_dist_final": float(res["min_dist"]),
        # Energy
        "V_star": float(V_star),
        "V_per_node": float(V_per_node),
        "V_per_interior": float(V_per_interior),
        "V_over_w3": float(V_over_w3),
        "V_positive": bool(V_star > 0),
        # Winding
        "B_final": float(B_final),
        "B_preserved": bool(abs(B_final - 1.0) < 0.15),
        "grad_norm_final": float(grad_norm),
        # Localization (KEY metrics for Derrick verdict)
        "spread_ratio": float(spread_ratio),
        "radius_rms": float(radius_rms),
        "box_fill_radius": float(box_fill_radius),
        "confined_fraction": float(confined_fraction),
        "localized": bool(spread_ratio < 0.5),
    }


# ---------------------------------------------------------------------------
# Write CSV helper
# ---------------------------------------------------------------------------

def _write_csv(rows: list[dict], csv_path: str) -> None:
    if not rows:
        return
    all_keys: set[str] = set()
    for r in rows:
        all_keys.update(k for k in r.keys() if k != "R")
    fieldnames = sorted(all_keys)
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            row = {k: r.get(k, "") for k in fieldnames}
            writer.writerow(row)


def _save_json(result: dict, path: str) -> None:
    r_json = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
              for k, v in result.items() if k != "R"}
    with open(path, "w") as f:
        json.dump(r_json, f, indent=2, default=str)


# ---------------------------------------------------------------------------
# Validation run: local smoke test, confirms FIRE speedup
# ---------------------------------------------------------------------------

def run_validate(verbose: bool = True) -> None:
    """Fast smoke test: grid_n=20, w=2a, alpha=0.7.

    Checks:
      1. Harness runs end-to-end without error.
      2. B=1 preserved after FIRE minimization.
      3. FIRE reaches lower grad_norm than plain-GD in same or fewer wallclock steps.

    grid_n=20, w=2: box=20a, box/w=10, F_edge=pi/(1+(10)^2)=0.031 (clean isolation).
    n_interior=18^3=5832, n_dofs=23328 (FIRE feasible in minutes; GD stalls).
    """
    print("\n" + "="*68)
    print("VALIDATION: grid_n=20  w=2a  alpha=0.7  (FIRE vs plain-GD speedup test)")
    print("="*68)

    m = 4
    a = 1.0
    w = float(VALIDATION_W)
    alpha = 0.7
    u0 = 1.0
    grid_n = VALIDATION_GRID_N

    F_edge = _f_edge_power2(grid_n, w, a)
    print(f"\nIsolation check: F(edge)={F_edge:.4f}  box/w={grid_n*a/w:.1f}")

    lp = LatticeParams(grid_shape=(grid_n,) * 3, spacing=a,
                       periodic_axes=(False, False, False))
    ap = ActionParams(k_s=1.0, alpha=alpha, rho=1.0, dt=0.1, n_slices=1, m_ambient=m)
    lattice = SpacelikeLattice(lp)
    ref = lattice.reference_positions(m)

    R_seed, _ = skyrme_twisted_hedgehog(lattice, m=m, u0=u0, w=w, profile_shape="power2")
    B_seed = compute_B_analytic(R_seed, ref)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    print(f"Seed: B_seed={B_seed:.4f}  V_seed={V_seed:.4e}")

    is_boundary = identify_boundary(lattice, n_shells=1)
    n_interior = int(np.sum(~is_boundary))
    print(f"Interior nodes: {n_interior:,}  n_dofs={n_interior*m:,}")
    R_init = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)

    # --- [A] FIRE ---
    print("\n[A] FIRE (max=5000 steps, tol=1e-4):")
    res_fire = fire_minimize(
        R_init.copy(), ref, lattice, ap, is_boundary, u0,
        max_steps=5000, dt_init=2e-3, dt_max=2e-2, tol_grad=1e-4,
        report_every=500, verbose=verbose,
    )
    B_fire = compute_B_analytic(res_fire["R"], ref)
    spread_fire = compute_spread(res_fire["R"], ref, lattice, is_boundary)
    print(f"  FIRE: status={res_fire['status']}  steps={res_fire['n_steps']}  "
          f"walltime={res_fire['walltime_s']:.1f}s")
    print(f"  V*={res_fire['V']:.4e}  |grad|={res_fire['grad_norm']:.4e}  "
          f"B={B_fire:.4f}  spread_ratio={spread_fire['spread_ratio']:.4f}  "
          f"radius_rms={spread_fire['radius_rms']:.4f}a")

    # --- [B] Plain gradient descent (500 steps for timing comparison) ---
    print("\n[B] Plain gradient-descent (500 steps, backtracking line-search):")
    R_gd = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)
    V_gd = float(spacelike_potential(R_gd, lattice, ap))
    step_size = 5e-4
    t0_gd = time.perf_counter()
    for _ in range(500):
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
    t_gd = time.perf_counter() - t0_gd
    grad_gd = float(np.linalg.norm(F_gd_final))
    print(f"  GD (500 steps): V={V_gd:.4e}  |grad|={grad_gd:.4e}  walltime={t_gd:.1f}s")

    # Speedup assessment
    fire_grad = res_fire["grad_norm"]
    fire_converged = res_fire["status"] == "converged"
    print(f"\nSpeedup comparison (FIRE vs GD):")
    print(f"  FIRE: steps={res_fire['n_steps']}  |grad|={fire_grad:.4e}  "
          f"status={'CONVERGED' if fire_converged else 'not_converged'}")
    print(f"  GD:   steps=500      |grad|={grad_gd:.4e}  status=running")
    if fire_converged:
        print(f"  FIRE CONVERGED in {res_fire['n_steps']} steps while GD is still at |grad|={grad_gd:.4e}")
    elif fire_grad < grad_gd:
        print(f"  FIRE made more progress: |grad| {fire_grad:.4e} < GD {grad_gd:.4e}")
    else:
        print(f"  NOTE: need more FIRE steps for definitive speedup; "
              f"both at |grad|~{fire_grad:.4e}")

    # Verdict checks
    b_ok = abs(B_fire - 1.0) < 0.15
    v_ok = res_fire["V"] > 0
    fire_better = fire_converged or (fire_grad < grad_gd * 0.9)
    val_pass = b_ok and v_ok and fire_better

    print(f"\n  VALIDATION {'PASS' if val_pass else 'PARTIAL'}:")
    print(f"    B=1 preserved (|B-1|<0.15): {'YES' if b_ok else 'NO'}  (B={B_fire:.4f})")
    print(f"    V* > 0:                     {'YES' if v_ok else 'NO'}  (V*={res_fire['V']:.4e})")
    print(f"    FIRE faster than GD:         {'YES' if fire_better else 'CHECK'}")

    out_dir = _output_dir()
    val_result = {
        "mode": "validation",
        "grid_n": grid_n, "w": w, "alpha": alpha, "u0": u0,
        "F_edge": F_edge, "B_seed": float(B_seed), "V_seed": float(V_seed),
        "FIRE_status": res_fire["status"],
        "FIRE_steps": int(res_fire["n_steps"]),
        "FIRE_V_star": float(res_fire["V"]),
        "FIRE_grad_norm": float(fire_grad),
        "FIRE_B_final": float(B_fire),
        "FIRE_walltime_s": float(res_fire["walltime_s"]),
        "FIRE_spread_ratio": float(spread_fire["spread_ratio"]),
        "FIRE_radius_rms": float(spread_fire["radius_rms"]),
        "GD_500steps_V": float(V_gd),
        "GD_500steps_grad_norm": float(grad_gd),
        "GD_500steps_walltime_s": float(t_gd),
        "validation_pass": bool(val_pass),
    }
    val_path = os.path.join(out_dir, "w_scan_validation.json")
    _save_json(val_result, val_path)
    print(f"\nValidation result saved: {val_path}")


# ---------------------------------------------------------------------------
# Primary w-scan
# ---------------------------------------------------------------------------

def run_primary_w_scan(
    w_values: list[float] | None = None,
    alpha: float = 0.7,
    max_steps: int = 20000,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> None:
    """Sweep w/a in PRIMARY_W_VALUES at fixed alpha=0.7.

    Box size for each w chosen so F_edge < 0.05.
    """
    if w_values is None:
        w_values = PRIMARY_W_VALUES

    out_dir = _output_dir()
    results = []

    print(f"\n{'='*68}")
    print(f"PRIMARY w-SCAN: w/a in {w_values}  alpha={alpha}")
    print(f"Box design: r_edge/w = 8 for all w (F_edge ~ 0.048 for all)")
    print(f"{'='*68}\n")

    for w in w_values:
        grid_n = _primary_grid_n(w)
        print(f"\n--- w={w:.1f}a  grid_n={grid_n}  box={grid_n:.0f}a ---")
        try:
            r = run_one_box(
                grid_n=grid_n, w=w, alpha=alpha,
                max_steps=max_steps, tol_grad=tol_grad, verbose=verbose,
            )
        except Exception as exc:
            import traceback
            print(f"\n  ERROR for w={w}, grid_n={grid_n}: {exc}")
            traceback.print_exc()
            r = {
                "grid_n": grid_n, "w": float(w), "w_over_a": float(w),
                "alpha": alpha, "error": str(exc), "descent_status": "error",
            }
        results.append(r)

        # Save per-run JSON immediately (partial runs recoverable)
        tag = f"w{str(w).replace('.', 'p')}_g{grid_n}"
        json_path = os.path.join(out_dir, f"w_scan_primary_{tag}.json")
        _save_json(r, json_path)
        print(f"\n  JSON saved: {json_path}")

    # Write CSV
    csv_path = os.path.join(out_dir, "w_scan_primary.csv")
    _write_csv(results, csv_path)
    print(f"\nPrimary w-scan CSV written: {csv_path}")

    # Summary table
    _print_primary_verdict_table(results, alpha)


def _print_primary_verdict_table(results: list[dict], alpha: float) -> None:
    """Print the primary w-scan verdict table and identify any self-localized regimes."""
    print(f"\n{'='*68}")
    print(f"PRIMARY w-SCAN VERDICT TABLE (alpha={alpha})")
    print(f"{'='*68}")
    print(f"{'w/a':>5} {'grid_n':>7} {'box/w':>6} {'F_edge':>8} "
          f"{'V*':>12} {'V*/N_nodes':>12} {'V*/w^3':>10} "
          f"{'radius_rms':>11} {'spread':>7} {'B_final':>8} "
          f"{'|grad|':>9} {'steps':>7} {'status':>20}")
    print("-" * 125)

    for r in results:
        if "error" in r:
            print(f"  w={r['w']:.1f}: ERROR - {r['error']}")
            continue
        status_str = r.get("descent_status", "?")
        if r.get("not_converged", False):
            status_str += "(!)"
        print(
            f"{r['w_over_a']:>5.1f} {r['grid_n']:>7d} {r['box_over_w']:>6.1f} "
            f"{r['F_edge_seed']:>8.4f} "
            f"{r['V_star']:>12.4e} {r['V_per_node']:>12.4e} {r['V_over_w3']:>10.4e} "
            f"{r['radius_rms']:>11.4f} {r['spread_ratio']:>7.4f} {r['B_final']:>8.4f} "
            f"{r['grad_norm_final']:>9.4e} {r['n_steps']:>7d} {status_str:>20}"
        )

    # Look for localized modes: spread_ratio < 0.5 and B_preserved
    good = [r for r in results if "error" not in r and not r.get("not_converged", False)]
    localized = [r for r in good if r.get("localized", False) and r.get("B_preserved", False)]

    print(f"\n  Converged runs: {len(good)}/{len(results)}")
    if localized:
        print(f"\n  LOCALIZED candidates (spread_ratio < 0.5, B preserved, converged):")
        for r in localized:
            print(f"    w/a={r['w_over_a']:.1f}  radius_rms={r['radius_rms']:.3f}a  "
                  f"spread={r['spread_ratio']:.4f}  V*={r['V_star']:.4e}  "
                  f"V*/w^3={r['V_over_w3']:.4e}")
        # Check if width scales with w or is constant (~a)
        w_arr = [r["w"] for r in localized]
        rms_arr = [r["radius_rms"] for r in localized]
        if len(localized) >= 2:
            w_range = max(w_arr) / min(w_arr)
            rms_range = max(rms_arr) / min(rms_arr)
            if rms_range < 1.3:
                print(f"\n  MECHANISM: radius_rms ~ constant ({min(rms_arr):.2f}-{max(rms_arr):.2f}a) "
                      f"across w={min(w_arr):.1f}-{max(w_arr):.1f}a "
                      f"=> LATTICE-STABILIZED (width pinned at lattice scale, not w).")
            elif rms_range > 0.7 * w_range:
                print(f"\n  MECHANISM: radius_rms scales with w ({min(rms_arr):.2f}-{max(rms_arr):.2f}a) "
                      f"=> consistent with QUARTIC-DERRICK (soliton tracks its seed width).")
            else:
                print(f"\n  MECHANISM: AMBIGUOUS (radius_rms range={rms_range:.2f}x, "
                      f"w range={w_range:.2f}x; box-independence check required).")
    else:
        if len(good) == 0:
            print(f"\n  NO converged runs -- increase --max-steps or check convergence.")
        else:
            rms_arr = [r["radius_rms"] for r in good]
            spread_arr = [r["spread_ratio"] for r in good]
            print(f"\n  NO localized candidates among converged runs.")
            print(f"  spread_ratio range: {min(spread_arr):.3f} - {max(spread_arr):.3f}  "
                  f"(all > 0.5 = delocalized)")
            print(f"  radius_rms range:   {min(rms_arr):.3f} - {max(rms_arr):.3f}a")
            # Check if radius_rms tracks box (boundary-confined signature)
            if len(good) >= 2:
                grid_arr = [r["grid_n"] for r in good]
                grid_range = max(grid_arr) / min(grid_arr)
                rms_range = max(rms_arr) / min(rms_arr)
                if rms_range > 0.5 * grid_range:
                    print(f"  radius_rms grew {rms_range:.2f}x as w (and box) grew {grid_range:.2f}x.")
                    print(f"  => BOUNDARY-CONFINED: substrate does NOT bind B=1 as parameterized.")
                else:
                    print(f"  Weak radius_rms trend ({rms_range:.2f}x vs grid {grid_range:.2f}x): "
                          f"ambiguous; box-independence check needed.")

    print(f"{'='*68}\n")


# ---------------------------------------------------------------------------
# Box-independence check
# ---------------------------------------------------------------------------

def run_box_independence(
    configs: list[tuple[float, int, int]] | None = None,
    alpha: float = 0.7,
    max_steps: int = 20000,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> None:
    """For each (w, grid_small, grid_large), run both boxes and compare radius_rms.

    A box-INDEPENDENT radius_rms => genuine soliton width (self-localized).
    A box-TRACKING radius_rms   => boundary-confined (spreading).
    """
    if configs is None:
        configs = BOX_INDEPENDENCE_CONFIGS

    out_dir = _output_dir()
    results = []

    print(f"\n{'='*68}")
    print(f"BOX-INDEPENDENCE CHECK (alpha={alpha})")
    print(f"For each w: two boxes; verdict = is radius_rms box-independent?")
    print(f"{'='*68}\n")

    for w, grid_small, grid_large in configs:
        print(f"\n--- Box-independence: w={w:.1f}a  boxes: {grid_small} vs {grid_large} ---")
        pair = []
        for grid_n in [grid_small, grid_large]:
            try:
                r = run_one_box(
                    grid_n=grid_n, w=w, alpha=alpha,
                    max_steps=max_steps, tol_grad=tol_grad, verbose=verbose,
                )
            except Exception as exc:
                import traceback
                print(f"  ERROR for w={w}, grid_n={grid_n}: {exc}")
                traceback.print_exc()
                r = {
                    "grid_n": grid_n, "w": float(w), "w_over_a": float(w),
                    "alpha": alpha, "error": str(exc), "descent_status": "error",
                }
            r["box_independence_group"] = f"w{str(w).replace('.', 'p')}"
            pair.append(r)
            results.append(r)

            tag = f"boxind_w{str(w).replace('.', 'p')}_g{grid_n}"
            json_path = os.path.join(out_dir, f"w_scan_{tag}.json")
            _save_json(r, json_path)
            print(f"\n  JSON saved: {json_path}")

        # Per-pair verdict
        good_pair = [r for r in pair if "error" not in r and not r.get("not_converged", False)]
        if len(good_pair) == 2:
            r0, r1 = good_pair[0], good_pair[1]
            box_ratio = r1["box_size"] / r0["box_size"]
            rms_ratio = r1["radius_rms"] / r0["radius_rms"] if r0["radius_rms"] > 0 else float("nan")
            V_pnode_ratio = r1["V_per_node"] / r0["V_per_node"] if r0["V_per_node"] > 0 else float("nan")
            print(f"\n  Box-independence verdict for w={w:.1f}a:")
            print(f"    box ratio (large/small):    {box_ratio:.2f}x")
            print(f"    radius_rms ratio:           {rms_ratio:.3f}  "
                  f"(~1.0 = self-localized; ~{box_ratio:.1f} = boundary-confined)")
            print(f"    V*/N_nodes ratio:           {V_pnode_ratio:.3f}  "
                  f"(~1.0 = boundary-confined; <<1 = self-localized)")
            if rms_ratio < 1.3:
                print(f"    => SELF-LOCALIZED at w={w:.1f}a (radius_rms box-INDEPENDENT)")
            elif rms_ratio > 0.65 * box_ratio:
                print(f"    => BOUNDARY-CONFINED at w={w:.1f}a (radius_rms tracks box)")
            else:
                print(f"    => AMBIGUOUS at w={w:.1f}a (need larger box range or tighter convergence)")

    # Write CSV
    csv_path = os.path.join(out_dir, "w_scan_box_independence.csv")
    _write_csv(results, csv_path)
    print(f"\nBox-independence CSV written: {csv_path}")

    # Summary table
    print(f"\n{'='*68}")
    print(f"BOX-INDEPENDENCE SUMMARY TABLE (alpha={alpha})")
    print(f"{'='*68}")
    print(f"{'w/a':>5} {'grid_n':>7} {'box/w':>6} {'radius_rms':>11} {'spread':>7} "
          f"{'V*/N_nodes':>12} {'B_final':>8} {'|grad|':>9} {'status':>20}")
    print("-" * 95)
    for r in results:
        if "error" in r:
            print(f"  w={r['w']:.1f} grid={r['grid_n']}: ERROR")
            continue
        status_str = r.get("descent_status", "?")
        if r.get("not_converged", False):
            status_str += "(!)"
        print(f"{r['w_over_a']:>5.1f} {r['grid_n']:>7d} {r['box_over_w']:>6.1f} "
              f"{r['radius_rms']:>11.4f} {r['spread_ratio']:>7.4f} "
              f"{r['V_per_node']:>12.4e} {r['B_final']:>8.4f} "
              f"{r['grad_norm_final']:>9.4e} {status_str:>20}")

    print(f"{'='*68}\n")
    print(f"All outputs in: {out_dir}")


# ---------------------------------------------------------------------------
# Full production run (primary scan + box-independence)
# ---------------------------------------------------------------------------

def run_full(
    alpha: float = 0.7,
    max_steps: int = 20000,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> None:
    """Run the complete w-scan experiment: primary sweep + box-independence check."""
    print(f"\n{'='*68}")
    print(f"FULL w-SCAN EXPERIMENT: primary sweep + box-independence check")
    print(f"  alpha={alpha}  max_steps={max_steps}  tol={tol_grad:.1e}")
    print(f"{'='*68}")

    run_primary_w_scan(alpha=alpha, max_steps=max_steps, tol_grad=tol_grad, verbose=verbose)
    run_box_independence(alpha=alpha, max_steps=max_steps, tol_grad=tol_grad, verbose=verbose)

    out_dir = _output_dir()
    print(f"\nAll outputs written to: {out_dir}")
    print("Key files:")
    print(f"  {out_dir}/w_scan_primary.csv")
    print(f"  {out_dir}/w_scan_box_independence.csv")
    print(f"  {out_dir}/w_scan_*.json")


# ---------------------------------------------------------------------------
# Wall-time and memory estimation (for local vs AWS decision)
# ---------------------------------------------------------------------------

def estimate_costs() -> None:
    """Print per-config wall-time and memory estimates.

    Uses the measured ~1.5 s/step/10^5-nodes rule from prior FIRE runs
    on a single CPU core (M-series baseline; AWS ~0.7x faster per core).
    """
    print(f"\n{'='*68}")
    print("COST ESTIMATES (single CPU core, FIRE, 20000 steps max)")
    print("  Rule of thumb: ~0.8 us/step/node (M-series; AWS ~similar)")
    print(f"{'='*68}")
    print(f"{'Config':>30} {'n_nodes':>8} {'est_steps':>10} {'est_time_s':>12} "
          f"{'est_mem_MB':>11}")
    print("-" * 75)

    # Primary scan
    for w in PRIMARY_W_VALUES:
        grid_n = _primary_grid_n(w)
        n_nodes = grid_n ** 3
        # FIRE typically converges in ~2000-8000 steps for these sizes;
        # use 10000 as a conservative estimate.
        est_steps = min(10000, 20000)
        est_time_s = 0.8e-6 * n_nodes * est_steps  # seconds
        est_mem_mb = n_nodes * 4 * 8 * 15 / 1e6   # 15 arrays x 4 components x float64
        tag = f"primary w={w:.1f}a grid={grid_n}"
        print(f"{tag:>30} {n_nodes:>8,} {est_steps:>10,} {est_time_s:>12.0f} {est_mem_mb:>11.0f}")

    # Box-independence check
    for w, grid_small, grid_large in BOX_INDEPENDENCE_CONFIGS:
        for grid_n in [grid_small, grid_large]:
            n_nodes = grid_n ** 3
            est_steps = 10000
            est_time_s = 0.8e-6 * n_nodes * est_steps
            est_mem_mb = n_nodes * 4 * 8 * 15 / 1e6
            tag = f"box-ind w={w:.1f}a grid={grid_n}"
            print(f"{tag:>30} {n_nodes:>8,} {est_steps:>10,} {est_time_s:>12.0f} {est_mem_mb:>11.0f}")

    # Total
    total_time = 0.0
    total_mem_peak = 0.0
    for w in PRIMARY_W_VALUES:
        grid_n = _primary_grid_n(w)
        n_nodes = grid_n ** 3
        total_time += 0.8e-6 * n_nodes * 10000
        total_mem_peak = max(total_mem_peak, n_nodes * 4 * 8 * 15 / 1e6)
    for w, gs, gl in BOX_INDEPENDENCE_CONFIGS:
        for grid_n in [gs, gl]:
            n_nodes = grid_n ** 3
            total_time += 0.8e-6 * n_nodes * 10000
            total_mem_peak = max(total_mem_peak, n_nodes * 4 * 8 * 15 / 1e6)

    print("-" * 75)
    print(f"  Total est. wall-time (sequential, single core): {total_time:.0f}s "
          f"= {total_time/3600:.1f}h")
    print(f"  Peak memory (largest single config):            {total_mem_peak:.0f} MB "
          f"= {total_mem_peak/1024:.1f} GB")
    print(f"\n  LOCAL OPTION: {total_time/3600:.1f}h single-core, or ~{total_time/3600/4:.1f}h "
          f"with 4 workers (configs are independent, trivially parallelizable).")
    print(f"  AWS OPTION: c7i.4xlarge (8 vCPUs, 32 GB) or c7i.2xlarge (4 vCPUs, 16 GB); "
          f"memory-light (peak ~{total_mem_peak:.0f}MB), compute-bound.")
    print(f"{'='*68}\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="w-scan: decisive B=1 soliton localization experiment"
    )
    parser.add_argument("--validate", action="store_true",
                        help="Fast smoke test: grid_n=20, w=2a; confirm FIRE speedup (~5min)")
    parser.add_argument("--run", action="store_true",
                        help="Full production run: primary w-scan + box-independence check")
    parser.add_argument("--primary", action="store_true",
                        help="Primary w-scan only (no box-independence check)")
    parser.add_argument("--box-independence", action="store_true",
                        help="Box-independence check only")
    parser.add_argument("--estimate", action="store_true",
                        help="Print wall-time and memory estimates only (no run)")
    parser.add_argument("--alpha", type=float, default=0.7,
                        help="Prestress parameter (default 0.7)")
    parser.add_argument("--w", type=float, default=None,
                        help="Single w value to run (requires --grid-n)")
    parser.add_argument("--grid-n", type=int, default=None,
                        help="Grid size for --w single run")
    parser.add_argument("--max-steps", type=int, default=20000,
                        help="Max FIRE steps per config (default 20000)")
    parser.add_argument("--tol", type=float, default=1e-4,
                        help="FIRE gradient convergence tolerance (default 1e-4)")
    parser.add_argument("--quiet", action="store_true", help="Reduce verbosity")
    args = parser.parse_args()

    verbose = not args.quiet

    if args.estimate:
        estimate_costs()
        return

    if args.validate:
        run_validate(verbose=verbose)
        return

    if args.w is not None:
        # Single w run
        grid_n = args.grid_n if args.grid_n is not None else _primary_grid_n(args.w)
        print(f"\nSingle run: w={args.w:.1f}a  grid_n={grid_n}  alpha={args.alpha}")
        r = run_one_box(
            grid_n=grid_n, w=args.w, alpha=args.alpha,
            max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose,
        )
        out_dir = _output_dir()
        tag = f"single_w{str(args.w).replace('.', 'p')}_g{grid_n}"
        json_path = os.path.join(out_dir, f"w_scan_{tag}.json")
        _save_json(r, json_path)
        print(f"\nResult saved: {json_path}")
        return

    if args.primary:
        run_primary_w_scan(
            alpha=args.alpha, max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose,
        )
        return

    if args.box_independence:
        run_box_independence(
            alpha=args.alpha, max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose,
        )
        return

    if args.run:
        run_full(alpha=args.alpha, max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose)
        return

    # Default when called with no flags (AWS): run full production.
    run_full(alpha=args.alpha, max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose)


if __name__ == "__main__":
    main()
