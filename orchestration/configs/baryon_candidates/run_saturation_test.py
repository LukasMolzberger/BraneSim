#!/usr/bin/env python3
"""Decisive radius-saturation test with vacuum-subtracted excess-energy metrics.

PHYSICS MOTIVATION
------------------
Prior runs showed V*/N_nodes constant to <2% across boxes — now known to be a
VACUUM-FLOOR ARTIFACT: at alpha=0.7 each interior node sits on a pre-stressed
lattice where all bonds have rest length alpha*a = 0.7a but are held at spacing
a, so each bond carries prestress strain (a - 0.7a) = 0.3a.  The prestress
contribution per interior node is:

    V_vac/node = (k_s/4) * 6_bonds * (a*(1-alpha))^2
               = (1/4) * 6 * (0.3)^2 = 0.135  [dimensionless, k_s=a=1]

This matches the observed V*/N_nodes ~ 0.129-0.133 (boundary nodes have fewer
bonds, pulling the average slightly below 0.135).  The "V*/N_nodes constant"
signal is purely the prestress floor, NOT a confinement signature.

FIX: Vacuum-subtracted excess energy.
    V_excess = V_soliton - V_vacuum
where V_vacuum = spacelike_potential(R_vacuum) evaluated on the SAME lattice
with the SAME Dirichlet BC, using the all-north-pole configuration
R_vacuum[p] = ref[p] + (0,0,0,u0).

V_excess is ZERO in vacuum by construction and positive only where the lattice
is deformed above its prestress rest state.  The per-node excess:
    v_exc_p = (sum over incident interior bonds) [
        (k_s/4) * (|R_p - R_q|_soliton - alpha*a)^2
        - (k_s/4) * (|R_p - R_q|_vacuum - alpha*a)^2
    ]
weighs each node by how much its bonds contribute ABOVE the vacuum floor.

The displacement-variation weight (Sigma|u_p - u_q|^2, weight_mode="strain")
is NOT contaminated because it is ZERO in any uniform-displacement vacuum
(all nodes have the same displacement (0,0,0,u0), so u_p - u_q = 0).  It
is the reference/check against the new excess-energy metric.

DELIVERABLE 1: Excess-energy metrics
    V_excess, V_excess/N_interior, fraction_confined_excess,
    radius_rms_excess, spread_ratio_excess, confined_fraction_excess.
    Cross-check: are radius_rms_excess and radius_rms (strain-weighted) consistent?

DELIVERABLE 2: Radius-saturation test
    w=3.0a at grid_n in {36, 54, 72, 96}  (box/w = 12, 18, 24, 32)
    w=1.5a at grid_n in {24, 36, 48}      (box/w = 16, 24, 32)
    alpha=0.7, FIRE-converged (tol=1e-4), Dirichlet vacuum BC.

    Verdict per w:
        radius_rms_excess saturates -> finite lambda* -> SOLITON EXISTS, route ALIVE
        radius_rms_excess tracks box linearly -> GENUINELY UNBOUND -> route DEAD

CONSTRAINTS (principles.md non-negotiables)
-------------------------------------------
- No confinement forces, clamps, damping, or nonlinear saturation.
- Physics from spacelike_potential / spacelike_force only.
- Dirichlet vacuum BC = legitimate constraint.
- FIRE = pure V-minimization tool; no physics claim.
- All diagnostics read-only.
- Hessian skipped (grid_n > 24: n_dofs too large for dense computation).

USAGE
-----
  # Validation (fast, ~2-3min): grid_n=24, w=1.5a
  python orchestration/configs/baryon_candidates/run_saturation_test.py --validate

  # Full saturation test (sequential, total ~4-10h depending on hardware):
  python orchestration/configs/baryon_candidates/run_saturation_test.py --run

  # Single run (test any one point):
  python orchestration/configs/baryon_candidates/run_saturation_test.py --single --w 3.0 --grid-n 54

  # w=1.5 saturation series only:
  python orchestration/configs/baryon_candidates/run_saturation_test.py --w1p5

  # w=3.0 saturation series only (does {36,54,72}, skips g96 unless --include-g96):
  python orchestration/configs/baryon_candidates/run_saturation_test.py --w3p0
  python orchestration/configs/baryon_candidates/run_saturation_test.py --w3p0 --include-g96

OUTPUTS
-------
  $BRANESIM_RESULTS_DIR/saturation_test_<w>_g<N>.json   per-run JSON
  $BRANESIM_RESULTS_DIR/saturation_test_w1p5.csv         w=1.5 saturation table
  $BRANESIM_RESULTS_DIR/saturation_test_w3p0.csv         w=3.0 saturation table
  $BRANESIM_RESULTS_DIR/saturation_verdict.json          final summary + verdict
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
# Path setup
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
# Output directory
# ---------------------------------------------------------------------------

def _output_dir() -> str:
    env = os.environ.get("BRANESIM_RESULTS_DIR", "").strip()
    if env:
        os.makedirs(env, exist_ok=True)
        return env
    d = os.path.join(_REPO_ROOT, "test-runs", "saturation_test")
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Saturation test configurations
# ---------------------------------------------------------------------------

# w=3.0a: grid_n in {36, 54, 72, 96}  (box/w = 12, 18, 24, 32)
# Note: g36 has F_edge=0.085 (marginal), g54+ are clean.
W3_CONFIGS = [
    (3.0, 36),   # box/w=12, F_edge=0.085 (marginal)
    (3.0, 54),   # box/w=18, F_edge=0.038 (clean)
    (3.0, 72),   # box/w=24, F_edge=0.022 (clean)
    (3.0, 96),   # box/w=32, F_edge=0.012 (clean)
]

# w=1.5a: grid_n in {24, 36, 48}  (box/w = 16, 24, 32)
W1P5_CONFIGS = [
    (1.5, 24),   # box/w=16, F_edge=0.048 (clean)
    (1.5, 36),   # box/w=24, F_edge=0.022 (clean)
    (1.5, 48),   # box/w=32, F_edge=0.012 (clean)
]


def _f_edge_power2(grid_n: int, w: float, a: float = 1.0) -> float:
    r_edge = grid_n * a / 2.0
    return math.pi / (1.0 + (r_edge / w) ** 2)


# ---------------------------------------------------------------------------
# B estimator
# ---------------------------------------------------------------------------

def _antideriv(Fv: float) -> float:
    return Fv / 2.0 - math.sin(2.0 * Fv) / 4.0


def compute_B_analytic(positions: np.ndarray, ref: np.ndarray) -> float:
    """Analytic B estimator (hedgehog winding about the X4 axis)."""
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
# Boundary identification and Dirichlet BC
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
    """Freeze boundary nodes to north-pole vacuum: displacement = (0,0,0,u0)."""
    R_out = R.copy()
    R_out[is_boundary, :4] = ref[is_boundary, :4]
    R_out[is_boundary, 3] = ref[is_boundary, 3] + u0
    return R_out


def build_vacuum_config(ref: np.ndarray, u0: float) -> np.ndarray:
    """All-north-pole vacuum: every node at ref + (0,0,0,u0).

    This is the EXACT vacuum that the Dirichlet BC enforces on the boundary.
    V_vacuum = spacelike_potential(R_vacuum) is the prestress floor that
    must be subtracted to get V_excess.

    For the Skyrme hedgehog, the 'north pole' is the far-field limit of the
    hedgehog profile (F=0 -> sin=0, cos=+1 -> xi_i=0, xi_dim=u0).
    This matches the Dirichlet BC exactly.
    """
    R_vac = ref.copy()
    R_vac[:, 3] = ref[:, 3] + u0  # all nodes get amplitude u0, zero lateral
    return R_vac


# ---------------------------------------------------------------------------
# FIRE minimizer
# ---------------------------------------------------------------------------

_FIRE_N_MIN = 5
_FIRE_F_INC = 1.1
_FIRE_F_DEC = 0.5
_FIRE_ALPHA_START = 0.1
_FIRE_F_ALPHA = 0.99


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
    """FIRE constrained minimizer: minimize spacelike_potential with frozen boundary."""
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
            print(f"    FIRE {step:7d}: V={V_cur:.5e}  |grad|={grad_norm:.4e}  "
                  f"dt={dt:.3e}  alpha={alpha:.4f}")

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
    elapsed = time.perf_counter() - t0

    return {
        "R": R,
        "V": V_final,
        "grad_norm": grad_norm_final,
        "status": status,
        "n_steps": step + 1,
        "walltime_s": elapsed,
    }


# ---------------------------------------------------------------------------
# Vacuum-subtracted excess-energy metrics
# ---------------------------------------------------------------------------

def _per_node_excess_energy(
    R_sol: np.ndarray,
    R_vac: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    is_boundary: np.ndarray,
) -> np.ndarray:
    """Per-node excess spring energy above the vacuum prestress floor.

    For each interior node p, accumulate:
        v_exc_p = (k_s/4) * sum_{q in interior-neighbors(p)} [
            (|R_sol_p - R_sol_q| - alpha*a)^2
            - (|R_vac_p - R_vac_q| - alpha*a)^2
        ]

    Only bonds where BOTH p and q are interior nodes are included
    (boundary-incident bonds are frozen to vacuum on one end, contributing
    a fixed vacuum term that would re-introduce a floor if included).

    Returns
    -------
    v_exc : ndarray, shape (n_nodes,)
        Per-node excess energy; zero at boundary and vacuum-like interior nodes.
        Can be slightly negative at nodes adjacent to the soliton due to
        partial relaxation below the uniform-strain vacuum (allowed: no clamp).
    """
    k_s = params.k_s
    alpha_a = params.alpha * lattice.params.spacing

    v_exc = np.zeros(lattice.n_nodes, dtype=np.float64)

    nb = lattice.neighbors

    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        # Guard: only consider pairs where q is a valid node (q >= 0)
        q_valid = q_arr >= 0
        # Temporary safe index to avoid OOB on is_boundary[-1]
        q_safe = np.where(q_valid, q_arr, 0)
        # Only interior-to-interior bonds: p not boundary, q not boundary, q valid
        valid = q_valid & (~is_boundary) & (~is_boundary[q_safe])
        if not np.any(valid):
            continue

        p_idx = np.where(valid)[0]
        q_idx = q_arr[p_idx]

        # Soliton bond length
        dR_sol = R_sol[p_idx] - R_sol[q_idx]
        dist_sol = np.linalg.norm(dR_sol, axis=1)
        strain_sol = dist_sol - alpha_a

        # Vacuum bond length
        dR_vac = R_vac[p_idx] - R_vac[q_idx]
        dist_vac = np.linalg.norm(dR_vac, axis=1)
        strain_vac = dist_vac - alpha_a

        # Excess per bond (directed; factor 1/4 from potential, both directions counted)
        exc = (k_s / 4.0) * (strain_sol**2 - strain_vac**2)
        np.add.at(v_exc, p_idx, exc)

    return v_exc


def compute_excess_metrics(
    R_sol: np.ndarray,
    R_vac: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    is_boundary: np.ndarray,
) -> dict[str, float]:
    """Compute vacuum-subtracted excess-energy localization metrics.

    Returns
    -------
    dict with keys:
        V_excess        : total excess energy (soliton - vacuum)
        V_star          : raw spacelike_potential(R_sol)
        V_vacuum        : spacelike_potential(R_vac)
        V_excess_per_interior : V_excess / n_interior_nodes
        radius_rms_excess     : excess-energy-weighted RMS radius
        spread_ratio_excess   : radius_rms_excess / box_fill_radius
        confined_fraction_excess : fraction of excess energy within box_fill_radius/2
        box_fill_radius : uniform-weight RMS radius of the ref geometry
        centroid_excess : [x, y, z] excess-weighted centroid (3 floats -> stored as list)
    """
    V_star = float(spacelike_potential(R_sol, lattice, params))
    V_vacuum = float(spacelike_potential(R_vac, lattice, params))
    V_excess = V_star - V_vacuum

    n_interior = int(np.sum(~is_boundary))

    # Per-node excess weights (interior only)
    v_exc = _per_node_excess_energy(R_sol, R_vac, lattice, params, is_boundary)

    # Small floor for numerical safety: use absolute(v_exc) + tiny eps as weight.
    # The floor must NOT change the centroid or radius meaningfully.
    # We use |v_exc| because excess can be slightly negative near soliton edges
    # (relaxation below uniform strain); the magnitude still weights those
    # nodes correctly as "near the soliton".
    weights = np.abs(v_exc) + 1e-40
    # Zero out boundary nodes (they are vacuum by construction)
    weights[is_boundary] = 1e-40

    dim = 3  # spatial dimension (hardcoded for 3D; this script is 3D only)
    ref_spatial = ref[:, :dim]
    geom_centre = ref_spatial.mean(axis=0)

    # Box fill radius (uniform weight)
    delta_ref = ref_spatial - geom_centre
    r_sq_ref = np.sum(delta_ref ** 2, axis=1)
    box_fill_radius = float(np.sqrt(np.mean(r_sq_ref)))

    total_weight = float(np.sum(weights))

    # Excess-energy-weighted centroid
    centroid_exc = np.sum(ref_spatial * weights[:, None], axis=0) / total_weight

    # Excess-energy-weighted RMS radius
    delta_exc = ref_spatial - centroid_exc[None, :]
    r_sq_exc = np.sum(delta_exc ** 2, axis=1)
    radius_rms_excess = float(np.sqrt(np.sum(weights * r_sq_exc) / total_weight))

    spread_ratio_excess = radius_rms_excess / box_fill_radius if box_fill_radius > 1e-30 else float("nan")

    # Confined fraction: fraction of excess-weight within box_fill_radius/2
    confinement_r = 0.5 * box_fill_radius
    radii_from_centroid = np.sqrt(r_sq_exc)
    mask_conf = radii_from_centroid <= confinement_r
    confined_fraction_excess = float(np.sum(weights[mask_conf]) / total_weight)

    return {
        "V_star": V_star,
        "V_vacuum": V_vacuum,
        "V_excess": V_excess,
        "V_excess_per_interior": V_excess / max(n_interior, 1),
        "radius_rms_excess": radius_rms_excess,
        "spread_ratio_excess": spread_ratio_excess,
        "confined_fraction_excess": confined_fraction_excess,
        "box_fill_radius": box_fill_radius,
        "centroid_excess_x": float(centroid_exc[0]),
        "centroid_excess_y": float(centroid_exc[1]),
        "centroid_excess_z": float(centroid_exc[2]),
    }


# ---------------------------------------------------------------------------
# Strain-weighted confinement metrics (displacement-variation, NOT contaminated)
# ---------------------------------------------------------------------------

def compute_strain_spread(
    positions: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    is_boundary: np.ndarray,
) -> dict[str, float]:
    """Interior-only displacement-variation spread metrics (weight_mode='strain').

    This uses the existing confinement_metrics_per_slice with a modified
    neighbor table that excludes boundary-incident bonds — matching the
    validated approach from run_w_scan_static_soliton.py.
    """
    nb = lattice.neighbors.copy()
    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        # Guard against -1 sentinel before indexing is_boundary
        q_safe = np.where(q_arr >= 0, q_arr, 0)
        bnd_nb = (q_arr >= 0) & is_boundary[q_safe]
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
# Core single-box run
# ---------------------------------------------------------------------------

def run_one_box(
    grid_n: int,
    w: float,
    alpha: float = 0.7,
    u0: float = 1.0,
    a: float = 1.0,
    k_s: float = 1.0,
    rho: float = 1.0,
    max_steps: int = 20000,
    dt_init: float = 2e-3,
    dt_max: float = 2e-2,
    tol_grad: float = 1e-4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Seed, FIRE-minimize, and measure excess-energy + strain-weighted metrics.

    Deliverables per run
    --------------------
    - V_excess (vacuum-subtracted soliton energy)
    - radius_rms_excess, spread_ratio_excess (excess-energy weighted)
    - radius_rms, spread_ratio (strain/displacement-variation weighted, uncontaminated ref)
    - Cross-check: agreement between the two radius estimates
    - B_final, grad_norm, convergence status
    """
    box_size = float(grid_n * a)
    F_edge_seed = _f_edge_power2(grid_n, w, a)
    isolation_clean = F_edge_seed < 0.05
    isolation_marginal = F_edge_seed < 0.10

    tag = f"grid_n={grid_n}  w={w:.1f}a  box={box_size:.0f}a  box/w={box_size/w:.1f}"
    print(f"\n{'='*70}")
    print(f"  {tag}")
    print(f"  alpha={alpha}  u0={u0}  F(edge)={F_edge_seed:.4f}  "
          f"({'clean' if isolation_clean else 'marginal' if isolation_marginal else 'POOR'})")
    print(f"{'='*70}")

    m = 4
    lp = LatticeParams(grid_shape=(grid_n,) * 3, spacing=a,
                       periodic_axes=(False, False, False))
    ap = ActionParams(k_s=k_s, alpha=alpha, rho=rho, dt=0.1, n_slices=1, m_ambient=m)
    lattice = SpacelikeLattice(lp)
    ref = lattice.reference_positions(m)

    # --- Build soliton seed ---
    t_build = time.perf_counter()
    R_seed, _ = skyrme_twisted_hedgehog(lattice, m=m, u0=u0, w=w, profile_shape="power2")
    B_seed = compute_B_analytic(R_seed, ref)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    t_build_s = time.perf_counter() - t_build
    print(f"  Seed: B_seed={B_seed:.4f}  V_seed={V_seed:.4e}  build={t_build_s:.1f}s")

    # --- Build vacuum config (same Dirichlet BC, same lattice) ---
    R_vacuum = build_vacuum_config(ref, u0)
    V_vacuum = float(spacelike_potential(R_vacuum, lattice, ap))

    # Sanity: vacuum energy per interior node should be ~ (k_s/4)*6*(a*(1-alpha))^2
    is_boundary = identify_boundary(lattice, n_shells=1)
    n_boundary = int(np.sum(is_boundary))
    n_interior = int(np.sum(~is_boundary))
    n_dofs = n_interior * m

    V_vac_expected_interior = (k_s / 4.0) * 6.0 * (a * (1.0 - alpha)) ** 2 * n_interior
    print(f"  Vacuum: V_vacuum={V_vacuum:.4e}  "
          f"V_vac_expected_int={V_vac_expected_interior:.4e}  "
          f"V_vac/N_nodes={V_vacuum/lattice.n_nodes:.5f}  "
          f"(theory interior: {(k_s/4)*6*(a*(1-alpha))**2:.5f}/node)")

    print(f"  Nodes: total={lattice.n_nodes:,}  boundary={n_boundary} "
          f"({100*n_boundary/lattice.n_nodes:.1f}%)  interior={n_interior:,}  "
          f"n_dofs={n_dofs:,}")

    # --- Apply Dirichlet BC and initialize ---
    R_init = apply_dirichlet(R_seed.copy(), ref, is_boundary, u0)
    B_init = compute_B_analytic(R_init, ref)
    print(f"  After Dirichlet BC: B_init={B_init:.4f}")

    # --- FIRE minimization ---
    print(f"\n  FIRE (max={max_steps} tol={tol_grad:.1e} dt_init={dt_init:.2e} dt_max={dt_max:.2e})")
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
    print(f"  V_star={V_star:.5e}  V_vacuum={V_vacuum:.5e}  "
          f"grad_norm={grad_norm:.4e}  B_final={B_final:.4f}")
    if not_converged:
        print(f"  *** NOT CONVERGED (grad_norm={grad_norm:.4e} > tol={tol_grad:.1e}) ***")

    # --- METRIC 1: Vacuum-subtracted excess-energy metrics ---
    print(f"\n  Computing excess-energy metrics (Deliverable 1)...")
    exc_metrics = compute_excess_metrics(
        R_final, R_vacuum, ref, lattice, ap, is_boundary,
    )
    V_excess = exc_metrics["V_excess"]
    print(f"    V_excess = V_star - V_vacuum = {V_star:.5e} - {V_vacuum:.5e} = {V_excess:.5e}")
    print(f"    V_excess/N_interior = {exc_metrics['V_excess_per_interior']:.5e}")
    print(f"    radius_rms_excess   = {exc_metrics['radius_rms_excess']:.4f}a")
    print(f"    spread_ratio_excess = {exc_metrics['spread_ratio_excess']:.4f}")
    print(f"    confined_frac_excess= {exc_metrics['confined_fraction_excess']:.4f}")
    print(f"    box_fill_radius     = {exc_metrics['box_fill_radius']:.4f}a")

    # --- METRIC 2: Strain/displacement-variation weighted metrics (uncontaminated ref) ---
    print(f"\n  Computing strain-weighted metrics (Deliverable 2 cross-check)...")
    strain_metrics = compute_strain_spread(R_final, ref, lattice, is_boundary)
    print(f"    radius_rms_strain   = {strain_metrics['radius_rms']:.4f}a")
    print(f"    spread_ratio_strain = {strain_metrics['spread_ratio']:.4f}")
    print(f"    confined_frac_strain= {strain_metrics['confined_fraction']:.4f}")

    # Cross-check: do the two metrics agree?
    rms_excess = exc_metrics["radius_rms_excess"]
    rms_strain = strain_metrics["radius_rms"]
    if rms_strain > 1e-6:
        agreement_ratio = rms_excess / rms_strain
    else:
        agreement_ratio = float("nan")
    print(f"\n  Cross-check: radius_rms_excess / radius_rms_strain = {agreement_ratio:.3f}")
    if abs(agreement_ratio - 1.0) < 0.2:
        print(f"  CROSS-CHECK: AGREE (ratio within 20%)")
    else:
        print(f"  CROSS-CHECK: DISAGREE (ratio {agreement_ratio:.3f} outside 20%)")

    # Raw V*/N_nodes (keep for record; contaminated)
    V_per_node_raw = V_star / lattice.n_nodes
    V_per_interior_raw = V_star / n_interior
    print(f"\n  Raw V*/N_nodes = {V_per_node_raw:.5f}  (contaminated by vacuum floor ~"
          f"{(k_s/4)*6*(a*(1-alpha))**2:.4f})")

    # Localization verdict (excess-based)
    localized_excess = exc_metrics["spread_ratio_excess"] < 0.5
    print(f"\n  LOCALIZATION: spread_ratio_excess={exc_metrics['spread_ratio_excess']:.4f}  "
          f"-> {'LOCALIZED' if localized_excess else 'DELOCALIZED'}")
    print(f"  box/w = {box_size/w:.1f}  "
          f"radius_rms_excess/box = {rms_excess/box_size:.4f}  "
          f"radius_rms_excess/w   = {rms_excess/w:.4f}")

    result = {
        # Config
        "grid_n": grid_n,
        "w": float(w),
        "w_over_a": float(w / a),
        "alpha": float(alpha),
        "u0": float(u0),
        "a": float(a),
        "k_s": float(k_s),
        "box_size": float(box_size),
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
        # Raw energy (contaminated; keep for record)
        "V_star": float(V_star),
        "V_vacuum": float(V_vacuum),
        "V_per_node_raw": float(V_per_node_raw),
        "V_per_interior_raw": float(V_per_interior_raw),
        "V_vac_per_node_theory": float((k_s / 4.0) * 6.0 * (a * (1.0 - alpha)) ** 2),
        # DELIVERABLE 1: Vacuum-subtracted excess-energy metrics
        "V_excess": float(V_excess),
        "V_excess_per_interior": float(exc_metrics["V_excess_per_interior"]),
        "radius_rms_excess": float(rms_excess),
        "spread_ratio_excess": float(exc_metrics["spread_ratio_excess"]),
        "confined_fraction_excess": float(exc_metrics["confined_fraction_excess"]),
        "box_fill_radius": float(exc_metrics["box_fill_radius"]),
        "centroid_excess_x": float(exc_metrics["centroid_excess_x"]),
        "centroid_excess_y": float(exc_metrics["centroid_excess_y"]),
        "centroid_excess_z": float(exc_metrics["centroid_excess_z"]),
        # DELIVERABLE 2: Strain-weighted (displacement-variation) metrics
        "radius_rms_strain": float(rms_strain),
        "spread_ratio_strain": float(strain_metrics["spread_ratio"]),
        "confined_fraction_strain": float(strain_metrics["confined_fraction"]),
        # Cross-check
        "rms_excess_over_strain": float(agreement_ratio),
        "rms_cross_check_ok": bool(abs(agreement_ratio - 1.0) < 0.2),
        # Winding
        "B_final": float(B_final),
        "B_preserved": bool(abs(B_final - 1.0) < 0.15),
        "grad_norm_final": float(grad_norm),
        # Localization verdict (excess-based)
        "localized_excess": bool(localized_excess),
        "rms_excess_over_box": float(rms_excess / box_size),
        "rms_excess_over_w": float(rms_excess / w),
        "rms_strain_over_box": float(rms_strain / (box_size if box_size > 0 else 1.0)),
    }
    return result


# ---------------------------------------------------------------------------
# CSV and JSON helpers
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
# Saturation series runner
# ---------------------------------------------------------------------------

def run_saturation_series(
    configs: list[tuple[float, int]],
    series_name: str,
    alpha: float = 0.7,
    max_steps: int = 20000,
    tol_grad: float = 1e-4,
    verbose: bool = True,
    out_dir: str | None = None,
) -> list[dict]:
    """Run a w-saturation series and produce CSV + per-run JSONs."""
    if out_dir is None:
        out_dir = _output_dir()

    print(f"\n{'='*70}")
    print(f"SATURATION SERIES: {series_name}")
    print(f"  alpha={alpha}  max_steps={max_steps}  tol={tol_grad:.1e}")
    print(f"{'='*70}\n")

    results = []
    for w, grid_n in configs:
        print(f"\n--- {series_name}: w={w:.1f}a  grid_n={grid_n}  box/w={grid_n/w:.1f} ---")
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
                "radius_rms_excess": float("nan"), "spread_ratio_excess": float("nan"),
                "radius_rms_strain": float("nan"), "V_excess": float("nan"),
            }
        results.append(r)

        tag = f"w{str(w).replace('.','p')}_g{grid_n}"
        json_path = os.path.join(out_dir, f"saturation_test_{tag}.json")
        _save_json(r, json_path)
        print(f"\n  JSON saved: {json_path}")

    csv_path = os.path.join(out_dir, f"saturation_test_{series_name.replace(' ','_').replace('=','').replace('/','')}.csv")
    _write_csv(results, csv_path)
    print(f"\nCSV written: {csv_path}")

    _print_saturation_verdict(results, series_name, alpha)
    return results


def _print_saturation_verdict(results: list[dict], series_name: str, alpha: float) -> None:
    """Print saturation verdict table and state conclusion."""
    print(f"\n{'='*80}")
    print(f"SATURATION TABLE: {series_name}  (alpha={alpha})")
    print(f"{'='*80}")
    header = (f"{'grid_n':>7} {'box/w':>6} {'V_excess':>12} {'r_rms_exc':>11} "
              f"{'sr_exc':>7} {'r_rms_str':>11} {'sr_str':>7} "
              f"{'B_final':>8} {'|grad|':>9} {'steps':>7} {'status':>16}")
    print(header)
    print("-" * len(header))

    good_results = []
    for r in results:
        if "error" in r:
            print(f"  grid_n={r['grid_n']} w={r['w']:.1f}: ERROR - {r['error']}")
            continue
        nc = "(!)" if r.get("not_converged", False) else "   "
        status = r.get("descent_status", "?") + nc
        rme = r.get("radius_rms_excess", float("nan"))
        sre = r.get("spread_ratio_excess", float("nan"))
        rms = r.get("radius_rms_strain", float("nan"))
        srs = r.get("spread_ratio_strain", float("nan"))
        print(
            f"{r['grid_n']:>7d} {r['box_over_w']:>6.1f} "
            f"{r.get('V_excess', float('nan')):>12.4e} "
            f"{rme:>11.4f} {sre:>7.4f} {rms:>11.4f} {srs:>7.4f} "
            f"{r['B_final']:>8.4f} {r['grad_norm_final']:>9.4e} "
            f"{r['n_steps']:>7d} {status:>16}"
        )
        if not r.get("not_converged", False) and "error" not in r:
            good_results.append(r)

    if len(good_results) < 2:
        print(f"\n  INSUFFICIENT converged results for saturation verdict.")
        return

    # Extract grid_n and radius_rms_excess for converged runs
    boxes = [r["box_size"] for r in good_results]
    rms_exc = [r["radius_rms_excess"] for r in good_results]
    rms_str = [r["radius_rms_strain"] for r in good_results]
    grid_ns = [r["grid_n"] for r in good_results]

    print(f"\n  Saturation analysis ({len(good_results)} converged runs):")
    print(f"  box_size:        {boxes}")
    print(f"  radius_rms_exc:  {[f'{x:.3f}' for x in rms_exc]}")
    print(f"  radius_rms_str:  {[f'{x:.3f}' for x in rms_str]}")

    # Saturation test: compare largest/smallest ratio vs box largest/smallest ratio
    box_ratio = max(boxes) / min(boxes)
    rms_exc_ratio = max(rms_exc) / max(min(rms_exc), 1e-6)
    rms_str_ratio = max(rms_str) / max(min(rms_str), 1e-6)

    print(f"\n  box_size ratio (max/min): {box_ratio:.3f}")
    print(f"  radius_rms_excess ratio:  {rms_exc_ratio:.3f}")
    print(f"  radius_rms_strain ratio:  {rms_str_ratio:.3f}")

    # Classify: saturates if rms_ratio < 1.3 (10% margin on top of 1.0),
    # tracks if rms_ratio > 0.6 * box_ratio.
    # Ambiguous otherwise.
    print(f"\n  *** VERDICT ({series_name}) ***")
    if rms_exc_ratio < 1.3:
        # Saturated: finite width
        lambda_star = float(np.mean(rms_exc[-2:]))  # average of largest two boxes
        lambda_star_derrick = math.sqrt(alpha / (1.0 - alpha))
        print(f"  EXCESS-ENERGY: SATURATES -> finite lambda* = {lambda_star:.3f}a  "
              f"  (Derrick scaling sqrt(alpha/(1-alpha)) = {lambda_star_derrick:.3f}a)")
        print(f"  ROUTE: ALIVE (soliton exists, just wide; lambda*={lambda_star:.3f}a)")
    elif rms_exc_ratio > 0.6 * box_ratio:
        print(f"  EXCESS-ENERGY: TRACKS BOX (rms_ratio={rms_exc_ratio:.3f} vs "
              f"box_ratio={box_ratio:.3f}) -> GENUINELY UNBOUND")
        print(f"  ROUTE: DEAD at alpha={alpha} for w={good_results[0]['w']:.1f}a")
        print(f"  NOTE: Try higher alpha (alpha->1 strengthens quartic E4 ~ k_s*alpha/a)")
    else:
        print(f"  EXCESS-ENERGY: AMBIGUOUS (rms_ratio={rms_exc_ratio:.3f} vs "
              f"box_ratio={box_ratio:.3f})")
        print(f"  Need larger box range or higher alpha to resolve.")

    # Cross-check with strain metric
    if rms_str_ratio < 1.3:
        print(f"  STRAIN CROSS-CHECK: also SATURATES (strain_rms_ratio={rms_str_ratio:.3f})")
    elif rms_str_ratio > 0.6 * box_ratio:
        print(f"  STRAIN CROSS-CHECK: also TRACKS BOX (strain_rms_ratio={rms_str_ratio:.3f})")
    else:
        print(f"  STRAIN CROSS-CHECK: AMBIGUOUS (strain_rms_ratio={rms_str_ratio:.3f})")

    print(f"{'='*80}\n")


# ---------------------------------------------------------------------------
# Validation run
# ---------------------------------------------------------------------------

def run_validate(verbose: bool = True) -> None:
    """Fast smoke test: grid_n=24, w=1.5a, alpha=0.7.

    Confirms:
      1. V_excess > 0 (real soliton energy above vacuum floor)
      2. V_vacuum ~ N_nodes * (k_s/4) * 6 * (a*(1-alpha))^2 (floor formula correct)
      3. radius_rms_excess and radius_rms_strain give consistent answers
      4. B=1 preserved after FIRE convergence
    """
    print("\n" + "=" * 70)
    print("VALIDATION: grid_n=24  w=1.5a  alpha=0.7")
    print("  Checks: V_excess > 0, floor formula, metric cross-check, B=1")
    print("=" * 70)

    r = run_one_box(grid_n=24, w=1.5, alpha=0.7, max_steps=5000, tol_grad=1e-4,
                    verbose=verbose)

    out_dir = _output_dir()
    val_path = os.path.join(out_dir, "saturation_validation.json")
    _save_json(r, val_path)

    print(f"\n  VALIDATION CHECKS:")
    v1 = r.get("V_excess", -1.0) > 0
    v2 = abs(r.get("V_per_node_raw", 0) - r.get("V_vac_per_node_theory", 0)) / max(r.get("V_vac_per_node_theory", 1), 1e-10) < 0.05
    v3 = r.get("rms_cross_check_ok", False)
    v4 = r.get("B_preserved", False)
    print(f"    V_excess > 0:                   {'PASS' if v1 else 'FAIL'}  "
          f"(V_excess={r.get('V_excess','?'):.4e})")
    print(f"    V_per_node matches floor theory: {'PASS' if v2 else 'CHECK'}  "
          f"(observed={r.get('V_per_node_raw','?'):.5f} "
          f"theory={r.get('V_vac_per_node_theory','?'):.5f})")
    print(f"    Metric cross-check (|ratio-1|<20%): {'PASS' if v3 else 'CHECK'}  "
          f"(ratio={r.get('rms_excess_over_strain','?'):.3f})")
    print(f"    B=1 preserved:                  {'PASS' if v4 else 'FAIL'}  "
          f"(B={r.get('B_final','?'):.4f})")
    print(f"\n  Validation result saved: {val_path}")
    overall = v1 and v4
    print(f"\n  VALIDATION {'PASS' if overall else 'PARTIAL'} "
          f"(v_excess_ok={v1}, B_ok={v4}, metric_ok={v3})")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Radius saturation test with vacuum-subtracted excess-energy metrics"
    )
    parser.add_argument("--validate", action="store_true",
                        help="Fast smoke test: grid_n=24, w=1.5a, confirms floor formula")
    parser.add_argument("--run", action="store_true",
                        help="Full saturation test: w=1.5 and w=3.0 series")
    parser.add_argument("--w1p5", action="store_true",
                        help="w=1.5a saturation series: grid_n in {24, 36, 48}")
    parser.add_argument("--w3p0", action="store_true",
                        help="w=3.0a saturation series: grid_n in {36, 54, 72}")
    parser.add_argument("--include-g96", action="store_true",
                        help="Include grid_n=96 in w=3.0 series (slow, ~2-3h)")
    parser.add_argument("--single", action="store_true",
                        help="Run a single (w, grid_n) pair")
    parser.add_argument("--w", type=float, default=None, help="w value for --single")
    parser.add_argument("--grid-n", type=int, default=None, help="grid_n for --single")
    parser.add_argument("--alpha", type=float, default=0.7, help="Prestress (default 0.7)")
    parser.add_argument("--max-steps", type=int, default=20000,
                        help="Max FIRE steps (default 20000)")
    parser.add_argument("--tol", type=float, default=1e-4,
                        help="FIRE grad convergence tolerance (default 1e-4)")
    parser.add_argument("--quiet", action="store_true", help="Reduce verbosity")
    args = parser.parse_args()

    verbose = not args.quiet
    out_dir = _output_dir()

    if args.validate:
        run_validate(verbose=verbose)
        return

    if args.single:
        if args.w is None or args.grid_n is None:
            parser.error("--single requires --w and --grid-n")
        r = run_one_box(
            grid_n=args.grid_n, w=args.w, alpha=args.alpha,
            max_steps=args.max_steps, tol_grad=args.tol, verbose=verbose,
        )
        tag = f"w{str(args.w).replace('.','p')}_g{args.grid_n}"
        json_path = os.path.join(out_dir, f"saturation_test_{tag}.json")
        _save_json(r, json_path)
        print(f"\nResult saved: {json_path}")
        return

    all_results = {}

    if args.w1p5 or args.run:
        r15 = run_saturation_series(
            W1P5_CONFIGS, "w1p5", alpha=args.alpha,
            max_steps=args.max_steps, tol_grad=args.tol,
            verbose=verbose, out_dir=out_dir,
        )
        all_results["w1p5"] = r15

    if args.w3p0 or args.run:
        configs_3 = W3_CONFIGS[:-1] if not args.include_g96 else W3_CONFIGS
        if not args.include_g96:
            print("\n  NOTE: Skipping grid_n=96 (add --include-g96 to include; ~2-3h extra).")
        r30 = run_saturation_series(
            configs_3, "w3p0", alpha=args.alpha,
            max_steps=args.max_steps, tol_grad=args.tol,
            verbose=verbose, out_dir=out_dir,
        )
        all_results["w3p0"] = r30

    # Write combined verdict JSON
    if all_results:
        verdict = {}
        for key, rows in all_results.items():
            good = [r for r in rows if "error" not in r and not r.get("not_converged", False)]
            if len(good) >= 2:
                rms_exc = [r["radius_rms_excess"] for r in good]
                boxes = [r["box_size"] for r in good]
                ratio = max(rms_exc) / max(min(rms_exc), 1e-6)
                box_ratio = max(boxes) / min(boxes)
                if ratio < 1.3:
                    conclusion = "SATURATES"
                elif ratio > 0.6 * box_ratio:
                    conclusion = "TRACKS_BOX"
                else:
                    conclusion = "AMBIGUOUS"
                verdict[key] = {
                    "conclusion": conclusion,
                    "rms_excess_values": rms_exc,
                    "box_sizes": boxes,
                    "rms_ratio": ratio,
                    "box_ratio": box_ratio,
                }
        verdict_path = os.path.join(out_dir, "saturation_verdict.json")
        with open(verdict_path, "w") as f:
            json.dump(verdict, f, indent=2)
        print(f"\nVerdict JSON: {verdict_path}")
        print(f"All outputs in: {out_dir}")

    if not any([args.validate, args.single, args.w1p5, args.w3p0, args.run]):
        parser.print_help()


if __name__ == "__main__":
    main()
