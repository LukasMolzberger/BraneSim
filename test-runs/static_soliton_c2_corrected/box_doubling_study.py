"""Box-doubling study for the C2 static B=1 soliton: self-localized vs boundary-confined?

Decisive test: hold (alpha=0.7, w=2a) fixed, vary ONLY grid_n in {13, 20, 26}.
Box sizes: 13a (~baseline), 20a (~1.5x), 26a (~2x = doubled.
Answers the Derrick question: does the soliton width track the box, or sit at a fixed
physical scale?

Key metrics per box:
  - B_final  (must stay ~1 throughout)
  - V*       (potential energy at converged minimum)
  - spread_ratio(strain, interior-only bonds)
  - radius_rms (ABSOLUTE, in lattice units -- the physical width proxy)
  - min|delta_R| (collapse watch)
  - N_neg + lambda_min from static Hessian (at grid_n=13 and 26; skip 20 if too large)

Verdict rule:
  SELF-LOCALIZED: V* and radius_rms are box-INDEPENDENT (flat with box size),
                  spread_ratio decreases as box grows (soliton stays small relative to box).
  BOUNDARY-CONFINED: radius_rms grows with box size, spread_ratio stays ~ const or
                     increases -- the soliton expands to fill the box.

Usage:
  python box_doubling_study.py            # runs all three boxes
  python box_doubling_study.py --quick    # skip Hessian at grid_n=26 if too slow

Outputs:
  test-runs/static_soliton_c2_corrected/box_doubling_sweep.csv
  test-runs/static_soliton_c2_corrected/box_doubling_report.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np

REPO_ROOT = "/Users/lukasmolzberger/PycharmProjects/BraneSim"
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from branesim.core.action import spacelike_potential, spacelike_force
from branesim.core.conventions import LatticeParams, ActionParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.confinement import confinement_metrics_per_slice
from branesim.initialization.seeds import skyrme_twisted_hedgehog

OUTPUT_DIR = os.path.join(REPO_ROOT, "test-runs/static_soliton_c2_corrected")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ============================================================
# Winding number (analytic radial integral -- authoritative)
# ============================================================

def compute_winding_analytic(
    positions: np.ndarray,
    ref: np.ndarray,
    w: float = 0.0,
) -> float:
    """B = -(2/pi) * [antideriv(F_edge) - antideriv(F_inner)].

    F_inner = pi (the continuum BC at r=0 for the Skyrme hedgehog).
    F_edge  = F at the outermost node (read from the field).

    This is the AUTHORITATIVE estimator per the harness spec.
    Not the Jacobian estimator (which is an even-grid artifact).
    """
    m_ambient = positions.shape[1]
    phi_4 = positions[:, :4] - ref[:, :4]

    norms = np.linalg.norm(phi_4, axis=1)
    eps = 1e-10 * (float(np.max(norms)) + 1e-30)
    valid = norms > eps

    cos_F = np.where(valid, phi_4[:, 3] / np.where(valid, norms, 1.0), 1.0)
    cos_F = np.clip(cos_F, -1.0, 1.0)
    F_nodes = np.arccos(cos_F)

    coords = ref[:, :3]
    centre = coords.mean(axis=0)
    dx = coords - centre
    r = np.linalg.norm(dx, axis=1)

    i_edge = int(np.argmax(r))
    F_edge = float(F_nodes[i_edge])

    def antideriv(Fv: float) -> float:
        return Fv / 2.0 - math.sin(2.0 * Fv) / 4.0

    F_inner = math.pi  # continuum BC: F(r=0)=pi
    return -(2.0 / math.pi) * (antideriv(F_edge) - antideriv(F_inner))


# ============================================================
# Boundary identification
# ============================================================

def identify_boundary_nodes(lattice: SpacelikeLattice, n_shells: int = 1) -> np.ndarray:
    nb = lattice.neighbors
    is_boundary = np.any(nb < 0, axis=1)
    for _ in range(n_shells - 1):
        nb_of_boundary = nb[is_boundary]
        valid_nb = nb_of_boundary[nb_of_boundary >= 0]
        is_boundary[valid_nb] = True
    return is_boundary


def apply_dirichlet_vacuum(
    R: np.ndarray,
    ref: np.ndarray,
    is_boundary: np.ndarray,
    u0: float,
) -> np.ndarray:
    R_out = R.copy()
    m_ambient = R.shape[1]
    R_out[is_boundary, :m_ambient] = ref[is_boundary, :m_ambient]
    R_out[is_boundary, 3] = ref[is_boundary, 3] + u0
    return R_out


# ============================================================
# Constrained gradient descent
# ============================================================

def constrained_gradient_descent(
    R_seed: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    is_boundary: np.ndarray,
    u0: float,
    *,
    max_steps: int = 8000,
    step_size: float = 5e-4,
    tol_grad: float = 1e-5,
    report_every: int = 2000,
    verbose: bool = True,
) -> dict[str, Any]:
    R = apply_dirichlet_vacuum(R_seed.copy(), ref, is_boundary, u0)
    a = lattice.params.spacing

    V_prev = float(spacelike_potential(R, lattice, params))
    step = step_size
    status = "max_steps_reached"
    t0 = time.perf_counter()
    min_dist_final = float("inf")

    for i in range(max_steps):
        F = spacelike_force(R, lattice, params)
        F[is_boundary] = 0.0

        grad_norm = float(np.linalg.norm(F))

        if i % report_every == 0 and verbose:
            V_cur = float(spacelike_potential(R, lattice, params))
            B_cur = compute_winding_analytic(R, ref)
            print(f"    step {i:5d}: V={V_cur:.5e}, |grad|={grad_norm:.4e}, B={B_cur:.4f}")

        if grad_norm < tol_grad:
            status = "converged"
            if verbose:
                print(f"    Converged at step {i}: |grad|={grad_norm:.4e}")
            break

        R_new = R + step * F
        V_new = float(spacelike_potential(R_new, lattice, params))

        if V_new < V_prev:
            R = R_new
            V_prev = V_new
            step = min(step * 1.05, step_size * 10)
        else:
            step = max(step * 0.5, step_size * 1e-3)
            R_try = R + step * F
            V_try = float(spacelike_potential(R_try, lattice, params))
            if V_try < V_prev:
                R = R_try
                V_prev = V_try

        R = apply_dirichlet_vacuum(R, ref, is_boundary, u0)

    # Final bond-distance check
    nb = lattice.neighbors
    min_dist = float("inf")
    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        valid = q_arr >= 0
        p_idx = np.where(valid)[0]
        q_idx = q_arr[p_idx]
        diffs = R[p_idx] - R[q_idx]
        dists = np.linalg.norm(diffs, axis=1)
        if len(dists) > 0:
            min_dist = min(min_dist, float(np.min(dists)))

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
        "n_steps": i + 1,
        "walltime_s": elapsed,
        "min_dist": min_dist,
    }


# ============================================================
# Strain-weighted spread_ratio (interior bonds only)
# ============================================================

def compute_strain_spread_ratio(
    positions: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    is_boundary: np.ndarray,
) -> dict[str, float]:
    """Interior-only strain-weighted confinement metrics.

    Excludes bonds incident to boundary nodes so the frozen vacuum shell
    doesn't inflate the spread estimate.
    """
    nb = lattice.neighbors.copy()
    # Mask boundary-incident bonds
    for nb_idx in range(nb.shape[1]):
        q_arr = nb[:, nb_idx]
        boundary_nb = (q_arr >= 0) & is_boundary[q_arr]
        nb[boundary_nb, nb_idx] = -1
    nb[is_boundary] = -1

    return confinement_metrics_per_slice(
        positions, ref,
        dim=3,
        confinement_radius_factor=0.5,
        weight_mode="strain",
        _neighbor_table=nb,
    )


# ============================================================
# Static Hessian (N_neg, lambda_min)
# ============================================================

def compute_static_hessian(
    R: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    is_boundary: np.ndarray,
    mass: float,
    fd_eps: float = 1e-5,
    max_n_dofs: int = 10000,
) -> dict[str, Any]:
    m_ambient = R.shape[1]
    n_nodes = lattice.n_nodes
    interior_nodes = np.where(~is_boundary)[0]
    n_interior = len(interior_nodes)
    n_dofs = n_interior * m_ambient

    if n_dofs > max_n_dofs:
        return {
            "skipped": True,
            "reason": f"n_dofs={n_dofs} > max_n_dofs={max_n_dofs}",
            "n_neg": None, "n_neg_genuine": None, "stable": None,
            "eig_min": None, "n_dofs": n_dofs,
        }

    shape = R.shape
    R_flat = R.ravel()

    interior_full_indices = np.array([
        interior_nodes[k // m_ambient] * m_ambient + k % m_ambient
        for k in range(n_dofs)
    ], dtype=np.int64)

    J = np.zeros((n_dofs, n_dofs), dtype=np.float64)
    t0 = time.perf_counter()

    for j in range(n_dofs):
        dRf = np.zeros(n_nodes * m_ambient)
        dRf[interior_full_indices[j]] = fd_eps
        Fp = spacelike_force((R_flat + dRf).reshape(shape), lattice, params).ravel()
        Fm = spacelike_force((R_flat - dRf).reshape(shape), lattice, params).ravel()
        dF = (Fp - Fm) / (2.0 * fd_eps)
        J[:, j] = dF[interior_full_indices]

    K = -0.5 * (J + J.T)
    eig = np.linalg.eigvalsh(K)

    elapsed = time.perf_counter() - t0

    pos_eig = eig[eig > 0]
    scale = float(np.median(pos_eig)) if len(pos_eig) > 0 else 1.0
    tol_zero = max(1e-6 * scale, 1e-8)
    tol_genuine_neg = 10 * tol_zero

    n_neg = int(np.sum(eig < -tol_zero))
    n_zero = int(np.sum(np.abs(eig) <= tol_zero))
    n_pos = int(np.sum(eig > tol_zero))
    n_neg_genuine = int(np.sum(eig < -tol_genuine_neg))

    return {
        "skipped": False,
        "n_neg": n_neg,
        "n_neg_genuine": n_neg_genuine,
        "n_zero": n_zero,
        "n_pos": n_pos,
        "n_dofs": n_dofs,
        "eig_min": float(eig[0]),
        "eig_max": float(eig[-1]),
        "eigenvalues_bottom10": eig[:10].tolist(),
        "stable": n_neg_genuine == 0,
        "hessian_walltime_s": elapsed,
        "tol_zero": tol_zero,
        "tol_genuine_neg": tol_genuine_neg,
    }


# ============================================================
# Single-box run
# ============================================================

def run_box(
    grid_n: int,
    alpha: float = 0.7,
    u0: float = 1.0,
    w: float = 2.0,
    a: float = 1.0,
    k_s: float = 1.0,
    rho: float = 1.0,
    n_boundary_shells: int = 1,
    max_steps: int = 8000,
    step_size: float = 5e-4,
    tol_grad: float = 1e-5,
    do_hessian: bool = True,
    hessian_max_dofs: int = 10000,
    verbose: bool = True,
) -> dict[str, Any]:
    box_size = grid_n * a
    box_half = box_size / 2.0
    F_at_edge_power2 = math.pi / (1.0 + (box_half / w) ** 2)

    print(f"\n{'='*65}")
    print(f"Box study: grid_n={grid_n}, box={box_size:.0f}a, box/w={box_size/w:.1f}")
    print(f"  alpha={alpha}, w={w}a, F(edge,power2)={F_at_edge_power2:.4f}")
    print(f"{'='*65}")

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
    mass = ap.mass(lp)
    ref = lattice.reference_positions(m_ambient)

    # --- Seed ---
    R_seed, _ = skyrme_twisted_hedgehog(
        lattice, m=m_ambient, u0=u0, w=w, profile_shape="power2"
    )
    B_seed = compute_winding_analytic(R_seed, ref, w=w)
    V_seed = float(spacelike_potential(R_seed, lattice, ap))
    print(f"  Seed: B_seed(analytic)={B_seed:.4f}, V_seed={V_seed:.4e}")

    # --- Boundary ---
    is_boundary = identify_boundary_nodes(lattice, n_shells=n_boundary_shells)
    n_boundary = int(np.sum(is_boundary))
    n_interior = int(np.sum(~is_boundary))
    print(f"  Boundary nodes: {n_boundary} ({100.0*n_boundary/lattice.n_nodes:.1f}%)")
    print(f"  Interior nodes: {n_interior} (n_dofs={n_interior*m_ambient})")

    R_init = apply_dirichlet_vacuum(R_seed.copy(), ref, is_boundary, u0)
    B_init = compute_winding_analytic(R_init, ref)
    print(f"  After Dirichlet BC: B_init={B_init:.4f}")

    # --- Gradient descent ---
    print(f"\n  Constrained gradient descent (max {max_steps} steps, tol={tol_grad:.0e})...")
    desc = constrained_gradient_descent(
        R_init, ref, lattice, ap, is_boundary, u0,
        max_steps=max_steps,
        step_size=step_size,
        tol_grad=tol_grad,
        report_every=max(max_steps // 4, 500),
        verbose=verbose,
    )
    R_final = desc["R"]
    V_star = desc["V"]
    B_final = compute_winding_analytic(R_final, ref)
    grad_norm_final = desc["grad_norm"]
    print(f"  Descent: status={desc['status']}, V*={V_star:.5e}, B={B_final:.4f}, "
          f"|grad|={grad_norm_final:.4e}, n_steps={desc['n_steps']}, "
          f"min_dist={desc['min_dist']:.4f}, time={desc['walltime_s']:.1f}s")

    # Check convergence warning
    not_converged = desc["status"] != "converged"
    if not_converged:
        print(f"  WARNING: NOT CONVERGED (grad_norm={grad_norm_final:.4e}). "
              f"V* and spread_ratio may be unreliable.")

    # --- Strain spread_ratio (interior bonds only) ---
    spread = compute_strain_spread_ratio(R_final, ref, lattice, is_boundary)
    spread_ratio = spread["spread_ratio"]
    radius_rms = spread["radius_rms"]       # absolute, in lattice units
    box_fill_radius = spread["box_fill_radius"]  # scales with grid_n
    confined_fraction = spread["confined_fraction"]

    print(f"\n  Localization:")
    print(f"    spread_ratio (strain, interior)  = {spread_ratio:.4f}  "
          f"(want << 0.5 for localized)")
    print(f"    radius_rms (absolute, lattice a) = {radius_rms:.4f}  "
          f"(box-INDEPENDENT if self-localized)")
    print(f"    box_fill_radius                  = {box_fill_radius:.4f}  "
          f"(reference: RMS radius of all nodes)")
    print(f"    confined_fraction (at 0.5*bfr)   = {confined_fraction:.4f}")

    # --- Hessian ---
    hessian = {"skipped": True, "n_neg": None, "n_neg_genuine": None,
               "eig_min": None, "stable": None, "n_dofs": n_interior * m_ambient}
    if do_hessian:
        n_dofs_est = n_interior * m_ambient
        print(f"\n  Static Hessian (n_dofs={n_dofs_est}, max={hessian_max_dofs})...")
        hessian = compute_static_hessian(
            R_final, lattice, ap, is_boundary, mass,
            max_n_dofs=hessian_max_dofs,
        )
        if not hessian.get("skipped", True):
            print(f"    N_neg={hessian['n_neg']}, N_neg_genuine={hessian['n_neg_genuine']}, "
                  f"N_zero={hessian['n_zero']}, lambda_min={hessian['eig_min']:.4e}, "
                  f"stable={hessian['stable']}, time={hessian['hessian_walltime_s']:.1f}s")
            print(f"    Bottom 10 eigenvalues: "
                  f"{[f'{e:.3e}' for e in hessian['eigenvalues_bottom10']]}")
        else:
            print(f"    Hessian skipped: {hessian.get('reason', '?')}")

    # --- Summary ---
    B_preserved = abs(B_final - 1.0) < 0.15
    V_positive = V_star > 0.0
    localized = spread_ratio < 0.5

    print(f"\n  SUMMARY grid_n={grid_n}:")
    print(f"    B_final={B_final:.4f}  ({'OK' if B_preserved else 'FAIL'})")
    print(f"    V*={V_star:.4e}  ({'V>0 OK' if V_positive else 'FAIL V<=0'})")
    print(f"    spread_ratio={spread_ratio:.4f}  "
          f"({'LOCALIZED' if localized else 'DELOCALIZED'})")
    print(f"    radius_rms={radius_rms:.4f}a  (absolute soliton width proxy)")
    print(f"    min|delta_R|={desc['min_dist']:.4f}a")
    n_neg = hessian.get("n_neg_genuine", None)
    eig_min = hessian.get("eig_min", None)
    print(f"    N_neg_genuine={n_neg}, lambda_min={eig_min}  "
          f"({'stable' if hessian.get('stable') else 'UNSTABLE' if hessian.get('stable') is False else 'skipped'})")
    if not_converged:
        print(f"    *** NOT CONVERGED -- results unreliable, interpret with caution ***")

    return {
        "grid_n": grid_n,
        "box_size": box_size,
        "box_half": box_half,
        "box_over_w": box_size / w,
        "alpha": alpha,
        "u0": u0,
        "w": w,
        "F_at_edge_power2": F_at_edge_power2,
        "n_boundary": n_boundary,
        "n_interior": n_interior,
        "n_dofs": n_interior * m_ambient,
        # winding
        "B_seed": B_seed,
        "B_final": B_final,
        "B_preserved": B_preserved,
        # energy
        "V_seed": V_seed,
        "V_star": V_star,
        "V_positive": V_positive,
        # descent
        "descent_status": desc["status"],
        "grad_norm_final": grad_norm_final,
        "n_descent_steps": desc["n_steps"],
        "min_dist_final": desc["min_dist"],
        "descent_walltime_s": desc["walltime_s"],
        "not_converged": not_converged,
        # localization
        "spread_ratio": spread_ratio,
        "radius_rms": radius_rms,       # ABSOLUTE -- key for Derrick verdict
        "box_fill_radius": box_fill_radius,
        "confined_fraction": confined_fraction,
        "localized": localized,
        # hessian
        "hessian_skipped": hessian.get("skipped", True),
        "N_neg": hessian.get("n_neg", None),
        "N_neg_genuine": hessian.get("n_neg_genuine", None),
        "lambda_min": hessian.get("eig_min", None),
        "hessian_stable": hessian.get("stable", None),
        "hessian_walltime_s": hessian.get("hessian_walltime_s", None),
        "hessian_bottom10": hessian.get("eigenvalues_bottom10", None),
    }


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Box-doubling study: self-localized vs boundary-confined?"
    )
    parser.add_argument("--quick", action="store_true",
                        help="Skip Hessian at grid_n=26 (too slow) -- "
                             "spread_ratio + V* still answer the localization question.")
    parser.add_argument("--max-steps", type=int, default=8000,
                        help="Max gradient-flow steps per box (default 8000).")
    parser.add_argument("--tol", type=float, default=1e-5,
                        help="Gradient convergence tolerance (default 1e-5).")
    args = parser.parse_args()

    # Grid sizes for the box-doubling study.
    # grid_n=13: baseline (existing result)
    # grid_n=20: intermediate (~1.5x)
    # grid_n=26: doubled (2x = 26/13)
    # Hessian feasibility: n_interior = (grid_n - 2)^3, n_dofs = 4 * n_interior
    #   grid_n=13: interior=11^3=1331, dofs=5324   (tractable ~15s)
    #   grid_n=20: interior=18^3=5832, dofs=23328  (too large for dense Hessian -- skip)
    #   grid_n=26: interior=24^3=13824, dofs=55296 (too large -- skip unless --quick not set
    #              and user has patience; with max_n_dofs guard it will auto-skip)
    # So: Hessian at 13 and 26 (will auto-skip if > max_n_dofs=10000).
    HESSIAN_MAX_DOFS = 10000  # auto-skip if n_dofs > this

    grid_sizes = [13, 20, 26]
    results = []

    for grid_n in grid_sizes:
        do_hessian = True  # always attempt; auto-skips if too large
        if args.quick and grid_n == 26:
            do_hessian = False

        r = run_box(
            grid_n=grid_n,
            alpha=0.7,
            u0=1.0,
            w=2.0,
            a=1.0,
            k_s=1.0,
            rho=1.0,
            n_boundary_shells=1,
            max_steps=args.max_steps,
            step_size=5e-4,
            tol_grad=args.tol,
            do_hessian=do_hessian,
            hessian_max_dofs=HESSIAN_MAX_DOFS,
            verbose=True,
        )
        results.append(r)

    # -----------------------------------------------------------------------
    # Write CSV
    # -----------------------------------------------------------------------
    csv_path = os.path.join(OUTPUT_DIR, "box_doubling_sweep.csv")
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

    print(f"\nCSV written to: {csv_path}")

    # -----------------------------------------------------------------------
    # Verdict: self-localized or boundary-confined?
    # -----------------------------------------------------------------------
    print(f"\n{'='*65}")
    print("BOX-DOUBLING VERDICT TABLE")
    print(f"{'='*65}")
    print(f"{'grid_n':>8} {'box':>6} {'V*':>12} {'radius_rms':>12} "
          f"{'spread_ratio':>13} {'confined_frac':>14} {'B_final':>8} "
          f"{'N_neg':>6} {'lam_min':>10} {'status':>15}")
    print("-" * 115)
    for r in results:
        status_str = r.get("descent_status", "?")
        if r.get("not_converged", False):
            status_str += "(!)"
        n_neg_str = str(r.get("N_neg_genuine", "skip"))
        lam_str = f"{r['lambda_min']:.3e}" if r.get("lambda_min") is not None else "skip"
        print(f"{r['grid_n']:>8} {r['box_size']:>6.0f}a {r['V_star']:>12.4e} "
              f"{r['radius_rms']:>12.4f} {r['spread_ratio']:>13.4f} "
              f"{r['confined_fraction']:>14.4f} {r['B_final']:>8.4f} "
              f"{n_neg_str:>6} {lam_str:>10} {status_str:>15}")

    # Compute trends
    converged_results = [r for r in results if not r.get("not_converged", False)]
    if len(converged_results) >= 2:
        r0 = converged_results[0]
        r1 = converged_results[-1]
        box_ratio = r1["box_size"] / r0["box_size"]
        radius_ratio = r1["radius_rms"] / r0["radius_rms"] if r0["radius_rms"] > 0 else float("nan")
        V_ratio = r1["V_star"] / r0["V_star"] if r0["V_star"] > 0 else float("nan")
        spread_delta = r1["spread_ratio"] - r0["spread_ratio"]

        print(f"\n  Box growth factor (largest / smallest converged): {box_ratio:.2f}x")
        print(f"  radius_rms growth factor:  {radius_ratio:.3f}  "
              f"(~1.0 = self-localized; ~{box_ratio:.1f} = boundary-confined)")
        print(f"  V* ratio (largest/smallest box): {V_ratio:.3f}  "
              f"(~1.0 = box-independent min)")
        print(f"  spread_ratio delta: {spread_delta:+.4f}  "
              f"(negative = shrinks relative to box = self-localized)")
    else:
        print("\n  WARNING: < 2 converged runs -- cannot compute trend ratios.")
        print("  Consider increasing --max-steps and re-running.")

    # Derrick verdict
    print(f"\n{'='*65}")
    print("DERRICK VERDICT:")
    if len(converged_results) < 2:
        print("  INCONCLUSIVE: not enough converged runs for a trend.")
    else:
        r_small = converged_results[0]
        r_large = converged_results[-1]
        radius_ratio = r_large["radius_rms"] / r_small["radius_rms"]
        box_ratio = r_large["box_size"] / r_small["box_size"]

        # Self-localized: radius_rms roughly constant (ratio << box_ratio)
        # Boundary-confined: radius_rms grows with box (ratio ~ box_ratio)
        # Threshold: if radius_ratio > 0.7 * box_ratio, likely boundary-confined
        if radius_ratio < 1.3:
            verdict = "SELF-LOCALIZED"
            detail = (f"radius_rms changed by factor {radius_ratio:.3f} while "
                      f"box grew by {box_ratio:.2f}x -- soliton width is box-independent.")
        elif radius_ratio > 0.7 * box_ratio:
            verdict = "BOUNDARY-CONFINED (Derrick wins)"
            detail = (f"radius_rms grew by factor {radius_ratio:.3f} tracking "
                      f"box growth of {box_ratio:.2f}x -- width scales with box.")
        else:
            verdict = "AMBIGUOUS"
            detail = (f"radius_rms grew by {radius_ratio:.3f}x; box grew {box_ratio:.2f}x. "
                      f"More steps or larger box needed.")

        print(f"  {verdict}")
        print(f"  {detail}")

        # V* independence check
        if abs(V_ratio - 1.0) < 0.3:
            print(f"  V* ratio={V_ratio:.3f}: energy minimum is approximately box-independent.")
        else:
            print(f"  V* ratio={V_ratio:.3f}: energy minimum depends on box size -- "
                  f"suggests boundary terms contribute to V*.")

    print(f"{'='*65}\n")

    # -----------------------------------------------------------------------
    # Per-box JSON
    # -----------------------------------------------------------------------
    for r in results:
        json_path = os.path.join(
            OUTPUT_DIR, f"box_doubling_g{r['grid_n']}.json"
        )
        r_json = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                  for k, v in r.items() if k != "hessian_bottom10"}
        r_json["hessian_bottom10"] = r.get("hessian_bottom10")
        with open(json_path, "w") as f:
            json.dump(r_json, f, indent=2, default=str)
        print(f"JSON saved: {json_path}")

    print(f"\nAll outputs in: {OUTPUT_DIR}")
    print(f"CSV:  {csv_path}")


if __name__ == "__main__":
    main()
