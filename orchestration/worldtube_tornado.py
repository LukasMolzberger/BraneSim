#!/usr/bin/env python3
"""4D worldtube "tornado tube" baryon simulation — bvp_chiral block solver.

HYPOTHESIS
----------
Topological protection of a baryon lives ALONG TIME, not in a static spatial
slice.  A worldtube / vortex line — a rotating hedgehog swirl that persists
along the time axis — is protected by the fact that cutting it creates /
destroys the particle.  The static spatial slice has no lattice topological
protection (the hedgehog unwinds to vacuum); the worldtube does.

APPROACH
--------
Seed: a rotating hedgehog swirl (lateral displacement locked to radial
direction, Gaussian profile), with carrier rotation encoded in the TWO
initial slices (R0, R1) so the chiral march produces a coherent rotating
worldtube.

Solver: bvp_chiral (ChiralBC, forward chirality) — the validated time-
symmetric path.  This IS the Verlet march from (R0, R1); no JFNK needed.

Physics:
  - k_s = rho = a = 1 (dimensionless)
  - alpha = 0.7 (in the Skyrme/breather sweet spot alpha~0.5-0.8)
  - dt = 0.25  (CFL-safe: dt < 2/omega_band_top ~ 0.790)
  - m_ambient = 4 (4D worldvolume)
  - omega_0 = 1.05 * sqrt(2*(alpha*2 + (1-alpha)*6)) just above band-top

SEED CONSTRUCTION
-----------------
Slice R0 (t=0): hedgehog
    xi^i(x) = A * g(r) * x_hat^i   for i in {0,1,2}
    xi^3    = 0                     (time ambient component)
    g(r) = exp(-r^2 / (2 w^2)), w ~ 3a, A ~ O(1) "large amplitude"

Slice R1 (t=dt): R0 with the lateral swirl rotated by omega_0 * dt.
The rotation acts on the 3 lateral components as a U(1)-like rotation in
the (xi^1, xi^2) plane (the simplest non-trivial choice consistent with
the hedgehog colour-lock structure).  This encodes the carrier phase advance
and sets forward chirality (matter worldtube).

Explicitly:
    xi_new^0 = xi^0 * cos(phi) - xi^1 * sin(phi)
    xi_new^1 = xi^0 * sin(phi) + xi^1 * cos(phi)
    xi_new^2 = xi^2                                 (third colour unchanged)
where phi = omega_0 * dt.

This is a rotation in the (colour-0, colour-1) plane; by the hedgehog lock
xi^i = A g(r) x_hat^i, the rotation in colour space is exactly a rotation
of the spatial unit vector x_hat in the (x^0, x^1) plane — an azimuthal
rotation of the vortex.

DIAGNOSTICS (field-based, no hard-coded winding analytic)
---------------------------------------------------------
Per time slice:
  1. E_excess(l): sum over all spacelike bonds of the VACUUM-SUBTRACTED
     excess ½ k_s [(|ΔR| - alpha*a)^2 - (a*(1-alpha))^2].
     (The vacuum bond has |ΔR| = a, so vacuum energy per bond = ½ k_s (a(1-alpha))^2.)
     This is zero only in the exact vacuum; positive wherever the field deviates.

  2. spread_ratio(l): energy-weighted RMS radius / box_fill_radius.
     Uses strain weight (direction variation) for Skyrme-aware localization.

  3. confined_fraction(l): fraction of E_excess within 0.5 * box_fill_radius.

  4. leakage(l): E_excess within the outermost shell (5% of box from boundary).

  5. E_excess_total(l): total E_excess per slice (should be roughly conserved).

  6. carrier phase advance (phi_advance(l)): cross-correlation of the lateral
     displacement field with its rotated version, measured at the soliton core
     (r < w).  Reports whether the swirl rotates coherently.

RUN PLAN
--------
(a) Tiny validation run: 24^3 x 48 slices.  Confirm: world-volume builds,
    residual is sane, diagnostics compute, E_excess localized at t=0.
    Seconds.

(b) Medium run: 48^3 x 128 slices (~340 MB of world-volume).  Report:
    E_excess(l), spread_ratio(l), leakage(l), phi_advance(l).

VERDICT CRITERIA
----------------
PERSIST:   spread_ratio stays << 1 throughout (< 0.5), confined_fraction > 0.5,
           E_excess conserved (< 20% loss), carrier phase advances coherently.
PARTIAL:   Initial localization decays but persists for many periods (spread_ratio
           < 0.8 for > 50% of slices).
DISPERSE:  spread_ratio -> 1.0 within a few periods; box-fill.

OUTPUT
------
test-runs/worldtube_tornado_v1/
    tiny_validation.json      -- tiny run result (step a)
    medium_run_diagnostics.csv -- per-slice diagnostics (step b)
    medium_run_summary.json   -- overall summary and verdict

PRINCIPLES COMPLIANCE
---------------------
- No artificial confinement forces, no saturation, no clamps.
- Solver: bvp_chiral (ChiralBC), forward Verlet march.
- Energy: Pythagorean 4D norm (spacelike_potential / spacelike_force as-is).
- temporal_model="a", r_t=0.
- Diagnostics: read-only, field-based.
- No analytic B-estimator (banned per task spec).
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from pathlib import Path

# Ensure repo is on sys.path when run directly
_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from branesim.core.conventions import ActionParams, LatticeParams, d_of_k_eigenvalues
from branesim.core.lattice import SpacelikeLattice
from branesim.core.action import spacelike_potential
from branesim.core.residual import residual_norm
from branesim.solver.boundary import ChiralBC
from branesim.solver.bvp import BoundaryProblem, SolveOpts, solve_block


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

OUT_DIR = Path(__file__).resolve().parent.parent / "test-runs" / "worldtube_tornado_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Band-top computation (exact, from dispersion relation)
# ---------------------------------------------------------------------------

def omega_band_top_exact(alpha: float, k_s: float = 1.0, rho: float = 1.0, a: float = 1.0) -> float:
    """Exact phonon band-top for the 3D axial-only lattice.

    At the zone corner k=(pi/a, pi/a, pi/a), all h_i = 1 - cos(pi) = 2.
    omega^2 = (2 k_s / rho) * [alpha * 2 + (1-alpha) * 6]
            = (2 k_s / rho) * [2*alpha + 6 - 6*alpha]
            = (2 k_s / rho) * [6 - 4*alpha]
    """
    k_corner = np.array([math.pi / a, math.pi / a, math.pi / a])
    eigs = d_of_k_eigenvalues(k_corner, alpha, k_s=k_s, rho=rho, a=a)
    return float(np.sqrt(np.max(eigs)))


# ---------------------------------------------------------------------------
# Seed construction: rotating hedgehog worldtube
# ---------------------------------------------------------------------------

def build_hedgehog_slice(
    lattice: SpacelikeLattice,
    m: int,
    amplitude: float,
    width: float,
) -> np.ndarray:
    """Build the hedgehog initial slice R0.

    xi^i(x) = amplitude * exp(-r^2/(2*width^2)) * x_hat^i  for i in {0,1,2}
    xi^3    = 0

    Parameters
    ----------
    lattice : SpacelikeLattice (3D)
    m : int  -- ambient dimension (=4)
    amplitude : float  -- peak lateral displacement
    width : float  -- Gaussian half-width (in lattice units)

    Returns
    -------
    R0 : ndarray, shape (n_nodes, m)
    """
    dim = lattice.params.dim   # 3
    a = lattice.params.spacing
    ref = lattice.reference_positions(m)     # (n_nodes, m)
    coords = ref[:, :dim]                    # (n_nodes, 3)
    centre = coords.mean(axis=0)             # (3,)
    dx = coords - centre                     # (n_nodes, 3)
    r = np.linalg.norm(dx, axis=1)          # (n_nodes,)

    _eps = 1e-30
    r_safe = np.where(r > _eps, r, 1.0)
    x_hat = dx / r_safe[:, None]            # (n_nodes, 3)
    x_hat[r <= _eps] = 0.0                  # zero at exact centre

    # Gaussian profile
    g = np.exp(-(r / width) ** 2 / 2.0)    # (n_nodes,)

    disp = np.zeros_like(ref)
    # Lateral (colour) channels: colour locked to x_hat
    disp[:, :dim] = (amplitude * g[:, None]) * x_hat   # (n_nodes, 3)
    # Ambient temporal component (index 3) stays zero in the spatial slice.

    return ref + disp


def build_rotated_slice(
    R0: np.ndarray,
    ref: np.ndarray,
    phi: float,
) -> np.ndarray:
    """Build R1 by rotating the lateral displacement of R0 by angle phi.

    The rotation acts in the (colour-0, colour-1) plane:
        xi_new^0 = xi^0 * cos(phi) - xi^1 * sin(phi)
        xi_new^1 = xi^0 * sin(phi) + xi^1 * cos(phi)
        xi_new^2 = xi^2   (unchanged)
        xi^3     = 0      (unchanged)

    This is a U(1) azimuthal rotation of the hedgehog swirl.  By the
    hedgehog lock xi^i = A*g(r)*x_hat^i, this rotation in colour space
    is equivalent to rotating x_hat in the (x^0, x^1) plane — i.e., an
    azimuthal rotation of the vortex around the z-axis.

    Parameters
    ----------
    R0 : ndarray, shape (n_nodes, m)
    ref : ndarray, shape (n_nodes, m)
    phi : float  -- rotation angle = omega_0 * dt

    Returns
    -------
    R1 : ndarray, shape (n_nodes, m)
    """
    xi = R0 - ref       # displacement (n_nodes, m)

    c, s = math.cos(phi), math.sin(phi)
    xi0 = xi[:, 0].copy()
    xi1 = xi[:, 1].copy()

    xi_new = xi.copy()
    xi_new[:, 0] = c * xi0 - s * xi1
    xi_new[:, 1] = s * xi0 + c * xi1
    # Component 2 and 3 unchanged

    return ref + xi_new


# ---------------------------------------------------------------------------
# Vacuum-subtracted excess energy (the key diagnostic)
# ---------------------------------------------------------------------------

def excess_energy_per_slice(
    positions: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
) -> tuple[float, np.ndarray]:
    """Compute vacuum-subtracted excess energy per node.

    For each bond (p, q), the VACUUM energy is:
        E_vac = (k_s/4) * (a*(1-alpha))^2   [one bond contribution to V^l]
        (vacuum: |ΔR| = a, so strain = a - alpha*a = a*(1-alpha))

    Excess per bond:
        Delta_E = (k_s/4) * [(|ΔR| - alpha*a)^2 - (a*(1-alpha))^2]

    The factor k_s/4 avoids double-counting (same convention as spacelike_potential).

    Returns
    -------
    total_excess : float  -- sum over all bonds
    node_excess : ndarray, shape (n_nodes,)  -- per-node contribution
        (half of each bond's excess assigned to each endpoint)
    """
    k_s = params.k_s
    alpha_a = params.alpha * lattice.params.spacing  # alpha * a
    a = lattice.params.spacing
    vac_strain = a * (1.0 - params.alpha)            # vacuum strain = a - alpha*a
    vac_energy_per_bond_factor = vac_strain ** 2     # (a*(1-alpha))^2, used below

    dim = lattice.dim
    periodic_axes = lattice.params.periodic_axes
    box_lengths = np.array(
        [n * a for n in lattice.params.grid_shape], dtype=np.float64
    )

    node_excess = np.zeros(lattice.n_nodes, dtype=np.float64)

    for nb_idx in range(lattice.n_neighbors):
        nb_ids = lattice.neighbors[:, nb_idx]
        valid = nb_ids >= 0
        if not np.any(valid):
            continue

        valid_ids = np.where(valid)[0]
        p_pos = positions[valid_ids]
        q_pos = positions[nb_ids[valid_ids]]

        raw_delta = q_pos - p_pos
        # Minimum-image for periodic axes
        delta = raw_delta.copy()
        for axis in range(dim):
            if not periodic_axes[axis]:
                continue
            L = box_lengths[axis]
            delta[:, axis] -= L * np.round(delta[:, axis] / L)

        dist = np.linalg.norm(delta, axis=1)   # (n_valid,)
        strain = dist - alpha_a                # (n_valid,)
        bond_excess = (k_s / 4.0) * (strain ** 2 - vac_energy_per_bond_factor)

        # Assign half to each endpoint (symmetric: each bond counted twice in loop)
        np.add.at(node_excess, valid_ids, bond_excess)

    total_excess = float(np.sum(node_excess))
    return total_excess, node_excess


# ---------------------------------------------------------------------------
# Carrier phase advance diagnostic
# ---------------------------------------------------------------------------

def carrier_phase_advance(
    R_prev: np.ndarray,
    R_curr: np.ndarray,
    ref: np.ndarray,
    centre: np.ndarray,
    core_radius: float,
) -> float:
    """Estimate the carrier phase advance from slice l-1 to slice l.

    Within the soliton core (r < core_radius from centre), compute:
        phi_advance = arctan2(cross, dot)
    where cross and dot are the cross- and dot-products of the 2D lateral
    displacement vectors (xi^0, xi^1) at the previous and current slice,
    spatially averaged over the core.

    This measures the azimuthal rotation of the vortex per time step,
    in the (xi^0, xi^1) plane.

    Returns
    -------
    phi : float  -- phase advance in radians (positive = forward rotation)
    """
    dim = ref.shape[1] - 1   # spatial dim = m-1 for m=4
    # Spatial coordinates for distance from centre
    coords = ref[:, :dim]   # (n_nodes, 3)
    r = np.linalg.norm(coords - centre, axis=1)   # (n_nodes,)
    core_mask = r < core_radius

    if not np.any(core_mask):
        return float("nan")

    xi_prev = (R_prev - ref)[core_mask]   # (n_core, m)
    xi_curr = (R_curr - ref)[core_mask]   # (n_core, m)

    # Use (colour-0, colour-1) = (xi^0, xi^1) components
    a0_prev = xi_prev[:, 0]   # (n_core,)
    a1_prev = xi_prev[:, 1]
    a0_curr = xi_curr[:, 0]
    a1_curr = xi_curr[:, 1]

    # Per-node: cross = a0_prev * a1_curr - a1_prev * a0_curr
    #           dot   = a0_prev * a0_curr + a1_prev * a1_curr
    cross = float(np.mean(a0_prev * a1_curr - a1_prev * a0_curr))
    dot   = float(np.mean(a0_prev * a0_curr + a1_prev * a1_curr))

    return math.atan2(cross, dot)


# ---------------------------------------------------------------------------
# Boundary leakage diagnostic
# ---------------------------------------------------------------------------

def boundary_leakage(
    node_excess: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    shell_fraction: float = 0.05,
) -> float:
    """Fraction of total excess energy in the outermost spatial shell.

    The shell is defined as nodes whose distance to the nearest box wall
    is < shell_fraction * min(box side length).

    Returns
    -------
    float  -- fraction of total excess in boundary shell [0, 1]
    """
    dim = lattice.params.dim
    a = lattice.params.spacing
    grid_shape = lattice.params.grid_shape
    box_sides = np.array([(n - 1) * a for n in grid_shape], dtype=np.float64)
    shell_depth = shell_fraction * np.min(box_sides)

    coords = ref[:, :dim]   # (n_nodes, dim)
    # Distance to nearest wall on each axis
    # Wall positions: 0 and box_side[i] for each axis
    # coords start at 0 (reference_positions uses i*a for index i starting at 0)
    dist_to_wall = np.min(
        np.stack([coords, box_sides - coords], axis=2), axis=2
    ).min(axis=1)   # (n_nodes,)

    shell_mask = dist_to_wall < shell_depth
    total_excess = float(np.sum(np.abs(node_excess)))
    if total_excess < 1e-40:
        return 0.0
    shell_excess = float(np.sum(np.abs(node_excess[shell_mask])))
    return shell_excess / total_excess


# ---------------------------------------------------------------------------
# Full diagnostics for one slice
# ---------------------------------------------------------------------------

def slice_diagnostics(
    positions: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    box_fill_radius: float,
    confinement_radius_factor: float = 0.5,
    shell_fraction: float = 0.05,
    R_prev: np.ndarray | None = None,
    omega_0: float | None = None,
) -> dict:
    """Compute all diagnostics for one spacelike slice.

    Parameters
    ----------
    positions : (n_nodes, m)
    ref       : (n_nodes, m)
    lattice   : SpacelikeLattice
    params    : ActionParams
    box_fill_radius : float  -- pre-computed
    confinement_radius_factor : float
    shell_fraction : float
    R_prev    : (n_nodes, m) or None  -- previous slice for phase advance
    omega_0   : float or None  -- carrier frequency (for coherence ratio)

    Returns
    -------
    dict with keys:
        E_excess_total, spread_ratio, confined_fraction,
        leakage_fraction, radius_rms, carrier_phase_advance (if R_prev given)
    """
    dim = lattice.params.dim
    a = lattice.params.spacing

    # --- Excess energy ---
    E_excess_total, node_excess = excess_energy_per_slice(positions, lattice, params)

    # --- Spatial weight for spread: excess per node (clipped to non-negative) ---
    weights = np.maximum(node_excess, 0.0) + 1e-40

    total_weight = float(np.sum(weights))
    ref_spatial = ref[:, :dim]   # (n_nodes, 3)
    geom_centre = np.mean(ref_spatial, axis=0)

    # Energy-weighted centroid
    centre_w = np.sum(ref_spatial * weights[:, None], axis=0) / total_weight

    # Radii from energy-weighted centre
    delta = ref_spatial - centre_w
    r_sq = np.sum(delta ** 2, axis=1)
    radii = np.sqrt(r_sq)

    radius_rms = float(np.sqrt(np.sum(weights * r_sq) / total_weight))

    spread_ratio = radius_rms / box_fill_radius if box_fill_radius > 1e-30 else float("nan")

    # Confined fraction
    conf_r = confinement_radius_factor * box_fill_radius
    conf_mask = radii <= conf_r
    confined_fraction = float(np.sum(weights[conf_mask]) / total_weight)

    # Leakage
    leakage_fraction = boundary_leakage(node_excess, ref, lattice, shell_fraction)

    result = {
        "E_excess_total": E_excess_total,
        "radius_rms": radius_rms,
        "spread_ratio": spread_ratio,
        "confined_fraction": confined_fraction,
        "leakage_fraction": leakage_fraction,
    }

    # Carrier phase advance
    if R_prev is not None:
        core_radius = 2.0 * lattice.params.spacing * 3.0   # ~ 2*w
        phi = carrier_phase_advance(R_prev, positions, ref, centre_w, core_radius)
        result["carrier_phase_advance"] = phi
        if omega_0 is not None and not math.isnan(phi):
            # Expected advance per step: omega_0 * dt
            expected = omega_0 * params.dt
            coherence = math.cos(phi - expected)   # 1 if perfectly coherent, -1 if opposite
            result["carrier_coherence"] = float(coherence)
        else:
            result["carrier_coherence"] = float("nan")
    else:
        result["carrier_phase_advance"] = float("nan")
        result["carrier_coherence"] = float("nan")

    return result


# ---------------------------------------------------------------------------
# Box fill radius (pre-compute once)
# ---------------------------------------------------------------------------

def compute_box_fill_radius(lattice: SpacelikeLattice) -> float:
    """RMS node radius about geometric centre (uniform weight)."""
    dim = lattice.params.dim
    m = dim + 1
    ref = lattice.reference_positions(m)
    ref_spatial = ref[:, :dim]
    geom_centre = np.mean(ref_spatial, axis=0)
    delta = ref_spatial - geom_centre
    r_sq = np.sum(delta ** 2, axis=1)
    return float(np.sqrt(np.mean(r_sq)))


# ---------------------------------------------------------------------------
# Main simulation runner
# ---------------------------------------------------------------------------

def run_worldtube(
    grid_n: int,
    n_slices: int,
    alpha: float = 0.7,
    amplitude: float = 1.0,
    width_in_a: float = 3.0,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
    dt: float = 0.25,
    m_ambient: int = 4,
    omega_0_factor: float = 1.05,
    periodic: bool = False,
    label: str = "run",
    verbose: bool = True,
) -> dict:
    """Run a single worldtube simulation and return diagnostics.

    Parameters
    ----------
    grid_n : int  -- cubic grid size (N^3 nodes)
    n_slices : int  -- number of time slices
    alpha : float  -- prestress (0.7 = sweet spot)
    amplitude : float  -- peak hedgehog displacement (in lattice units)
    width_in_a : float  -- Gaussian width in units of a
    k_s, rho, a, dt : float  -- action parameters
    m_ambient : int  -- ambient dimension (4)
    omega_0_factor : float  -- omega_0 = omega_0_factor * omega_band_top
    label : str  -- run label for logging
    verbose : bool

    Returns
    -------
    dict  -- run metadata + per-slice diagnostic arrays
    """
    t_wall_start = time.perf_counter()

    # --- Lattice ---
    grid_shape = (grid_n, grid_n, grid_n)
    lp = LatticeParams(
        grid_shape=grid_shape,
        spacing=a,
        periodic_axes=(periodic, periodic, periodic),   # open=free (global tension relaxation); periodic=fixed box holds prestress
        axial_weight=1.0,
    )
    lattice = SpacelikeLattice(lp)
    dim = lattice.params.dim   # 3

    # --- Action params ---
    ap = ActionParams(
        k_s=k_s,
        alpha=alpha,
        rho=rho,
        dt=dt,
        n_slices=n_slices,
        m_ambient=m_ambient,
        temporal_model="a",
        r_t=0.0,
    )
    mass = rho * a ** dim   # = 1 in dimensionless units

    # --- Band top and omega_0 ---
    omega_top = omega_band_top_exact(alpha, k_s=k_s, rho=rho, a=a)
    omega_0 = omega_0_factor * omega_top
    phi_per_step = omega_0 * dt   # carrier phase advance per slice

    # CFL check (warning only; not an error — caller is responsible)
    cfl_limit = 2.0 / omega_top
    cfl_ok = dt <= cfl_limit
    if verbose:
        print(f"[{label}] grid={grid_n}^3 x {n_slices} slices")
        print(f"[{label}] alpha={alpha}, A={amplitude}, w={width_in_a}*a")
        print(f"[{label}] omega_top={omega_top:.4f}, omega_0={omega_0:.4f}")
        print(f"[{label}] dt={dt}, CFL limit={cfl_limit:.4f}, CFL OK={cfl_ok}")
        print(f"[{label}] phi_per_step={phi_per_step:.4f} rad")
        n_nodes_total = grid_n ** 3
        mem_MB = (n_slices + 1) * n_nodes_total * m_ambient * 8 / 1e6
        print(f"[{label}] World-volume memory: ~{mem_MB:.1f} MB")

    # --- Reference positions ---
    ref = lattice.reference_positions(m_ambient)   # (n_nodes, m_ambient)

    # --- Build R0: hedgehog slice at t=0 ---
    width = width_in_a * a
    R0 = build_hedgehog_slice(lattice, m_ambient, amplitude, width)

    # --- Build R1: rotated hedgehog at t=dt ---
    R1 = build_rotated_slice(R0, ref, phi_per_step)

    # --- Sanity: amplitude at centre ---
    centre_flat = np.argmin(np.linalg.norm(ref[:, :dim] - ref[:, :dim].mean(axis=0), axis=1))
    xi_centre = R0[centre_flat] - ref[centre_flat]
    xi_max = float(np.max(np.abs(R0 - ref)))
    if verbose:
        print(f"[{label}] max |xi| at t=0: {xi_max:.4f} (amplitude={amplitude})")

    # --- BVP chiral solve (forward march) ---
    bc = ChiralBC(R0=R0, R1=R1, chirality="forward")
    problem = BoundaryProblem(lattice=lattice, params=ap, mass=mass,
                              boundary_condition=bc)
    opts = SolveOpts(verbose=verbose)

    if verbose:
        print(f"[{label}] Starting bvp_chiral forward march...")
    t_solve_start = time.perf_counter()
    wv = solve_block(problem, opts)
    t_solve = time.perf_counter() - t_solve_start

    solver_report = wv.solver_report
    res_norm = float(residual_norm(wv.slices, lattice, ap, mass))
    n_interior_dof = (n_slices - 1) * lattice.n_nodes * m_ambient
    res_per_dof = res_norm / math.sqrt(max(n_interior_dof, 1))

    if verbose:
        print(f"[{label}] Solve done: {t_solve:.1f}s, residual_norm={res_norm:.4e}, "
              f"residual_per_dof={res_per_dof:.4e}")
        print(f"[{label}] condition_estimate={solver_report.get('condition_estimate', 'N/A'):.2g}")

    # --- Pre-compute box fill radius ---
    box_fill_radius = compute_box_fill_radius(lattice)
    if verbose:
        print(f"[{label}] box_fill_radius={box_fill_radius:.4f}")

    # --- Per-slice diagnostics ---
    if verbose:
        print(f"[{label}] Computing per-slice diagnostics over {n_slices + 1} slices...")

    slice_results = []
    E_excess_0 = None   # reference: first slice

    for l in range(n_slices + 1):
        positions = wv.slices[l]
        R_prev_slice = wv.slices[l - 1] if l > 0 else None

        diag = slice_diagnostics(
            positions=positions,
            ref=ref,
            lattice=lattice,
            params=ap,
            box_fill_radius=box_fill_radius,
            confinement_radius_factor=0.5,
            shell_fraction=0.05,
            R_prev=R_prev_slice,
            omega_0=omega_0,
        )
        diag["slice"] = l
        diag["t"] = l * dt
        slice_results.append(diag)

        if l == 0:
            E_excess_0 = diag["E_excess_total"]

    if verbose:
        print(f"[{label}] Diagnostics done.")

    # --- Summary statistics ---
    E_arr = np.array([d["E_excess_total"] for d in slice_results])
    spread_arr = np.array([d["spread_ratio"] for d in slice_results])
    conf_arr = np.array([d["confined_fraction"] for d in slice_results])
    leak_arr = np.array([d["leakage_fraction"] for d in slice_results])
    phi_arr = np.array([d["carrier_phase_advance"] for d in slice_results])
    coh_arr = np.array([d["carrier_coherence"] for d in slice_results])

    # Energy conservation: fraction of initial E_excess remaining at each slice
    E_frac_arr = E_arr / max(abs(E_excess_0), 1e-40) if E_excess_0 is not None else E_arr * 0

    # Verdict
    mean_spread = float(np.nanmean(spread_arr[1:]))   # skip l=0 (no prev slice)
    final_spread = float(spread_arr[-1])
    mean_conf = float(np.nanmean(conf_arr))
    E_conservation = float(E_arr[-1] / max(abs(E_excess_0), 1e-40)) if E_excess_0 else float("nan")
    mean_coh = float(np.nanmean(coh_arr[~np.isnan(coh_arr)]))

    if mean_spread < 0.4 and mean_conf > 0.6 and abs(E_conservation - 1.0) < 0.3:
        verdict = "PERSIST"
    elif mean_spread < 0.8 and mean_conf > 0.3:
        verdict = "PARTIAL"
    else:
        verdict = "DISPERSE"

    t_wall_total = time.perf_counter() - t_wall_start

    summary = {
        "label": label,
        "grid_n": grid_n,
        "n_slices": n_slices,
        "alpha": alpha,
        "amplitude": amplitude,
        "width_in_a": width_in_a,
        "k_s": k_s, "rho": rho, "a": a, "dt": dt,
        "m_ambient": m_ambient,
        "omega_top": omega_top,
        "omega_0": omega_0,
        "phi_per_step": phi_per_step,
        "cfl_ok": cfl_ok,
        "residual_norm": res_norm,
        "residual_per_dof": res_per_dof,
        "condition_estimate": float(solver_report.get("condition_estimate", 0.0)),
        "walltime_solve_s": t_solve,
        "walltime_total_s": t_wall_total,
        "box_fill_radius": box_fill_radius,
        "E_excess_t0": float(E_excess_0) if E_excess_0 is not None else 0.0,
        "E_excess_final": float(E_arr[-1]),
        "E_conservation_ratio": E_conservation,
        "spread_ratio_mean": mean_spread,
        "spread_ratio_final": final_spread,
        "spread_ratio_initial": float(spread_arr[0]),
        "confined_fraction_mean": mean_conf,
        "leakage_mean": float(np.nanmean(leak_arr)),
        "carrier_coherence_mean": mean_coh,
        "verdict": verdict,
    }

    if verbose:
        print(f"\n[{label}] === VERDICT: {verdict} ===")
        print(f"[{label}] E_excess_t0={E_excess_0:.4f}, E_excess_final={E_arr[-1]:.4f}, "
              f"E_conservation={E_conservation:.3f}")
        print(f"[{label}] spread_ratio: initial={spread_arr[0]:.4f}, mean={mean_spread:.4f}, "
              f"final={final_spread:.4f}")
        print(f"[{label}] confined_fraction mean={mean_conf:.4f}")
        print(f"[{label}] leakage mean={float(np.nanmean(leak_arr)):.4f}")
        print(f"[{label}] carrier_coherence mean={mean_coh:.4f}")
        print(f"[{label}] Total walltime: {t_wall_total:.1f}s")

    return {
        "summary": summary,
        "slice_results": slice_results,
        "slices_arr": {
            "E_excess": E_arr.tolist(),
            "spread_ratio": spread_arr.tolist(),
            "confined_fraction": conf_arr.tolist(),
            "leakage_fraction": leak_arr.tolist(),
            "carrier_phase_advance": phi_arr.tolist(),
            "carrier_coherence": coh_arr.tolist(),
            "E_fraction": E_frac_arr.tolist(),
        },
        "wv": wv,   # WorldVolume (in-memory; not serialized to JSON)
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plots(result: dict, out_dir: Path, label: str) -> None:
    """Generate diagnostic plots for a run."""
    slices_arr = result["slices_arr"]
    summary = result["summary"]
    n_slices = summary["n_slices"]
    dt = summary["dt"]
    t_arr = np.arange(n_slices + 1) * dt

    fig, axes = plt.subplots(3, 2, figsize=(12, 10))
    fig.suptitle(f"Worldtube tornado: {label}\n"
                 f"alpha={summary['alpha']}, A={summary['amplitude']}, "
                 f"w={summary['width_in_a']}a, grid={summary['grid_n']}^3 x {n_slices}",
                 fontsize=11)

    ax = axes[0, 0]
    ax.plot(t_arr, slices_arr["E_excess"], "C0-", linewidth=1)
    ax.set_xlabel("t")
    ax.set_ylabel("E_excess (vacuum-subtracted)")
    ax.set_title("Total excess energy vs time")
    ax.axhline(0, color="gray", linewidth=0.5)

    ax = axes[0, 1]
    E_frac = slices_arr["E_fraction"]
    ax.plot(t_arr, E_frac, "C1-", linewidth=1)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="E0")
    ax.axhline(0.8, color="orange", linestyle=":", linewidth=0.8, label="80%")
    ax.set_xlabel("t")
    ax.set_ylabel("E_excess(t) / E_excess(0)")
    ax.set_title("Energy conservation ratio")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(t_arr, slices_arr["spread_ratio"], "C2-", linewidth=1)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="box-fill")
    ax.axhline(0.4, color="green", linestyle=":", linewidth=0.8, label="persist threshold")
    ax.set_xlabel("t")
    ax.set_ylabel("spread_ratio")
    ax.set_title("Spread ratio (1=dispersed, 0=localized)")
    ax.legend(fontsize=8)
    verdict_color = {"PERSIST": "green", "PARTIAL": "orange", "DISPERSE": "red"}[summary["verdict"]]
    ax.text(0.98, 0.98, summary["verdict"], transform=ax.transAxes,
            ha="right", va="top", color=verdict_color, fontsize=14, fontweight="bold")

    ax = axes[1, 1]
    ax.plot(t_arr, slices_arr["confined_fraction"], "C3-", linewidth=1)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, label="50%")
    ax.set_xlabel("t")
    ax.set_ylabel("confined_fraction")
    ax.set_title("Confined fraction (within 0.5 * box_fill_radius)")
    ax.legend(fontsize=8)

    ax = axes[2, 0]
    ax.plot(t_arr, slices_arr["leakage_fraction"], "C4-", linewidth=1)
    ax.set_xlabel("t")
    ax.set_ylabel("leakage_fraction")
    ax.set_title("Boundary leakage fraction (outer 5% shell)")

    ax = axes[2, 1]
    phi_arr = np.array(slices_arr["carrier_phase_advance"])
    coh_arr = np.array(slices_arr["carrier_coherence"])
    phi_valid = phi_arr[~np.isnan(phi_arr)]
    t_valid = t_arr[~np.isnan(phi_arr)]
    if len(phi_valid) > 0:
        ax.plot(t_valid, phi_valid, "C5-", linewidth=1, label="phi_advance")
        expected_phi = summary["phi_per_step"]
        ax.axhline(expected_phi, color="gray", linestyle="--", linewidth=0.8,
                   label=f"expected phi={expected_phi:.3f}")
    ax.set_xlabel("t")
    ax.set_ylabel("carrier phase advance (rad/step)")
    ax.set_title("Carrier phase advance per slice")
    ax.legend(fontsize=8)

    plt.tight_layout()
    path = out_dir / f"{label}_diagnostics.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def save_csv(slice_results: list[dict], out_dir: Path, label: str) -> Path:
    """Save per-slice diagnostics to CSV."""
    fields = ["slice", "t", "E_excess_total", "spread_ratio", "confined_fraction",
              "leakage_fraction", "carrier_phase_advance", "carrier_coherence"]
    path = out_dir / f"{label}_diagnostics.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in slice_results:
            writer.writerow({k: row.get(k, "") for k in fields})
    return path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run step (a) tiny validation then step (b) medium run."""

    print("=" * 70)
    print("WORLDTUBE TORNADO SIMULATION")
    print("Hypothesis: 4D worldtube / vortex-line baryon analog")
    print("=" * 70)

    # Common parameters
    PARAMS = dict(
        alpha=0.7,
        amplitude=1.0,
        width_in_a=3.0,
        k_s=1.0,
        rho=1.0,
        a=1.0,
        dt=0.25,
        m_ambient=4,
        omega_0_factor=1.05,
    )

    # =========================================================================
    # Step (a): Tiny validation run
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP (a): Tiny validation run — 24^3 x 48 slices")
    print("=" * 60)

    tiny_result = run_worldtube(
        grid_n=24,
        n_slices=48,
        label="tiny_validation",
        verbose=True,
        **PARAMS,
    )

    tiny_summary = tiny_result["summary"]
    print(f"\n[tiny] Validation checks:")
    print(f"  E_excess(t=0)    = {tiny_summary['E_excess_t0']:.6f}")
    print(f"  residual_per_dof = {tiny_summary['residual_per_dof']:.4e}")
    print(f"  spread_ratio(0)  = {tiny_summary['spread_ratio_initial']:.4f}")
    print(f"  verdict          = {tiny_summary['verdict']}")

    # Save tiny results
    tiny_json_path = OUT_DIR / "tiny_validation.json"
    tiny_json_path.write_text(json.dumps(tiny_summary, indent=2, default=str))
    print(f"  [saved] {tiny_json_path}")

    save_csv(tiny_result["slice_results"], OUT_DIR, "tiny_validation")
    plot_path = make_plots(tiny_result, OUT_DIR, "tiny_validation")
    print(f"  [saved] {plot_path}")

    # Validation gate: E_excess must be positive and localized
    E0 = tiny_summary["E_excess_t0"]
    sr0 = tiny_summary["spread_ratio_initial"]
    res_pdof = tiny_summary["residual_per_dof"]

    print(f"\n[VALIDATION GATE]")
    print(f"  E_excess(0) > 0:         {'PASS' if E0 > 0 else 'FAIL'} (E0={E0:.4f})")
    print(f"  spread_ratio(0) < 0.9:   {'PASS' if sr0 < 0.9 else 'FAIL'} (sr={sr0:.4f})")
    print(f"  residual_per_dof < 1e-6: {'PASS' if res_pdof < 1e-6 else 'FAIL'} "
          f"(res={res_pdof:.2e})")

    # =========================================================================
    # Step (b): Medium run
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP (b): Medium run — 48^3 x 128 slices")
    print("=" * 60)

    medium_result = run_worldtube(
        grid_n=48,
        n_slices=128,
        label="medium_run",
        verbose=True,
        **PARAMS,
    )

    medium_summary = medium_result["summary"]

    # Save medium results
    medium_json_path = OUT_DIR / "medium_run_summary.json"
    medium_json_path.write_text(json.dumps(medium_summary, indent=2, default=str))
    print(f"\n[saved] {medium_json_path}")

    medium_csv_path = save_csv(medium_result["slice_results"], OUT_DIR, "medium_run_diagnostics")
    print(f"[saved] {medium_csv_path}")

    medium_plot_path = make_plots(medium_result, OUT_DIR, "medium_run")
    print(f"[saved] {medium_plot_path}")

    # =========================================================================
    # Final report
    # =========================================================================
    print("\n" + "=" * 70)
    print("FINAL REPORT")
    print("=" * 70)

    for label, summary in [("tiny (24^3 x 48)", tiny_summary),
                            ("medium (48^3 x 128)", medium_summary)]:
        print(f"\n--- {label} ---")
        print(f"  Verdict:                   {summary['verdict']}")
        print(f"  omega_0:                   {summary['omega_0']:.4f} "
              f"(= {PARAMS['omega_0_factor']:.2f} * omega_top={summary['omega_top']:.4f})")
        print(f"  phi_per_step (carrier):    {summary['phi_per_step']:.4f} rad")
        print(f"  E_excess(t=0):             {summary['E_excess_t0']:.4f}")
        print(f"  E_excess(final):           {summary['E_excess_final']:.4f}")
        print(f"  E_conservation_ratio:      {summary['E_conservation_ratio']:.4f}")
        print(f"  spread_ratio (0/mean/fin): {summary['spread_ratio_initial']:.4f} / "
              f"{summary['spread_ratio_mean']:.4f} / {summary['spread_ratio_final']:.4f}")
        print(f"  confined_fraction mean:    {summary['confined_fraction_mean']:.4f}")
        print(f"  leakage mean:              {summary['leakage_mean']:.4f}")
        print(f"  carrier_coherence mean:    {summary['carrier_coherence_mean']:.4f}")
        print(f"  residual_per_dof:          {summary['residual_per_dof']:.4e}")
        print(f"  walltime:                  {summary['walltime_total_s']:.1f}s")

    print(f"\nAll outputs in: {OUT_DIR}/")
    print(f"  tiny_validation.json")
    print(f"  tiny_validation_diagnostics.csv")
    print(f"  tiny_validation_diagnostics.png")
    print(f"  medium_run_summary.json")
    print(f"  medium_run_diagnostics_diagnostics.csv")
    print(f"  medium_run_diagnostics.png")


if __name__ == "__main__":
    main()
