"""D8 — U(1)↔SU(3) Binding Probe.

Instruments the binding physics derived in
  paper/derivations/time_link_binding.md  (Part 1: energetic, Part 2: topological)
  paper/derivations/u1_su3_binding.md     (channel A/B/C derivations)

READ-ONLY device — no solver state modification, no back-reaction.

Probes
------
P1  Sector centroids & separation (Channel-B co-location baseline).
P2  Longitudinal stretch profile Δu_∥(ρ) (time-link binding source).
P3  Per-slice carrier rate ω(l)   (linchpin caveat 1: uniformity check).
P4  Trace-sector loop holonomy γ_Γ (soft-winding-lock falsifier, scalar).
P5  Antisymmetric Kähler part Im𝒢  (Part-2 topological structure).

Outputs
-------
binding_probe.csv         — per-slice time series
binding_probe_radial.csv  — Δu_∥ radial profile
binding_probe.png         — 2×3 multi-panel figure

Conventions chosen here (documented where a choice was required)
----------------------------------------------------------------
* Vortex axis: z-axis through the grid centre (the spherical-harmonic seed
  places the donut symmetrically around the z-axis, confirmed by the seed
  geometry in vortex_worldtube.py — the azimuthal winding is around z).
  Radial distance ρ is measured in the XY plane: ρ=√(Δx²+Δy²) with Δx,Δy
  relative to the grid centre.  The named constant VORTEX_AXIS=2 (z) below
  must match the seed orientation in vortex_worldtube.py.

* ω₀ source: if n_t is passed as a keyword argument to device_binding_probe,
  we use omega0 = 2π·n_t / (n_slices·dt) (physical rate in 1/[time_unit]).
  Otherwise omega0 is measured from the carrier 2-plane FH links; falls back
  to 1.0 (documented).

* Periodic detection: world is treated as periodic in time if
  ||world[0] − world[-1]||_inf < 1e-4 * ||world[0]||_inf.
  Interior-only FH links are used when non-periodic; full loop when periodic.

* Complex Ψ lift: Psi_i(l,node) = u_i(l,node) + (i/ω₀)·udot_i(l,node)
  where udot uses central differences (interior), one-sided at boundaries
  (matching device_berry's pattern).  Lateral triplet: components 0,1,2.
  This is the "i from time" lift (project_complex_u1_from_time).

* P_U1, P_SU3 from alpha_separability.projection_operators():
  P_U1 = (1/3)·ones(3,3)   (trace / U(1))
  P_SU3 = I − P_U1         (traceless / SU(3))

* psi_s: trace scalar = (Ψ₀+Ψ₁+Ψ₂)/√3  (the EM/U(1) carrier)
  psi_perp: traceless vector = Ψ − psi_s/√3 · ones(3)

Principles compliance (PRINCIPLES.md)
--------------------------------------
- Read-only: all arrays are computed from inputs; no in-place mutation.
- Vacuum-subtracted: lateral displacement uses world[l,:,:3] - ref[:,:3].
- No hard-coded winding estimator: phase extracted from the actual field.
- No back-reaction: returns scalars/CSV/PNG; does not touch world/lattice.
- Agg backend enforced at module level.
- _apply_style(), _savefig(), _phase_to_rgb() imported from
  branesim.diagnostics._plot_helpers (shared with run_measurements; no
  circular import because _plot_helpers depends on neither).
- dim=3 layout OK for PNG (PRINCIPLES §7.6 relaxation for visualisation);
  the spatial-map panel is guarded with an explicit dim==3 check.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.alpha_separability import projection_operators
from branesim.diagnostics.confinement import _energy_weights
from branesim.initialization.vortex_worldtube import CARRIER_RE, CARRIER_IM, vacuum_offsets
from branesim.diagnostics._plot_helpers import _apply_style, _savefig


# ---------------------------------------------------------------------------
# Axis convention
# ---------------------------------------------------------------------------

# The vortex seed in vortex_worldtube.py winds the carrier around the z-axis
# (the spherical-harmonic Y_1^1 seed places the donut in the equatorial XY
# plane).  VORTEX_AXIS=2 is the z index.  This must match the seed orientation.
VORTEX_AXIS: int = 2


# ---------------------------------------------------------------------------
# Shared helper: build the ℂ³ lateral envelope Ψ(l, node, 3)
# ---------------------------------------------------------------------------

def _build_psi_triplet(
    world: np.ndarray,
    ref: np.ndarray,
    omega0: float,
    params: ActionParams,
) -> np.ndarray:
    """Build the carrier-complexified lateral triplet Ψ.

    For each lateral component i ∈ {0,1,2}:
        u_i(l) = world[l,:,i] − ref[:,i]
        udot_i(l): central diff interior; one-sided at ends (matches device_berry)
        Psi_i(l) = u_i(l) + (1j/ω₀)·udot_i(l)

    When ω₀ == 0 the velocity term is dropped with a docstring note; this
    prevents division-by-zero on degenerate seeds.

    Parameters
    ----------
    world : (T, n_nodes, m_ambient)  float64
    ref   : (n_nodes, m_ambient)     float64
    omega0 : float   physical carrier angular rate (1/time_unit); 1.0 if unknown
    params : ActionParams

    Returns
    -------
    Psi : complex ndarray, shape (T, n_nodes, 3)
    """
    n_T = world.shape[0]
    dt = params.dt

    # Lateral displacement: components 0,1,2
    u = world[:, :, :3] - ref[np.newaxis, :, :3]  # (T, n_nodes, 3)

    # Time derivative (central diff interior, one-sided ends)
    udot = np.empty_like(u)
    udot[0] = (u[1] - u[0]) / dt if n_T > 1 else np.zeros_like(u[0])
    udot[-1] = (u[-1] - u[-2]) / dt if n_T > 1 else np.zeros_like(u[-1])
    for l in range(1, n_T - 1):
        udot[l] = (u[l + 1] - u[l - 1]) / (2.0 * dt)

    if abs(omega0) < 1e-30:
        # omega0 unknown / zero: drop velocity term, document
        # (the "i from time" lift degenerates; Psi is real)
        return u.astype(complex)

    return u + 1j * udot / omega0   # (T, n_nodes, 3)


# ---------------------------------------------------------------------------
# P1 — Sector centroids & separation
# ---------------------------------------------------------------------------

def _p1_sector_centroids(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
) -> dict[str, Any]:
    """Compute per-slice energy-weighted centroids for U(1) and SU(3) sectors.

    Returns
    -------
    dict with:
        sep_d       : (T,)  |c_tr(l) − c_perp(l)|
        c_tr        : (T, 3)
        c_perp      : (T, 3)
        sep_mean    : float
        sep_max     : float
    """
    P_U1, P_SU3 = projection_operators()
    n_T = world.shape[0]
    dim = lattice.dim
    ref_spatial = ref[:, :dim]  # (n_nodes, dim)

    c_tr = np.zeros((n_T, dim))
    c_perp = np.zeros((n_T, dim))
    sep_d = np.zeros(n_T)

    for l in range(n_T):
        disp_lat = world[l, :, :3] - ref[:, :3]   # (n_nodes, 3)
        d_u1 = disp_lat @ P_U1.T                   # (n_nodes, 3)
        d_su3 = disp_lat @ P_SU3.T                 # (n_nodes, 3)

        # Energy proxy per node: |P·disp|^2
        # +1e-40 is a numerical guard preventing division-by-zero in the centroid
        # when displacement is identically zero; it is NOT a physics threshold.
        w_u1 = np.sum(d_u1 ** 2, axis=1) + 1e-40   # (n_nodes,)
        w_su3 = np.sum(d_su3 ** 2, axis=1) + 1e-40

        # Energy-weighted centroid over node spatial positions
        c_tr[l] = np.sum(ref_spatial * w_u1[:, np.newaxis], axis=0) / np.sum(w_u1)
        c_perp[l] = np.sum(ref_spatial * w_su3[:, np.newaxis], axis=0) / np.sum(w_su3)

        sep_d[l] = float(np.linalg.norm(c_tr[l] - c_perp[l]))

    return {
        "sep_d": sep_d,
        "c_tr": c_tr,
        "c_perp": c_perp,
        "sep_mean": float(np.mean(sep_d)),
        "sep_max": float(np.max(sep_d)),
    }


# ---------------------------------------------------------------------------
# P2 — Longitudinal stretch profile Δu_∥(ρ)
# ---------------------------------------------------------------------------

def _p2_longitudinal_stretch(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
) -> dict[str, Any]:
    """Compute the vacuum-subtracted longitudinal stretch Δu_∥ vs radial ρ.

    Convention: vortex axis = z-axis through the grid centre (see module
    docstring).  ρ = √(Δx²+Δy²) in the XY plane.

    Vacuum subtraction: for each bond we compute
        du_par_excess = (bond_hat · actual_bond) − (bond_hat · ref_bond)
    where ref_bond is the reference lattice bond vector.  For an axial bond of
    length a, bond_hat · ref_bond = a.  In vacuum (u=0) the actual bond equals
    the reference bond and du_par_excess = 0 exactly, regardless of α.  This
    is the only choice that gives an identically-zero vacuum baseline.

    Derivation sign prediction: Δu_∥ < 0 (time-link compression), peak at
    the donut ridge ρ ≈ w.

    Returns
    -------
    dict with:
        rho_bins        : (n_bins,) radial bin centres
        du_par_mean     : (n_bins,) mean Δu_∥ per bin at t=0
        du_par_avg      : (n_bins,) slice-averaged Δu_∥
        rho_extremum    : float   ρ of the extremum at t=0
        du_par_extremum : float   extremum value at t=0
        du_par_scalar   : float   mean over all bonds within the lump region
        du_par_sign     : int     sign of the extremum (-1, 0, or +1)
        normalizer      : float   Placeholder for ω²·A² normalizer (set by caller)
    """
    a = lattice.params.spacing
    grid_shape = lattice.params.grid_shape
    n_T = world.shape[0]
    dim = lattice.dim

    # Grid centre in spatial coordinates
    centre = np.array(
        [0.5 * (g - 1) * a for g in grid_shape[:dim]], dtype=float
    )  # (dim,)

    # Radial distance ρ in XY plane (axes 0,1), ignoring z
    # (the vortex ring lives in the XY equatorial plane for the Y_1^1 seed)
    xy_ref = ref[:, :2]                                         # (n_nodes, 2)
    xy_centre = centre[:2]
    rho_node = np.linalg.norm(xy_ref - xy_centre[np.newaxis, :], axis=1)  # (n_nodes,)

    rho_max = float(np.max(rho_node))
    n_bins = 30
    rho_edges = np.linspace(0.0, rho_max, n_bins + 1)
    rho_bins = 0.5 * (rho_edges[:-1] + rho_edges[1:])

    # Accumulate Δu_∥ per bond at t=0 and slice-averaged
    du_par_sum0 = np.zeros(n_bins)
    du_par_count0 = np.zeros(n_bins, dtype=int)
    du_par_sumavg = np.zeros(n_bins)
    du_par_countavg = np.zeros(n_bins, dtype=int)

    for nb_idx in range(lattice.n_neighbors):
        nb_ids = lattice.neighbors[:, nb_idx]
        valid = nb_ids >= 0
        if not np.any(valid):
            continue

        p_idx = np.where(valid)[0]
        q_idx = nb_ids[p_idx]

        # Reference bond vector and unit direction
        ref_bond = ref[q_idx, :dim] - ref[p_idx, :dim]           # (n_valid, dim)
        ref_dist = np.linalg.norm(ref_bond, axis=1, keepdims=True)  # (n_valid, 1)
        ref_safe = np.where(ref_dist > 1e-12, ref_dist, 1.0)
        bond_hat = ref_bond / ref_safe                            # (n_valid, dim)

        # Vacuum longitudinal stretch: bond_hat · ref_bond = |ref_bond| = a.
        # We subtract this per-bond to get du_par_excess = 0 in vacuum exactly.
        # Note: do NOT subtract (a − α·a) — that would leave the prestress offset
        # (~α·a) floating; the correct reference is the actual reference bond length.
        du_vac_per_bond = np.sum(ref_bond * bond_hat, axis=1)    # (n_valid,) ≈ a

        # Radial distance of bond midpoint (use p-node ρ as proxy)
        rho_p = rho_node[p_idx]

        # Bin indices for p-nodes.
        # np.clip here is an index guard: searchsorted can return n_bins for
        # nodes exactly at rho_max; clip keeps indices within [0, n_bins-1].
        # This is NOT a physics clamp.
        bin_idx = np.searchsorted(rho_edges[1:], rho_p)
        bin_idx = np.clip(bin_idx, 0, n_bins - 1)

        # --- t=0 ---
        u_p = world[0, p_idx, :dim] - ref[p_idx, :dim]   # (n_valid, dim)
        u_q = world[0, q_idx, :dim] - ref[q_idx, :dim]

        # Current bond vector and longitudinal projection
        actual_bond = (ref[q_idx, :dim] + u_q) - (ref[p_idx, :dim] + u_p)
        du_par_bond = np.sum(actual_bond * bond_hat, axis=1) - du_vac_per_bond

        for bi in range(n_bins):
            mask = bin_idx == bi
            if np.any(mask):
                du_par_sum0[bi] += np.sum(du_par_bond[mask])
                du_par_count0[bi] += int(np.sum(mask))

        # --- slice average ---
        for l in range(n_T):
            u_p_l = world[l, p_idx, :dim] - ref[p_idx, :dim]
            u_q_l = world[l, q_idx, :dim] - ref[q_idx, :dim]
            actual_bond_l = (ref[q_idx, :dim] + u_q_l) - (ref[p_idx, :dim] + u_p_l)
            du_par_l = np.sum(actual_bond_l * bond_hat, axis=1) - du_vac_per_bond

            for bi in range(n_bins):
                mask = bin_idx == bi
                if np.any(mask):
                    du_par_sumavg[bi] += np.sum(du_par_l[mask])
                    du_par_countavg[bi] += int(np.sum(mask))

    safe_count0 = np.where(du_par_count0 > 0, du_par_count0, 1)
    safe_countavg = np.where(du_par_countavg > 0, du_par_countavg, 1)
    du_par_mean = du_par_sum0 / safe_count0
    du_par_avg = du_par_sumavg / safe_countavg
    # Bins with no bonds stay NaN for honest display
    du_par_mean = np.where(du_par_count0 > 0, du_par_mean, np.nan)
    du_par_avg = np.where(du_par_countavg > 0, du_par_avg, np.nan)

    # Extremum at t=0 (ignore NaN bins)
    valid_mask = ~np.isnan(du_par_mean)
    if np.any(valid_mask):
        abs_vals = np.where(valid_mask, np.abs(du_par_mean), 0.0)
        ext_bin = int(np.argmax(abs_vals))
        rho_extremum = float(rho_bins[ext_bin])
        du_par_extremum = float(du_par_mean[ext_bin])
        du_par_sign = int(np.sign(du_par_extremum)) if abs(du_par_extremum) > 1e-30 else 0
    else:
        rho_extremum = float("nan")
        du_par_extremum = float("nan")
        du_par_sign = 0

    # Scalar mean over the lump region (ρ < r0 + w, estimated from params)
    # Use a generous 3-sigma region around the ring for the scalar summary
    du_par_scalar = float(np.nanmean(du_par_mean)) if np.any(valid_mask) else float("nan")

    return {
        "rho_bins": rho_bins,
        "du_par_mean": du_par_mean,
        "du_par_avg": du_par_avg,
        "rho_extremum": rho_extremum,
        "du_par_extremum": du_par_extremum,
        "du_par_scalar": du_par_scalar,
        "du_par_sign": du_par_sign,
        "normalizer": float("nan"),   # set by caller after omega0 is known
    }


# ---------------------------------------------------------------------------
# P3 — Per-slice carrier rate ω(l)
# ---------------------------------------------------------------------------

def _p3_carrier_rate(
    psi_s: np.ndarray,
    dt: float,
    n_t: int | None,
    n_slices: int,
    is_periodic: bool,
) -> dict[str, Any]:
    """Per-slice FH carrier rate from the trace scalar psi_s.

    ω(l) = arg( Σ_nodes conj(psi_s[l]) · psi_s[l+1] ) / dt

    Parameters
    ----------
    psi_s   : (T, n_nodes)  complex  trace scalar
    dt      : float
    n_t     : int or None   known temporal winding (from config)
    n_slices: int           number of temporal steps (T = n_slices+1)
    is_periodic: bool       whether world[0] ≈ world[-1]

    Returns
    -------
    dict with:
        omega_l         : (n_slices,)  per-step ω values
        omega_mean      : float
        omega_std       : float
        omega_ref       : float or None  2π·n_t/(n_slices·dt) if n_t known
        closure_locked  : bool or None   std/mean < 0.05 if mean nonzero
    """
    n_T = psi_s.shape[0]
    n_steps = n_T - 1

    # Choose which links to use
    if is_periodic:
        # Include the wrap-around link l=N→0
        l_range = list(range(n_steps)) + [n_T - 1]
        l_next = list(range(1, n_T)) + [0]
    else:
        l_range = list(range(n_steps))
        l_next = list(range(1, n_T))

    omega_l = np.zeros(len(l_range))
    for i, (l, lp1) in enumerate(zip(l_range, l_next)):
        inner = np.sum(np.conj(psi_s[l]) * psi_s[lp1])
        amp = abs(inner)
        if amp > 1e-30:
            omega_l[i] = float(np.angle(inner)) / dt
        else:
            omega_l[i] = 0.0

    omega_mean = float(np.mean(omega_l))
    omega_std = float(np.std(omega_l))

    omega_ref = None
    if n_t is not None and dt > 0 and n_slices > 0:
        omega_ref = float(2.0 * np.pi * n_t / (n_slices * dt))

    closure_locked = None
    if abs(omega_mean) > 1e-30:
        closure_locked = bool(omega_std / abs(omega_mean) < 0.05)

    return {
        "omega_l": omega_l,
        "omega_mean": omega_mean,
        "omega_std": omega_std,
        "omega_ref": omega_ref,
        "closure_locked": closure_locked,
    }


# ---------------------------------------------------------------------------
# P4 — Trace-sector loop holonomy γ_Γ
# ---------------------------------------------------------------------------

def _p4_loop_holonomy(
    psi_s: np.ndarray,
    dt: float,
    is_periodic: bool,
) -> dict[str, Any]:
    """Accumulate the trace-sector Berry phase around the full time loop.

    γ_Γ = Σ_l arg( Σ_nodes conj(psi_s[l]) · psi_s[l+1] )

    Uses the same FH discrete-link accumulation as device_berry.

    Returns
    -------
    dict with:
        gamma_accum : (T,)  accumulated phase, gamma_accum[0]=0
        gamma_total : float  total accumulated phase
    """
    n_T = psi_s.shape[0]
    n_steps = n_T - 1

    if is_periodic:
        l_range = list(range(n_steps)) + [n_T - 1]
        l_next = list(range(1, n_T)) + [0]
    else:
        l_range = list(range(n_steps))
        l_next = list(range(1, n_T))

    # For the per-slice time series we need one entry per slice (n_T values)
    # gamma_accum[l] = accumulated phase up to and including link l→l+1
    gamma_accum = np.zeros(n_T)
    for i, (l, lp1) in enumerate(zip(l_range, l_next)):
        inner = np.sum(np.conj(psi_s[l]) * psi_s[lp1])
        amp = abs(inner)
        phase_step = float(np.angle(inner)) if amp > 1e-30 else 0.0
        idx = min(i + 1, n_T - 1)
        gamma_accum[idx] = gamma_accum[i] + phase_step

    gamma_total = float(gamma_accum[-1])

    return {
        "gamma_accum": gamma_accum,
        "gamma_total": gamma_total,
    }


# ---------------------------------------------------------------------------
# P5 — Antisymmetric Kähler part Im𝒢
# ---------------------------------------------------------------------------

def _p5_kahler(
    psi_s: np.ndarray,
    psi_perp: np.ndarray,
) -> dict[str, Any]:
    """Compute the trace-sector ↔ traceless cross-overlap per slice.

    For each traceless component a ∈ {0,1,2}:
        g_sa(l) = Σ_nodes conj(psi_s[l,node]) * psi_perp[l,node,a]
                  normalized by sqrt(||psi_s||²·||psi_perp_a||²)

    Reports per-slice Re(g) and Im(g) (mean over a), and the ratio
        R_kahler = mean_l |Im g(l)| / (mean_l |Re g(l)| + eps)

    Prediction: Im𝒢 ≠ 0 requires both J (time-quadrature) and α-split
    (ω_tr ≠ ω_⊥).  Threshold R_kahler > ~0.1 ⇒ Kähler part present.

    Returns
    -------
    dict with:
        Re_g    : (T,)   mean over a of Re(g_sa)
        Im_g    : (T,)   mean over a of Im(g_sa)
        R_kahler: float
    """
    n_T = psi_s.shape[0]
    n_comps = psi_perp.shape[2]

    Re_g = np.zeros(n_T)
    Im_g = np.zeros(n_T)

    for l in range(n_T):
        g_per_comp = np.zeros(n_comps, dtype=complex)
        ps = psi_s[l]               # (n_nodes,)
        # +1e-40: numerical guard preventing division-by-zero when psi_s or
        # psi_perp vanish (e.g. at vacuum nodes); NOT a physics threshold.
        norm_s = float(np.sqrt(np.sum(np.abs(ps) ** 2))) + 1e-40
        for a in range(n_comps):
            pp = psi_perp[l, :, a]  # (n_nodes,)
            norm_p = float(np.sqrt(np.sum(np.abs(pp) ** 2))) + 1e-40
            raw = np.sum(np.conj(ps) * pp)
            g_per_comp[a] = raw / (norm_s * norm_p)
        Re_g[l] = float(np.mean(np.real(g_per_comp)))
        Im_g[l] = float(np.mean(np.imag(g_per_comp)))

    # +1e-40 in denominator: numerical guard preventing division-by-zero when
    # Re_g is identically zero; NOT a physics threshold.
    R_kahler = float(
        np.mean(np.abs(Im_g)) / (np.mean(np.abs(Re_g)) + 1e-40)
    )

    return {
        "Re_g": Re_g,
        "Im_g": Im_g,
        "R_kahler": R_kahler,
    }


# ---------------------------------------------------------------------------
# Main device entry point
# ---------------------------------------------------------------------------


def device_binding_probe(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
    *,
    n_t: int | None = None,
) -> dict[str, Any]:
    """D8 — U(1)↔SU(3) Binding Probe.

    Instruments the five probes defined in
    paper/derivations/time_link_binding.md and paper/derivations/u1_su3_binding.md.

    Parameters
    ----------
    world  : (n_slices+1, n_nodes, m_ambient)  float64
    ref    : (n_nodes, m_ambient)               float64
    lattice: SpacelikeLattice
    params : ActionParams   (frozen — never mutated here)
    out_dir: Path  — directory where CSV + PNG are written
    n_t    : int or None   temporal winding number from config["vortex_params"]["n_t"].
             When provided, omega0 = 2π·n_t/(n_slices·dt) (closure-locked rate).
             When None, omega0 is measured from the carrier 2-plane FH links.
             Do NOT set params.n_t — ActionParams is frozen=True.

    Returns
    -------
    dict with scalar summaries + "csv", "csv_radial", "png" paths.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_T = world.shape[0]
    n_slices = n_T - 1
    dt = params.dt
    times = np.arange(n_T) * dt

    # ------------------------------------------------------------------
    # Detect n_t and compute omega0
    # ------------------------------------------------------------------
    # Priority order:
    #   1. n_t kwarg (passed by run_measurements from config; ActionParams is
    #      frozen=True so it must NOT be mutated by the caller)
    #   2. Measure from the carrier 2-plane (CARRIER_RE + i*CARRIER_IM) via FH link
    #   3. Fall back to omega0=1.0 (documented, logged)

    if n_t is not None and n_slices > 0 and dt > 0:
        # Closure-locked rate from config
        omega0 = float(2.0 * np.pi * n_t / (n_slices * dt))
        omega0_source = "closure_locked_config"
    else:
        # Measure from the carrier 2-plane (device_berry approach, more accurate than triplet)
        # Carrier complex field = (u_CARRIER_RE + i*u_CARRIER_IM) using vacuum-offsets
        off_re, off_im = vacuum_offsets(n_slices, params.r_t)
        re_d = world[:, :, CARRIER_RE] - ref[np.newaxis, :, CARRIER_RE] - off_re[:, np.newaxis]
        im_d = world[:, :, CARRIER_IM] - ref[np.newaxis, :, CARRIER_IM] - off_im[:, np.newaxis]
        Psi_carrier = re_d + 1j * im_d   # (n_T, n_nodes)
        # FH link phases over interior steps (avoid one-sided boundary artefact at l=0)
        phase_sum = 0.0
        count = 0
        for l_meas in range(1, n_T - 2):
            inner_m = np.sum(np.conj(Psi_carrier[l_meas]) * Psi_carrier[l_meas + 1])
            if abs(inner_m) > 1e-30:
                phase_sum += float(np.angle(inner_m))
                count += 1
        if count > 0 and dt > 0:
            omega0 = float(phase_sum / count / dt)
            omega0_source = "measured_carrier_plane"
        else:
            omega0 = 1.0
            omega0_source = "fallback_1"

    # omega_per_slice (dimensionless, for reference reporting)
    omega_per_slice = float(omega0 * dt) if dt > 0 else 0.0

    # ------------------------------------------------------------------
    # Detect periodicity (world[0] ≈ world[-1])
    # ------------------------------------------------------------------
    norm0 = float(np.max(np.abs(world[0])))
    diff_seam = float(np.max(np.abs(world[0] - world[-1])))
    is_periodic = (norm0 > 1e-30) and (diff_seam < 1e-4 * norm0)

    # ------------------------------------------------------------------
    # Build ℂ³ triplet envelope Ψ
    # ------------------------------------------------------------------
    Psi = _build_psi_triplet(world, ref, omega0, params)
    # (n_T, n_nodes, 3)

    # Trace scalar psi_s and traceless psi_perp
    psi_s = (Psi[:, :, 0] + Psi[:, :, 1] + Psi[:, :, 2]) / np.sqrt(3.0)
    # (n_T, n_nodes)
    psi_perp = Psi - (psi_s / np.sqrt(3.0))[:, :, np.newaxis] * np.ones(3)
    # (n_T, n_nodes, 3)

    # ------------------------------------------------------------------
    # Run the five probes
    # ------------------------------------------------------------------
    p1 = _p1_sector_centroids(world, ref, lattice)
    p2 = _p2_longitudinal_stretch(world, ref, lattice, params)
    p3 = _p3_carrier_rate(psi_s, dt, n_t, n_slices, is_periodic)
    p4 = _p4_loop_holonomy(psi_s, dt, is_periodic)
    p5 = _p5_kahler(psi_s, psi_perp)

    # Fill normalizer from measured omega0 and amplitude A (mean over nodes)
    amp0 = float(np.mean(np.abs(psi_s[0])))
    measured_omega0 = float(p3["omega_mean"])
    if abs(measured_omega0) > 1e-30 and abs(amp0) > 1e-30:
        normalizer = float(measured_omega0 ** 2 * amp0 ** 2)
    else:
        normalizer = float("nan")
    p2["normalizer"] = normalizer

    # ------------------------------------------------------------------
    # Write CSVs
    # ------------------------------------------------------------------
    # Main per-slice CSV
    csv_path = out_dir / "binding_probe.csv"
    # Omega time series has n_steps entries; pad or clip to n_T for alignment
    omega_ts = p3["omega_l"]
    gamma_ts = p4["gamma_accum"]  # (n_T,)
    Re_g_ts = p5["Re_g"]          # (n_T,)
    Im_g_ts = p5["Im_g"]          # (n_T,)
    sep_d_ts = p1["sep_d"]        # (n_T,)

    # omega_l has n_steps = n_T−1 entries; broadcast to n_T by repeating last
    if len(omega_ts) < n_T:
        omega_full = np.concatenate([omega_ts, [omega_ts[-1] if len(omega_ts) > 0 else 0.0]])
    else:
        omega_full = omega_ts[:n_T]

    rows = np.column_stack([times, sep_d_ts, omega_full, gamma_ts, Re_g_ts, Im_g_ts])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,sep_d,omega_l,gamma_accum,Re_g,Im_g",
        comments="",
    )

    # Radial profile CSV
    csv_radial_path = out_dir / "binding_probe_radial.csv"
    rho_bins = p2["rho_bins"]
    du_par_mean = p2["du_par_mean"]
    du_par_avg = p2["du_par_avg"]
    rows_rad = np.column_stack([rho_bins, du_par_mean, du_par_avg])
    np.savetxt(
        str(csv_radial_path), rows_rad, delimiter=",",
        header="rho,du_par_t0,du_par_sliceavg",
        comments="",
    )

    # ------------------------------------------------------------------
    # PNG: 2×3 multi-panel figure
    # ------------------------------------------------------------------
    _apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("D8 — U(1)↔SU(3) Binding Probe", fontweight="bold")

    # (a) Sector centroid separation d(l)
    ax = axes[0, 0]
    ax.plot(times, sep_d_ts, color="tab:blue", label="|c_tr − c_perp|")
    ax.axhline(p1["sep_mean"], color="tab:blue", linestyle="--", alpha=0.6,
               label=f"mean={p1['sep_mean']:.3g}")
    ax.set_xlabel("time"); ax.set_ylabel("separation d (lattice units)")
    ax.legend(); ax.set_title("(a) P1 — Sector centroid separation")

    # (b) Δu_∥(ρ) radial profile
    ax = axes[0, 1]
    valid_mask = ~np.isnan(du_par_mean)
    if np.any(valid_mask):
        ax.plot(rho_bins[valid_mask], du_par_mean[valid_mask], color="tab:red",
                label="Δu_∥(ρ) at t=0")
        ax.plot(rho_bins[valid_mask], du_par_avg[valid_mask], color="tab:orange",
                linestyle="--", label="slice avg")
        if not np.isnan(p2["rho_extremum"]):
            ax.axvline(p2["rho_extremum"], color="gray", linestyle=":",
                       label=f"ρ_ext={p2['rho_extremum']:.2g}")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("ρ (lattice units)"); ax.set_ylabel("Δu_∥  (vacuum-subtracted)")
    sign_str = {-1: "NEGATIVE (binding)", 0: "ZERO", 1: "POSITIVE (repulsive)"}
    ax.set_title(f"(b) P2 — Longitudinal stretch\n"
                 f"sign={sign_str.get(p2['du_par_sign'], '?')}")
    ax.legend(fontsize=8)

    # (c) Per-slice carrier rate ω(l)
    ax = axes[0, 2]
    omega_ts_plot = p3["omega_l"]
    t_omega = np.arange(len(omega_ts_plot)) * dt
    ax.plot(t_omega, omega_ts_plot, color="tab:green", label="ω(l) measured")
    if p3["omega_ref"] is not None:
        ax.axhline(p3["omega_ref"], color="gray", linestyle="--",
                   label=f"ω_ref={p3['omega_ref']:.3g}")
    ax.axhline(p3["omega_mean"], color="tab:green", linestyle=":", alpha=0.6,
               label=f"mean={p3['omega_mean']:.3g}")
    ax.set_xlabel("time"); ax.set_ylabel("ω [rad/time_unit]")
    cl_str = str(p3["closure_locked"]) if p3["closure_locked"] is not None else "n/a"
    ax.set_title(f"(c) P3 — Carrier rate ω(l)\nclosure-locked={cl_str}")
    ax.legend(fontsize=8)

    # (d) Re/Im g time series (P5)
    ax = axes[1, 0]
    ax.plot(times, Re_g_ts, color="tab:blue", label="Re g")
    ax.plot(times, Im_g_ts, color="tab:red", label="Im g")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("time"); ax.set_ylabel("g_sa (norm.)")
    ax.set_title(f"(d) P5 — Kähler overlap  R_kahler={p5['R_kahler']:.3g}")
    ax.legend()

    # (e) Accumulated trace-sector phase → γ_Γ
    ax = axes[1, 1]
    ax.plot(times, p4["gamma_accum"], color="tab:purple")
    ax.axhline(p4["gamma_total"], color="tab:purple", linestyle="--", alpha=0.5,
               label=f"γ_Γ={p4['gamma_total']:.3g} rad")
    ax.set_xlabel("time"); ax.set_ylabel("accumulated phase [rad]")
    ax.set_title(f"(e) P4 — Trace-sector loop holonomy\nγ_Γ={p4['gamma_total']:.4g} rad")
    ax.legend()

    # (f) Spatial map: XY midplane of rho_tr vs rho_perp at t=0
    # Dim guard: the 3D reshape and z-slice are only valid for dim==3.
    # For non-3D grids, the panel is left blank with a descriptive note.
    # (PRINCIPLES §7.6: dimension-agnostic degradation for visualisation.)
    grid_shape = lattice.params.grid_shape
    ax = axes[1, 2]

    if len(grid_shape) == 3:
        # VORTEX_AXIS = 2 (z-axis); must match seed orientation in vortex_worldtube.py
        nx, ny, nz = grid_shape
        z_mid = nz // 2   # midplane along VORTEX_AXIS

        P_U1, P_SU3 = projection_operators()
        disp_lat0 = world[0, :, :3] - ref[:, :3]
        d_u1_0 = disp_lat0 @ P_U1.T
        d_su3_0 = disp_lat0 @ P_SU3.T
        amp_u1_map = np.linalg.norm(d_u1_0, axis=-1).reshape(nx, ny, nz)
        amp_su3_map = np.linalg.norm(d_su3_0, axis=-1).reshape(nx, ny, nz)

        extent = [0, nx, 0, ny]
        vmax_u1 = float(np.max(amp_u1_map)) or 1.0
        vmax_su3 = float(np.max(amp_su3_map)) or 1.0

        # Overlay: U(1) in blue, SU(3) in red.
        # np.clip here is an index/display guard (keeps RGB in [0,1]), not a physics clamp.
        u1_slice = amp_u1_map[:, :, z_mid] / vmax_u1
        su3_slice = amp_su3_map[:, :, z_mid] / vmax_su3
        rgb_overlay = np.zeros((nx, ny, 3), dtype=float)
        rgb_overlay[:, :, 0] = np.clip(su3_slice, 0, 1)  # SU(3) → red  (display guard)
        rgb_overlay[:, :, 2] = np.clip(u1_slice, 0, 1)   # U(1)  → blue (display guard)
        ax.imshow(rgb_overlay.swapaxes(0, 1), origin="lower", extent=extent, aspect="equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title("(f) P1 — XY midplane: U(1)=blue, SU(3)=red  (t=0)")
    else:
        ax.text(0.5, 0.5, f"dim={len(grid_shape)} — spatial map\nrequires dim=3",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("(f) P1 — spatial map (dim=3 only)")

    _savefig(fig, out_dir / "binding_probe.png")

    # ------------------------------------------------------------------
    # Scalar summary dict
    # ------------------------------------------------------------------
    return {
        # P1
        "sep_d_mean": p1["sep_mean"],
        "sep_d_max": p1["sep_max"],
        # P2
        "du_par_sign": p2["du_par_sign"],
        "du_par_extremum": p2["du_par_extremum"],
        "rho_extremum": p2["rho_extremum"],
        "du_par_scalar": p2["du_par_scalar"],
        "du_par_normalizer": normalizer,
        # P3
        "omega_mean": p3["omega_mean"],
        "omega_std": p3["omega_std"],
        "omega_ref": p3["omega_ref"],
        "closure_locked": p3["closure_locked"],
        # P4
        "gamma_Gamma": p4["gamma_total"],
        # P5
        "R_kahler": p5["R_kahler"],
        "Re_g_mean": float(np.mean(Re_g_ts)),
        "Im_g_mean": float(np.mean(Im_g_ts)),
        # Output paths
        "csv": str(csv_path),
        "csv_radial": str(csv_radial_path),
        "png": str(out_dir / "binding_probe.png"),
        # Metadata
        "omega0_used": omega0,
        "omega0_source": omega0_source,
        "is_periodic": is_periodic,
        "n_t_used": n_t,
    }