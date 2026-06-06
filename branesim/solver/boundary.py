"""Boundary-condition schemes for the two-time BVP.

ARCHITECTURE.md §2 design decision D2 (resolved 2026-05-30).

Two schemes are provided:

1. **Dirichlet two-time** — fix slices l=0 and l=N to prescribed data.
   Known to be ill-posed at resonant N (determinant 2i·sin(Nθ(k)) vanishes).
   Condition estimate: max_k 1/|sin(Nθ(k))|.  Used as a negative control in
   the validation suite.

2. **Chiral Cauchy BC** (κ bounded, N-independent) — the well-posed scheme.
   Verdict (a): the correct chiral BC is *two adjacent past slices* (R0, R1)
   marched forward via the Störmer–Verlet stencil.  No FFT, no characteristic
   projection, no future condition.

   Physics:
     - A real field requires a+(k) = conj(a+(-k)); a single past slice A^0
       cannot determine the temporal characteristic split without A^1.
     - The Verlet recurrence   R^{l+1} = 2·cosθ·R^l − R^{l−1}   is the
       discrete forward march from Cauchy data (R0, R1); it is real by
       construction (real stencil + real data).
     - Forward chirality: march l=0 → N   (matter / retarded worldtube).
     - Backward chirality: reverse stencil R^{l-1} = 2·cosθ·R^l − R^{l+1}
       applied from (RN := R0, R_{N-1} := R1), producing earlier slices
       (antimatter worldtube).
     - DC modes (θ → 0) handled automatically: the recurrence with cosθ=1
       becomes a linear ramp, which is the exact solution.
     - Condition estimate: the per-mode Cauchy propagator maps two past slices
       to the full column; its 2-norm condition number is O(1) and N-independent
       (bounded by 1/sin(θ) for the 2×2 companion matrix, which is O(1) for
       propagating modes θ ∈ (0,π)).  Derivation + numerical certificate:
       OPEN_PROBLEMS.md A2 (verdict a) and ARCHITECTURE.md §2 D2.

   The old FFT-based ``apply_chiral`` (which set a₊ := A0_k per spatial bin and
   multiplied by e^{−iθN}) is deleted per §7.5; it violated reality by conflating
   spatial-FFT amplitude with the temporal characteristic amplitude.

This module does not modify solver state; it returns modified copies of the
world-volume's boundary rows and reports condition estimates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from branesim.core.conventions import ActionParams, LatticeParams, d_of_k_eigenvalues
from branesim.core.action import spacelike_force


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

THETA_TOL: float = 1e-6   # kept for diagnostic use only; not a branch point


# ---------------------------------------------------------------------------
# BoundaryCondition dataclasses
# ---------------------------------------------------------------------------


@dataclass
class DirichletBC:
    """Fix slices l=0 and l=N to prescribed arrays.

    Attributes
    ----------
    R0 : ndarray, shape (n_nodes, m_ambient)
        Prescribed configuration on the past boundary slice (l=0).
    RN : ndarray, shape (n_nodes, m_ambient)
        Prescribed configuration on the future boundary slice (l=N).
    """

    R0: np.ndarray
    RN: np.ndarray

    def condition_estimate(
        self,
        lattice_params: LatticeParams,
        action_params: ActionParams,
    ) -> float:
        """Estimate the operator condition number for Dirichlet two-time BVP.

        The modal operator determinant is 2i·sin(Nθ(k)).  The condition number
        scales as max_k 1/|sin(Nθ(k))|, which diverges at resonances Nθ=mπ.

        Returns
        -------
        float
            max_k  1 / max(|sin(N·θ(k))|, eps)
            A finite large value is returned if a mode is exactly resonant.
        """
        theta_vals = _theta_for_all_modes(lattice_params, action_params)
        N = action_params.n_slices
        sin_vals = np.abs(np.sin(N * theta_vals))
        eps = 1e-14
        return float(np.max(1.0 / np.maximum(sin_vals, eps)))


@dataclass
class ChiralBC:
    """Chiral Cauchy BC — two adjacent past slices, marched forward (verdict a).

    The chiral "solve" is the forward (or backward) Verlet march from the two
    past slices.  No FFT, no characteristic projection, no future condition is
    needed.  Reality is automatic: real stencil + real data → real world-volume.

    Attributes
    ----------
    R0 : ndarray, shape (n_nodes, m_ambient)
        Configuration on past slice l=0.
    R1 : ndarray, shape (n_nodes, m_ambient)
        Configuration on past slice l=1.
    chirality : {"forward", "backward"}
        ``"forward"`` produces the matter worldtube (march l=0..N).
        ``"backward"`` produces the antimatter worldtube: the reverse stencil
        R^{l-1} = 2·cosθ·R^l − R^{l+1} is used from (RN:=R0, R_{N-1}:=R1).
    theta_tol : float
        Diagnostic-only: used by ``condition_estimate`` to skip near-DC modes
        (θ ≤ theta_tol) when reporting the worst-mode condition number. It does
        NOT affect the chiral solve (the march), which handles θ→0 exactly.
    """

    R0: np.ndarray
    R1: np.ndarray
    chirality: Literal["forward", "backward"] = "forward"
    theta_tol: float = THETA_TOL

    def condition_estimate(
        self,
        lattice_params: LatticeParams,
        action_params: ActionParams,
    ) -> float:
        """Condition estimate for the chiral Cauchy scheme.

        The per-mode Cauchy propagator maps (R0_k, R1_k) to all N+1 slices
        via a 2×2 companion-matrix recurrence.  Its 2-norm condition is
        bounded by a constant that grows at most as 1/sin(θ_min) for the
        worst propagating mode — O(1/θ) for small θ, N-independent.  For
        the default lattice regime (θ ∈ (0.01, π)) the bound is < 100.

        Returns
        -------
        float
            A bounded N-independent condition estimate computed from the
            per-mode companion-matrix 2-norm.  Always finite.
        """
        theta_vals = _theta_for_all_modes(lattice_params, action_params)
        # Per-mode bound: cond of 2x2 companion ~ 1/sin(theta) for theta in (0,pi).
        # Cap at a large but finite value for DC/near-DC modes.
        prop = theta_vals > self.theta_tol
        if not np.any(prop):
            return 1.0
        sin_vals = np.sin(theta_vals[prop])
        # 1e-10 is a divide-by-zero floor for this reporting value only (not a
        # physics threshold); near-DC modes are already excluded by theta_tol.
        return float(np.max(1.0 / np.maximum(np.abs(sin_vals), 1e-10)))


# ---------------------------------------------------------------------------
# Internal helpers: θ(k) for all Fourier modes of a periodic grid
# ---------------------------------------------------------------------------


def _theta_for_all_modes(
    lattice_params: LatticeParams,
    action_params: ActionParams,
) -> np.ndarray:
    """Compute θ(k) = arccos(1 − Δt²·ω²(k)/2) for every Fourier mode.

    For a periodic lattice with grid_shape and spacing ``a``, the Fourier
    wavevectors are k_i = 2π·n_i / (N_i·a) for n_i ∈ {0,…,N_i−1}.

    The eigenvectors are Cartesian unit vectors (backbone #15), so each
    polarization axis α gives one independent mode family.  For each
    combination of (polarization axis α, wavevector k) the eigenfrequency
    ω_α²(k) = D_α(k) is the α-th eigenvalue from d_of_k_eigenvalues.

    Returns θ clipped to [0, π].
    """
    grid_shape = lattice_params.grid_shape
    dim = lattice_params.dim
    dt = action_params.dt
    a = lattice_params.spacing

    # Build the full Fourier wavevector grid (C-order ravel):
    # shape (n_nodes, dim) with one row per (n_0,...,n_{d-1}) combination.
    ranges = [np.arange(n) for n in grid_shape]
    grids = np.meshgrid(*ranges, indexing="ij")
    multi = np.stack([g.ravel() for g in grids], axis=1)  # (n_nodes, dim)

    # Physical wavevectors
    k_grid = 2.0 * np.pi * multi / (np.array(grid_shape) * a)  # (n_nodes, dim)

    # Gather ω²(k) for all polarization axes × wavevectors
    # Each row k gives dim eigenvalues (one per polarization axis).
    # Total number of mode scalars: n_nodes * dim.
    theta_vals = []
    for idx in range(len(k_grid)):
        k = k_grid[idx]
        omega_sq = d_of_k_eigenvalues(k, action_params.alpha,
                                      action_params.k_s, action_params.rho, a)
        # omega_sq has shape (dim,); one θ per polarization axis
        arg = 1.0 - 0.5 * dt * dt * omega_sq
        # Clip to [-1, 1] for arccos (floating-point safety only — not a threshold)
        arg_clipped = np.clip(arg, -1.0, 1.0)
        theta_vals.append(np.arccos(arg_clipped))

    return np.concatenate(theta_vals)  # shape (n_nodes * dim,)


# ---------------------------------------------------------------------------
# Public API: apply BCs and compute constraint arrays
# ---------------------------------------------------------------------------


def apply_dirichlet(
    world: np.ndarray,
    bc: DirichletBC,
) -> np.ndarray:
    """Return a copy of ``world`` with boundary slices pinned to bc values.

    Parameters
    ----------
    world : ndarray, shape (N+1, n_nodes, m_ambient)
    bc : DirichletBC

    Returns
    -------
    ndarray, shape (N+1, n_nodes, m_ambient)
        Copy with world[0] = bc.R0 and world[-1] = bc.RN.
    """
    out = world.copy()
    out[0] = bc.R0
    out[-1] = bc.RN
    return out


def apply_chiral(
    bc: ChiralBC,
    lattice: "SpacelikeLattice",  # noqa: F821  (forward ref; imported in bvp)
    params: ActionParams,
    mass: float,
) -> tuple[np.ndarray, float]:
    """Produce the full world-volume by marching from the two past slices.

    This is the correct chiral Cauchy BC (verdict a, 2026-05-30).  The
    world-volume is built by the Störmer–Verlet stencil:

        Forward: R^{l+1} = 2 R^l - R^{l-1} + (dt²/m) F(R^l)   l=1..N-1

        Backward (reverse stencil, antimatter worldtube):
            Treat (R0, R1) as (R^N, R^{N-1}); march in −l:
            R^{l-1} = 2 R^l - R^{l+1} + (dt²/m) F(R^l)

    Reality is automatic: real stencil + real data → real world-volume.
    No FFT, no characteristic projection, no future boundary condition needed.

    Parameters
    ----------
    bc : ChiralBC
        Holds R0 (slice l=0), R1 (slice l=1), and chirality flag.
    lattice : SpacelikeLattice
        Spacelike neighbor topology and lattice parameters.
    params : ActionParams
        Action parameters (dt, n_slices N, k_s, alpha, rho).
    mass : float
        Node mass m = rho * a^dim.

    Returns
    -------
    world : ndarray, shape (N+1, n_nodes, m_ambient)
        Full world-volume produced by the march.  world[0] = bc.R0,
        world[1] = bc.R1 (forward) or world[N] = bc.R0, world[N-1] = bc.R1
        (backward), with all other slices filled by the recurrence.
    condition_estimate : float
        N-independent bounded condition estimate from bc.condition_estimate().
    """
    N = params.n_slices
    dt2_over_m = (params.dt * params.dt) / mass
    n_nodes, m_ambient = bc.R0.shape

    world = np.empty((N + 1, n_nodes, m_ambient), dtype=np.float64)

    if bc.chirality == "forward":
        # Forward march: matter worldtube l=0..N
        world[0] = bc.R0
        world[1] = bc.R1
        R_prev = world[0]
        R_curr = world[1]
        for l in range(1, N):
            F = spacelike_force(R_curr, lattice, params)
            R_next = 2.0 * R_curr - R_prev + dt2_over_m * F
            world[l + 1] = R_next
            R_prev = R_curr
            R_curr = R_next
    else:
        # Backward march: antimatter worldtube.
        # (R0, R1) are the two most-future slices: R^N := R0, R^{N-1} := R1.
        # Reverse stencil: R^{l-1} = 2 R^l - R^{l+1} + (dt²/m) F(R^l)
        world[N] = bc.R0
        world[N - 1] = bc.R1
        R_next = world[N]
        R_curr = world[N - 1]
        for l in range(N - 1, 0, -1):
            F = spacelike_force(R_curr, lattice, params)
            R_prev = 2.0 * R_curr - R_next + dt2_over_m * F
            world[l - 1] = R_prev
            R_next = R_curr
            R_curr = R_prev

    cond = bc.condition_estimate(lattice.params, params)
    return world, cond


# ---------------------------------------------------------------------------
# Condition-number diagnostic: Dirichlet vs chiral at a given N and lattice
# ---------------------------------------------------------------------------


def dirichlet_condition_estimate(
    lattice_params: LatticeParams,
    action_params: ActionParams,
) -> float:
    """Compute the Dirichlet two-time condition estimate for a given N and lattice.

    Returns max_k 1/|sin(N·θ(k))|; large when any mode is near resonance.
    """
    theta_vals = _theta_for_all_modes(lattice_params, action_params)
    N = action_params.n_slices
    sin_vals = np.abs(np.sin(N * theta_vals))
    eps = 1e-14
    return float(np.max(1.0 / np.maximum(sin_vals, eps)))


def find_resonant_N(
    lattice_params: LatticeParams,
    action_params_template: ActionParams,
    N_range: range | list[int],
    condition_threshold: float = 1e3,
) -> list[tuple[int, float]]:
    """Scan a range of N values and return those near resonance.

    A resonant N has condition estimate > condition_threshold.

    Returns list of (N, condition) for resonant extents.
    """
    results = []
    for N in N_range:
        ap = ActionParams(
            k_s=action_params_template.k_s,
            alpha=action_params_template.alpha,
            rho=action_params_template.rho,
            dt=action_params_template.dt,
            n_slices=N,
            m_ambient=action_params_template.m_ambient,
            r_t=action_params_template.r_t,
        )
        cond = dirichlet_condition_estimate(lattice_params, ap)
        if cond > condition_threshold:
            results.append((N, cond))
    return results
