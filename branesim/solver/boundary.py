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
       papers/core/derivations/status.md A2 (verdict a) and ARCHITECTURE.md §2 D2.

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
class PeriodicBC:
    """Rotating-frame-periodic BC — a closed (cyclic) time loop, no fixed slices.

    The worldtube is periodic over the loop: R^P ≡ R^0.  ALL P = n_slices slices
    are unknowns; none is pinned (pinning a whole slice would re-create the
    Dirichlet two-time resonance, κ~1/sin Pθ → 1e14).  The carrier winding n_t
    is a topological integer carried by the seed (the JFNK initial guess) and
    preserved by the smooth Newton flow unless the amplitude collapses.

    "Rotating frame": the carrier phase advances 2π·n_t over the loop, so the
    physical field is periodic up to a global U(1) carrier rotation — which, for
    integer n_t, is the identity, hence the lab-frame cyclic condition R^P ≡ R^0.
    The carrier phase / rigid ambient symmetries (E(4) + the carrier SO(2)) are a
    low-dimensional null space of the linearization; lgmres tolerates it from a
    seed start, and an optional point anchor (``gauge="anchor"``) pins one node's
    slice-0 position to remove the translational part.

    Attributes
    ----------
    R0 : ndarray, shape (n_nodes, m_ambient)
        Slice-0 of the seed — used for shape, and (if ``gauge="anchor"``) as the
        pinned-anchor reference.  The full initial guess is passed to
        ``solve_block`` via ``initial_world``.
    gauge : {"none", "anchor"}
        ``"none"`` (default): rely on lgmres tolerating the symmetry null space
        from the seed start.  ``"anchor"``: additionally pin node ``gauge_node``
        on slice 0 to its R0 value (removes the 4 ambient translations).
    gauge_node : int
        Node index to anchor when ``gauge="anchor"``.
    """

    R0: np.ndarray
    gauge: Literal["none", "anchor"] = "none"
    gauge_node: int = 0

    def condition_estimate(
        self,
        lattice_params: LatticeParams,
        action_params: ActionParams,
    ) -> float:
        """Periodic-operator condition estimate (contrast with Dirichlet's 1/sin Pθ).

        The closed-loop linear operator at temporal harmonic j and spatial mode
        k has eigenvalue ``|ω²(k) − ω_t,j²|`` with the temporal harmonic
        ``ω_t,j² = (2/dt²)(1 − cos 2πj/P)``.  The condition number is
        max/min over the NONZERO eigenvalues — large only at a genuine
        standing-wave resonance ω(k) = ω_t,j, NOT at every Dirichlet node, so it
        grows polynomially in P rather than diverging as 1/sin Pθ.
        """
        P = action_params.n_slices
        dt = action_params.dt
        omega2 = _omega2_for_all_modes(lattice_params, action_params)   # (n_modes,)
        j = np.arange(P)
        omega_t2 = (2.0 / (dt * dt)) * (1.0 - np.cos(2.0 * np.pi * j / P))  # (P,)
        eig = np.abs(omega2[:, None] - omega_t2[None, :]).ravel()
        nonzero = eig[eig > 1e-12]
        if nonzero.size == 0:
            return 1.0
        return float(np.max(nonzero) / np.min(nonzero))


@dataclass
class ChiralBC:
    """Chiral Cauchy BC — two adjacent past slices, marched forward (verdict a).

    The chiral "solve" is the forward (or backward) Verlet march from the two
    past slices.  No FFT, no characteristic projection, no future condition is
    needed.  Reality is automatic: real stencil + real data → real world-volume.
    This is the **r_t=0 linear/Verlet limit** (free-wave Cauchy march); a
    prestressed r_t>0 worldtube is implicit and time-periodic — ``apply_chiral``
    raises for r_t>0 (use ``solve_block`` Dirichlet or ``solver/breather.py``).

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


def _omega2_for_all_modes(
    lattice_params: LatticeParams,
    action_params: ActionParams,
) -> np.ndarray:
    """ω²(k) = D_α(k) for every (polarization axis α, Fourier wavevector k).

    Same mode enumeration as :func:`_theta_for_all_modes`, but returns the raw
    squared eigenfrequencies (no arccos), as needed by the periodic-operator
    condition estimate.  Shape (n_nodes * dim,).
    """
    grid_shape = lattice_params.grid_shape
    a = lattice_params.spacing
    ranges = [np.arange(n) for n in grid_shape]
    grids = np.meshgrid(*ranges, indexing="ij")
    multi = np.stack([g.ravel() for g in grids], axis=1)
    k_grid = 2.0 * np.pi * multi / (np.array(grid_shape) * a)
    out = []
    for idx in range(len(k_grid)):
        out.append(d_of_k_eigenvalues(k_grid[idx], action_params.alpha,
                                      action_params.k_s, action_params.rho, a))
    return np.concatenate(out)


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
    # The chiral solve below is the explicit Störmer–Verlet Cauchy march, which
    # equals the discrete Euler–Lagrange equation ONLY in the r_t=0 linear/Verlet
    # limit (free-wave propagation).  For r_t>0 the temporal central-force spring
    # makes the forward step implicit, so this march would report a residual-zero
    # claim on a non-stationary world-volume.  A prestressed (r_t>0) worldtube is
    # not a two-past-slice Cauchy march: it is a bound, time-periodic eigenstate —
    # solve it with the block root-find (Dirichlet solve_block) or the periodic
    # eigen-BVP in solver/breather.py.  (ARCHITECTURE.md §2 D2 / A4.)
    if params.r_t > 0.0:
        raise NotImplementedError(
            f"apply_chiral() is the r_t=0 linear/Verlet Cauchy march, but "
            f"params.r_t={params.r_t} > 0.  The prestressed worldtube is implicit "
            "and time-periodic; use solve_block() or solver/breather.py instead."
        )

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
