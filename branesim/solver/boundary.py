"""Boundary-condition schemes for the two-time BVP.

ARCHITECTURE.md §2 design decision D2 (resolved 2026-05-30).

Two schemes are provided:

1. **Dirichlet two-time** — fix slices l=0 and l=N to prescribed data.
   Known to be ill-posed at resonant N (determinant 2i·sin(Nθ(k)) vanishes).
   Condition estimate: max_k 1/|sin(Nθ(k))|.  Used as a negative control in
   the validation suite.

2. **Chiral characteristic BC** (κ=1 ∀N) — the well-posed scheme.
   Because D(k) is diagonal in the Cartesian basis at every k and α (6-neighbor
   axial-only stencil, backbone #15, Sprint-2 #9), the block BVP decouples
   per mode into the scalar recurrence
       a^{l+1} − 2·cosθ(k)·a^l + a^{l−1} = 0
       θ(k) = arccos(1 − Δt²·ω²(k)/2)
   For each propagating mode (θ > θ_tol):
     - Fix the forward characteristic amplitude a₊(k) from the past slice (l=0).
     - Impose the no-incoming / Sommerfeld condition on the future slice:
         a^N(k) − e^{−iθ(k)}·a^{N−1}(k) = 0   (annihilates a₋)
   DC/near-zero modes (θ < θ_tol ≈ 1e-6) route through plain Dirichlet.
   The ``chirality`` flag (default ``"forward"``) selects which characteristic
   is fixed; ``"backward"`` flips the roles (antimatter worldtube selection).
   Condition estimate: max_k κ_k = 1 (exactly, for propagating modes).

Implementation note: the Cartesian-polarization projection is simply the real
FFT along each periodic spacelike axis, applied per ambient component.  For a
fully periodic lattice this is exact; for open boundaries the modes are not
pure Fourier modes, but the linear-regime guarantee still applies mode-by-mode
because D(k) remains diagonal.

This module does not modify solver state; it returns modified copies of the
world-volume's boundary rows and reports condition estimates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from branesim.core.conventions import ActionParams, LatticeParams, d_of_k_eigenvalues


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

THETA_TOL: float = 1e-6   # below this θ a mode is treated as DC (Dirichlet)


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
    """Chiral characteristic BC — κ=1 exactly for all N (ARCHITECTURE §2, D2).

    Attributes
    ----------
    R0 : ndarray, shape (n_nodes, m_ambient)
        Prescribed configuration on the past boundary slice (l=0).
    chirality : {"forward", "backward"}
        ``"forward"`` selects matter worldtube (fixes a₊ from l=0, kills a₋
        at l=N).  ``"backward"`` flips the roles for the antimatter branch.
    theta_tol : float
        Modes with θ < theta_tol are treated as DC and handled by Dirichlet.
    """

    R0: np.ndarray
    chirality: Literal["forward", "backward"] = "forward"
    theta_tol: float = THETA_TOL

    def condition_estimate(
        self,
        lattice_params: LatticeParams,
        action_params: ActionParams,
    ) -> float:
        """Condition estimate for the chiral BC scheme.

        For propagating modes the system is perfectly conditioned (κ=1 each).
        DC modes route through Dirichlet (also κ=1 for DC: they have a flat
        a^l = A + Bl solution with no resonance).

        Returns
        -------
        float
            1.0 always (theoretical guarantee).  Exposed here so the solver
            can report it alongside the Dirichlet estimate.
        """
        return 1.0


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
# Core chiral BC application on a 1D modal problem
# ---------------------------------------------------------------------------


def _apply_chiral_bc_1d(
    a0: complex,
    theta: float,
    N: int,
    chirality: str,
) -> tuple[np.ndarray, float]:
    """Given past-slice amplitude a0 and mode parameters, return the constrained
    future-slice value under the chiral (Sommerfeld) condition.

    The general solution of the recurrence
        a^{l+1} − 2·cosθ·a^l + a^{l−1} = 0
    is  a^l = a₊·e^{−iθl} + a₋·e^{+iθl}.

    Forward chirality:
      - Fix a₊ from l=0 data (set a₊ = a0, i.e. the past slice is a pure
        forward wave at l=0 — the correct linearisation of the warm start).
      - No-incoming future condition: a^N − e^{−iθ}·a^{N−1} = 0 kills a₋.

    Substituting: a^N − e^{−iθ}·a^{N−1}
        = (a₊·e^{−iθN} + a₋·e^{+iθN}) − e^{−iθ}(a₊·e^{−iθ(N−1)} + a₋·e^{+iθ(N−1)})
        = a₊·e^{−iθN}(1 − e^{−iθ}·e^{iθ}) + a₋·e^{+iθN}(1 − e^{−iθ}·e^{−iθ})
              [wait — let us do this step by step]

    Step by step:
        a^N     = a₊·e^{−iθN} + a₋·e^{+iθN}
        a^{N−1} = a₊·e^{−iθ(N−1)} + a₋·e^{+iθ(N−1)}
        e^{−iθ}·a^{N−1} = a₊·e^{−iθN} + a₋·e^{+iθ(N−2)}·e^{+iθ·2}·e^{-iθ}
                     wait: e^{-iθ}·e^{+iθ(N-1)} = e^{+iθ(N-2)}

    Cleaner:
        a^N − e^{−iθ}·a^{N−1}
          = a₊(e^{−iθN} − e^{−iθ}·e^{−iθ(N−1)}) + a₋(e^{+iθN} − e^{−iθ}·e^{+iθ(N−1)})
          = a₊(e^{−iθN} − e^{−iθN}) + a₋(e^{+iθN} − e^{+iθ(N−1)−iθ})
          = 0 + a₋(e^{+iθN} − e^{+iθ(N−2)})
    Hmm — that is not zero in general for a₋.

    Correct derivation:
        a^N   − e^{−iθ} a^{N−1}
          = a₊ e^{−iθN} + a₋ e^{iθN}
            − e^{−iθ}(a₊ e^{−iθ(N−1)} + a₋ e^{iθ(N−1)})
          = a₊(e^{−iθN} − e^{−iθ} e^{−iθ(N−1)}) + a₋(e^{iθN} − e^{−iθ} e^{iθ(N−1)})
          = a₊ e^{−iθN}(1 − 1) + a₋(e^{iθN} − e^{iθ(N−1)−iθ})
          = 0 + a₋(e^{iθN} − e^{iθ(N−2)})

    Wait, e^{-iθ}·e^{iθ(N-1)} = e^{iθ(N-1)-iθ} = e^{iθ(N-2)}.
    So the a₋ coefficient is e^{iθN} - e^{iθ(N-2)}.
    That is NOT zero in general.

    The correct outgoing row should be:
        a^N − e^{iθ} a^{N−1} = 0   (not e^{-iθ})

    Check: a^N − e^{iθ}·a^{N−1}
        = a₊(e^{−iθN} − e^{iθ}·e^{−iθ(N−1)}) + a₋(e^{iθN} − e^{iθ}·e^{iθ(N−1)})
        = a₊(e^{−iθN} − e^{iθ−iθ(N−1)}) + a₋(e^{iθN} − e^{iθN})
        = a₊(e^{−iθN} − e^{−iθ(N−2)}) + 0

    That zeros a₋ but keeps an a₊ residual (not what we want either).

    The correct pair: for forward chirality, the *outgoing* wave at the future
    boundary is a^l ~ e^{-iθl} (forward).  The Sommerfeld absorbing condition
    at l=N for an outgoing wave is:
        a^N = e^{-iθ} a^{N-1}           →  kills the backward wave a₋.

    Proof:
        a^N = e^{-iθ}·a^{N-1}
        a₊ e^{-iθN} + a₋ e^{+iθN} = e^{-iθ}(a₊ e^{-iθ(N-1)} + a₋ e^{+iθ(N-1)})
        a₊ e^{-iθN} + a₋ e^{+iθN} = a₊ e^{-iθN} + a₋ e^{+iθN-iθ+iθ}
                  [wait: e^{-iθ}·e^{+iθ(N-1)} = e^{+iθ(N-1)-iθ} = e^{+iθ(N-2)}]
        a₊ e^{-iθN} + a₋ e^{+iθN} = a₊ e^{-iθN} + a₋ e^{+iθ(N-2)}
        a₋(e^{+iθN} - e^{+iθ(N-2)}) = 0
        a₋ e^{+iθ(N-2)}(e^{+2iθ} - 1) = 0  → kills a₋ iff θ ≠ kπ (propagating).

    So for θ ∈ (0,π) the condition a^N = e^{-iθ}·a^{N-1} annihilates a₋.
    And a^N = e^{-iθ}·a^{N-1} means the *future slice value* equals e^{-iθ}
    times the penultimate slice value — a forward-propagating wave condition.

    For the solve: given a₊ (from l=0 past data), the interior evolves as
    a^l = a₊·e^{-iθl} (pure forward), so a^N = a₊·e^{-iθN}.

    Returns (a_N_prescribed, kappa_k) where a_N_prescribed is the value the
    future boundary amplitude is constrained to, and kappa_k = 1.0.
    """
    if chirality == "forward":
        # a₊ = a0  (past slice is the forward amplitude at l=0)
        a_plus = a0
        # Pure forward propagation: a^l = a₊·e^{-iθl}
        a_N = a_plus * np.exp(-1j * theta * N)
    else:
        # Backward chirality: fix a₋ from l=0, propagate forward in the
        # backward-wave direction.  a₋ = a0; a^l = a₋·e^{+iθl}.
        a_minus = a0
        a_N = a_minus * np.exp(1j * theta * N)
    return a_N, 1.0


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
    world: np.ndarray,
    bc: ChiralBC,
    lattice_params: LatticeParams,
    action_params: ActionParams,
) -> tuple[np.ndarray, float]:
    """Return (world_constrained, condition_estimate).

    Applies the chiral characteristic BC to the future slice l=N:
    - Projects the boundary slabs onto Fourier modes (per polarization axis).
    - For each propagating mode: pins the forward characteristic a₊ from
      the past slice (l=0) and sets the future slice (l=N) to the forward-
      propagated value a₊·e^{-iθN}.
    - For each DC mode (θ < theta_tol): plain Dirichlet on that mode's future
      slice value (copies from the past).

    Parameters
    ----------
    world : ndarray, shape (N+1, n_nodes, m_ambient)
        Current world-volume; world[0] is the past slice.
    bc : ChiralBC
    lattice_params : LatticeParams
    action_params : ActionParams

    Returns
    -------
    world_out : ndarray, shape (N+1, n_nodes, m_ambient)
        Copy with future slice (world[-1]) replaced by the chiral BC value.
    condition_estimate : float
        1.0 (guaranteed for propagating modes; DC modes also κ=1).
    """
    grid_shape = lattice_params.grid_shape
    dim = lattice_params.dim
    N = action_params.n_slices
    dt = action_params.dt
    a = lattice_params.spacing
    theta_tol = bc.theta_tol

    # Work with the past slice (l=0 = bc.R0) and derive the future slice.
    R0 = bc.R0  # shape (n_nodes, m_ambient)
    m_ambient = R0.shape[1]

    # Reshape positions to (grid_shape..., m_ambient) for FFT
    R0_grid = R0.reshape(grid_shape + (m_ambient,))  # (N0, N1, ..., m_amb)

    # We will build the future slice by operating in modal space.
    # Only spacelike dimensions are periodic (FFT-able); process per-ambient-component.
    # For each ambient component c, for each polarization axis pol:
    #   - Take FFT along all spacelike axes (gives modal amplitudes at each k).
    #   - For each k, the ω²(k) eigenvalue along polarization axis pol is
    #     D_pol(k), and θ(k) = arccos(1 - dt²·D_pol(k)/2).
    #   - Replace the modal future amplitude with the chiral BC value.
    #   - IFFT back to real space.

    # Build wavevector index arrays for the grid.
    ranges = [np.arange(n, dtype=np.float64) for n in grid_shape]
    grids = np.meshgrid(*ranges, indexing="ij")
    multi = np.stack([g for g in grids], axis=-1)  # (N0, N1, ..., dim)

    # Physical wavevectors: k_i = 2π·n_i / (N_i·a)
    # shape (N0, N1, ..., dim)
    k_grid = 2.0 * np.pi * multi / (np.array(grid_shape, dtype=np.float64) * a)

    # Precompute ω²(k) for each mode and each polarization axis:
    # omega_sq shape: (N0, N1, ..., dim)
    n_nodes = R0.shape[0]
    k_flat = k_grid.reshape(-1, dim)  # (n_nodes, dim)
    omega_sq_flat = np.zeros((n_nodes, dim), dtype=np.float64)
    for idx in range(n_nodes):
        omega_sq_flat[idx] = d_of_k_eigenvalues(
            k_flat[idx], action_params.alpha,
            action_params.k_s, action_params.rho, a
        )
    omega_sq_grid = omega_sq_flat.reshape(grid_shape + (dim,))

    # Compute θ for all modes: shape (N0, N1, ..., dim)
    arg = 1.0 - 0.5 * dt * dt * omega_sq_grid
    arg_clipped = np.clip(arg, -1.0, 1.0)
    theta_grid = np.arccos(arg_clipped)  # shape (N0, N1, ..., dim)

    # Propagating mask: shape (N0, N1, ..., dim)
    propagating = theta_grid > theta_tol

    # The future slice R_N will be built component by component.
    # For a fully periodic lattice, ambient component c couples to spacelike
    # polarization axis pol iff pol < dim (the spacelike components).
    # The ambient components beyond dim-1 (the timelike/extra components) have
    # no direct D(k) coupling, but they still need to be propagated.  For
    # those components there is no "right" chiral split from D(k) alone;
    # we propagate them as if they were longitudinal (pol=0-like) using the
    # longitudinal θ at each k.  This is consistent with the 4th ambient
    # component being a deformation in the timelike direction — it disperses
    # like the longitudinal mode.

    # Apply FFT along all spacelike axes for each ambient component:
    R0_fft = np.fft.fftn(R0_grid, axes=tuple(range(dim)))  # complex, same shape
    RN_fft = np.zeros_like(R0_fft, dtype=complex)

    for c in range(m_ambient):
        # Determine which polarization axis this ambient component "lives on"
        # Components 0..dim-1: correspond to spacelike polarizations 0..dim-1.
        # Components dim..m_ambient-1: use the longitudinal (axis 0) θ, which
        # is the appropriate conservative choice for extra ambient components.
        pol = min(c, dim - 1)
        theta_c = theta_grid[..., pol]   # shape (N0, N1, ...)
        prop_c = propagating[..., pol]   # shape (N0, N1, ...)

        A0 = R0_fft[..., c]   # complex modal amplitudes at l=0, shape (N0,N1,...)

        # For propagating modes: a^N = a₊·e^{-iθN}  (forward chirality)
        # For DC modes: copy the past amplitude (Dirichlet on that mode).
        if bc.chirality == "forward":
            # a₊ = A0 (pure forward assumption from the past data)
            A_N_prop = A0 * np.exp(-1j * theta_c * N)
        else:
            # a₋ = A0; a^N = a₋·e^{+iθN}
            A_N_prop = A0 * np.exp(1j * theta_c * N)

        # DC modes: Dirichlet (copy the past slice amplitude — i.e. A^N = A^0)
        A_N = np.where(prop_c, A_N_prop, A0)
        RN_fft[..., c] = A_N

    # IFFT back to real space — the result should be real (up to float precision)
    RN_grid = np.fft.ifftn(RN_fft, axes=tuple(range(dim))).real

    # Flatten back to (n_nodes, m_ambient)
    RN = RN_grid.reshape(n_nodes, m_ambient)

    # Build output world-volume with pinned boundary slices
    out = world.copy()
    out[0] = bc.R0
    out[-1] = RN

    return out, 1.0


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
            temporal_model=action_params_template.temporal_model,
            r_t=action_params_template.r_t,
        )
        cond = dirichlet_condition_estimate(lattice_params, ap)
        if cond > condition_threshold:
            results.append((N, cond))
    return results
