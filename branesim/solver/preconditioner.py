"""FFT (frozen-vacuum-linearized) preconditioner for the PeriodicBC JFNK solve.

=============================================================================
REVERTIBILITY (read first)
=============================================================================
This whole module is **opt-in and isolated** so it can be removed cleanly:

  * It does nothing unless ``SolveOpts.preconditioner == "fft_linear"`` (default
    ``"none"`` → identical behaviour to before this module existed).
  * To revert entirely: delete this file and remove the three call sites that
    reference it — (1) the ``inner_M=...`` wiring in ``branesim/solver/bvp.py``
    (search "preconditioner"), (2) the ``preconditioner`` / ``precond_floor``
    fields on ``SolveOpts``, (3) ``preconditioner="fft_linear"`` in
    ``branesim/experiments/vortex_seed_render.py`` ``_relax``.  Nothing else
    imports it.

=============================================================================
WHY (the problem it solves)
=============================================================================
The rotating-frame-periodic JFNK solve floors at ‖R‖≈1.5 by outer-iter ~8 on
the 96³ vortex (``condition_estimate ≈ 4.4e7``).  The inner ``lgmres`` is
UNPRECONDITIONED, so at that conditioning it removes the well-conditioned error
fast and then stalls on the stiff subspace it cannot resolve within
``inner_maxiter`` cycles.  A preconditioner ``M⁻¹ ≈ J⁻¹`` is the principled fix.

=============================================================================
HOW (the operator — exact, then the one approximation)
=============================================================================
The linearized periodic residual Jacobian is DIAGONAL in the joint Fourier
basis (3D spatial FFT × 1D temporal FFT) and the ambient-component basis,
because (``core.conventions.d_of_k_eigenvalues``):

    the axial-stencil dynamical matrix D(k) is diagonal in the Cartesian basis
    at every k and α — its eigenvectors are the ambient unit vectors ê_a.

So per spatial wavevector k, ambient component a, and temporal harmonic j, the
closed-loop linear operator has the scalar eigenvalue (the same one
``PeriodicBC.condition_estimate`` uses):

    λ(k, a, j) = ω²_a(k) − ω²_t,j

    ω²_a(k) = (2 k_s/ρ) [ α·h_a + (1−α)·Σ_b h_b ],   h_i = 1 − cos(k_i a)
              (a = 0,1,2 spatial; for the TIMELIKE component a = 3 there is no
               longitudinal spatial bond, so h_a → 0:  ω²_3 = (2k_s/ρ)(1−α)Σh)
    ω²_t,j  = (2/dt²)(1 − cos 2π j/P)                 (periodic temporal 2nd-diff)

The preconditioner applies the inverse of |this operator| via FFTs:

    M⁻¹ v = IFFT_{t,x,y,z}[ FFT_{t,x,y,z}(v) / max(|λ|, ε) ]   (per component a)

THE ONE APPROXIMATION: λ is the operator linearized about the *prestressed
vacuum* (frozen coefficients) — it ignores (i) the nonlinear vortex-core
correction (the very stiffness that drives cond up — left for the outer Newton
to mop up) and (ii) the r_t time-link spring's own stiffness/prestress (the
temporal term is the kinematic 2nd-difference, matching ``condition_estimate``).
A preconditioner need not be exact; it only must capture the dominant spectral
spread, which the spatial-Laplacian + temporal-harmonic structure does.

ε = ``precond_floor`` · max|λ| floors the near-null modes (the k=0,j=0
translational/zero modes and the ω²_a(k)=ω²_t,j resonances) so M⁻¹ does NOT
blow them up; using ``1/max(|λ|,ε)`` makes M⁻¹ a symmetric positive-definite
≈|J|⁻¹, the standard choice for an indefinite wave operator under GMRES/lgmres.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator

from branesim.core.conventions import ActionParams, LatticeParams


def build_fft_preconditioner(
    lattice_params: LatticeParams,
    action_params: ActionParams,
    n_slices: int,
    m_ambient: int,
    precond_floor: float = 1e-2,
) -> LinearOperator:
    """Build the frozen-vacuum FFT preconditioner ``M⁻¹`` for the PeriodicBC solve.

    Parameters
    ----------
    lattice_params, action_params : the solve's parameters (α, k_s, ρ, dt, …).
    n_slices : int
        Period ``P`` of the time loop (the solver flattens slices 0..P-1).
    m_ambient : int
        Ambient dimension (4): components 0,1,2 spatial, 3 timelike.
    precond_floor : float
        ε / max|λ| — the relative floor on |λ| that caps amplification of the
        near-null (zero-mode / resonant) subspace.  Larger = gentler / safer,
        smaller = more aggressive.

    Returns
    -------
    scipy.sparse.linalg.LinearOperator
        ``M`` with ``M.matvec(v) = M⁻¹ v``, shape (P·N·m, P·N·m), to pass as
        ``inner_M`` to ``scipy.optimize.newton_krylov``.  Real-valued in/out.
    """
    # Dimension-agnostic (PRINCIPLES §7.6: solver core must not hard-code dim).
    grid = tuple(int(n) for n in lattice_params.grid_shape)
    dim = len(grid)
    P = int(n_slices)
    n_nodes = int(np.prod(grid))

    # --- spatial ω²_a(k) on the FFT grid, per ambient component a -------------
    # h_i(k_i) = 1 − cos(k_i a) with k_i = 2π n_i/(N_i a)  ⇒  h_i = 1 − cos(2π n_i/N_i)
    # (a cancels; matches d_of_k_eigenvalues' mode enumeration).  Symmetric in
    # n_i → −n_i, so the np.fft frequency ordering needs no care here.
    h_axes = [1.0 - np.cos(2.0 * np.pi * np.arange(n) / n) for n in grid]
    h_grids = np.meshgrid(*h_axes, indexing="ij")   # list of dim arrays, shape `grid`
    Hsum = sum(h_grids)                              # Σ_b h_b, shape `grid`
    pref = 2.0 * action_params.k_s / action_params.rho
    alpha = action_params.alpha
    # Per-component spatial eigenvalue ω²_a(k): a longitudinal axis (a < dim) uses
    # h_a; components without a longitudinal bond (a >= dim, e.g. the timelike
    # channel) have h_a → 0 → transverse-only (1−α)Σh.
    comp_h = [h_grids[c] if c < dim else np.zeros_like(Hsum) for c in range(m_ambient)]
    omega2_spatial = np.stack(
        [pref * (alpha * comp_h[c] + (1.0 - alpha) * Hsum) for c in range(m_ambient)],
        axis=-1,
    )  # shape (grid..., m)

    # --- temporal ω²_t,j (periodic 2nd-difference harmonics) ------------------
    dt = action_params.dt
    j = np.arange(P)
    omega2_t = (2.0 / (dt * dt)) * (1.0 - np.cos(2.0 * np.pi * j / P))   # (P,)

    # --- λ(j,k,a) = ω²_a(k) − ω²_t,j, then the floored inverse multiplier -----
    t_shape = (P,) + (1,) * dim + (1,)               # broadcast ω²_t over grid+comp
    lam = omega2_spatial[None, ...] - omega2_t.reshape(t_shape)   # (P, grid..., m)
    lam_abs = np.abs(lam)
    eps = precond_floor * max(float(lam_abs.max()), 1e-30)  # absolute floor guards
    inv = 1.0 / np.maximum(lam_abs, eps)             # (P, grid..., m); SPD ≈ |J|⁻¹

    shape_full = (P,) + grid + (m_ambient,)
    n_dof = P * n_nodes * m_ambient
    fft_axes = tuple(range(dim + 1))  # 0 = temporal, 1..dim = spatial (comp axis left)

    def _matvec(v: np.ndarray) -> np.ndarray:
        V = np.asarray(v, dtype=np.float64).reshape(shape_full)
        Vk = np.fft.fftn(V, axes=fft_axes)
        Vk *= inv
        out = np.fft.ifftn(Vk, axes=fft_axes).real
        return out.reshape(n_dof)

    return LinearOperator((n_dof, n_dof), matvec=_matvec, rmatvec=_matvec, dtype=np.float64)
