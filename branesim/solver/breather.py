"""Time-periodic discrete-breather eigen-solver.

Step-3 particle-search solver.  Finds time-periodic, spatially localized
solutions of the discrete brane equation of motion (the Pythagorean
central-force energy) by root-finding the cyclic time-collocation residual.

Physics context
---------------
The Pythagorean central-force link norm hardens for a staggered (q=π)
transverse carrier.  In a 1D transverse chain the exact nonlinear frequency
of a spatially-uniform staggered mode satisfies

    ω²(A) = (4 k_s / m) · (1 − α / √(1 + 4 A² / a²))

At small A this reproduces the Duffing law

    ω²(A) ≈ ω_max² + (3Q / 4m) A²,
    ω_max² = 4 k_s (1 − α) / m,   Q = 8 k_s α / a²

The hardening pushes ω above the acoustic band top ω_max, making the
mode above-band and therefore non-radiating → a discrete breather.
The solver finds the SPATIALLY LOCALIZED version by coupling the uniform
nonlinear frequency to the near-zero linear evanescent tail.

Method: time-collocation Newton (JFNK, matrix-free)
----------------------------------------------------
Represent one temporal period by P slices R⁰..R^{P-1}, each shape
(n_nodes, m_ambient).  The period T is a continuous unknown.
Effective timestep: dt_eff = T / P.

Cyclic residual at slice l:
    ℛ^l = m (R^{l+1} − 2R^l + R^{l−1}) / dt_eff² − F(R^l)
where F = spacelike_force (same Pythagorean central-force energy as the IVP
solver).  Cyclic: R^P ≡ R^0, R^{−1} ≡ R^{P−1}.

P must be even (required for the half-period anti-amplitude constraint).

Unknowns: {R^0..R^{P-1}} (P·n_nodes·m) + T (1)  = P·n·m + 1 total.
Equations (P·n·m + 1):
  · P·n·m − 2  cyclic residual equations (all except (l=0,peak,lat) and
               (l=P/2,peak,lat), which are replaced by amplitude conditions)
  · 1          amplitude condition at l=0:    R^0[peak,lat]·sign = +A
  · 1          anti-amplitude condition at l=P/2: R^{P/2}[peak,lat]·sign = −A
               Together these select the non-trivial oscillating branch (A at
               l=0, −A at half-period) and rule out the static equilibrium
               solution (which would give +A at both l=0 and l=P/2).
  · 1          gauge fix (time-translation zero mode): the peak node's
               lateral finite-difference velocity = 0 at l=0:
               R^1[peak, lat] − R^{P-1}[peak, lat] = 0
               (turning-point condition: carrier at its maximum)

Solved with scipy.optimize.newton_krylov (LGMRES inner, matrix-free JFNK).

SADDLE DISCIPLINE: the solver targets ‖ℛ‖ = 0 (root-find).
The brane action S is Lorentzian (saddle, unbounded below) — minimising S
would diverge or solve the wrong Euclidean problem.  This module enforces
the root-find discipline: OBJECTIVE == "residual_norm" (not "action").
(ARCHITECTURE.md §1.3, PRINCIPLES.md §1.2, papers/core/derivations/status.md A1.)

API
---
    solve_breather(lattice, params_like, mass, *, P, amplitude, seed, opts)
        → dict with keys: slices, T, omega, residual_norm, converged

    continue_breather(lattice, params_like, mass, *, P, amplitudes, ...)
        → list[dict]

Dimension-agnostic: uses len(grid_shape) for all loops; no hard-coded 3D.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.optimize import newton_krylov, NoConvergence

from branesim.core.action import spacelike_force
from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.seeds import skyrme_twisted_hedgehog


# ---------------------------------------------------------------------------
# Module-level saddle-discipline assertion
# ---------------------------------------------------------------------------

OBJECTIVE: str = "residual_norm"
assert OBJECTIVE == "residual_norm", (
    "breather.py: solver targets residual_norm (root-find), NOT the Lorentzian "
    "action S.  S is a saddle, unbounded below — minimising it would diverge or "
    "solve the wrong Euclidean problem.  (ARCHITECTURE.md §1.3, PRINCIPLES.md §1.2)"
)


# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BreatherOpts:
    """Options for the breather Newton-Krylov solve.

    Parameters
    ----------
    tol : float
        Convergence tolerance for ‖G‖ (default 1e-8).
    max_iter : int
        Maximum Newton-Krylov outer iterations (default 500).
    inner_maxiter : int
        Maximum LGMRES inner iterations (default 2000).  Larger values are
        needed because the system matrix is poorly conditioned near the band
        edge; 2000 gives reliable convergence for chains up to ~31 nodes.
    method : str
        Inner Krylov method: 'lgmres', 'gmres', 'bicgstab' (default 'lgmres').
    verbose : bool
        Print convergence progress to stdout (default False).
    """

    tol: float = 1e-8
    max_iter: int = 1000
    inner_maxiter: int = 2000
    method: str = "lgmres"
    verbose: bool = False


# ---------------------------------------------------------------------------
# Physical helpers
# ---------------------------------------------------------------------------


def omega_max(k_s: float, alpha: float, mass: float, a: float = 1.0) -> float:
    """Band-top frequency ω_max = sqrt(4 k_s (1−α) / m).

    This is the staggered (q=π) linear frequency of the transverse chain.
    The breather hardening law starts here: ω(A=0) = ω_max.
    """
    return math.sqrt(4.0 * k_s * (1.0 - alpha) / mass)


def omega_exact(
    k_s: float, alpha: float, mass: float, amplitude: float, a: float = 1.0
) -> float:
    """Exact nonlinear frequency of the spatially-uniform staggered mode.

    Derived from F_y = −ω²·m·y for the two-neighbor force:
        ω²(A) = (4 k_s / m) · (1 − α / √(1 + 4 A² / a²))

    This is ABOVE ω_max for all A > 0 (hardening / above-band).
    At A = 0 it reduces to ω_max exactly.
    """
    return math.sqrt(4.0 * k_s / mass * (1.0 - alpha / math.sqrt(1.0 + 4.0 * amplitude ** 2 / a ** 2)))


def omega_duffing(
    k_s: float, alpha: float, mass: float, amplitude: float, a: float = 1.0
) -> float:
    """Small-A Duffing approximation: ω² ≈ ω_max² + (3Q / 4m) A².

    Q = 8 k_s α / a² (per spec definition).  Accurate for A/a ≤ 0.3.
    """
    Q = 8.0 * k_s * alpha / a ** 2
    om_max2 = 4.0 * k_s * (1.0 - alpha) / mass
    # max(0.0, ...) guards only the degenerate edge α=1, A=0 where the band
    # bottom is exactly 0; for all α < 1 or A > 0 the argument is strictly positive.
    return math.sqrt(max(0.0, om_max2 + (3.0 * Q / (4.0 * mass)) * amplitude ** 2))


def omega_longitudinal_top(k_s: float, mass: float, a: float = 1.0) -> float:
    """Closed-form LONGITUDINAL (compression) band-top frequency.

    The transverse staggered mode has band-top ω_max = √(4 k_s (1−α)/m)
    (see :func:`omega_max`), set by the prestress *tension* k_s(1−α).
    The longitudinal branch is restored by the FULL Hookean stiffness k_s
    (∂²/∂r² of ½k_s(r−αa)² = k_s, independent of α), so its band-top is

        ω_L,max = √(4 k_s / m)   (q = π, 1D nearest-neighbor).

    This is the dangerous one: the transverse breather frequency obeys
        ω_max  ≤  ω(A)  <  ω_L,max          for all A,
    because ω²(A) = (4k_s/m)(1 − α/√(1+4A²/a²)) → 4k_s/m as A → ∞.
    So the transverse breather is ABOVE its own (transverse) band but always
    INSIDE the longitudinal continuum.  Radiation into longitudinal phonons is
    therefore only forbidden if the relevant HARMONIC nω clears ω_L,max — which
    is what :func:`harmonic_resonance_check` tests.  (1D NN closed form; use
    :func:`phonon_band_top` for the exact finite-lattice / higher-D value.)
    """
    return math.sqrt(4.0 * k_s / mass)


def omega_band_top_analytic(
    k_s: float,
    mass: float,
    dim: int,
    a: float = 1.0,
) -> float:
    """Dim-D analytic upper bound for the phonon band top of a staggered cubic lattice.

    Returns
    -------
    float
        ω_band_top = √(4 · dim · k_s / m).

    Derivation
    ----------
    In the dim-D hypercubic lattice with nearest-neighbour central-force springs
    the nonlinear breather frequency approaches

        lim_{A→∞} ω(A) = √(4·dim·k_s/m)

    because each of the 2·dim bonds contributes a stiffness k_s in the staggered
    (q=π) limit, and the transverse tension term α/√(1+4A²/a²) → 0.  This is
    therefore a conservative OVER-ESTIMATE of the true phonon band top: the actual
    finite-lattice spectrum maximum is below this value for α > 0 (verified against
    the numeric phonon_band_top at 16³; see NOTES.md §analytic-band-top-validation).

    The key consequence for :func:`harmonic_resonance_check`: because this is an
    OVER-ESTIMATE, passing band_top=omega_band_top_analytic(...) makes the
    radiationless verdict strictly MORE CONSERVATIVE — a harmonic is only called
    "radiationless" if it clears this higher threshold.  The false-negative rate
    (calling a radiating breather "radiationless") is zero under this estimate.

    Dim-agnostic: for dim=1 this equals :func:`omega_longitudinal_top` exactly.
    For dim=3 and the smoke-run parameters (k_s=1, m=1) it returns √12 ≈ 3.464,
    which is ~38% above the numeric finite-lattice value of ~2.518 (see NOTES.md).

    Motivation: the dense phonon eigensystem of :func:`phonon_band_top` costs
    O(N³) for N = n_nodes · m_ambient — infeasible at 32³+ grids (313 s at 16³).
    This closed form costs O(1) and unblocks large-grid sweeps.  Always validate
    the analytic vs numeric values on a small grid before relying on it for a
    physics verdict (the NOTES.md validation provides this certificate).

    Parameters
    ----------
    k_s : float    — spring constant (same units as ActionParams.k_s).
    mass : float   — node mass m (= rho * a^dim for uniform lattice).
    dim : int      — spatial dimension of the brane (e.g. 3 for 3D).
    a : float      — lattice spacing (default 1.0; cancels in the final formula).
    """
    return math.sqrt(4.0 * dim * k_s / mass)


# ---------------------------------------------------------------------------
# Force Jacobian  J = dF/dR  (matrix-free directional derivative + dense build)
# ---------------------------------------------------------------------------


def _jacobian_matvec(
    R: np.ndarray,
    v: np.ndarray,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    fd_eps: float = 1e-6,
) -> np.ndarray:
    """Directional derivative  J(R)·v = (dF/dR)·v  by central finite difference.

    Matrix-free: never forms the N×N Jacobian.  Two force evaluations.
    The step h is chosen so that ‖h·v‖ ≈ fd_eps (scale-invariant in ‖v‖).

    Parameters
    ----------
    R : ndarray, shape (n_nodes, m_ambient)  — base configuration.
    v : ndarray, shape (n_nodes, m_ambient)  — perturbation direction.

    Returns
    -------
    ndarray, shape (n_nodes, m_ambient)  — J(R)·v.
    """
    nrm = float(np.linalg.norm(v))
    h = fd_eps / (nrm + 1e-300)
    Fp = spacelike_force(R + h * v, lattice, params_like)
    Fm = spacelike_force(R - h * v, lattice, params_like)
    return (Fp - Fm) / (2.0 * h)


def _force_jacobian(
    R: np.ndarray,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    fd_eps: float = 1e-6,
) -> np.ndarray:
    """Dense force Jacobian  J[i,j] = dF_i/dR_j  by central finite difference.

    Returns an (N, N) array with N = n_nodes·m_ambient (C-order flatten).
    Used for the phonon spectrum (small one-off cost) and, optionally, for the
    dense monodromy assembly on small lattices.  O(N) force evaluations — fine
    for 1D / small 3D; for large 3D use the matrix-free path in
    :func:`floquet_multipliers`.
    """
    shape = R.shape
    N = R.size
    Rf = R.ravel()
    J = np.empty((N, N), dtype=np.float64)
    for j in range(N):
        dRf = np.zeros(N)
        dRf[j] = fd_eps
        Fp = spacelike_force((Rf + dRf).reshape(shape), lattice, params_like).ravel()
        Fm = spacelike_force((Rf - dRf).reshape(shape), lattice, params_like).ravel()
        J[:, j] = (Fp - Fm) / (2.0 * fd_eps)
    return J


# ---------------------------------------------------------------------------
# Phonon spectrum / band top  (the radiation continuum)
# ---------------------------------------------------------------------------


def phonon_spectrum(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    fd_eps: float = 1e-6,
) -> dict[str, Any]:
    """Linear phonon spectrum on the reference (prestressed) configuration.

    Linearise  m δR̈ = F(R)  about the reference:  m δR̈ = −K δR  with the
    stiffness  K = −dF/dR|_ref = d²V/dR²|_ref.  The eigenvalues of K/m are the
    squared phonon frequencies ω² of every mode on this finite lattice
    (all branches — transverse AND longitudinal — and all directions).  The
    largest is the band top: the top of the radiation continuum [0, ω_band_top].

    Dimension-agnostic: works for any lattice.dim.  K is symmetrised to remove
    finite-difference asymmetry.  Genuinely negative eigenvalues (unstable or
    zero-mode) are returned raw in omega2 — they are physical information.

    Returns
    -------
    dict with keys:
        omegas       : ndarray — phonon frequencies (sqrt of positive eigenvalues
                                 only; negative eigenvalues → 0.0 in this array,
                                 but see omega2 for the true signed values).
        band_top     : float   — max phonon frequency (radiation continuum top).
        omega2       : ndarray — raw eigenvalues of K/m (may include negatives for
                                 unstable / zero modes).  Ascending order.
        n_modes      : int
    """
    m_ambient = params_like.m_ambient or (lattice.dim + 1)
    R_ref = lattice.reference_positions(m_ambient)
    J = _force_jacobian(R_ref, lattice, params_like, fd_eps=fd_eps)
    K = -0.5 * (J + J.T)  # stiffness = −dF/dR, symmetrised
    eig = np.linalg.eigvalsh(K)
    omega2 = eig / mass  # raw: may contain near-zero or negative values for rigid modes / FD noise
    # sqrt domain guard: negative eigenvalues are real instabilities / zero modes;
    # we report 0.0 in omegas for display but preserve the signed value in omega2.
    omegas = np.sqrt(np.maximum(omega2, 0.0))
    return {
        "omegas": omegas,
        "omega2": omega2,
        "band_top": float(omegas[-1]),
        "n_modes": int(omegas.size),
    }


def phonon_band_top(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    fd_eps: float = 1e-6,
) -> float:
    """Top of the phonon radiation continuum (max linear phonon frequency).

    Convenience wrapper over :func:`phonon_spectrum`.  This is the frequency a
    breather harmonic must CLEAR to be radiationless.
    """
    return phonon_spectrum(lattice, params_like, mass, fd_eps=fd_eps)["band_top"]


def harmonic_resonance_check(
    omega: float,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    n_max: int = 6,
    band_top: float | None = None,
    rel_tol: float = 1e-3,
    fd_eps: float = 1e-6,
) -> dict[str, Any]:
    """Test whether the breather fundamental and its harmonics avoid the phonon band.

    An above-band discrete breather is radiationless only if NO harmonic nω
    (n = 1, 2, …) coincides with the linear phonon continuum [0, ω_band_top]:
    a harmonic landing in the band opens a resonant radiation channel
    (the lattice analog of the sine-Gordon "breather sits in a spectral gap"
    stability condition).

    Two band tops matter and are reported separately:

      · transverse_top  = ω_max = √(4k_s(1−α)/m)   — the breather's OWN branch.
        ω > transverse_top is the original above-band (non-self-radiating) claim.
      · band_top        = max phonon frequency (numeric; longitudinal-dominated).
        Because ω(A) < ω_L,max always, the FUNDAMENTAL (n=1) is generically
        INSIDE the full continuum — its coupling to the longitudinal branch is
        parity-suppressed for a symmetric staggered breather (leading radiation
        at 2ω), so the operationally decisive quantity is the lowest n ≥ 2 with
        nω inside the band.

    Parameters
    ----------
    omega : float        — breather fundamental angular frequency (result["omega"]).
    n_max : int          — highest harmonic to test (default 6).
    band_top : float, optional — precomputed full band top; computed if None.
    rel_tol : float      — fractional margin; nω counts as in-band if
                           nω ≤ band_top·(1+rel_tol).

    Returns
    -------
    dict with keys:
        harmonics            : list of dicts {n, freq, in_band, margin}
                               where margin = (nω − band_top)/band_top.
        band_top             : float
        transverse_top       : float
        fundamental_above_transverse : bool  — ω > ω_max (own-branch above-band).
        lowest_in_band_harmonic : int | None — smallest n with nω in band
                               (None if all harmonics clear the continuum).
        lowest_in_band_n_ge_2 : int | None  — same but restricted to n ≥ 2
                               (the symmetry-allowed radiation channel).
        radiationless        : bool  — True iff no n ≥ 2 harmonic is in band.
    """
    if band_top is None:
        band_top = phonon_band_top(lattice, params_like, mass, fd_eps=fd_eps)

    k_s = params_like.k_s
    alpha = params_like.alpha
    a = lattice.params.spacing
    transverse_top = omega_max(k_s, alpha, mass, a)

    thresh = band_top * (1.0 + rel_tol)
    harmonics: list[dict[str, Any]] = []
    lowest_in_band: int | None = None
    lowest_in_band_n_ge_2: int | None = None
    for n in range(1, n_max + 1):
        freq = n * omega
        in_band = freq <= thresh
        harmonics.append({
            "n": n,
            "freq": float(freq),
            "in_band": bool(in_band),
            "margin": float((freq - band_top) / band_top),
        })
        if in_band and lowest_in_band is None:
            lowest_in_band = n
        if in_band and n >= 2 and lowest_in_band_n_ge_2 is None:
            lowest_in_band_n_ge_2 = n

    return {
        "harmonics": harmonics,
        "band_top": float(band_top),
        "transverse_top": float(transverse_top),
        "fundamental_above_transverse": bool(omega > transverse_top),
        "lowest_in_band_harmonic": lowest_in_band,
        "lowest_in_band_n_ge_2": lowest_in_band_n_ge_2,
        "radiationless": lowest_in_band_n_ge_2 is None,
    }


# ---------------------------------------------------------------------------
# Seed construction
# ---------------------------------------------------------------------------


def _build_seed(
    lattice: SpacelikeLattice,
    k_s: float,
    alpha: float,
    mass: float,
    a: float,
    P: int,
    amplitude: float,
    m_ambient: int,
) -> tuple[np.ndarray, float]:
    """Build a linear seed for the Newton iteration.

    Returns (slices, T_seed) where slices has shape (P, n_nodes, m_ambient).

    The seed is: R^l_p = ref_p + f(p) · (−1)^{p_center} · e_lat · cos(2π l / P)
    where:
      · f(p) = amplitude · exp(−((p − center) / sigma)²) is a Gaussian envelope
      · p_center is the peak node index (center of chain)
      · e_lat is the lateral unit vector (component m_ambient−1 in 1D, or more
        generally the last ambient component)
      · T_seed = 2π / ω_max (small-amplitude period, above-band starting point)

    The staggered phase (−1)^{node_index_along_axis} is applied uniformly:
    phase_p = (−1)^{multi_index_sum(p) mod 2} so the seed generalises to
    higher dimensions automatically.

    For the gauge-fix condition (peak velocity = 0 at l=0) to be satisfied
    automatically by the seed we need:
        R^1[peak, lat] = R^{P-1}[peak, lat]
    which holds when cos(2π/P) = cos(−2π/P), i.e., always for cosine seeds.
    """
    n_nodes = lattice.n_nodes
    ref = lattice.reference_positions(m_ambient)  # (n_nodes, m_ambient)

    dim = lattice.dim
    mi = lattice.multi_indices  # (n_nodes, dim)

    # Center node: multi-index at center of grid
    center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
    # Spatial distance to center (in lattice units)
    d_to_center = np.linalg.norm((mi - center_mi) * a, axis=1)  # (n_nodes,)

    # Dimension-aware band-top and nonlinear frequency.
    # The staggered (q=π in all axes) mode has 2*dim neighbors each contributing
    # 2*(1-α)*k_s/m transverse stiffness (factor 2 from staggered sign): so
    #   ω²_max(ndim) = 2*dim × 2*(1-α)*k_s/m = 4*dim*(1-α)*k_s/m
    # and the exact nonlinear frequency generalises as:
    #   ω²(A, ndim) = 4*dim*k_s/m × (1 - α/√(1 + 4A²/a²))
    # For dim=1 these reduce to the 1D omega_max / omega_exact formulas.
    om_max_nd = math.sqrt(4.0 * dim * k_s * (1.0 - alpha) / mass)
    om_nd2 = 4.0 * dim * k_s / mass * (1.0 - alpha / math.sqrt(1.0 + 4.0 * amplitude ** 2 / a ** 2))
    om_nd = math.sqrt(max(om_nd2, om_max_nd ** 2))  # at least band top (avoids sub-band seed)

    # Evanescent wavenumber: along each decay axis the per-axis evanescent stiffness
    # is 2*(1-α)*k_s/m (same as 1D — it's a per-bond quantity, not a sum over bonds).
    # cosh(κa) = 1 + (ω² − ω²_max) / (2*(1-α)*k_s/m)
    delta_om2 = om_nd ** 2 - om_max_nd ** 2
    prefac = 2.0 * (1.0 - alpha) * k_s / mass
    cosh_arg = 1.0 + delta_om2 / prefac if prefac > 1e-12 else 2.0
    cosh_arg = max(cosh_arg, 1.0 + 1e-6)  # must be > 1 for arccosh
    kappa_inv = a / math.acosh(cosh_arg)
    sigma = max(kappa_inv, 0.5 * a)  # at least half a lattice spacing

    # Gaussian envelope (Gaussian approximates the sech tail at small A)
    envelope = amplitude * np.exp(-(d_to_center / sigma) ** 2)  # (n_nodes,)

    # Staggered sign: (-1)^{sum of multi-indices}
    parity = (-1) ** (mi.sum(axis=1) % 2)  # (n_nodes,) of +1/-1

    # Lateral component index: last ambient component.
    # The ambient has m_ambient = dim + 1 components; the last one (index m_ambient-1)
    # is the timelike/amplitude direction as seen by the inside observer (backbone #22).
    lat_comp = m_ambient - 1

    # Seed period: use the dimension-aware exact nonlinear frequency so the seed
    # is in the basin of attraction of the true above-band breather for any ndim.
    # Using T = 2π/ω_max would seed at the band edge and Newton may converge
    # to the linear mode instead of the nonlinear breather.
    T_seed = 2.0 * math.pi / om_nd

    # Build slices
    slices = np.zeros((P, n_nodes, m_ambient), dtype=np.float64)
    for l_idx in range(P):
        phase_time = math.cos(2.0 * math.pi * l_idx / P)
        slices[l_idx] = ref.copy()
        slices[l_idx, :, lat_comp] += envelope * parity * phase_time

    return slices, T_seed


# ---------------------------------------------------------------------------
# Cyclic residual (pure function — no state mutation)
# ---------------------------------------------------------------------------


def _cyclic_residual(
    slices: np.ndarray,
    T: float,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    P: int,
) -> np.ndarray:
    """Compute the cyclic time-collocation residual for all P slices.

    ℛ^l = m (R^{l+1} − 2R^l + R^{l−1}) / dt_eff² − F(R^l)

    with R^P ≡ R^0 and R^{−1} ≡ R^{P−1} (cyclic).

    Parameters
    ----------
    slices : ndarray, shape (P, n_nodes, m_ambient)
    T : float  — current period estimate
    lattice, params_like, mass, P : physics parameters

    Returns
    -------
    res : ndarray, shape (P, n_nodes, m_ambient)
    """
    dt_eff = T / P
    dt_eff2 = dt_eff * dt_eff

    n_nodes, m_ambient = slices.shape[1], slices.shape[2]
    res = np.empty((P, n_nodes, m_ambient), dtype=np.float64)

    for l_idx in range(P):
        R_prev = slices[(l_idx - 1) % P]
        R_curr = slices[l_idx]
        R_next = slices[(l_idx + 1) % P]

        accel = (R_next - 2.0 * R_curr + R_prev) / dt_eff2
        F = spacelike_force(R_curr, lattice, params_like)
        res[l_idx] = mass * accel - F

    return res


# ---------------------------------------------------------------------------
# Pack / unpack for the Newton-Krylov vector
# ---------------------------------------------------------------------------


def _pack(slices: np.ndarray, u: float) -> np.ndarray:
    """Pack (P, n_nodes, m_ambient) slices + log-period scalar u = ln(T) into a flat vector.

    The Newton-Krylov unknown for the period is u = ln(T), not T directly.
    This ensures T = exp(u) > 0 for any real u without any hard clamp.
    """
    return np.concatenate([slices.ravel(), [u]])


def _unpack(
    z: np.ndarray, P: int, n_nodes: int, m_ambient: int
) -> tuple[np.ndarray, float]:
    """Unpack flat vector z → (slices, u) where u = ln(T) is the log-period unknown.

    T = exp(u) > 0 always.  The caller is responsible for converting u → T.
    """
    n_dof = P * n_nodes * m_ambient
    slices = z[:n_dof].reshape(P, n_nodes, m_ambient)
    u = float(z[n_dof])
    return slices, u


# ---------------------------------------------------------------------------
# Residual function factory for Newton-Krylov
# ---------------------------------------------------------------------------


def _make_G(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    P: int,
    amplitude: float,
    peak_node: int,
    lat_comp: int,
    peak_sign: float,
    opts: BreatherOpts,
) -> "callable":
    """Return the function G(z) = 0 that the JFNK solver drives to zero.

    The flat vector z = [flat(R^0..R^{P-1}), u] has shape P*n*m + 1,
    where u = ln(T) is the log-reparametrisation of the period.  This
    guarantees T = exp(u) > 0 for any real u without any hard clamp,
    keeping G smooth and continuous in all its arguments.

    G has P*n*m + 1 components:

      G[amp0_idx]   = R^0[peak,lat] * peak_sign − A          (amplitude at l=0)
      G[ampP2_idx]  = R^{P/2}[peak,lat] * peak_sign − (−A)  (anti-amplitude at l=P/2)
      G[P*n*m]      = R^1[peak,lat] − R^{P-1}[peak,lat]     (gauge: turning point at l=0)
      G[all other]  = cyclic residual ℛ^l at that (l,node,comp)

    The two amplitude constraints (l=0 and l=P/2) together:
      - Fix the amplitude of the oscillation (selecting the non-trivial branch).
      - Rule out the static equilibrium solution (where both would equal +A).
      - Break the half-period time-translation symmetry.
    The gauge fix breaks the remaining residual time-translation zero mode
    (the full-period shift by an arbitrary number of slices).

    P must be even so that l=P/2 is a valid integer slice index.
    """
    if P % 2 != 0:
        raise ValueError(f"P must be even for the half-period anti-amplitude constraint; got P={P}")

    n_nodes = lattice.n_nodes
    m_ambient = params_like.m_ambient or (lattice.dim + 1)

    n_per_slice = n_nodes * m_ambient

    # Flat index of the l=0 amplitude-constraint slot (replaces residual at l=0, peak, lat)
    amp0_idx = 0 * n_per_slice + peak_node * m_ambient + lat_comp

    # Flat index of the l=P/2 anti-amplitude slot (replaces residual at l=P/2, peak, lat)
    ampP2_idx = (P // 2) * n_per_slice + peak_node * m_ambient + lat_comp

    # Flat indices for gauge fix: R^1[peak, lat] and R^{P-1}[peak, lat]
    idx_R1_peak_lat = 1 * n_per_slice + peak_node * m_ambient + lat_comp
    idx_RP1_peak_lat = (P - 1) * n_per_slice + peak_node * m_ambient + lat_comp

    def G(z: np.ndarray) -> np.ndarray:
        slices, u = _unpack(z, P, n_nodes, m_ambient)

        # Log-reparametrisation: T = exp(u) is always positive for any real u.
        # No hard clamp needed — the smooth exponential map replaces np.clip.
        T = float(math.exp(u))

        # Compute cyclic residual (main equations)
        res = _cyclic_residual(slices, T, lattice, params_like, mass, P)
        g = res.ravel().copy()  # shape P*n*m

        # Replace l=0 residual at (peak, lat) with amplitude constraint:
        # R^0[peak, lat] * peak_sign = +amplitude
        g[amp0_idx] = slices[0, peak_node, lat_comp] * peak_sign - amplitude

        # Replace l=P/2 residual at (peak, lat) with anti-amplitude constraint:
        # R^{P/2}[peak, lat] * peak_sign = -amplitude
        # (At half period, the cosine carrier is at -A; this rules out static solutions.)
        g[ampP2_idx] = slices[P // 2, peak_node, lat_comp] * peak_sign + amplitude

        # Append gauge fix: R^1[peak,lat] - R^{P-1}[peak,lat] = 0
        # (Turning-point condition: zero velocity of the peak node at l=0.)
        g_gauge = float(z[idx_R1_peak_lat] - z[idx_RP1_peak_lat])

        return np.append(g, g_gauge)

    return G


# ---------------------------------------------------------------------------
# Topological-breather (Skyrmion-carrier) seed builder
# ---------------------------------------------------------------------------


def breather_seed_skyrmion(
    lattice: SpacelikeLattice,
    m: int,
    u0: float,
    w: float,
    P: int,
    profile: str = "power2",
) -> tuple[np.ndarray, float]:
    """Build a topological-breather seed worldtube for the 3D skyrmion carrier.

    The seed implements the "common-carrier breathing skyrmion" ansatz:

        R_p^l = ref_p + u0 * cos(2π l/P)
                          * [ sin F(r_p) * x̂_p^i  (components 0..dim-1)
                            ; cos F(r_p)            (component dim = X⁴) ]

    where F(r) is the Skyrme profile (power2: F(r) = π/(1+(r/w)²)),
    x̂_p = (coords_p - centre)/|coords_p - centre|, and the WHOLE
    (lateral, X⁴) 4-vector shares ONE scalar carrier  b(l) = cos(2π l/P).

    The period seed is T_seed = 2π/ω_neon where ω_neon is the
    dimension-aware band-top frequency (as in _build_seed).

    Physical properties of the seed:
    - At l=0: carrier = 1  →  reduces exactly to skyrme_twisted_hedgehog.
    - At l=P/2: carrier = -1 → the anti-amplitude slice.
    - B=1 winding preserved: degree depends only on the direction map
      (sinF * x̂, cosF) which is time-frozen; the carrier rescales S³
      radius but does not change the winding.
    - The X⁴ component at the peak (r=0) node at l=0 is
      u0 * cos(F(0)) = u0 * cos(π) = -u0  (nonzero → clean pin point).

    Parameters
    ----------
    lattice : SpacelikeLattice
        Must have dim=3 (the Skyrme winding requires S² → S² for B=1).
        Works at any grid size.
    m : int
        Ambient dimension.  Must be >= dim+1 = 4 for the X⁴ channel.
    u0 : float
        S³ radius (carrier amplitude); the peak displacement magnitude.
    w : float
        Skyrme profile half-width in the same units as the lattice spacing.
    P : int
        Number of temporal slices per period.  Must be even.
    profile : str
        Skyrme profile for F(r): ``"power2"`` or ``"tanh"``.

    Returns
    -------
    slices : ndarray, shape (P, n_nodes, m)
        Seed worldtube.
    T_seed : float
        Seed period estimate (2π / band-top frequency).
    """
    if P % 2 != 0:
        raise ValueError(f"P must be even; got P={P}")
    if m < lattice.dim + 1:
        raise ValueError(
            f"breather_seed_skyrmion requires m >= dim+1 = {lattice.dim + 1}; got m={m}"
        )

    n_nodes = lattice.n_nodes
    k_s = 1.0  # the seed T estimate only needs the band top; caller can override
    # Use dim-aware band-top as the seed frequency (same formula as _build_seed)
    # ω²_max(ndim) = 4*dim*(1-α)*k_s/m  but for the seed T we use
    # a conservative estimate at the band-top edge.
    # We delegate to the lattice params: if params were available we would use
    # them, but here we just compute the seed T from the number of slices and
    # a default band-top estimate, keeping the breather_seed_skyrmion API
    # self-contained.  The caller (solve_breather) can override T_seed.
    #
    # Conservative estimate: use the transverse band-top at α=0.5 in 3D
    # as a safe starting point.  The exact value doesn't matter much; Newton
    # will correct it.  We compute it analytically so no lattice-action
    # dependency creeps in.
    dim = lattice.dim
    # Approximate band-top for seed: use staggered 3D mode at alpha=0.5, k_s=1, m=1
    # ω²_max = 4*dim*(1-0.5)*1/1 for a rough estimate.
    # The real period is a free unknown; this just gives the seed T.
    om_approx = math.sqrt(4.0 * dim * 0.5)  # ≈ 2*sqrt(dim)
    T_seed = 2.0 * math.pi / om_approx

    # Build the static (l=0) Skyrme slice: this is exactly skyrme_twisted_hedgehog
    R0_static, _ = skyrme_twisted_hedgehog(
        lattice, m=m, u0=u0, w=w, profile_shape=profile
    )
    ref = lattice.reference_positions(m)  # (n_nodes, m)
    disp_static = R0_static - ref  # (n_nodes, m)  — the static displacement

    # Build P slices by scaling the static displacement with the carrier
    slices = np.zeros((P, n_nodes, m), dtype=np.float64)
    for l_idx in range(P):
        carrier = math.cos(2.0 * math.pi * l_idx / P)
        slices[l_idx] = ref + carrier * disp_static

    return slices, T_seed


# ---------------------------------------------------------------------------
# Topological-breather multi-component constraint factory
# ---------------------------------------------------------------------------


def _make_G_topological(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    P: int,
    u0: float,
    peak_node: int,
    x4_comp: int,
    x4_pin_value: float,
    opts: BreatherOpts,
) -> "callable":
    """Return G(z)=0 for the 4-component topological breather solve.

    Constraint scheme (makes the system square):
    - 1 amplitude pin:     R^0[peak, x4_comp] = x4_pin_value   (X⁴ at peak, l=0)
                           x4_pin_value = u0 * cos(F(0)) = u0 * cos(π) = -u0
                           (nonzero and sign-definite → clean, well-conditioned pin)
    - 1 anti-amplitude:    R^{P/2}[peak, x4_comp] = -x4_pin_value   (= +u0 at half-period)
                           This rules out static solutions and selects the breathing branch.
    - 1 gauge fix:         R^1[peak, x4_comp] - R^{P-1}[peak, x4_comp] = 0
                           (turning point: carrier velocity = 0 at l=0)

    The 3 lateral components (indices 0..dim-1) at (peak, l=0) and (peak, l=P/2)
    are NOT independently pinned.  They ride the SAME carrier as X⁴ via
    Newton convergence; no independent lateral pin is imposed.  This is the
    "common-carrier" discipline: 3 constraints total (1 pin + 1 anti-amp + 1 gauge),
    same as the scalar breather, so the system remains square.

    Parameters
    ----------
    peak_node : int
        Index of the centre node (r=0).
    x4_comp : int
        Ambient component index for X⁴ (= dim = m_ambient-1).
    x4_pin_value : float
        The target value for R^0[peak, x4_comp].  Should equal u0*cos(F(0)) = -u0.
    """
    if P % 2 != 0:
        raise ValueError(f"P must be even; got P={P}")

    n_nodes = lattice.n_nodes
    m_ambient = params_like.m_ambient or (lattice.dim + 1)
    n_per_slice = n_nodes * m_ambient

    # Flat index for X⁴ component at the peak node
    # l=0 slice: amplitude pin
    amp0_idx = 0 * n_per_slice + peak_node * m_ambient + x4_comp
    # l=P/2 slice: anti-amplitude pin
    ampP2_idx = (P // 2) * n_per_slice + peak_node * m_ambient + x4_comp

    # Gauge fix indices: l=1 and l=P-1 at (peak, x4_comp)
    idx_R1_peak_x4 = 1 * n_per_slice + peak_node * m_ambient + x4_comp
    idx_RP1_peak_x4 = (P - 1) * n_per_slice + peak_node * m_ambient + x4_comp

    def G(z: np.ndarray) -> np.ndarray:
        slices, u = _unpack(z, P, n_nodes, m_ambient)
        T = float(math.exp(u))

        # Cyclic residual (P*n*m equations)
        res = _cyclic_residual(slices, T, lattice, params_like, mass, P)
        g = res.ravel().copy()

        # Replace l=0 residual at (peak, X⁴) with the amplitude pin:
        # R^0[peak, X⁴] = x4_pin_value  (= -u0 for the Skyrme carrier)
        g[amp0_idx] = slices[0, peak_node, x4_comp] - x4_pin_value

        # Replace l=P/2 residual at (peak, X⁴) with anti-amplitude pin:
        # R^{P/2}[peak, X⁴] = -x4_pin_value  (= +u0 at half-period)
        g[ampP2_idx] = slices[P // 2, peak_node, x4_comp] + x4_pin_value

        # Append gauge fix: turning point at l=0 for the X⁴ carrier
        g_gauge = float(z[idx_R1_peak_x4] - z[idx_RP1_peak_x4])

        return np.append(g, g_gauge)

    return G


# ---------------------------------------------------------------------------
# Main solver
# ---------------------------------------------------------------------------


def solve_breather(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    P: int,
    amplitude: float,
    seed: np.ndarray | None = None,
    T_seed: float | None = None,
    opts: BreatherOpts | None = None,
    mode: str = "scalar",
    skyrme_profile: str = "power2",
    skyrme_w: float | None = None,
) -> dict[str, Any]:
    """Find a time-periodic discrete breather by JFNK root-finding.

    NEVER minimises the action S (Lorentzian saddle).
    Drives ‖G‖ → 0 where G is the cyclic time-collocation residual +
    amplitude constraint + gauge fix.  (OBJECTIVE == "residual_norm".)

    Parameters
    ----------
    lattice : SpacelikeLattice
        Spacelike topology (any dimension d).
    params_like : ActionParams
        Carries k_s, alpha, rho, a, m_ambient.  The dt and n_slices fields
        are ignored; the period T is the continuous unknown.
    mass : float
        Node mass m (= rho * a^dim for uniform lattice).
    P : int
        Number of temporal slices per period.  Must be even (required for the
        half-period anti-amplitude constraint that rules out static solutions);
        typically 16..64.
    amplitude : float
        Peak displacement amplitude at l=0.
        In ``"scalar"`` mode: peak lateral displacement A.
        In ``"topological"`` mode: S³ radius u0 (same as skyrme_twisted_hedgehog u0).
    seed : ndarray, shape (P, n_nodes, m_ambient), optional
        Initial guess for the slices.  If None, built automatically:
        - ``"scalar"``: staggered Gaussian envelope × cos carrier.
        - ``"topological"``: Skyrme-twisted hedgehog × cos carrier (via
          :func:`breather_seed_skyrmion`).
    T_seed : float, optional
        Initial period guess.  Used when seed is not None.  If None, defaults
        to 2π / ω_max (scalar) or 2π / band-top estimate (topological).
    opts : BreatherOpts, optional
        Solver options (tolerances, iteration counts, etc.).
    mode : str
        Breather mode.  One of:

        ``"scalar"`` (default)
            Single-component staggered transverse breather.  Pins the LAST
            ambient component (X⁴) at the peak node.  Works for any dimension.

        ``"topological"``
            Common-carrier Skyrme skyrmion breather.  Pins the X⁴ component
            (amplitude channel, component dim) at the peak (r=0) node via the
            Skyrme-profile value cos(F(0)) = cos(π) = -1.  Requires
            m_ambient >= dim+1 (canonical: dim=3, m=4).  The 3 lateral
            components ride the SAME carrier implicitly; only X⁴ is pinned.
    skyrme_profile : str
        Profile shape for the topological seed's F(r) function.
        ``"power2"`` (default) or ``"tanh"``.  Ignored in ``"scalar"`` mode.
    skyrme_w : float, optional
        Skyrme profile half-width for the topological seed.  Defaults to
        ``1.5 * lattice.params.spacing`` if None.  Ignored in ``"scalar"`` mode.

    Returns
    -------
    dict with keys:
        slices       : ndarray, shape (P, n_nodes, m_ambient) — real float64
        T            : float — converged period
        omega        : float — converged angular frequency 2π/T
        residual_norm : float — ‖G‖ at the solution
        converged    : bool
        objective    : str — always "residual_norm" (saddle discipline)
        walltime_s   : float
        peak_node    : int
        lat_comp     : int
        mode         : str — echoed mode
    """
    if opts is None:
        opts = BreatherOpts()

    if P % 2 != 0:
        raise ValueError(f"P must be even (required for half-period constraint); got P={P}")
    if mode not in ("scalar", "topological"):
        raise ValueError(f"mode must be 'scalar' or 'topological'; got {mode!r}")

    k_s = params_like.k_s
    alpha = params_like.alpha
    a = lattice.params.spacing
    m_ambient = params_like.m_ambient or (lattice.dim + 1)
    n_nodes = lattice.n_nodes

    # Identify peak node (center of the grid, r=0)
    mi = lattice.multi_indices  # (n_nodes, dim)
    center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
    d_to_center = np.linalg.norm((mi - center_mi) * a, axis=1)
    peak_node = int(np.argmin(d_to_center))

    # Last ambient component (timelike / X⁴ direction, backbone #22)
    lat_comp = m_ambient - 1

    if mode == "scalar":
        # ---- SCALAR PATH (existing 1D/2D/3D staggered breather) ----
        # Sign of the staggered carrier at the peak node
        peak_parity = int(mi[peak_node].sum() % 2)
        peak_sign = float((-1) ** peak_parity)

        # Build seed
        if seed is None:
            slices_seed, T_seed_auto = _build_seed(
                lattice, k_s, alpha, mass, a, P, amplitude, m_ambient
            )
            if T_seed is None:
                T_seed = T_seed_auto
        else:
            slices_seed = seed.copy()
            if T_seed is None:
                T_seed = 2.0 * math.pi / omega_max(k_s, alpha, mass, a)

        G = _make_G(
            lattice, params_like, mass, P, amplitude,
            peak_node, lat_comp, peak_sign, opts,
        )

    else:
        # ---- TOPOLOGICAL PATH (Skyrme common-carrier breather) ----
        if m_ambient < lattice.dim + 1:
            raise ValueError(
                f"Topological mode requires m_ambient >= dim+1 = {lattice.dim + 1}; "
                f"got m_ambient={m_ambient}"
            )

        # X⁴ component index (the amplitude channel)
        x4_comp = lattice.dim  # = m_ambient - 1 for canonical m = dim+1

        # Skyrme profile width: default 1.5*a
        if skyrme_w is None:
            skyrme_w = 1.5 * a

        # Build the topological seed (Skyrme carrier)
        if seed is None:
            slices_seed, T_seed_auto = breather_seed_skyrmion(
                lattice, m=m_ambient, u0=amplitude,
                w=skyrme_w, P=P, profile=skyrme_profile,
            )
            if T_seed is None:
                T_seed = T_seed_auto
        else:
            slices_seed = seed.copy()
            if T_seed is None:
                # Rough estimate: use approximate band-top
                dim = lattice.dim
                T_seed = 2.0 * math.pi / math.sqrt(4.0 * dim * 0.5)

        # The X⁴ pin value: u0 * cos(F(0)) = amplitude * cos(π) = -amplitude
        # F(0) = π for both power2 and tanh profiles (Skyrme boundary cond.)
        x4_pin_value = -amplitude  # = u0 * cos(π) = -u0

        G = _make_G_topological(
            lattice, params_like, mass, P, amplitude,
            peak_node, x4_comp, x4_pin_value, opts,
        )

        # For output consistency: lat_comp is still the last ambient component
        lat_comp = x4_comp

    # ---- Common Newton-Krylov solve path ----

    # Log-reparametrisation: T = exp(u) > 0 for any real u — no hard clamp.
    u0_log = math.log(T_seed)
    z0 = _pack(slices_seed, u0_log)

    # Measure initial residual
    g0 = G(z0)
    res_init = float(np.linalg.norm(g0))

    t0 = time.perf_counter()
    converged = False
    z_sol = z0.copy()
    res_final = res_init

    try:
        z_sol = newton_krylov(
            G,
            z0,
            method=opts.method,
            verbose=opts.verbose,
            f_tol=opts.tol,
            iter=opts.max_iter,
            inner_maxiter=opts.inner_maxiter,
        )
        g_final = G(z_sol)
        res_final = float(np.linalg.norm(g_final))
        converged = res_final <= opts.tol
    except NoConvergence as exc:
        if exc.args:
            z_sol = exc.args[0]
        g_sol = G(z_sol)
        res_final = float(np.linalg.norm(g_sol))
        converged = False

    elapsed = time.perf_counter() - t0

    # Recover T from the log-reparametrised unknown u = ln(T).
    slices_sol, u_sol = _unpack(z_sol, P, n_nodes, m_ambient)
    T_sol = math.exp(u_sol)
    omega_sol = 2.0 * math.pi / T_sol

    return {
        "slices": slices_sol.astype(np.float64),
        "T": float(T_sol),
        "omega": float(omega_sol),
        "residual_norm": res_final,
        "converged": converged,
        "objective": OBJECTIVE,  # always "residual_norm" — saddle discipline
        "walltime_s": elapsed,
        "peak_node": peak_node,
        "lat_comp": lat_comp,
        "residual_initial": res_init,
        "mode": mode,
    }


# ---------------------------------------------------------------------------
# Amplitude continuation
# ---------------------------------------------------------------------------


def continue_breather(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    P: int,
    amplitudes: list[float] | np.ndarray,
    opts: BreatherOpts | None = None,
) -> list[dict[str, Any]]:
    """Amplitude continuation of the discrete breather.

    Solves ``solve_breather`` for each amplitude in ``amplitudes`` in order,
    using the previous solution as a warm-start for the next amplitude.
    The continuation proceeds even if a step fails to converge (the failed
    result is included with ``converged=False``).

    Each amplitude step uses the Gaussian-envelope seed (not scaled from the
    previous solution), but the seed period T_seed is always set to
    2π / ω_exact(A) so that Newton starts in the basin of attraction of the
    above-band nonlinear breather rather than at the linear band edge.
    Empirically, the scaled-previous-solution warm-start gives a worse initial
    residual than the Gaussian seed, so each step is solved independently.

    Parameters
    ----------
    lattice, params_like, mass, P : see solve_breather.
    amplitudes : sequence of floats
        Amplitude values to solve at, in the order provided.
        Typically an increasing sequence starting from a small value.
    opts : BreatherOpts, optional

    Returns
    -------
    list of dicts, one per amplitude, each as returned by solve_breather
    with an additional key ``amplitude``.
    """
    if opts is None:
        opts = BreatherOpts()

    results = []

    for A in amplitudes:
        A_float = float(A)

        # Each step uses the Gaussian-envelope seed (fresh) and T seeded from
        # omega_exact(A).  The scaled-previous-solution warm-start gives a
        # higher initial residual than the Gaussian seed for this problem.
        result = solve_breather(
            lattice, params_like, mass,
            P=P,
            amplitude=A_float,
            seed=None,   # always use fresh Gaussian seed
            T_seed=None,  # _build_seed computes T from omega_exact(A)
            opts=opts,
        )
        result["amplitude"] = A_float
        results.append(result)

    return results


# ---------------------------------------------------------------------------
# Floquet stability of the periodic orbit  (existence ≠ stability)
# ---------------------------------------------------------------------------


def _monodromy_matvec(
    z: np.ndarray,
    slices: np.ndarray,
    T: float,
    mass: float,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    fd_eps: float,
) -> np.ndarray:
    """Apply the one-period monodromy M to a state perturbation z = [δR^0; δR^{-1}].

    The breather orbit satisfies the discrete (leapfrog) equation of motion
        m (R^{l+1} − 2R^l + R^{l−1}) / dt_eff² = F(R^l),   dt_eff = T/P.
    Its linearisation about the orbit is the variational recurrence
        δR^{l+1} = 2δR^l − δR^{l−1} + (dt_eff²/m) · J(R^l) · δR^l
    with J = dF/dR.  Marching this for one full period (l = 0 … P−1, cyclic)
    maps the 2N-state z^0 = [δR^0; δR^{−1}] → z^P = M z^0.  J(R^l)·δR^l is
    formed matrix-free (two force evaluations per slice).

    This is the variational stability of the SAME time-collocation
    discretisation that produced the orbit (a symplectic leapfrog map), so its
    multipliers measure the stability of the actual computed orbit and converge
    to the continuum Floquet problem as P → ∞.
    """
    P, n_nodes, m_ambient = slices.shape
    N = n_nodes * m_ambient
    dt2_over_m = (T / P) ** 2 / mass

    dR = z[:N].reshape(n_nodes, m_ambient).copy()
    dR_prev = z[N:].reshape(n_nodes, m_ambient).copy()

    for l_idx in range(P):
        JdR = _jacobian_matvec(slices[l_idx], dR, lattice, params_like, fd_eps=fd_eps)
        dR_next = 2.0 * dR - dR_prev + dt2_over_m * JdR
        dR_prev = dR
        dR = dR_next

    return np.concatenate([dR.ravel(), dR_prev.ravel()])


def _spectral_radius_power(
    M_apply: "callable",
    n_state: int,
    *,
    n_iter: int = 200,
    burn_frac: float = 0.5,
    n_restarts: int = 2,
    seed: int = 0,
) -> float:
    """Spectral radius of the monodromy by normalised power iteration.

    Robust, matrix-free estimate of max|ρ| = the per-period growth factor — the
    quantity the stability verdict needs.  Never fails (unlike ARPACK on a
    symplectic spectrum clustered on the unit circle).

    The per-step growth ‖M z_k‖ converges to max|ρ| only AFTER the iterate aligns
    with the dominant invariant subspace.  Because the monodromy is non-normal,
    the early steps overshoot (transient/pseudospectral amplification), so we
    discard a burn-in fraction and geometric-mean only the tail — this removes
    the few-percent bias a full-sequence average would carry.  A complex
    dominant pair makes the tail oscillate around max|ρ|; the geometric mean
    still recovers it.  A couple of random restarts guard against an unlucky
    start vector orthogonal to the dominant eigenspace.  Deterministic seed.
    """
    rng = np.random.default_rng(seed)
    n_burn = int(burn_frac * n_iter)
    best = 0.0
    for _ in range(n_restarts):
        z = rng.standard_normal(n_state)
        z /= np.linalg.norm(z)
        logs: list[float] = []
        for _ in range(n_iter):
            z = M_apply(z)
            nrm = float(np.linalg.norm(z))
            if nrm == 0.0 or not math.isfinite(nrm):
                break
            logs.append(math.log(nrm))
            z /= nrm
        tail = logs[n_burn:] if len(logs) > n_burn else logs
        if tail:
            best = max(best, math.exp(sum(tail) / len(tail)))
    return best


def floquet_multipliers(
    slices: np.ndarray,
    T: float,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    n_multipliers: int = 6,
    dense_threshold: int = 400,
    want_multipliers: bool = False,
    power_iter: int = 200,
    fd_eps: float = 1e-6,
    stability_tol: float = 1e-2,
) -> dict[str, Any]:
    """Floquet (monodromy) stability of a converged breather orbit.

    Builds the one-period monodromy M of the variational leapfrog map and
    returns its eigenvalues (Floquet multipliers ρ).  The orbit is linearly
    stable iff every multiplier lies on the unit circle; an instability shows
    up as a multiplier with |ρ| > 1 (the per-period growth factor).

    For a Hamiltonian (symplectic) orbit the multipliers come in reciprocal
    conjugate quartets {ρ, 1/ρ, ρ̄, 1/ρ̄}, and there is always a marginal pair
    at ρ = 1 (time-translation / the amplitude-family direction); these have
    |ρ| = 1 and do NOT count as instabilities.  ``stability_tol`` absorbs the
    finite-difference + discretisation drift of that marginal pair.

    Cost & method: 2N = 2·n_nodes·m_ambient state dimension.  If 2N ≤
    ``dense_threshold`` the full spectrum is assembled (2N matvecs) and the
    radius is exact.  Otherwise (large 3D) the radius is obtained matrix-free by
    normalised power iteration (``_spectral_radius_power``), which never fails;
    the individual dominant multipliers are an OPT-IN bonus (``want_multipliers``
    → ARPACK, off by default because it is slow/unreliable on this clustered
    symplectic spectrum).  The returned ``method`` records which path ran.

    ACCURACY CAVEAT (matrix-free path): a genuinely stable orbit has ALL |ρ| = 1
    (no spectral gap), so power iteration converges only algebraically — its
    worst-case positive bias is ~1% at the default ``power_iter`` and halves per
    doubling of it.  Consequences: the gross bare-breather instability (ρ ≈
    1.5–2) and the predicted topological "2× drop" are detected robustly, but a
    DEFINITIVE stable/unstable call at the |ρ|−1 ~ 1% level needs a larger
    ``power_iter`` (≈240–480) or the dense path on a small enough lattice.

    Parameters
    ----------
    slices : ndarray, shape (P, n_nodes, m_ambient)  — the converged orbit.
    T : float            — converged period.
    n_multipliers : int  — # dominant multipliers requested when want_multipliers.
    dense_threshold : int — assemble the full 2N×2N M when 2N ≤ this.
    want_multipliers : bool — also attempt ARPACK for the multiplier list
                              (matrix-free path only; default False).
    power_iter : int     — power-iteration steps for the matrix-free radius
                           (default 200; raise to tighten a borderline verdict).
    stability_tol : float — |ρ| − 1 ≤ stability_tol counts as on the unit circle.

    Returns
    -------
    dict with keys:
        multipliers      : ndarray (complex) — multipliers, |ρ| descending
                           (EMPTY on the matrix-free path unless want_multipliers).
        spectral_radius  : float  — max |ρ|.
        growth_per_period: float  — == spectral_radius (factor per period T).
        growth_rate      : float  — ln(spectral_radius)/T (per unit time).
        stable           : bool   — spectral_radius ≤ 1 + stability_tol.
        n_unstable       : int    — # multipliers with |ρ| > 1 + stability_tol
                           (inferred from the radius if no multipliers computed).
        dense            : bool    — whether the full spectrum was computed.
        method           : str     — "dense" | "power" | "arnoldi+power".
        n_state          : int     — 2N.
    """
    P, n_nodes, m_ambient = slices.shape
    N = n_nodes * m_ambient
    n_state = 2 * N

    def M_apply(z: np.ndarray) -> np.ndarray:
        return _monodromy_matvec(z, slices, T, mass, lattice, params_like, fd_eps)

    method = "dense"
    if n_state <= dense_threshold:
        # Assemble the full monodromy by applying M to each basis vector.
        M = np.empty((n_state, n_state), dtype=np.float64)
        e = np.zeros(n_state)
        for j in range(n_state):
            e[j] = 1.0
            M[:, j] = M_apply(e)
            e[j] = 0.0
        multipliers = np.linalg.eigvals(M)
        dense = True
        order = np.argsort(-np.abs(multipliers))
        multipliers = multipliers[order]
        spectral_radius = float(np.abs(multipliers[0]))
    else:
        dense = False
        # The spectral radius (per-period growth factor) is the verdict driver;
        # compute it robustly by normalised power iteration.  This NEVER fails
        # and directly measures max|ρ|, independent of ARPACK's trouble with a
        # symplectic spectrum clustered on the unit circle (which can converge
        # ZERO eigenvalues and grind for thousands of iterations).
        spectral_radius = _spectral_radius_power(M_apply, n_state, n_iter=power_iter)
        method = "power"
        multipliers = np.empty(0, dtype=complex)

        # ARPACK is an OPT-IN bonus (want_multipliers): it yields the individual
        # dominant multipliers when it converges, but is slow/unreliable on this
        # clustered symplectic spectrum, so it is OFF by default — the verdict
        # already stands on the power-iteration radius.  A fail-fast budget keeps
        # it from grinding when it cannot converge.
        if want_multipliers:
            from scipy.sparse.linalg import (
                LinearOperator, eigs, ArpackNoConvergence,
            )

            op = LinearOperator((n_state, n_state), matvec=M_apply, dtype=np.float64)
            k = min(n_multipliers, n_state - 2)
            ncv = min(n_state - 1, max(2 * k + 1, 20))
            try:
                multipliers = eigs(
                    op, k=k, which="LM", ncv=ncv,
                    maxiter=300, return_eigenvectors=False,
                )
                method = "arnoldi+power"
            except ArpackNoConvergence as exc:
                multipliers = exc.eigenvalues
            if multipliers.size:
                order = np.argsort(-np.abs(multipliers))
                multipliers = multipliers[order]
                # Power iteration is the authority on the radius; ARPACK can miss
                # the dominant when it only partially converges, so take the max.
                spectral_radius = max(spectral_radius, float(np.abs(multipliers[0])))

    # Count clearly-unstable multipliers when we have them; otherwise infer from
    # the always-available spectral radius.
    if multipliers.size:
        n_unstable = int(np.sum(np.abs(multipliers) > 1.0 + stability_tol))
    else:
        n_unstable = 1 if spectral_radius > 1.0 + stability_tol else 0

    growth_rate = (
        float(math.log(spectral_radius) / T)
        if (spectral_radius > 0 and T > 0)
        else float("nan")
    )

    return {
        "multipliers": multipliers,
        "spectral_radius": spectral_radius,
        "growth_per_period": spectral_radius,
        "growth_rate": growth_rate,
        "stable": spectral_radius <= 1.0 + stability_tol,
        "n_unstable": n_unstable,
        "dense": dense,
        "method": method,
        "n_state": n_state,
    }


def analyze_breather(
    result: dict[str, Any],
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    n_max_harmonics: int = 6,
    n_multipliers: int = 6,
    dense_threshold: int = 400,
    fd_eps: float = 1e-6,
    stability_tol: float = 1e-2,
) -> dict[str, Any]:
    """Run both post-solve diagnostics on a converged breather result.

    Convenience wrapper: takes the dict returned by :func:`solve_breather`
    and attaches a harmonic-resonance verdict (radiation continuum) and a
    Floquet-stability verdict (linear stability of the periodic orbit).

    Returns
    -------
    dict with keys ``resonance`` (see :func:`harmonic_resonance_check`) and
    ``floquet`` (see :func:`floquet_multipliers`), plus a top-level boolean
    ``radiationless_and_stable``.
    """
    resonance = harmonic_resonance_check(
        result["omega"], lattice, params_like, mass,
        n_max=n_max_harmonics, fd_eps=fd_eps,
    )
    floquet = floquet_multipliers(
        result["slices"], result["T"], lattice, params_like, mass,
        n_multipliers=n_multipliers, dense_threshold=dense_threshold,
        fd_eps=fd_eps, stability_tol=stability_tol,
    )
    return {
        "resonance": resonance,
        "floquet": floquet,
        "radiationless_and_stable": bool(
            resonance["radiationless"] and floquet["stable"]
        ),
    }


# ---------------------------------------------------------------------------
# Screening gate: solve → analyze → single verdict  (baryon-search front end)
# ---------------------------------------------------------------------------


def screen_breather(
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    P: int,
    amplitude: float,
    seed: np.ndarray | None = None,
    T_seed: float | None = None,
    opts: BreatherOpts | None = None,
    n_max_harmonics: int = 6,
    n_multipliers: int = 6,
    dense_threshold: int = 400,
    fd_eps: float = 1e-6,
    stability_tol: float = 1e-2,
) -> dict[str, Any]:
    """Solve a breather candidate and screen it on BOTH stability gates.

    The front end for a soliton/baryon search: a candidate is only interesting
    if it (a) converges to an exact periodic orbit, (b) is radiationless (no
    n ≥ 2 harmonic in the phonon continuum), and (c) is Floquet-stable.  This
    runs :func:`solve_breather` then, if converged, :func:`analyze_breather`,
    and collapses the result to a single ``verdict`` string so a sweep can rank
    candidates without re-deriving the logic at every call site.

    Channel-agnostic: pass any ``seed`` (e.g. a multi-channel hedgehog IC).  Note
    the current :func:`solve_breather` pins a SINGLE (peak, last-ambient)
    amplitude/gauge constraint, so genuinely multi-channel (in-brane hedgehog)
    candidates need the Phase-2 multi-component constraint generalisation before
    this gate is meaningful for them; for the validated single transverse carrier
    it is ready to use.

    Returns
    -------
    dict with keys:
        verdict      : str — one of
            "NO_CONVERGENCE" — solver did not reach tol.
            "RADIATING"      — converged but an n ≥ 2 harmonic is in band.
            "UNSTABLE"       — converged, radiationless, but Floquet-unstable.
            "PASS"           — converged, radiationless AND Floquet-stable.
        converged, residual_norm, omega, T : echoed solve metrics.
        radiationless, lowest_in_band_n_ge_2 : resonance summary (None if no conv).
        spectral_radius, stable, n_unstable  : Floquet summary (None if no conv).
        solve        : full :func:`solve_breather` result dict.
        analysis     : full :func:`analyze_breather` result (None if no conv).
    """
    res = solve_breather(
        lattice, params_like, mass,
        P=P, amplitude=amplitude, seed=seed, T_seed=T_seed, opts=opts,
    )

    out: dict[str, Any] = {
        "converged": res["converged"],
        "residual_norm": res["residual_norm"],
        "omega": res["omega"],
        "T": res["T"],
        "radiationless": None,
        "lowest_in_band_n_ge_2": None,
        "spectral_radius": None,
        "stable": None,
        "n_unstable": None,
        "solve": res,
        "analysis": None,
    }

    if not res["converged"]:
        out["verdict"] = "NO_CONVERGENCE"
        return out

    report = analyze_breather(
        res, lattice, params_like, mass,
        n_max_harmonics=n_max_harmonics, n_multipliers=n_multipliers,
        dense_threshold=dense_threshold, fd_eps=fd_eps, stability_tol=stability_tol,
    )
    resonance, floquet = report["resonance"], report["floquet"]
    out.update({
        "analysis": report,
        "radiationless": resonance["radiationless"],
        "lowest_in_band_n_ge_2": resonance["lowest_in_band_n_ge_2"],
        "spectral_radius": floquet["spectral_radius"],
        "stable": floquet["stable"],
        "n_unstable": floquet["n_unstable"],
    })

    if not resonance["radiationless"]:
        out["verdict"] = "RADIATING"
    elif not floquet["stable"]:
        out["verdict"] = "UNSTABLE"
    else:
        out["verdict"] = "PASS"
    return out
