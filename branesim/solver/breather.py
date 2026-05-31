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
(ARCHITECTURE.md §1.3, principles.md §1.2, OPEN_PROBLEMS.md A1.)

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


# ---------------------------------------------------------------------------
# Module-level saddle-discipline assertion
# ---------------------------------------------------------------------------

OBJECTIVE: str = "residual_norm"
assert OBJECTIVE == "residual_norm", (
    "breather.py: solver targets residual_norm (root-find), NOT the Lorentzian "
    "action S.  S is a saddle, unbounded below — minimising it would diverge or "
    "solve the wrong Euclidean problem.  (ARCHITECTURE.md §1.3, principles.md §1.2)"
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
        Peak lateral displacement A at l=0.  Fixes the non-trivial branch.
    seed : ndarray, shape (P, n_nodes, m_ambient), optional
        Initial guess for the slices.  If None, built from the small-amplitude
        linear mode (staggered Gaussian envelope × cos carrier).
    T_seed : float, optional
        Initial period guess.  Used when seed is not None.  If None, defaults
        to 2π / ω_max.
    opts : BreatherOpts, optional
        Solver options (tolerances, iteration counts, etc.).

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
    """
    if opts is None:
        opts = BreatherOpts()

    if P % 2 != 0:
        raise ValueError(f"P must be even (required for half-period constraint); got P={P}")

    k_s = params_like.k_s
    alpha = params_like.alpha
    a = lattice.params.spacing
    m_ambient = params_like.m_ambient or (lattice.dim + 1)
    n_nodes = lattice.n_nodes

    # Identify peak node (center of the chain) and lateral component
    mi = lattice.multi_indices  # (n_nodes, dim)
    center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
    d_to_center = np.linalg.norm((mi - center_mi) * a, axis=1)
    peak_node = int(np.argmin(d_to_center))

    # Lateral component: last ambient component (index m_ambient-1).
    # The ambient has m_ambient = dim + 1 components; the last one is the
    # timelike/amplitude direction seen by the inside observer (backbone #22).
    lat_comp = m_ambient - 1

    # Sign of the staggered carrier at the peak node
    peak_parity = int(mi[peak_node].sum() % 2)
    peak_sign = float((-1) ** peak_parity)

    # Build seed if not provided
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

    # Build the residual function
    G = _make_G(
        lattice, params_like, mass, P, amplitude,
        peak_node, lat_comp, peak_sign, opts,
    )

    # Log-reparametrisation: pack u0 = ln(T_seed) so that the Newton unknown
    # is u = ln(T) and T = exp(u) > 0 always — no hard clamp ever needed.
    u0 = math.log(T_seed)
    z0 = _pack(slices_seed, u0)

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


def floquet_multipliers(
    slices: np.ndarray,
    T: float,
    lattice: SpacelikeLattice,
    params_like: ActionParams,
    mass: float,
    *,
    n_multipliers: int = 6,
    dense_threshold: int = 400,
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

    Cost: 2N = 2·n_nodes·m_ambient state dimension.  If 2N ≤ ``dense_threshold``
    the full spectrum is assembled (2N matvecs) and returned.  Otherwise only
    the ``n_multipliers`` largest-magnitude multipliers are found matrix-free
    via Arnoldi (``scipy.sparse.linalg.eigs``) — enough for the spectral radius,
    hence the stability verdict, without forming M.

    Parameters
    ----------
    slices : ndarray, shape (P, n_nodes, m_ambient)  — the converged orbit.
    T : float            — converged period.
    n_multipliers : int  — # dominant multipliers in the matrix-free path.
    dense_threshold : int — assemble the full 2N×2N M when 2N ≤ this.
    stability_tol : float — |ρ| − 1 ≤ stability_tol counts as on the unit circle.

    Returns
    -------
    dict with keys:
        multipliers      : ndarray (complex) — multipliers, |ρ| descending.
        spectral_radius  : float  — max |ρ|.
        growth_per_period: float  — == spectral_radius (factor per period T).
        growth_rate      : float  — ln(spectral_radius)/T (per unit time).
        stable           : bool   — spectral_radius ≤ 1 + stability_tol.
        n_unstable       : int    — # multipliers with |ρ| > 1 + stability_tol.
        dense            : bool    — whether the full spectrum was computed.
        n_state          : int     — 2N.
    """
    P, n_nodes, m_ambient = slices.shape
    N = n_nodes * m_ambient
    n_state = 2 * N

    def M_apply(z: np.ndarray) -> np.ndarray:
        return _monodromy_matvec(z, slices, T, mass, lattice, params_like, fd_eps)

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
    else:
        from scipy.sparse.linalg import LinearOperator, eigs, ArpackNoConvergence

        op = LinearOperator((n_state, n_state), matvec=M_apply, dtype=np.float64)
        k = min(n_multipliers, n_state - 2)
        # The monodromy is symplectic: many multipliers cluster on the unit
        # circle (degenerate in magnitude), which makes ARPACK's "LM" search
        # slow to fully converge.  Use a generous Krylov subspace + iteration
        # budget, and on partial convergence keep the converged subset — for
        # which="LM" those ARE the largest-magnitude multipliers, so the
        # spectral radius (hence the stability verdict) is still reliable.
        ncv = min(n_state - 1, max(2 * k + 1, 20))
        try:
            multipliers = eigs(
                op, k=k, which="LM", ncv=ncv,
                maxiter=10 * n_state, return_eigenvectors=False,
            )
        except ArpackNoConvergence as exc:
            multipliers = exc.eigenvalues
            if multipliers.size == 0:
                raise
        dense = False

    # Sort by magnitude, descending.
    order = np.argsort(-np.abs(multipliers))
    multipliers = multipliers[order]
    mags = np.abs(multipliers)
    spectral_radius = float(mags[0])
    n_unstable = int(np.sum(mags > 1.0 + stability_tol))
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
