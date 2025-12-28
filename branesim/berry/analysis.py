"""branesim.berry.analysis

Dimension-agnostic Berry analysis for brane simulations.

This module provides a *pragmatic* U(1) (rank-1) Berry analysis mounted on top
of the brane's displacement/velocity time series.

At each brane lattice point p, choose E embedding components (E=2,3,4...). Let
u_p(t) ∈ R^E be the displacement and v_p(t) ∈ R^E the velocity.

Hardcoded carrier
-----------------
We define a hardcoded carrier angular frequency ω0 from the Compton wavelength
λ_C (as used in the photon experiments):

    ω0 = 2π c / λ_C

Complex "polarization" vector
-----------------------------
Given ω0, we define a complex polarization vector

    a(t) = √ω0 u(t) + i v(t)/√ω0      ∈ C^E

This is the multi-component analogue of the standard harmonic-oscillator
quadrature combination. A global real scaling of a(t) does not affect the
normalized ray.

Normalized ray
--------------

    ψ(t) = a(t) / ||a(t)||           ∈ C^E

where it is defined (||a|| above a threshold).

Pancharatnam phase increment (time)
----------------------------------
Between adjacent frames we compute, per point:

    Δγ_k = arg( <ψ_k | ψ_{k+1}> )

and the cumulative Berry-like phase (relative to the first frame):

    γ_{k+1} = γ_k + Δγ_k

We also report a discrete-time Berry connection proxy:

    A_t,k ≈ Δγ_k / Δt_k

Notes
-----
- This is **not** a full adiabatic-eigenspace Berry calculation; it is a
  convenient "Berry machinery" on top of a fixed complex structure (hardcoded
  ω0).
- Undefined regions are handled via an amplitude mask; consumers should use the
  provided alpha mask for visualization.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class BerryTimeSeries:
    """Outputs of the U(1) Berry analysis.

    Shapes
    ------
    Let T = number of frames, N = number of lattice points.

    - phase: (T, N) cumulative phase γ(t) [rad]
    - connection: (T, N) discrete connection proxy A_t(t) [rad/s]
      (last frame repeats the last computed value)
    - delta_phase: (T, N) per-step increments Δγ [rad]
      (last frame is 0)
    - amplitude: (T, N) ||a(t)|| used for masking
    - alpha: (T, N) recommended alpha ∈ [0,1] for "defined" visualization
    - overlap_abs: (T, N) |<ψ_k|ψ_{k+1}>| for k<T-1 (last frame repeats)
    - omega0: scalar ω0 used
    - times_s: (T,) physical times in seconds

    The normalized rays ψ(t) can be returned optionally via compute_berry_time_series(...).
    """

    phase: np.ndarray
    connection: np.ndarray
    delta_phase: np.ndarray
    amplitude: np.ndarray
    alpha: np.ndarray
    overlap_abs: np.ndarray
    omega0: float
    times_s: np.ndarray
    psi: Optional[np.ndarray] = None  # (T, N, E)


def compton_omega0() -> float:
    """Return the hardcoded carrier angular frequency ω0 from physical_constants.

    Uses:
        omega0 = 2*pi*c / lambda_C

    Returns:
        omega0 in rad/s

    Raises:
        ImportError if branesim.config.physical_constants is unavailable.
    """
    from branesim.config.physical_constants import PhysicalConstants

    constants = PhysicalConstants()
    return 2.0 * np.pi * float(constants.c) / float(constants.lambda_C)


def _as_numpy_frames(frames: Sequence[np.ndarray]) -> List[np.ndarray]:
    out: List[np.ndarray] = []
    for f in frames:
        a = np.asarray(f)
        if a.ndim != 2:
            raise ValueError(f"Each frame must be a 2D array (N,E). Got shape {a.shape}.")
        out.append(a)
    return out


def complex_amplitude_from_u_v(u: np.ndarray, v: np.ndarray, omega0: float) -> np.ndarray:
    """Compute complex amplitude a = √ω0 u + i v/√ω0.

    Parameters
    ----------
    u, v:
        Arrays of shape (N,E) in consistent physical units (e.g. meters and m/s).
    omega0:
        Carrier angular frequency in rad/s.

    Returns
    -------
    a:
        Complex array of shape (N,E).
    """
    u = np.asarray(u)
    v = np.asarray(v)
    if u.shape != v.shape:
        raise ValueError(f"u and v must have the same shape. Got {u.shape} vs {v.shape}.")
    if u.ndim != 2:
        raise ValueError(f"u and v must be 2D arrays (N,E). Got ndim={u.ndim}.")

    w = float(omega0)
    if w <= 0:
        raise ValueError(f"omega0 must be positive. Got {omega0}.")

    sw = np.sqrt(w)
    return sw * u.astype(np.float64) + 1j * (v.astype(np.float64) / sw)


def _norm_rows(x: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    return np.sqrt(np.maximum(np.sum(np.abs(x) ** 2, axis=-1), 0.0) + eps * 0.0)


def compute_berry_time_series(
    frames_u: Sequence[np.ndarray],
    frames_v: Sequence[np.ndarray],
    times_s: Sequence[float],
    omega0: Optional[float] = None,
    amp_eps_rel: float = 1e-4,
    overlap_eps_rel: float = 1e-6,
    alpha_gamma: float = 1.0,
    alpha_scale: float = 0.95,
    return_psi: bool = False,
) -> BerryTimeSeries:
    """Compute a per-point U(1) Berry-like phase time series.

    Parameters
    ----------
    frames_u, frames_v:
        Lists of frames (T elements), each frame is (N,E). Units should be
        consistent with times_s (typically meters and m/s, times in seconds).
    times_s:
        Physical times of each frame (T elements).
    omega0:
        If None, uses :func:`compton_omega0`.
    amp_eps_rel:
        Relative amplitude threshold for "defined" polarization. Actual cutoff
        is amp_eps_rel * max(amplitude).
    overlap_eps_rel:
        Relative overlap magnitude threshold. If |<ψ_k|ψ_{k+1}>| is below this
        relative cutoff (w.r.t. 1), alpha is forced to 0 to avoid noisy phase.
    alpha_gamma, alpha_scale:
        Map amplitude → alpha = alpha_scale * (amp/amp_max)^alpha_gamma, then
        masked by the overlap threshold.
    return_psi:
        If True, returns ψ(t) as part of the dataclass.

    Returns
    -------
    BerryTimeSeries
    """
    U = _as_numpy_frames(frames_u)
    V = _as_numpy_frames(frames_v)

    if len(U) != len(V):
        raise ValueError(f"frames_u and frames_v must have same length. Got {len(U)} vs {len(V)}.")

    T = len(U)
    if T < 2:
        raise ValueError("Need at least 2 frames to compute Berry phase increments.")

    times = np.asarray(times_s, dtype=np.float64)
    if times.shape != (T,):
        raise ValueError(f"times_s must have shape ({T},). Got {times.shape}.")

    # Ensure shape consistency across frames
    N, E = U[0].shape
    for k in range(T):
        if U[k].shape != (N, E):
            raise ValueError(f"frames_u[{k}] shape {U[k].shape} != {(N, E)}")
        if V[k].shape != (N, E):
            raise ValueError(f"frames_v[{k}] shape {V[k].shape} != {(N, E)}")

    w0 = float(omega0) if omega0 is not None else float(compton_omega0())

    # Build complex amplitudes
    A = np.empty((T, N, E), dtype=np.complex128)
    amp = np.empty((T, N), dtype=np.float64)
    for k in range(T):
        A[k] = complex_amplitude_from_u_v(U[k], V[k], w0)
        amp[k] = _norm_rows(A[k])

    amp_max = float(np.max(amp)) if np.max(amp) > 0 else 1.0
    amp_eps = float(amp_eps_rel) * amp_max

    # Normalize to rays
    psi: Optional[np.ndarray] = None
    if return_psi:
        psi = np.empty((T, N, E), dtype=np.complex128)

    # We normalize for overlap computation either way; store in local variable
    PSI = np.empty((T, N, E), dtype=np.complex128)
    for k in range(T):
        denom = np.maximum(amp[k], amp_eps)
        PSI[k] = A[k] / denom[:, None]
        if return_psi:
            psi[k] = PSI[k]

    # Overlaps and increments
    delta = np.zeros((T, N), dtype=np.float64)
    conn = np.zeros((T, N), dtype=np.float64)
    overlap_abs = np.zeros((T, N), dtype=np.float64)

    for k in range(T - 1):
        dt = float(times[k + 1] - times[k])
        if dt <= 0:
            raise ValueError(f"times_s must be strictly increasing. Got dt={dt} at k={k}.")

        ov = np.sum(np.conj(PSI[k]) * PSI[k + 1], axis=-1)  # (N,)
        overlap_abs[k] = np.abs(ov)
        delta_k = np.angle(ov)
        delta[k] = delta_k
        conn[k] = delta_k / dt

    # Repeat last for display-friendly per-frame arrays
    overlap_abs[T - 1] = overlap_abs[T - 2]
    conn[T - 1] = conn[T - 2]
    delta[T - 1] = 0.0

    # Cumulative phase
    phase = np.cumsum(delta, axis=0)

    # Recommended alpha: amplitude-based, plus overlap sanity
    # (phase only meaningful where amplitude is non-negligible and overlap is non-degenerate)
    amp_norm = np.clip(amp / max(amp_max, 1e-30), 0.0, 1.0)
    alpha = alpha_scale * (amp_norm ** float(alpha_gamma))

    # Mask regions with weak amplitude
    alpha = np.where(amp >= amp_eps, alpha, 0.0)

    # Mask regions where overlap is too small (avoid random phase from orthogonality)
    overlap_cut = float(overlap_eps_rel)
    alpha = np.where(overlap_abs >= overlap_cut, alpha, 0.0)

    # Clamp
    alpha = np.clip(alpha, 0.0, 1.0)

    return BerryTimeSeries(
        phase=phase,
        connection=conn,
        delta_phase=delta,
        amplitude=amp,
        alpha=alpha,
        overlap_abs=overlap_abs,
        omega0=w0,
        times_s=times,
        psi=psi,
    )
