from dataclasses import dataclass
from typing import Tuple, Optional

import numpy as np


@dataclass
class PhotonModeParameters:
    """
    Parameters for a circularly polarized photon mode in tubular coordinates.

    The mode lives on a closed curve C(z) with local orthonormal frame (t, n, b).
    The electric field E is circularly polarized in the (n, b) plane.

    Along the path we introduce a longitudinal wave amplitude that varies like
        A_long(z) ~ cos(phi_wave(z)),
    while the transverse profile is a Gaussian in the (n, b) directions.

    The internal phase phi_internal(z) runs from 0 to total_phase (default 4π)
    once around the closed loop. This drives the polarization orientation.
    A higher-frequency phase phi_wave(z) gives multiple crests along the loop.
    """
    peak_amplitude: float = 1.0

    # Gaussian widths in local transverse coordinates
    sigma_n: float = 0.12
    sigma_b: float = 0.12

    # How far we sample the Gaussian envelope (in units of σ)
    extent_sigma: float = 3.5

    # Polar sampling resolution of the cross section
    num_radial_samples: int = 8
    num_angular_samples: int = 32

    # Total internal phase advance along the loop (for polarization orientation)
    total_phase: float = 4.0 * np.pi

    # Optional global phase offset
    phase_offset: float = 0.0

    # Number of longitudinal wave cycles along one closed loop
    wave_cycles: int = 8

    # Relative magnitude of B compared to E (for visualization scaling)
    B_over_E: float = 1.0


def _normalize(v: np.ndarray, axis: int = -1, eps: float = 1e-9) -> np.ndarray:
    """Normalize vectors along a given axis."""
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v, axis=axis, keepdims=True)
    norm = np.maximum(norm, eps)
    return v / norm


def compute_circular_polarization_EB(
    tangents: np.ndarray,
    normals: np.ndarray,
    binormals: np.ndarray,
    params: PhotonModeParameters,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute circularly polarized electric and magnetic field vectors
    along the closed centerline C(z).

    Internal polarization phase (spinor-like):
        phi_internal(z) ∈ [0, total_phase),
        e_pol(z) = cos(phi_internal) * n(z) + sin(phi_internal) * b(z).

    Longitudinal wave amplitude:
        phi_wave(z) = 0.5 * wave_cycles * phi_internal(z),
        A_long(z) = cos(phi_wave(z)),

    so we get 'wave_cycles' full cos-oscillations along the loop, while
    ensuring A_long(z) matches at the ends because phi_internal runs
    from 0 to 4π and cos is 2π-periodic.

    Electric field:
        E(z) = peak_amplitude * A_long(z) * e_pol(z).

    Magnetic field:
        B(z) ∝ t(z) × E_hat(z),

    so that E, B and t are mutually orthogonal.

    Parameters
    ----------
    tangents, normals, binormals : (N, 3) np.ndarray
        Orthonormal frame vectors (t, n, b) along C(z).
    params : PhotonModeParameters
        Photon mode parameters.

    Returns
    -------
    E : (N, 3) np.ndarray
        Electric field vectors at each centerline point.
    B : (N, 3) np.ndarray
        Magnetic field vectors at each centerline point.
    phase_internal : (N,) np.ndarray
        Internal polarization phase φ_internal(z) ∈ [0, total_phase).
    A_long : (N,) np.ndarray
        Longitudinal wave amplitude A_long(z).
    """
    t = _normalize(tangents, axis=-1)
    n = _normalize(normals, axis=-1)
    b = _normalize(binormals, axis=-1)

    assert t.shape == n.shape == b.shape
    assert t.shape[1] == 3

    N = t.shape[0]

    # Internal phase along the loop: [0, total_phase)
    phase_internal = np.linspace(
        0.0,
        params.total_phase,
        N,
        endpoint=False,
        dtype=float,
    )
    phase_internal = phase_internal + params.phase_offset

    # Polarization orientation in (n, b) plane (driven by 4π phase)
    cos_phi = np.cos(phase_internal)
    sin_phi = np.sin(phase_internal)
    pol_dir = n * cos_phi[:, np.newaxis] + b * sin_phi[:, np.newaxis]

    # Longitudinal amplitude: higher-frequency wave along the loop
    phase_wave = 0.5 * params.wave_cycles * phase_internal
    A_long = np.cos(phase_wave)  # 'wave_cycles' cos periods along the loop

    E = params.peak_amplitude * A_long[:, np.newaxis] * pol_dir

    # Magnetic field orthogonal to both t and E
    E_hat = _normalize(E, axis=-1)
    B_dir = _normalize(np.cross(t, E_hat), axis=-1)
    B = params.B_over_E * params.peak_amplitude * A_long[:, np.newaxis] * B_dir

    return E, B, phase_internal, A_long


def sample_gaussian_envelope(
    centerline: np.ndarray,
    normals: np.ndarray,
    binormals: np.ndarray,
    params: PhotonModeParameters,
    longitudinal_modulation: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sample a Gaussian field envelope in the transverse (n, b) directions
    around the centerline, optionally modulated along the path.

    Transverse profile:
        A_trans(r) = exp(-0.5 * r^2),
    where r is an elliptical radial coordinate in (n, b) units of σ.

    Longitudinal modulation:
        A_long(z) = longitudinal_modulation[z]  (if provided),
        otherwise A_long(z) = 1.

    Total amplitude:
        A_total(z, r) = peak_amplitude * A_long(z) * A_trans(r).
    """
    pts = np.asarray(centerline, dtype=float)
    n = _normalize(normals, axis=-1)
    b = _normalize(binormals, axis=-1)

    assert pts.shape == n.shape == b.shape
    assert pts.shape[1] == 3

    N = pts.shape[0]
    Nr = params.num_radial_samples
    Ntheta = params.num_angular_samples

    # Longitudinal modulation A_long(z)
    if longitudinal_modulation is None:
        A_long = np.ones(N, dtype=float)
    else:
        A_long = np.asarray(longitudinal_modulation, dtype=float)
        if A_long.shape != (N,):
            raise ValueError(
                f"longitudinal_modulation must have shape ({N},), got {A_long.shape}"
            )

    # Dimensionless radial coordinate r (in units of σ)
    r_max = params.extent_sigma
    r = np.linspace(0.0, r_max, Nr, dtype=float)
    theta = np.linspace(0.0, 2.0 * np.pi, Ntheta, endpoint=False, dtype=float)
    R, Theta = np.meshgrid(r, theta, indexing="ij")  # (Nr, Ntheta)

    # Elliptical polar coordinates in (n, b) plane
    Xn = R * np.cos(Theta) * params.sigma_n  # (Nr, Ntheta)
    Xb = R * np.sin(Theta) * params.sigma_b  # (Nr, Ntheta)

    # Transverse Gaussian: in these coordinates, r^2 = (Xn/σ_n)^2 + (Xb/σ_b)^2
    A_trans = np.exp(-0.5 * R ** 2)  # (Nr, Ntheta)

    # Combine longitudinal and transverse amplitude
    A_long_expanded = A_long[:, np.newaxis, np.newaxis]  # (N, 1, 1)
    amplitudes = (
        params.peak_amplitude * A_long_expanded * A_trans[np.newaxis, :, :]
    )  # (N, Nr, Nθ)

    # Broadcast centerline + frame + transverse offsets to 3D positions
    pts_expanded = pts[:, np.newaxis, np.newaxis, :]    # (N, 1, 1, 3)
    n_expanded = n[:, np.newaxis, np.newaxis, :]        # (N, 1, 1, 3)
    b_expanded = b[:, np.newaxis, np.newaxis, :]        # (N, 1, 1, 3)

    Xn3 = Xn[np.newaxis, :, :, np.newaxis]              # (1, Nr, Nθ, 1)
    Xb3 = Xb[np.newaxis, :, :, np.newaxis]              # (1, Nr, Nθ, 1)

    field_points = pts_expanded + Xn3 * n_expanded + Xb3 * b_expanded  # (N, Nr, Nθ, 3)

    # Flatten for plotting
    field_points_flat = field_points.reshape(-1, 3)
    amplitudes_flat = amplitudes.reshape(-1)

    return field_points_flat, amplitudes_flat