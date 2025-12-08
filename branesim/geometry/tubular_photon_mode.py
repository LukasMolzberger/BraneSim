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
        A_long(z) ~ cos(phi(z)),
    while the transverse profile is a Gaussian in the (n, b) directions.

    The internal phase phi(z) runs from 0 to total_phase (default 4π) once
    around the closed loop. This guarantees that both amplitude and orientation
    match continuously at the "ends" of the strip.
    """
    peak_amplitude: float = 1.0

    # Gaussian widths in local transverse coordinates
    sigma_n: float = 0.08
    sigma_b: float = 0.08

    # How far we sample the Gaussian envelope (in units of σ)
    extent_sigma: float = 3.0

    # Polar sampling resolution of the cross section
    num_radial_samples: int = 6
    num_angular_samples: int = 24

    # Total internal phase advance along the loop
    total_phase: float = 4.0 * np.pi

    # Optional global phase offset
    phase_offset: float = 0.0

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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute circularly polarized electric and magnetic field vectors
    along the closed centerline C(z).

    We use a single internal phase φ(z) ∈ [0, total_phase) and define

        e_pol(z) = cos φ(z) * n(z) + sin φ(z) * b(z)   (unit polarization dir)
        A_long(z) = cos φ(z)                           (longitudinal amplitude)

    so that the *orientation* rotates in the (n, b) plane and the *magnitude*
    of the wave oscillates along the path.

    The electric field is
        E(z) = peak_amplitude * A_long(z) * e_pol(z),

    and the magnetic field is constructed as
        B(z) ∝ t(z) × E_hat(z),

    which ensures E, B and t are mutually orthogonal.

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
    phase : (N,) np.ndarray
        Internal phase values φ(z) ∈ [0, total_phase).
        The longitudinal amplitude is A_long(z) = cos(phase[z]).
    """
    t = _normalize(tangents, axis=-1)
    n = _normalize(normals, axis=-1)
    b = _normalize(binormals, axis=-1)

    assert t.shape == n.shape == b.shape
    assert t.shape[1] == 3

    N = t.shape[0]

    # Internal phase along the loop: [0, total_phase)
    phase = np.linspace(
        0.0,
        params.total_phase,
        N,
        endpoint=False,
        dtype=float,
    )
    phase = phase + params.phase_offset

    cos_phi = np.cos(phase)
    sin_phi = np.sin(phase)

    # Polarization direction in (n, b) plane
    pol_dir = n * cos_phi[:, np.newaxis] + b * sin_phi[:, np.newaxis]

    # Longitudinal amplitude: actual wave crests/troughs along the path
    A_long = np.cos(phase)  # same phase → 2 full wavelengths for total_phase = 4π

    E = params.peak_amplitude * A_long[:, np.newaxis] * pol_dir

    # Magnetic field orthogonal to both t and E
    E_hat = _normalize(E, axis=-1)
    B_dir = _normalize(np.cross(t, E_hat), axis=-1)
    B = params.B_over_E * params.peak_amplitude * A_long[:, np.newaxis] * B_dir

    return E, B, phase


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

    The total amplitude is then
        A_total(z, r) = peak_amplitude * A_long(z) * A_trans(r).

    Parameters
    ----------
    centerline : (N, 3) np.ndarray
        Points C(z_i) along the closed curve in R^3.
    normals, binormals : (N, 3) np.ndarray
        Local orthonormal frame vectors at each centerline point.
    params : PhotonModeParameters
        Controls peak amplitude, Gaussian widths and sampling resolution.
    longitudinal_modulation : (N,) np.ndarray or None
        Optional longitudinal amplitude A_long(z). If None, a constant
        factor A_long(z) = 1 is used for all z.

    Returns
    -------
    field_points : (N * Nr * Nθ, 3) np.ndarray
        3D sample points in space.
    amplitudes : (N * Nr * Nθ,) np.ndarray
        Field amplitude A_total at each sample point.
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