"""
Envelope and phase functions for carrier compilation.

This module provides functions for building amplitude envelopes and
carrier phases, which are then combined with polarization bases to
produce complex carrier fields.
"""

import torch
import numpy as np


def gaussian_envelope(
    coords: torch.Tensor,      # [N, d]
    center: torch.Tensor,      # [d]
    sigma: float,
    amplitude: float,
) -> torch.Tensor:
    """
    Compute Gaussian envelope centered at a point.

    Returns:
        envelope: [N] amplitude values A(x) = amplitude * exp(-r²/(2σ²))
    """
    # Compute squared distance from center
    r_sq = torch.sum((coords - center) ** 2, dim=1)

    # Gaussian envelope
    envelope = amplitude * torch.exp(-r_sq / (2 * sigma ** 2))

    return envelope


def plane_wave_phase(
    coords: torch.Tensor,      # [N, d]
    k_vector: torch.Tensor,    # [d]
) -> torch.Tensor:
    """
    Compute plane wave phase φ(x) = k · x.

    Returns:
        phase: [N] phase values
    """
    phase = torch.matmul(coords, k_vector)
    return phase


def tubular_envelope_electron(
    coords: torch.Tensor,          # [N, 3] lattice coordinates
    centerline_points: torch.Tensor,  # [M, 3] centerline samples
    tube_sigma: float,             # radial width
    twist_winding: int = 1,
    device=None,
    dtype=None,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Compute tubular envelope for electron initialization.

    Maps each lattice point to tubular coordinates (s, ρ, θ) relative
    to the centerline, then computes:
    - Radial Gaussian envelope
    - Double-ridge angular profile with twist
    - Arclength coordinate along centerline
    - Half-angle twist for spinorial transport

    Args:
        coords: [N, 3] lattice point coordinates
        centerline_points: [M, 3] sampled centerline positions
        tube_sigma: Gaussian width for radial falloff
        twist_winding: Number of times ridges rotate per loop (ℓ)
        device: torch device
        dtype: torch dtype

    Returns:
        envelope: [N] combined amplitude envelope
        arclength: [N] arclength coordinate s along centerline
        alpha: [N] twist angle α(s)
        alpha_half: [N] half-angle α/2 for spinorial polarization
    """
    if device is None:
        device = coords.device
    if dtype is None:
        dtype = coords.dtype

    N = coords.shape[0]
    M = centerline_points.shape[0]

    # Compute distances from each lattice point to all centerline points
    # coords: [N, 1, 3], centerline: [1, M, 3]
    coords_exp = coords.unsqueeze(1)  # [N, 1, 3]
    centerline_exp = centerline_points.unsqueeze(0)  # [1, M, 3]
    dists = torch.linalg.norm(coords_exp - centerline_exp, dim=2)  # [N, M]

    # Find closest centerline point for each lattice point
    closest_indices = torch.argmin(dists, dim=1)  # [N]
    rho = dists[torch.arange(N), closest_indices]  # [N] radial distance

    # Compute arclength (approximate as proportional to index)
    # For a closed curve with M samples: s ≈ (2πR) * (index / M)
    # We'll normalize to [0, 2π] for convenience
    s_normalized = (closest_indices.float() / M) * 2 * np.pi

    # Radial envelope: Gaussian in ρ
    G_rho = torch.exp(-rho ** 2 / (2 * tube_sigma ** 2))

    # Compute angular coordinate θ around the tube
    # This requires computing the local Frenet frame, but for simplicity
    # we'll use a simpler approach: project onto radial direction
    closest_points = centerline_points[closest_indices]  # [N, 3]
    radial_vec = coords - closest_points  # [N, 3]
    radial_dist = torch.linalg.norm(radial_vec, dim=1, keepdim=True).clamp(min=1e-10)
    radial_unit = radial_vec / radial_dist  # [N, 3] normalized

    # Simple angular coordinate using atan2 on first two components
    # (this is approximate but sufficient for envelope purposes)
    theta = torch.atan2(radial_unit[:, 1], radial_unit[:, 0])  # [N]

    # Twist angle: α(s) = ℓ * s
    alpha = twist_winding * s_normalized  # [N]

    # Double-ridge profile: two peaks at θ = α and θ = α + π
    sigma_theta = 0.5  # angular width in radians
    ridge1 = torch.exp(-((theta - alpha) % (2*np.pi))**2 / (2 * sigma_theta**2))
    ridge2 = torch.exp(-((theta - (alpha + np.pi)) % (2*np.pi))**2 / (2 * sigma_theta**2))

    # Combined envelope
    envelope = G_rho * (ridge1 + ridge2)

    # Half-angle for spinorial transport
    alpha_half = alpha / 2.0

    return envelope, s_normalized, alpha, alpha_half