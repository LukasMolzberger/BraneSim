"""
Berry phase computation for 1D chains.

Computes the discrete Berry phase profile along a 1D lattice using the
gauge-invariant phase increment between neighboring states:

    Δφ_i = arg⟨u_i|u_{i+1}⟩

The cumulative Berry phase is then:

    γ_0 = 0,  γ_{i+1} = γ_i + Δφ_i

This gives one Berry phase value per lattice point.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
import torch


@dataclass(frozen=True)
class BerryPhase1DConfig:
    """
    Configuration for 1D Berry phase computation.

    Attributes
    ----------
    spacing : float
        Grid spacing in simulation units (usually h_sim = 1.0)
    amplitude_threshold : float
        Minimum amplitude to include a point in Berry phase calculation.
        Points with |ψ| < threshold are masked out to avoid noise.
    eps : float
        Small constant for numerical stability
    unwrap : bool
        If True, accumulate raw phase increments (may grow beyond [-π, π]).
        Currently always uses cumulative sum; "unwrapping" in the phase-unwrap
        sense would require additional logic.
    force_cpu_on_mps : bool
        If True and device is MPS, move tensors to CPU for complex operations.
        MPS complex support can be spotty, so this is recommended.
    """
    spacing: float
    amplitude_threshold: float = 1e-8
    eps: float = 1e-12
    unwrap: bool = True
    force_cpu_on_mps: bool = True


def berry_phase_profile_along_x(
    u_hat: torch.Tensor,
    amp: torch.Tensor,
    cfg: BerryPhase1DConfig
) -> dict[str, torch.Tensor]:
    """
    Compute Berry phase profile γ(x) along a 1D chain.

    Uses the discrete Berry connection:
        A_x[i] ≈ (1/h) · arg⟨u_i|u_{i+1}⟩

    and integrates to get the Berry phase:
        γ[i+1] = γ[i] + h·A_x[i] = γ[i] + arg⟨u_i|u_{i+1}⟩

    Parameters
    ----------
    u_hat : torch.Tensor
        Normalized complex state, shape [N] or [N, C]
        Should have |u_hat[i]| ≈ 1 at each point (from pointwise_normalize)
    amp : torch.Tensor
        Real amplitude before normalization, shape [N]
        Used for thresholding/masking low-amplitude regions
    cfg : BerryPhase1DConfig
        Configuration parameters

    Returns
    -------
    dict
        Dictionary containing:
        - gamma_unwrapped [N]: cumulative Berry phase (may grow beyond [-π, π])
        - gamma_wrapped [N]: wrapped to [-π, π] at each point
        - dphi [N-1]: phase increment between neighbors, Δφ_i = arg⟨u_i|u_{i+1}⟩
        - A_x [N-1]: Berry connection along x, A_x[i] = Δφ_i / h
        - mask [N]: boolean mask indicating which points have sufficient amplitude
        - valid_edge [N-1]: boolean mask for edges between valid points

    Notes
    -----
    - Invalid edges (where either endpoint is below threshold) have dphi = 0
      to avoid contaminating the cumulative phase with noise.
    - The Berry connection A_x has units [rad / sim-length].
    - For physical interpretation, convert sim-length to meters using the mapper.

    Examples
    --------
    >>> from branesim.analysis import (
    ...     complex_band_state_from_quadrature,
    ...     pointwise_normalize,
    ...     berry_phase_profile_along_x,
    ...     BerryPhase1DConfig
    ... )
    >>>
    >>> # Construct complex state from real fields
    >>> xi = state.positions[:, 3]
    >>> xi_dot = state.velocities[:, 3]
    >>> omega_sim = 2*pi * c_wave_sim / wavelength_sim
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>>
    >>> # Compute Berry phase profile
    >>> cfg = BerryPhase1DConfig(spacing=h_sim, amplitude_threshold=1e-6)
    >>> result = berry_phase_profile_along_x(psi_hat, amp, cfg)
    >>>
    >>> gamma = result['gamma_wrapped']  # Berry phase per point [N]
    >>> A_x = result['A_x']  # Berry connection [N-1]
    """
    # Handle MPS device limitations with complex numbers
    if cfg.force_cpu_on_mps and u_hat.device.type == "mps":
        u_hat = u_hat.to("cpu")
        amp = amp.to("cpu")

    # Mask low-amplitude regions
    mask = amp > cfg.amplitude_threshold
    valid_edge = mask[:-1] & mask[1:]  # Edge is valid if both endpoints are valid

    # Compute overlap ⟨u_i|u_{i+1}⟩
    if u_hat.ndim == 1:
        # Scalar case: simple conjugate multiply
        overlap = torch.conj(u_hat[:-1]) * u_hat[1:]
    else:
        # Vector case [N, C]: inner product along component dimension
        overlap = torch.sum(torch.conj(u_hat[:-1, :]) * u_hat[1:, :], dim=-1)

    # Set invalid edges to unity overlap (phase increment = 0)
    # This prevents noise from contaminating the cumulative phase
    overlap = torch.where(valid_edge, overlap, torch.ones_like(overlap))

    # Extract phase increment Δφ_i = arg⟨u_i|u_{i+1}⟩
    # torch.angle returns values in [-π, π]
    dphi = torch.angle(overlap)  # shape [N-1]

    # Cumulative Berry phase: γ[i+1] = γ[i] + Δφ_i
    gamma = torch.zeros(u_hat.shape[0], device=dphi.device, dtype=dphi.real.dtype)
    if cfg.unwrap:
        # Simple cumulative sum (may exceed [-π, π])
        gamma[1:] = torch.cumsum(dphi, dim=0)
    else:
        # Still use cumsum for now; true phase unwrapping would need more logic
        gamma[1:] = torch.cumsum(dphi, dim=0)

    # Wrapped version: fold back into [-π, π]
    two_pi = 2.0 * math.pi
    gamma_wrapped = (gamma + math.pi) % two_pi - math.pi

    # Berry connection: A_x[i] = Δφ_i / h
    # Units: [rad / sim-length]
    A_x = dphi / float(cfg.spacing)

    return {
        "gamma_unwrapped": gamma,
        "gamma_wrapped": gamma_wrapped,
        "dphi": dphi,
        "A_x": A_x,
        "mask": mask,
        "valid_edge": valid_edge,
    }