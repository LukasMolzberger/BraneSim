"""
Complex band state construction from real quadratures.

For a real wave field with position q and velocity v, constructs a complex
"positive-frequency" band state: ψ = q + i·v/ω

This is the standard analytic signal construction used in quantum mechanics
and wave optics to extract phase information from real fields.
"""

from __future__ import annotations
import torch


def _complex_dtype(real_dtype: torch.dtype) -> torch.dtype:
    """Map real dtype to corresponding complex dtype."""
    if real_dtype == torch.float64:
        return torch.complex128
    return torch.complex64  # float32 -> complex64, fallback


def complex_band_state_from_quadrature(
    q: torch.Tensor,
    v: torch.Tensor,
    omega: float,
    eps: float = 1e-12,
) -> torch.Tensor:
    """
    Build a complex 'positive-frequency' band state from real quadratures.

    Uses the standard analytic signal construction:
        ψ = q + i·v/ω

    This extracts the complex amplitude of a wave with carrier frequency ω,
    allowing computation of Berry phases and other phase-dependent quantities.

    Parameters
    ----------
    q : torch.Tensor
        Position/displacement field, shape [N] or [N, C], real-valued
    v : torch.Tensor
        Velocity field, shape [N] or [N, C], real-valued
    omega : float
        Carrier angular frequency (rad/s in physical units, or rad/sim-time in sim units)
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    torch.Tensor
        Complex band state, same shape as input, on same device

    Examples
    --------
    >>> # 1D wave with known frequency
    >>> xi = state.positions[:, 3]  # Amplitude field (ξ)
    >>> xi_dot = state.velocities[:, 3]
    >>> omega_sim = 2*pi * c_wave / wavelength
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    """
    cdtype = _complex_dtype(q.dtype)
    omega_t = torch.tensor(float(omega), device=q.device, dtype=q.dtype)
    denom = omega_t + eps
    return q.to(cdtype) + 1j * (v / denom).to(cdtype)


def pointwise_normalize(
    psi: torch.Tensor,
    eps: float = 1e-12,
):
    """
    Normalize complex state pointwise to unit magnitude.

    For Berry phase calculations, we need the normalized state |ψ̂⟩ = |ψ⟩/|ψ|
    at each point. This function computes both the normalized state and the
    amplitude (for masking/thresholding).

    Parameters
    ----------
    psi : torch.Tensor
        Complex state, shape [N] or [N, C]
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (psi_hat, amp) where:
        - psi_hat: normalized state, same shape as input
        - amp: pointwise amplitude, shape [N]

    Examples
    --------
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> # Now psi_hat has |psi_hat[i]| = 1 at each point i (where amp > eps)
    """
    if psi.ndim == 1:
        amp = torch.abs(psi)
        psi_hat = psi / (amp + eps)
        return psi_hat, amp

    # Multi-component case [N, C]
    amp = torch.linalg.norm(psi, dim=-1)  # [N]
    psi_hat = psi / (amp.unsqueeze(-1) + eps)
    return psi_hat, amp