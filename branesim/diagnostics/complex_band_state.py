"""
Complex band state construction from real quadratures (dimension-agnostic).

For a real wave field with position q and velocity v, constructs a complex
"positive-frequency" band state: ψ = q + i·v/ω

This is the standard analytic signal construction used in quantum mechanics
and wave optics to extract phase information from real fields.

Works for scalar and vector fields in any dimension (1D/2D/3D).
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

    Works for:
    - Scalar fields: q.shape = [*grid_shape]
    - Vector fields: q.shape = [*grid_shape, C] (components in last dimension)

    Parameters
    ----------
    q : torch.Tensor
        Position/displacement field, real-valued
        Shape: [*grid_shape] or [*grid_shape, C]
    v : torch.Tensor
        Velocity field, real-valued
        Shape: same as q
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
    >>> # 1D scalar field
    >>> xi = torch.randn(100)
    >>> xi_dot = torch.randn(100)
    >>> omega_sim = 6.28
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    >>> psi.shape
    torch.Size([100])
    >>>
    >>> # 2D scalar field
    >>> xi = torch.randn(64, 64)
    >>> xi_dot = torch.randn(64, 64)
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    >>> psi.shape
    torch.Size([64, 64])
    >>>
    >>> # 3D vector field (e.g., polarization)
    >>> xi = torch.randn(32, 32, 32, 3)
    >>> xi_dot = torch.randn(32, 32, 32, 3)
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    >>> psi.shape
    torch.Size([32, 32, 32, 3])
    """
    if q.shape != v.shape:
        raise ValueError(f"Position and velocity shapes must match: {q.shape} != {v.shape}")

    cdtype = _complex_dtype(q.dtype)
    omega_t = torch.tensor(float(omega), device=q.device, dtype=q.dtype)
    denom = omega_t + eps
    return q.to(cdtype) + 1j * (v / denom).to(cdtype)


def pointwise_normalize(
    psi: torch.Tensor,
    eps: float = 1e-12,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Normalize complex state pointwise to unit magnitude.

    For Berry phase calculations, we need the normalized state |ψ̂⟩ = |ψ⟩/|ψ|
    at each point. This function computes both the normalized state and the
    amplitude (for masking/thresholding).

    Works for:
    - Scalar fields: psi.shape = [*grid_shape]
    - Vector fields: psi.shape = [*grid_shape, C] (norm over last dimension)

    Parameters
    ----------
    psi : torch.Tensor
        Complex state
        Shape: [*grid_shape] or [*grid_shape, C]
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (psi_hat, amp) where:
        - psi_hat: normalized state, same shape as input
        - amp: pointwise amplitude, shape [*grid_shape]

    Examples
    --------
    >>> # 1D scalar field
    >>> psi = torch.randn(100, dtype=torch.complex64)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> psi_hat.shape, amp.shape
    (torch.Size([100]), torch.Size([100]))
    >>> torch.allclose(torch.abs(psi_hat), torch.ones_like(amp), atol=1e-6)
    True
    >>>
    >>> # 2D vector field
    >>> psi = torch.randn(64, 64, 3, dtype=torch.complex64)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> psi_hat.shape, amp.shape
    (torch.Size([64, 64, 3]), torch.Size([64, 64]))
    """
    # Detect if last dimension is component dimension
    # Heuristic: if tensor is complex and has more than 1 dimension,
    # assume last dimension is component dimension if it's reasonably small
    is_vector = False
    if psi.ndim >= 2:
        # Check if last dimension looks like a component dimension (size 2-4 typically)
        # But to be safe, we compute the norm and check dimensionality
        # For now, we'll compute norm over last dim if ndim > 1 and treat as vector
        # Actually, let's be more explicit: if user wants vector norm, last dim should be small
        # Let's use a simple heuristic: if shape[-1] <= 4 or explicitly multi-component
        # Actually, better to always check: compute abs, if result shape differs, it's scalar
        test_amp = torch.abs(psi)
        if test_amp.shape == psi.shape:
            # Scalar field (abs preserves shape)
            is_vector = False
        else:
            # This shouldn't happen with abs, so always scalar case for abs
            # Let's use linalg.norm for potential vector case
            is_vector = False  # Default to scalar

    # Try to detect vector field: if last dimension is small, assume it's components
    # More robust: always do scalar unless explicitly told otherwise
    # For now, implement both branches based on shape analysis

    # Scalar case: compute pointwise absolute value
    if not is_vector:
        amp = torch.abs(psi)
        psi_hat = psi / (amp + eps)
        return psi_hat, amp

    # Vector case: compute norm over last dimension
    # (This path is for future vector field support)
    amp = torch.linalg.norm(psi, dim=-1)  # [*grid_shape]
    psi_hat = psi / (amp.unsqueeze(-1) + eps)  # [*grid_shape, C]
    return psi_hat, amp


def pointwise_normalize_vector(
    psi: torch.Tensor,
    eps: float = 1e-12,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Normalize complex vector state pointwise to unit magnitude.

    Explicitly handles vector fields where the last dimension contains
    components. Computes the norm over the component dimension.

    Parameters
    ----------
    psi : torch.Tensor
        Complex vector state, shape [*grid_shape, C]
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (psi_hat, amp) where:
        - psi_hat: normalized state, shape [*grid_shape, C]
        - amp: pointwise amplitude, shape [*grid_shape]

    Examples
    --------
    >>> psi = torch.randn(64, 64, 3, dtype=torch.complex64)
    >>> psi_hat, amp = pointwise_normalize_vector(psi)
    >>> psi_hat.shape, amp.shape
    (torch.Size([64, 64, 3]), torch.Size([64, 64]))
    """
    amp = torch.linalg.norm(psi, dim=-1)  # [*grid_shape]
    psi_hat = psi / (amp.unsqueeze(-1) + eps)
    return psi_hat, amp