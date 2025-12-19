"""
ω-free analytic signal construction via positive-frequency projection.

This module provides tools to construct complex analytic signals from real fields
without requiring velocity or an explicit carrier frequency ω. The analytic signal
is obtained by taking the FFT along a chosen spatial axis and zeroing negative
frequencies (positive-frequency projection).

This is the correct approach for extracting phase information from a real field
when we want to diagnose Berry-like phase transport without assuming a single
carrier frequency or filtering to a specific band.

Works for scalar and vector fields in any dimension (1D/2D/3D).
"""

from __future__ import annotations
import torch

from .types import GridSpec


def _complex_dtype(real_dtype: torch.dtype) -> torch.dtype:
    """Map real dtype to corresponding complex dtype."""
    if real_dtype == torch.float64:
        return torch.complex128
    return torch.complex64


def _analytic_filter(N: int, device, dtype) -> torch.Tensor:
    """
    Frequency-domain multiplier h[k] such that ifft(fft(x)*h) is the analytic signal.

    For a real signal x, the analytic signal is obtained by:
    1. Take FFT: X[k] = FFT(x)
    2. Zero negative frequencies: X_a[k] = h[k] * X[k]
    3. Take IFFT: x_a = IFFT(X_a)

    The filter h[k] is:
    - h[0] = 1 (DC component, keep as-is)
    - h[k] = 2 for positive frequencies (double to preserve energy)
    - h[k] = 0 for negative frequencies (zero out)
    - h[N/2] = 1 for Nyquist frequency (even N only)

    Parameters
    ----------
    N : int
        Length of signal
    device : torch.device
        Device for output
    dtype : torch.dtype
        Data type for output (should be real dtype)

    Returns
    -------
    torch.Tensor
        Filter array of length N
    """
    h = torch.zeros(N, device=device, dtype=dtype)
    if N % 2 == 0:
        # Even length: [DC, pos freqs, Nyquist, neg freqs]
        h[0] = 1                    # DC
        h[N // 2] = 1               # Nyquist
        h[1:N // 2] = 2             # Positive frequencies
        # h[N//2+1:] = 0 stays zero (negative frequencies)
    else:
        # Odd length: [DC, pos freqs, neg freqs]
        h[0] = 1                    # DC
        h[1:(N + 1) // 2] = 2       # Positive frequencies
        # h[(N+1)//2:] = 0 stays zero (negative frequencies)
    return h


def analytic_signal_along_axis(q: torch.Tensor, axis: int, spatial_ndim: int) -> torch.Tensor:
    """
    ω-free analytic signal along a chosen spatial axis via positive-frequency projection.

    Constructs a complex analytic signal from a real field by:
    1. Taking FFT along the specified axis
    2. Zeroing negative frequencies (keeping DC and positive frequencies)
    3. Taking IFFT to get complex field

    This gives a complex field whose instantaneous phase can be used for Berry
    diagnostics without requiring velocity or a carrier frequency ω.

    **Important:** This approach does NOT isolate a single eigen-branch. If the
    field contains multiple superposed modes or degeneracies, the phase will mix
    contributions from all modes. For clean Berry diagnostics, the input field
    should ideally be a narrow-band wavepacket dominated by a single mode.

    Works for:
    - Scalar fields: q.shape == [*grid_shape]
    - Vector fields: q.shape == [*grid_shape, C] (components in last dimension)

    where spatial_ndim == len(grid_shape).

    Parameters
    ----------
    q : torch.Tensor
        Real-valued field
        Shape: [*grid_shape] or [*grid_shape, C]
    axis : int
        Spatial axis along which to compute analytic signal (0 to spatial_ndim-1)
    spatial_ndim : int
        Number of spatial dimensions (e.g., 1 for 1D, 2 for 2D, 3 for 3D)

    Returns
    -------
    torch.Tensor
        Complex analytic signal, same shape as input

    Examples
    --------
    >>> # 1D scalar field: cosine wave
    >>> import torch
    >>> x = torch.linspace(0, 10*torch.pi, 256)
    >>> q = torch.cos(x)
    >>> psi = analytic_signal_along_axis(q, axis=0, spatial_ndim=1)
    >>> # psi should be ≈ exp(i*x) (up to boundaries)
    >>> phase = torch.angle(psi)
    >>> # phase should be ≈ x (linear)
    >>>
    >>> # 2D scalar field: wave along x-axis
    >>> q = torch.cos(x[:, None]).expand(256, 64)
    >>> psi = analytic_signal_along_axis(q, axis=0, spatial_ndim=2)
    >>> psi.shape
    torch.Size([256, 64])
    >>>
    >>> # 2D vector field (e.g., polarization)
    >>> q = torch.randn(128, 128, 3)
    >>> psi = analytic_signal_along_axis(q, axis=0, spatial_ndim=2)
    >>> psi.shape
    torch.Size([128, 128, 3])
    """
    if axis < 0 or axis >= spatial_ndim:
        raise ValueError(f"axis must be in [0,{spatial_ndim-1}], got {axis}")
    if q.ndim not in (spatial_ndim, spatial_ndim + 1):
        raise ValueError(
            f"q.ndim={q.ndim} not compatible with spatial_ndim={spatial_ndim}. "
            f"Expected q.shape = [*grid_shape] or [*grid_shape, C]"
        )

    has_components = (q.ndim == spatial_ndim + 1)
    cdtype = _complex_dtype(q.dtype)

    # Strategy: move chosen axis to last spatial position for processing
    # Then unfold components if present, process all at once, fold back

    # Move chosen axis to position (spatial_ndim - 1)
    x = q.movedim(axis, spatial_ndim - 1)
    # Now x.shape = [*transverse, N] or [*transverse, N, C]

    if has_components:
        # x.shape = [*transverse, N, C]
        # Move components before N: [*transverse, C, N]
        x = x.movedim(-1, -2)

    # Now x.shape = [*transverse, N] or [*transverse, C, N]
    N = x.shape[-1]
    batch = int(x.numel() // N)
    x2d = x.reshape(batch, N)  # [batch, N]

    # Compute FFT along N
    X = torch.fft.fft(x2d.to(cdtype), dim=-1)  # [batch, N]

    # Apply analytic filter
    h = _analytic_filter(N, device=q.device, dtype=X.real.dtype)
    X = X * h[None, :]  # [batch, N]

    # Inverse FFT
    psi2d = torch.fft.ifft(X, dim=-1)  # [batch, N], complex

    # Reshape back
    psi = psi2d.reshape(*x.shape)  # [*transverse, N] or [*transverse, C, N]

    if has_components:
        # [*transverse, C, N] -> [*transverse, N, C]
        psi = psi.movedim(-2, -1)

    # Move axis back to original position
    psi = psi.movedim(spatial_ndim - 1, axis)

    return psi


def pointwise_normalize_from_grid(
    psi: torch.Tensor,
    grid: GridSpec,
    eps: float = 1e-12,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Normalize complex state pointwise to unit magnitude.

    Explicitly handles:
    - Scalar fields: psi.ndim == grid.D
      Uses pointwise absolute value
    - Vector fields: psi.ndim == grid.D + 1 (components in last dimension)
      Uses norm over component dimension

    For Berry phase calculations, we need the normalized state |ψ̂⟩ = |ψ⟩/|ψ|
    at each point. This function computes both the normalized state and the
    amplitude (for masking/thresholding).

    Parameters
    ----------
    psi : torch.Tensor
        Complex state
        Shape: [*grid.shape] or [*grid.shape, C]
    grid : GridSpec
        Grid specification (used to infer dimensionality)
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (psi_hat, amp) where:
        - psi_hat: normalized state, same shape as input
        - amp: pointwise amplitude, shape [*grid.shape]

    Examples
    --------
    >>> import torch
    >>> from branesim.diagnostics.types import GridSpec
    >>>
    >>> # 1D scalar field
    >>> grid = GridSpec(shape=(100,), spacing_sim=1.0)
    >>> psi = torch.randn(100, dtype=torch.complex64)
    >>> psi_hat, amp = pointwise_normalize_from_grid(psi, grid)
    >>> psi_hat.shape, amp.shape
    (torch.Size([100]), torch.Size([100]))
    >>> torch.allclose(torch.abs(psi_hat), torch.ones_like(amp), atol=1e-6)
    True
    >>>
    >>> # 2D vector field
    >>> grid = GridSpec(shape=(64, 64), spacing_sim=1.0)
    >>> psi = torch.randn(64, 64, 3, dtype=torch.complex64)
    >>> psi_hat, amp = pointwise_normalize_from_grid(psi, grid)
    >>> psi_hat.shape, amp.shape
    (torch.Size([64, 64, 3]), torch.Size([64, 64]))
    """
    if psi.ndim == grid.D:
        # Scalar field: pointwise absolute value
        amp = torch.abs(psi)
        psi_hat = psi / (amp + eps)
        return psi_hat, amp
    elif psi.ndim == grid.D + 1:
        # Vector field: norm over last dimension
        amp = torch.linalg.norm(psi, dim=-1)  # [*grid.shape]
        psi_hat = psi / (amp.unsqueeze(-1) + eps)  # [*grid.shape, C]
        return psi_hat, amp
    else:
        raise ValueError(
            f"psi.ndim={psi.ndim} incompatible with grid.D={grid.D}. "
            f"Expected psi.shape = {grid.shape} or {grid.shape} + (C,)"
        )


def pointwise_normalize_scalar(
    psi: torch.Tensor,
    eps: float = 1e-12,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Normalize complex scalar state pointwise to unit magnitude.

    Helper function for scalar fields where we don't need grid info.

    Parameters
    ----------
    psi : torch.Tensor
        Complex scalar state, any shape
    eps : float, optional
        Small constant to prevent division by zero

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (psi_hat, amp) where both have same shape as input
    """
    amp = torch.abs(psi)
    psi_hat = psi / (amp + eps)
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