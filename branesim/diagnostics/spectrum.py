"""
Local spectrum diagnostics via STFT (dimension-agnostic).

This module provides local/windowed spectrum analysis using Short-Time Fourier
Transform (STFT) over spatial coordinates. Unlike global FFT, STFT gives
position-resolved spectral information: S(x, k) instead of just S(k).

**Design principle:**
- All functions are dimension-agnostic (work for 1D/2D/3D without branching)
- Use tensor operations (movedim, reshape, unfold) to handle arbitrary dimensions
- No global FFT functions (no spatial_power_spectrum_1d/nd, no radial averaging)

Works for scalar and vector fields in any dimension.
"""

from __future__ import annotations
from typing import Literal
import math
import torch

from .types import GridSpec


def _get_window(
    win_len: int,
    window: Literal["none", "hann", "hamming", "blackman"] | None,
    device,
    dtype,
) -> torch.Tensor | None:
    """
    Create window function of given length.

    Parameters
    ----------
    win_len : int
        Window length
    window : str | None
        Window type
    device : torch.device
        Device for output
    dtype : torch.dtype
        Data type for output

    Returns
    -------
    torch.Tensor | None
        Window array of length win_len, or None if window="none"
    """
    if window is None or window == "none":
        return None

    if window == "hann":
        return torch.hann_window(win_len, device=device, dtype=dtype)
    elif window == "hamming":
        return torch.hamming_window(win_len, device=device, dtype=dtype)
    elif window == "blackman":
        return torch.blackman_window(win_len, device=device, dtype=dtype)
    else:
        raise ValueError(f"Unknown window type: {window}")


def local_power_spectrum_along_axis(
    field: torch.Tensor,
    grid: GridSpec,
    axis: int,
    win_len: int,
    hop: int,
    window: Literal["none", "hann", "hamming", "blackman"] | None = "hann",
    transverse_reduction: Literal["mean", "none"] = "mean",
    component_reduction: Literal["sum", "none"] = "sum",
    normalize: Literal["none", "win_len"] = "win_len",
    return_complex: bool = False,
) -> dict[str, torch.Tensor]:
    """
    Compute local power spectrum along a spatial axis using STFT.

    This function is **dimension-agnostic**: it works for 1D, 2D, 3D (and beyond)
    fields by using generic tensor operations.

    The STFT is computed by:
    1. Sliding a window along the chosen axis
    2. Computing FFT within each window
    3. Computing power spectrum |FFT|^2
    4. Optionally averaging over transverse directions and/or components

    Works for:
    - Scalar fields: field.shape == grid.shape
    - Vector fields: field.shape == (*grid.shape, C) (components last)

    Parameters
    ----------
    field : torch.Tensor
        Real or complex field
        Shape: [*grid.shape] or [*grid.shape, C]
    grid : GridSpec
        Grid specification
    axis : int
        Spatial axis along which to compute STFT (0 to grid.D-1)
    win_len : int
        Length of window (in grid points)
    hop : int
        Hop size between windows (stride in grid points)
    window : str | None, optional
        Window function: "none", "hann", "hamming", "blackman"
        Default: "hann"
    transverse_reduction : {"mean", "none"}, optional
        How to reduce over transverse directions:
        - "mean": return mean and std over transverse coords
        - "none": return full array with transverse dims preserved
        Default: "mean"
    component_reduction : {"sum", "none"}, optional
        How to reduce over vector components (for vector fields):
        - "sum": sum power over components before other reductions
        - "none": keep components separate
        Default: "sum"
    normalize : {"none", "win_len"}, optional
        Power normalization:
        - "none": raw |FFT|^2
        - "win_len": divide by win_len (energy conservation)
        Default: "win_len"

    Returns
    -------
    dict[str, torch.Tensor]
        Dictionary containing:
        - "x_centers": window center positions [n_win] (in sim coords)
        - "k_axis": wavenumber array [n_k] (rad / sim-length)

        If transverse_reduction == "mean":
        - "power_mean": mean power [n_win, n_k]
        - "power_std": std power [n_win, n_k]

        If transverse_reduction == "none":
        - "power_full": full power [*transverse, n_win, n_k]
          or [*transverse, n_win, n_k, C] if component_reduction="none"

    Examples
    --------
    >>> import torch
    >>> from branesim.diagnostics.types import GridSpec
    >>>
    >>> # 1D field: sinusoid
    >>> N = 512
    >>> x = torch.linspace(0, 10*torch.pi, N)
    >>> field = torch.sin(5 * x)  # k=5
    >>> grid = GridSpec(shape=(N,), spacing_sim=x[1].item() - x[0].item())
    >>> spec = local_power_spectrum_along_axis(
    ...     field, grid, axis=0, win_len=128, hop=32
    ... )
    >>> spec["x_centers"].shape
    torch.Size([...])  # number of windows
    >>> spec["k_axis"].shape
    torch.Size([65])  # rfft: win_len//2 + 1
    >>> spec["power_mean"].shape
    torch.Size([..., 65])
    >>>
    >>> # 2D field: wave along x, constant in y
    >>> field_2d = torch.sin(5 * x[:, None]).expand(N, 64)
    >>> grid_2d = GridSpec(shape=(N, 64), spacing_sim=1.0)
    >>> spec_2d = local_power_spectrum_along_axis(
    ...     field_2d, grid_2d, axis=0, win_len=128, hop=32
    ... )
    >>> # power_mean is averaged over y direction
    >>>
    >>> # 2D vector field (3 components)
    >>> field_vec = torch.randn(N, 64, 3)
    >>> spec_vec = local_power_spectrum_along_axis(
    ...     field_vec, grid_2d, axis=0, win_len=128, hop=32,
    ...     component_reduction="sum"
    ... )
    >>> # power_mean has components summed before averaging
    """
    # Validate inputs
    if axis < 0 or axis >= grid.D:
        raise ValueError(f"axis must be in [0, {grid.D-1}], got {axis}")

    has_components = (field.ndim == grid.D + 1)
    if field.ndim not in (grid.D, grid.D + 1):
        raise ValueError(
            f"field.ndim={field.ndim} incompatible with grid.D={grid.D}. "
            f"Expected field.shape = {grid.shape} or {grid.shape} + (C,)"
        )

    if not grid.is_compatible(field):
        raise ValueError(f"Field shape {field.shape} incompatible with grid {grid.shape}")

    if win_len > grid.shape[axis]:
        raise ValueError(
            f"win_len={win_len} larger than grid.shape[axis]={grid.shape[axis]}"
        )

    device = field.device
    dtype = field.dtype if not torch.is_complex(field) else field.real.dtype

    # Strategy:
    # 1. Move analyzed axis to last spatial position
    # 2. If vector field: move components before axis, fold into batch
    # 3. Flatten transverse dims into batch
    # 4. Unfold windows along axis
    # 5. Apply window, FFT, power
    # 6. Reshape back and reduce as requested

    spatial_ndim = grid.D

    # Step 1: Move axis to position (spatial_ndim - 1)
    x = field.movedim(axis, spatial_ndim - 1)
    # x.shape = [*transverse, N] or [*transverse, N, C]

    if has_components:
        # x.shape = [*transverse, N, C]
        # Move C before N: [*transverse, C, N]
        x = x.movedim(-1, -2)
        # Now x.shape = [*transverse, C, N]
        transverse_shape = x.shape[:spatial_ndim - 1]
        C = x.shape[spatial_ndim - 1]
        N = x.shape[spatial_ndim]
        # Flatten: [batch, N] where batch = prod(transverse) * C
        batch_total = int(x.numel() // N)
        x_flat = x.reshape(batch_total, N)
    else:
        # x.shape = [*transverse, N]
        transverse_shape = x.shape[:spatial_ndim - 1]
        N = x.shape[spatial_ndim - 1]
        C = None
        # Flatten: [batch, N]
        batch_total = int(x.numel() // N)
        x_flat = x.reshape(batch_total, N)

    # Step 2: Unfold windows
    # x_flat.shape = [batch, N]
    # unfold -> [batch, N_windows, win_len]
    x_unfolded = x_flat.unfold(dimension=1, size=win_len, step=hop)
    # x_unfolded.shape = [batch, n_win, win_len]
    n_win = x_unfolded.shape[1]

    # Step 3: Apply window function
    if window is not None and window != "none":
        win = _get_window(win_len, window, device, dtype)
        if win is not None:
            # Broadcast: [1, 1, win_len]
            x_unfolded = x_unfolded * win[None, None, :]

    # Step 4: FFT per window
    # x_unfolded.shape = [batch, n_win, win_len]
    # We want FFT along last dim
    if torch.is_complex(field):
        X = torch.fft.fft(x_unfolded, dim=-1)  # [batch, n_win, win_len]
        n_k = win_len
    else:
        # rfft expects real input, don't convert to complex dtype
        X = torch.fft.rfft(x_unfolded, dim=-1)  # [batch, n_win, win_len//2+1]
        n_k = win_len // 2 + 1

    # Step 5: Compute power
    power_flat = torch.abs(X) ** 2  # [batch, n_win, n_k]

    if normalize == "win_len":
        power_flat = power_flat / win_len

    # Step 6: Reshape back to [*transverse, C, n_win, n_k] or [*transverse, n_win, n_k]
    if has_components:
        # batch_total = prod(transverse) * C
        # Reshape to [*transverse, C, n_win, n_k]
        batch_transverse = int(batch_total // C)
        power_reshaped = power_flat.reshape(*transverse_shape, C, n_win, n_k)
    else:
        # Reshape to [*transverse, n_win, n_k]
        power_reshaped = power_flat.reshape(*transverse_shape, n_win, n_k)

    # Step 7: Reduce over components if requested (for vector fields)
    if has_components and component_reduction == "sum":
        # Sum over component dimension (axis = spatial_ndim - 1)
        power_reshaped = power_reshaped.sum(dim=spatial_ndim - 1)
        # Now power_reshaped.shape = [*transverse, n_win, n_k]

    # Step 8: Reduce over transverse dimensions if requested
    if transverse_reduction == "mean":
        if spatial_ndim == 1:
            # No transverse dimensions
            power_mean = power_reshaped  # shape: [n_win, n_k]
            power_std = torch.zeros_like(power_mean)
        else:
            # Reduce over all transverse dimensions
            reduce_dims = list(range(power_reshaped.ndim - 2))  # All except last two (n_win, n_k)
            if len(reduce_dims) > 0:
                power_mean = power_reshaped.mean(dim=reduce_dims)
                power_std = power_reshaped.std(dim=reduce_dims)
            else:
                # Edge case: no transverse dims to reduce
                power_mean = power_reshaped
                power_std = torch.zeros_like(power_mean)

        result_power = {
            "power_mean": power_mean,
            "power_std": power_std,
        }
    else:  # transverse_reduction == "none"
        result_power = {
            "power_full": power_reshaped,
        }

    # Step 9: Compute window center positions and wavenumbers
    # Window centers in grid indices
    win_indices = torch.arange(n_win, device=device, dtype=dtype) * hop + win_len / 2.0
    # Convert to sim coordinates
    x_centers = win_indices * grid.spacing_sim

    # Wavenumber array
    if torch.is_complex(field):
        freqs = torch.fft.fftfreq(win_len, d=float(grid.spacing_sim), device=device)
    else:
        freqs = torch.fft.rfftfreq(win_len, d=float(grid.spacing_sim), device=device)
    k_axis = 2.0 * math.pi * freqs

    return {
        "x_centers": x_centers,
        "k_axis": k_axis,
        **result_power,
    }


def _complex_dtype(real_dtype: torch.dtype) -> torch.dtype:
    """Map real dtype to corresponding complex dtype."""
    if real_dtype == torch.float64:
        return torch.complex128
    return torch.complex64