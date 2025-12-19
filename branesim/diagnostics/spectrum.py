"""
Dimension-agnostic spectrum diagnostics (FFT-based).

Provides tools for computing spatial power spectra in 1D, 2D, and 3D.
Supports:
- 1D FFT along a specific axis
- Full D-dimensional FFT
- Optional windowing to reduce spectral leakage
- Radial averaging for 2D/3D isotropic spectra
"""

from __future__ import annotations
from typing import Literal
import math
import torch
import numpy as np

from .types import GridSpec


def spatial_power_spectrum_1d(
    field: torch.Tensor,
    grid: GridSpec,
    axis: int,
    window: Literal["none", "hann", "hamming", "blackman"] | None = None
) -> dict[str, torch.Tensor]:
    """
    Compute 1D spatial power spectrum along a specified axis.

    For each transverse coordinate, computes the 1D FFT along the specified
    axis and returns the power spectrum. Statistics over transverse directions
    are returned (mean and std).

    Parameters
    ----------
    field : torch.Tensor
        Real or complex field, shape [*grid.shape]
    grid : GridSpec
        Grid specification
    axis : int
        Axis along which to compute FFT (0 to D-1)
    window : str | None
        Window function to apply before FFT:
        - "none" or None: no windowing (default)
        - "hann": Hann window
        - "hamming": Hamming window
        - "blackman": Blackman window

    Returns
    -------
    dict
        Dictionary containing:
        - k_axis: wavenumber array [k_len] (positive frequencies only)
        - power_mean: mean power over transverse directions [k_len]
        - power_std: standard deviation over transverse directions [k_len]
        - power_full: power for all transverse slices [*transverse_shape, k_len]

    Examples
    --------
    >>> # 1D field
    >>> field = torch.randn(256)
    >>> grid = GridSpec(shape=(256,), spacing_sim=1.0)
    >>> result = spatial_power_spectrum_1d(field, grid, axis=0)
    >>> k = result['k_axis']
    >>> power = result['power_mean']
    >>>
    >>> # 2D field: spectrum along x-axis, averaged over y
    >>> field = torch.randn(128, 128)
    >>> grid = GridSpec(shape=(128, 128), spacing_sim=1.0)
    >>> result = spatial_power_spectrum_1d(field, grid, axis=0)
    >>> k_x = result['k_axis']
    >>> power_x_mean = result['power_mean']
    >>> power_x_std = result['power_std']
    """
    if not grid.is_compatible(field):
        raise ValueError(f"Field shape {field.shape} incompatible with grid {grid.shape}")

    # Ensure field is on CPU for numpy interop (if needed for windowing)
    device = field.device
    dtype = field.dtype

    # Apply window if requested
    if window is not None and window != "none":
        # Create window along the specified axis
        axis_len = grid.shape[axis]
        if window == "hann":
            win = torch.hann_window(axis_len, device=device, dtype=torch.float32)
        elif window == "hamming":
            win = torch.hamming_window(axis_len, device=device, dtype=torch.float32)
        elif window == "blackman":
            win = torch.blackman_window(axis_len, device=device, dtype=torch.float32)
        else:
            raise ValueError(f"Unknown window: {window}")

        # Reshape window to broadcast along the axis
        win_shape = [1] * field.ndim
        win_shape[axis] = axis_len
        win = win.reshape(win_shape)

        field = field * win

    # Compute FFT along the specified axis
    # torch.fft.rfft for real input, torch.fft.fft for complex input
    if torch.is_complex(field):
        fft_result = torch.fft.fft(field, dim=axis)
    else:
        fft_result = torch.fft.rfft(field, dim=axis)

    # Compute power: |FFT|^2
    power = torch.abs(fft_result) ** 2

    # Normalize power (optional: by 1/N for energy conservation)
    N = grid.shape[axis]
    power = power / N

    # Compute wavenumber array
    # k = 2π * f / spacing, where f is frequency in cycles per point
    if torch.is_complex(field):
        freqs = torch.fft.fftfreq(N, d=float(grid.spacing_sim), device=device)
        k_axis = 2.0 * math.pi * freqs
    else:
        freqs = torch.fft.rfftfreq(N, d=float(grid.spacing_sim), device=device)
        k_axis = 2.0 * math.pi * freqs

    # Statistics over transverse directions
    if grid.D == 1:
        # No transverse directions
        power_mean = power
        power_std = torch.zeros_like(power)
        power_full = power
    else:
        # Compute mean and std over all dimensions except axis
        reduce_dims = [i for i in range(power.ndim) if i != axis]
        power_mean = power.mean(dim=reduce_dims)
        power_std = power.std(dim=reduce_dims)
        power_full = power

    return {
        "k_axis": k_axis,
        "power_mean": power_mean,
        "power_std": power_std,
        "power_full": power_full,
    }


def spatial_power_spectrum_nd(
    field: torch.Tensor,
    grid: GridSpec,
    window: Literal["none", "hann", "hamming", "blackman"] | None = None
) -> dict[str, torch.Tensor]:
    """
    Compute full N-dimensional spatial power spectrum.

    Parameters
    ----------
    field : torch.Tensor
        Real or complex field, shape [*grid.shape]
    grid : GridSpec
        Grid specification
    window : str | None
        Window function to apply before FFT (applied separably to each axis)

    Returns
    -------
    dict
        Dictionary containing:
        - k_grids: tuple of wavenumber arrays for each dimension
        - power: N-dimensional power spectrum, shape [*k_shape]

    Examples
    --------
    >>> # 2D field
    >>> field = torch.randn(128, 128)
    >>> grid = GridSpec(shape=(128, 128), spacing_sim=1.0)
    >>> result = spatial_power_spectrum_nd(field, grid)
    >>> k_x, k_y = result['k_grids']
    >>> power_2d = result['power']
    """
    if not grid.is_compatible(field):
        raise ValueError(f"Field shape {field.shape} incompatible with grid {grid.shape}")

    device = field.device

    # Apply window if requested (separable windows along each axis)
    if window is not None and window != "none":
        windowed = field
        for axis in range(grid.D):
            axis_len = grid.shape[axis]
            if window == "hann":
                win = torch.hann_window(axis_len, device=device, dtype=torch.float32)
            elif window == "hamming":
                win = torch.hamming_window(axis_len, device=device, dtype=torch.float32)
            elif window == "blackman":
                win = torch.blackman_window(axis_len, device=device, dtype=torch.float32)
            else:
                raise ValueError(f"Unknown window: {window}")

            # Reshape window to broadcast along the axis
            win_shape = [1] * windowed.ndim
            win_shape[axis] = axis_len
            win = win.reshape(win_shape)
            windowed = windowed * win
        field = windowed

    # Compute full N-D FFT
    if torch.is_complex(field):
        fft_result = torch.fft.fftn(field, dim=list(range(grid.D)))
    else:
        fft_result = torch.fft.rfftn(field, dim=list(range(grid.D)))

    # Compute power
    power = torch.abs(fft_result) ** 2

    # Normalize power
    N_total = grid.num_points
    power = power / N_total

    # Compute wavenumber grids
    k_grids = []
    for axis in range(grid.D):
        N = grid.shape[axis]
        # For rfftn, last axis uses rfftfreq
        if not torch.is_complex(field) and axis == grid.D - 1:
            freqs = torch.fft.rfftfreq(N, d=float(grid.spacing_sim), device=device)
        else:
            freqs = torch.fft.fftfreq(N, d=float(grid.spacing_sim), device=device)
        k = 2.0 * math.pi * freqs
        k_grids.append(k)

    return {
        "k_grids": tuple(k_grids),
        "power": power,
    }


def radial_average_spectrum_2d(
    power_2d: torch.Tensor,
    k_x: torch.Tensor,
    k_y: torch.Tensor,
    num_bins: int = 50
) -> dict[str, torch.Tensor]:
    """
    Compute radially averaged power spectrum from 2D spectrum.

    Useful for isotropic systems where we want P(|k|) instead of P(k_x, k_y).

    Parameters
    ----------
    power_2d : torch.Tensor
        2D power spectrum, shape [n_kx, n_ky]
    k_x : torch.Tensor
        Wavenumber array for x-axis, shape [n_kx]
    k_y : torch.Tensor
        Wavenumber array for y-axis, shape [n_ky]
    num_bins : int
        Number of radial bins

    Returns
    -------
    dict
        Dictionary containing:
        - k_radial: radial wavenumber array [num_bins]
        - power_radial: radially averaged power [num_bins]
        - counts: number of (k_x, k_y) points in each bin [num_bins]

    Examples
    --------
    >>> field = torch.randn(128, 128)
    >>> grid = GridSpec(shape=(128, 128), spacing_sim=1.0)
    >>> result_2d = spatial_power_spectrum_nd(field, grid)
    >>> k_x, k_y = result_2d['k_grids']
    >>> power_2d = result_2d['power']
    >>> result_radial = radial_average_spectrum_2d(power_2d, k_x, k_y)
    >>> k_r = result_radial['k_radial']
    >>> P_r = result_radial['power_radial']
    """
    # Create 2D grid of |k| values
    K_X, K_Y = torch.meshgrid(k_x, k_y, indexing='ij')
    K_mag = torch.sqrt(K_X**2 + K_Y**2)

    # Flatten
    k_mag_flat = K_mag.flatten()
    power_flat = power_2d.flatten()

    # Compute radial bins
    k_max = k_mag_flat.max().item()
    k_bins = torch.linspace(0, k_max, num_bins + 1, device=k_x.device)
    k_radial = 0.5 * (k_bins[:-1] + k_bins[1:])

    # Bin the power
    power_radial = torch.zeros(num_bins, device=k_x.device, dtype=power_2d.dtype)
    counts = torch.zeros(num_bins, device=k_x.device, dtype=torch.int64)

    for i in range(num_bins):
        mask = (k_mag_flat >= k_bins[i]) & (k_mag_flat < k_bins[i + 1])
        if mask.any():
            power_radial[i] = power_flat[mask].mean()
            counts[i] = mask.sum()

    return {
        "k_radial": k_radial,
        "power_radial": power_radial,
        "counts": counts,
    }


def radial_average_spectrum_3d(
    power_3d: torch.Tensor,
    k_x: torch.Tensor,
    k_y: torch.Tensor,
    k_z: torch.Tensor,
    num_bins: int = 50
) -> dict[str, torch.Tensor]:
    """
    Compute radially averaged power spectrum from 3D spectrum.

    Parameters
    ----------
    power_3d : torch.Tensor
        3D power spectrum, shape [n_kx, n_ky, n_kz]
    k_x, k_y, k_z : torch.Tensor
        Wavenumber arrays for each axis
    num_bins : int
        Number of radial bins

    Returns
    -------
    dict
        Dictionary containing:
        - k_radial: radial wavenumber array [num_bins]
        - power_radial: radially averaged power [num_bins]
        - counts: number of points in each bin [num_bins]
    """
    # Create 3D grid of |k| values
    K_X, K_Y, K_Z = torch.meshgrid(k_x, k_y, k_z, indexing='ij')
    K_mag = torch.sqrt(K_X**2 + K_Y**2 + K_Z**2)

    # Flatten
    k_mag_flat = K_mag.flatten()
    power_flat = power_3d.flatten()

    # Compute radial bins
    k_max = k_mag_flat.max().item()
    k_bins = torch.linspace(0, k_max, num_bins + 1, device=k_x.device)
    k_radial = 0.5 * (k_bins[:-1] + k_bins[1:])

    # Bin the power
    power_radial = torch.zeros(num_bins, device=k_x.device, dtype=power_3d.dtype)
    counts = torch.zeros(num_bins, device=k_x.device, dtype=torch.int64)

    for i in range(num_bins):
        mask = (k_mag_flat >= k_bins[i]) & (k_mag_flat < k_bins[i + 1])
        if mask.any():
            power_radial[i] = power_flat[mask].mean()
            counts[i] = mask.sum()

    return {
        "k_radial": k_radial,
        "power_radial": power_radial,
        "counts": counts,
    }