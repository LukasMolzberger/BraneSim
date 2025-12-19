"""
Visualization tools for spectrum diagnostics.

Provides plotting functions for:
- 1D power spectra
- 2D/3D radially averaged spectra
- Dispersion relations
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import torch


def plot_power_spectrum_1d(
    k: np.ndarray | torch.Tensor,
    power: np.ndarray | torch.Tensor,
    power_std: np.ndarray | torch.Tensor | None = None,
    k_label: str = "k",
    k_unit: str = "1/sim",
    title: str = "Power Spectrum",
    loglog: bool = True,
    figsize: tuple[float, float] = (10, 6),
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot 1D power spectrum.

    Parameters
    ----------
    k : np.ndarray | torch.Tensor
        Wavenumber array
    power : np.ndarray | torch.Tensor
        Power array (mean if from multiple transverse slices)
    power_std : np.ndarray | torch.Tensor | None
        Standard deviation of power (if multiple slices)
    k_label : str
        Label for wavenumber axis
    k_unit : str
        Unit for wavenumber
    title : str
        Plot title
    loglog : bool
        If True, use log-log scale
    figsize : tuple
        Figure size
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure

    Examples
    --------
    >>> from branesim.diagnostics.spectrum import spatial_power_spectrum_1d
    >>> result = spatial_power_spectrum_1d(field, grid, axis=0)
    >>> fig = plot_power_spectrum_1d(
    ...     result['k_axis'],
    ...     result['power_mean'],
    ...     result['power_std']
    ... )
    """
    # Convert to numpy if needed
    if isinstance(k, torch.Tensor):
        k = k.cpu().numpy()
    if isinstance(power, torch.Tensor):
        power = power.cpu().numpy()
    if power_std is not None and isinstance(power_std, torch.Tensor):
        power_std = power_std.cpu().numpy()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if loglog:
        # Filter out k=0 and negative values
        mask = k > 0
        k_plot = k[mask]
        power_plot = power[mask]
        if power_std is not None:
            power_std_plot = power_std[mask]
        else:
            power_std_plot = None

        ax.loglog(k_plot, power_plot, 'b-', linewidth=1.5)

        if power_std_plot is not None:
            ax.fill_between(k_plot,
                           power_plot - power_std_plot,
                           power_plot + power_std_plot,
                           alpha=0.3, color='blue', label='±1σ')
            ax.legend()
    else:
        ax.plot(k, power, 'b-', linewidth=1.5)

        if power_std is not None:
            ax.fill_between(k,
                           power - power_std,
                           power + power_std,
                           alpha=0.3, color='blue', label='±1σ')
            ax.legend()

    ax.set_xlabel(f"{k_label} ({k_unit})")
    ax.set_ylabel("Power")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_radial_spectrum(
    k_radial: np.ndarray | torch.Tensor,
    power_radial: np.ndarray | torch.Tensor,
    k_label: str = "|k|",
    k_unit: str = "1/sim",
    title: str = "Radially Averaged Power Spectrum",
    loglog: bool = True,
    figsize: tuple[float, float] = (10, 6),
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot radially averaged power spectrum (for 2D/3D).

    Parameters
    ----------
    k_radial : np.ndarray | torch.Tensor
        Radial wavenumber array
    power_radial : np.ndarray | torch.Tensor
        Radially averaged power
    k_label : str
        Label for wavenumber axis
    k_unit : str
        Unit for wavenumber
    title : str
        Plot title
    loglog : bool
        If True, use log-log scale
    figsize : tuple
        Figure size
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure

    Examples
    --------
    >>> from branesim.diagnostics.spectrum import (
    ...     spatial_power_spectrum_nd,
    ...     radial_average_spectrum_2d
    ... )
    >>> result_2d = spatial_power_spectrum_nd(field, grid)
    >>> k_x, k_y = result_2d['k_grids']
    >>> power_2d = result_2d['power']
    >>> result_radial = radial_average_spectrum_2d(power_2d, k_x, k_y)
    >>> fig = plot_radial_spectrum(
    ...     result_radial['k_radial'],
    ...     result_radial['power_radial']
    ... )
    """
    # Convert to numpy if needed
    if isinstance(k_radial, torch.Tensor):
        k_radial = k_radial.cpu().numpy()
    if isinstance(power_radial, torch.Tensor):
        power_radial = power_radial.cpu().numpy()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if loglog:
        # Filter out k=0 and negative/zero power
        mask = (k_radial > 0) & (power_radial > 0)
        k_plot = k_radial[mask]
        power_plot = power_radial[mask]

        ax.loglog(k_plot, power_plot, 'b-', linewidth=1.5, marker='o', markersize=3)
    else:
        ax.plot(k_radial, power_radial, 'b-', linewidth=1.5, marker='o', markersize=3)

    ax.set_xlabel(f"{k_label} ({k_unit})")
    ax.set_ylabel("Power")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_dispersion_relation(
    k: np.ndarray | torch.Tensor,
    omega_measured: np.ndarray | torch.Tensor,
    omega_expected: np.ndarray | torch.Tensor | None = None,
    k_label: str = "k",
    k_unit: str = "1/m",
    omega_label: str = "ω",
    omega_unit: str = "rad/s",
    title: str = "Dispersion Relation",
    figsize: tuple[float, float] = (10, 6),
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot dispersion relation ω(k).

    Compares measured dispersion to expected (if provided).

    Parameters
    ----------
    k : np.ndarray | torch.Tensor
        Wavenumber array
    omega_measured : np.ndarray | torch.Tensor
        Measured angular frequency array
    omega_expected : np.ndarray | torch.Tensor | None
        Expected angular frequency (e.g., from theory)
    k_label : str
        Label for wavenumber axis
    k_unit : str
        Unit for wavenumber
    omega_label : str
        Label for frequency axis
    omega_unit : str
        Unit for frequency
    title : str
        Plot title
    figsize : tuple
        Figure size
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    # Convert to numpy if needed
    if isinstance(k, torch.Tensor):
        k = k.cpu().numpy()
    if isinstance(omega_measured, torch.Tensor):
        omega_measured = omega_measured.cpu().numpy()
    if omega_expected is not None and isinstance(omega_expected, torch.Tensor):
        omega_expected = omega_expected.cpu().numpy()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax.plot(k, omega_measured, 'bo', markersize=4, label='Measured')

    if omega_expected is not None:
        ax.plot(k, omega_expected, 'r-', linewidth=2, label='Expected')

    ax.set_xlabel(f"{k_label} ({k_unit})")
    ax.set_ylabel(f"{omega_label} ({omega_unit})")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_power_spectrum_2d(
    k_x: np.ndarray | torch.Tensor,
    k_y: np.ndarray | torch.Tensor,
    power_2d: np.ndarray | torch.Tensor,
    k_unit: str = "1/sim",
    title: str = "2D Power Spectrum",
    cmap: str = "viridis",
    logscale: bool = True,
    figsize: tuple[float, float] = (10, 8),
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot 2D power spectrum as a heatmap.

    Parameters
    ----------
    k_x : np.ndarray | torch.Tensor
        Wavenumber array for x-axis
    k_y : np.ndarray | torch.Tensor
        Wavenumber array for y-axis
    power_2d : np.ndarray | torch.Tensor
        2D power spectrum, shape [n_kx, n_ky]
    k_unit : str
        Unit for wavenumber
    title : str
        Plot title
    cmap : str
        Colormap
    logscale : bool
        If True, plot log10(power)
    figsize : tuple
        Figure size
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    # Convert to numpy if needed
    if isinstance(k_x, torch.Tensor):
        k_x = k_x.cpu().numpy()
    if isinstance(k_y, torch.Tensor):
        k_y = k_y.cpu().numpy()
    if isinstance(power_2d, torch.Tensor):
        power_2d = power_2d.cpu().numpy()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    extent = [k_x.min(), k_x.max(), k_y.min(), k_y.max()]

    if logscale:
        # Avoid log(0)
        power_plot = np.log10(power_2d + 1e-20)
        label = "log10(Power)"
    else:
        power_plot = power_2d
        label = "Power"

    im = ax.imshow(power_plot.T, origin='lower', extent=extent,
                   cmap=cmap, aspect='auto')
    plt.colorbar(im, ax=ax, label=label)

    ax.set_xlabel(f"$k_x$ ({k_unit})")
    ax.set_ylabel(f"$k_y$ ({k_unit})")
    ax.set_title(title)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig