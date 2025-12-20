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

def plot_spectrogram_1d(
    x_centers: np.ndarray | torch.Tensor,
    k_axis: np.ndarray | torch.Tensor,
    power_mean: np.ndarray | torch.Tensor,
    x_label: str = "x",
    x_unit: str = "sim",
    k_label: str = "k",
    k_unit: str = "rad/sim",
    title: str = "Local Power Spectrum (STFT)",
    cmap: str = "viridis",
    logscale: bool = True,
    plot_peak_k: bool = True,
    figsize: tuple[float, float] = (12, 6),
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot 1D spectrogram (position-resolved power spectrum).

    Shows S(x, k) as a heatmap with optional overlay of dominant wavenumber k_peak(x).

    Parameters
    ----------
    x_centers : np.ndarray | torch.Tensor
        Window center positions [n_win]
    k_axis : np.ndarray | torch.Tensor
        Wavenumber array [n_k]
    power_mean : np.ndarray | torch.Tensor
        Power spectrum [n_win, n_k]
    x_label : str
        Label for position axis
    x_unit : str
        Unit for position
    k_label : str
        Label for wavenumber axis
    k_unit : str
        Unit for wavenumber
    title : str
        Plot title
    cmap : str
        Colormap
    logscale : bool
        If True, plot log10(power)
    plot_peak_k : bool
        If True, overlay curve showing k_peak(x)
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
    >>> from branesim.diagnostics.spectrum import local_power_spectrum_along_axis
    >>> spec = local_power_spectrum_along_axis(field, grid, axis=0, win_len=128, hop=32)
    >>> fig = plot_spectrogram_1d(
    ...     spec['x_centers'],
    ...     spec['k_axis'],
    ...     spec['power_mean'],
    ...     x_label="x", x_unit="nm"
    ... )
    """
    # Convert to numpy if needed
    if isinstance(x_centers, torch.Tensor):
        x_centers = x_centers.cpu().numpy()
    if isinstance(k_axis, torch.Tensor):
        k_axis = k_axis.cpu().numpy()
    if isinstance(power_mean, torch.Tensor):
        power_mean = power_mean.cpu().numpy()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Prepare power for plotting
    if logscale:
        # Avoid log(0)
        power_plot = np.log10(power_mean + 1e-20)
        clabel = "log₁₀(Power)"
    else:
        power_plot = power_mean
        clabel = "Power"

    # Plot spectrogram as heatmap
    # power_mean shape: [n_win, n_k]
    # We want x on horizontal axis, k on vertical
    extent = [x_centers.min(), x_centers.max(), k_axis.min(), k_axis.max()]

    im = ax.imshow(power_plot.T, origin='lower', extent=extent,
                   cmap=cmap, aspect='auto', interpolation='nearest')
    plt.colorbar(im, ax=ax, label=clabel)

    # Overlay peak wavenumber curve
    if plot_peak_k:
        # Find peak k for each window
        k_peak_idx = power_mean.argmax(axis=1)  # [n_win]
        k_peak = k_axis[k_peak_idx]

        ax.plot(x_centers, k_peak, 'r-', linewidth=2, label=f'{k_label}_peak')
        ax.legend(loc='upper right')

    ax.set_xlabel(f"{x_label} ({x_unit})")
    ax.set_ylabel(f"{k_label} ({k_unit})")
    ax.set_title(title)
    ax.grid(True, alpha=0.2, color='white')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_spectrogram_1d_multi_time(
    spectrogram_data: list[dict],
    x_label: str = "x",
    x_unit: str = "sim",
    k_label: str = "k",
    k_unit: str = "rad/sim",
    cmap: str = "viridis",
    logscale: bool = True,
    plot_peak_k: bool = True,
    figsize: tuple[float, float] | None = None,
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot 1D spectrograms for multiple timestamps in separate subplots.

    Parameters
    ----------
    spectrogram_data : list[dict]
        List of dictionaries, each containing:
        - 'x_centers': np.ndarray of window center positions
        - 'k_axis': np.ndarray of wavenumber values
        - 'power_mean': np.ndarray of power spectrum [n_win, n_k]
        - 't_phys_s': float, timestamp in seconds
    x_label : str
        Label for position axis
    x_unit : str
        Unit for position
    k_label : str
        Label for wavenumber axis
    k_unit : str
        Unit for wavenumber
    cmap : str
        Colormap
    logscale : bool
        If True, plot log10(power)
    plot_peak_k : bool
        If True, overlay curve showing k_peak(x)
    figsize : tuple | None
        Figure size (if None, auto-computed based on number of timestamps)
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    n_times = len(spectrogram_data)

    # Auto-compute figsize if not provided
    # Use vertical stacking (n rows, 1 column) matching standard brane plots
    if figsize is None:
        figsize = (14, max(4, 3 * n_times))

    fig, axes = plt.subplots(n_times, 1, figsize=figsize)
    if n_times == 1:
        axes = [axes]

    # Collect all power data to determine common colorbar range
    all_power = []
    for data in spectrogram_data:
        power_mean = data['power_mean']
        if isinstance(power_mean, torch.Tensor):
            power_mean = power_mean.cpu().numpy()
        all_power.append(power_mean)

    if logscale:
        vmin = np.log10(min(p.min() for p in all_power) + 1e-20)
        vmax = np.log10(max(p.max() for p in all_power) + 1e-20)
    else:
        vmin = min(p.min() for p in all_power)
        vmax = max(p.max() for p in all_power)

    for i, (ax, data) in enumerate(zip(axes, spectrogram_data)):
        # Extract data
        x_centers = data['x_centers']
        k_axis = data['k_axis']
        power_mean = data['power_mean']
        t_phys_s = data['t_phys_s']

        # Convert to numpy if needed
        if isinstance(x_centers, torch.Tensor):
            x_centers = x_centers.cpu().numpy()
        if isinstance(k_axis, torch.Tensor):
            k_axis = k_axis.cpu().numpy()
        if isinstance(power_mean, torch.Tensor):
            power_mean = power_mean.cpu().numpy()

        # Prepare power for plotting
        if logscale:
            power_plot = np.log10(power_mean + 1e-20)
            clabel = "log₁₀(Power)"
        else:
            power_plot = power_mean
            clabel = "Power"

        # Plot spectrogram as heatmap
        extent = [x_centers.min(), x_centers.max(), k_axis.min(), k_axis.max()]

        im = ax.imshow(power_plot.T, origin='lower', extent=extent,
                       cmap=cmap, aspect='auto', interpolation='nearest',
                       vmin=vmin, vmax=vmax)

        # Overlay peak wavenumber curve
        if plot_peak_k:
            k_peak_idx = power_mean.argmax(axis=1)
            k_peak = k_axis[k_peak_idx]
            ax.plot(x_centers, k_peak, 'r-', linewidth=2, alpha=0.8)

        ax.set_ylabel(f"{k_label} ({k_unit})", fontsize=11)
        ax.grid(True, alpha=0.2, color='white')

        # Time label with proper precision
        t_fs = t_phys_s * 1e15
        ax.text(0.02, 0.95, f't = {t_fs:.3f} fs',
                transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    axes[-1].set_xlabel(f"{x_label} ({x_unit})", fontsize=12)

    # Add colorbar on the right side of the figure
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(clabel, fontsize=11)

    fig.suptitle("Local Power Spectrum (STFT)", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig
