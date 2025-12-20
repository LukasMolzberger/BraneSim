"""
Visualization tools for Berry phase diagnostics.

Dimension-specific plotting functions that consume DiagnosticResult objects
and display Berry connection, phase profiles, curvature, and quality metrics.

Key principle: Always show quality metrics (amplitude, overlap, validity masks)
alongside phase data to make reliability transparent.
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import torch

from branesim.diagnostics.types import DiagnosticResult, GridSpec


def plot_berry_1d_profile(
    result: DiagnosticResult,
    grid: GridSpec,
    coords: np.ndarray | None = None,
    coord_label: str = "x",
    coord_unit: str = "sim",
    figsize: tuple[float, float] = (12, 8),
    save_path: Path | str | None = None,
    show_quality: bool = True,
) -> Figure:
    """
    Plot 1D Berry phase profile with quality metrics.

    Creates a multi-panel figure showing:
    - Amplitude envelope
    - Berry phase γ(x)
    - Berry connection A_x
    - Overlap magnitude (if available)
    - Valid edge mask

    Parameters
    ----------
    result : DiagnosticResult
        Berry phase diagnostic result
    grid : GridSpec
        Grid specification (must be 1D)
    coords : np.ndarray | None
        Coordinate array for x-axis (if None, use indices)
    coord_label : str
        Label for coordinate axis
    coord_unit : str
        Unit for coordinate axis
    figsize : tuple
        Figure size
    save_path : Path | str | None
        If provided, save figure to this path
    show_quality : bool
        If True, plot quality metrics (amplitude, overlap, etc.)

    Returns
    -------
    Figure
        Matplotlib figure

    Examples
    --------
    >>> result = berry_phase_profile_along_x(psi_hat, amp, cfg)
    >>> grid = GridSpec(shape=(100,), spacing_sim=1.0)
    >>> fig = plot_berry_1d_profile(result, grid, coords=x_nm, coord_label="x", coord_unit="nm")
    >>> plt.show()
    """
    if grid.D != 1:
        raise ValueError("plot_berry_1d_profile requires 1D grid")

    # Generate coordinates if not provided
    if coords is None:
        coords = np.arange(grid.shape[0]) * grid.spacing_sim

    # Extract data
    gamma = _get_field(result.data, "gamma_wrapped", "gamma_unwrapped")
    A_axis = _get_field(result.data, "A_axis", "A_x")
    dphi = _get_field(result.data, "dphi")

    # Extract quality fields
    amp = _get_field(result.quality, "amp", result.data, "amp")
    overlap_abs = _get_field(result.quality, "overlap_abs")
    valid_edge = _get_field(result.quality, "valid_edge")
    mask_point = _get_field(result.quality, "mask_point", "mask")

    # Create figure
    n_panels = 2 + (2 if show_quality else 0)
    fig, axes = plt.subplots(n_panels, 1, figsize=figsize, sharex=True)

    # Panel 0: Amplitude (if available)
    if show_quality and amp is not None:
        ax = axes[0]
        ax.plot(coords, amp, 'k-', linewidth=1)
        ax.set_ylabel(f"Amplitude")
        ax.set_title("Amplitude Envelope")
        ax.grid(True, alpha=0.3)

        # Show threshold if available
        if "amplitude_threshold" in result.meta:
            thresh = result.meta["amplitude_threshold"]
            ax.axhline(thresh, color='r', linestyle='--', alpha=0.5, label=f"Threshold: {thresh:.2e}")
            ax.legend()

        panel_offset = 1
    else:
        panel_offset = 0

    # Panel 1: Berry phase γ(x)
    if gamma is not None:
        ax = axes[panel_offset]
        ax.plot(coords, gamma, 'b-', linewidth=1.5)
        ax.set_ylabel("Berry Phase γ (rad)")
        ax.set_title("Berry Phase Profile")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linestyle='-', alpha=0.2)

        # Mark invalid regions if mask available
        if mask_point is not None:
            invalid = ~mask_point
            if invalid.any():
                ax.fill_between(coords, gamma.min(), gamma.max(),
                                where=invalid, alpha=0.2, color='red',
                                label='Low amplitude')
                ax.legend()

        panel_offset += 1

    # Panel 2: Berry connection A_x
    if A_axis is not None:
        ax = axes[panel_offset]
        # A_axis is defined on edges, so use edge coordinates
        coords_edge = 0.5 * (coords[:-1] + coords[1:])
        ax.plot(coords_edge, A_axis, 'g-', linewidth=1.5)
        ax.set_ylabel(f"Berry Connection A (rad/{coord_unit})")
        ax.set_title("Berry Connection")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linestyle='-', alpha=0.2)

        # Mark invalid edges
        if valid_edge is not None:
            invalid = ~valid_edge
            if invalid.any():
                ax.scatter(coords_edge[invalid], A_axis[invalid],
                          color='red', s=10, alpha=0.5, label='Invalid edge')
                ax.legend()

        panel_offset += 1

    # Panel 3: Overlap magnitude (if available)
    if show_quality and overlap_abs is not None:
        ax = axes[panel_offset]
        coords_edge = 0.5 * (coords[:-1] + coords[1:])
        ax.plot(coords_edge, overlap_abs, 'orange', linewidth=1)
        ax.set_ylabel("Overlap |⟨u|u'⟩|")
        ax.set_title("Wavefunction Overlap")
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1.1])

        # Show threshold
        if "overlap_threshold" in result.meta:
            thresh = result.meta["overlap_threshold"]
            ax.axhline(thresh, color='r', linestyle='--', alpha=0.5, label=f"Threshold: {thresh:.2e}")
            ax.legend()

    # X-axis label on bottom panel
    axes[-1].set_xlabel(f"{coord_label} ({coord_unit})")

    # Overall title
    title = f"Berry Phase Diagnostics: {result.name}"
    if result.t_phys_s is not None:
        title += f" (t = {result.t_phys_s:.2e} s)"
    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_berry_2d_map(
    result: DiagnosticResult,
    grid: GridSpec,
    field_name: str = "curvature",
    coords: tuple[np.ndarray, np.ndarray] | None = None,
    coord_labels: tuple[str, str] = ("x", "y"),
    coord_unit: str = "sim",
    figsize: tuple[float, float] = (10, 8),
    cmap: str = "RdBu_r",
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot 2D map of Berry curvature or other 2D field.

    Parameters
    ----------
    result : DiagnosticResult
        Berry diagnostic result
    grid : GridSpec
        Grid specification (must be 2D)
    field_name : str
        Field to plot ("curvature", "gamma_wrapped", etc.)
    coords : tuple[np.ndarray, np.ndarray] | None
        Coordinate arrays (x, y)
    coord_labels : tuple[str, str]
        Labels for axes
    coord_unit : str
        Unit for coordinates
    figsize : tuple
        Figure size
    cmap : str
        Colormap
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure

    Examples
    --------
    >>> result = berry_plaquette_curvature(psi_hat, amp, axes=(0, 1), cfg)
    >>> grid = GridSpec(shape=(64, 64), spacing_sim=1.0)
    >>> fig = plot_berry_2d_map(result, grid, field_name="curvature")
    """
    if grid.D != 2:
        raise ValueError("plot_berry_2d_map requires 2D grid")

    # Get field
    field = _get_field(result.data, field_name)
    if field is None:
        raise ValueError(f"Field '{field_name}' not found in result")

    # Generate coordinates if not provided
    if coords is None:
        x = np.arange(field.shape[0]) * grid.spacing_sim
        y = np.arange(field.shape[1]) * grid.spacing_sim
    else:
        x, y = coords

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Plot 2D map
    extent = [x.min(), x.max(), y.min(), y.max()]
    im = ax.imshow(field.T, origin='lower', extent=extent, cmap=cmap, aspect='auto')
    plt.colorbar(im, ax=ax, label=field_name)

    ax.set_xlabel(f"{coord_labels[0]} ({coord_unit})")
    ax.set_ylabel(f"{coord_labels[1]} ({coord_unit})")

    # Title
    title = f"{field_name}: {result.name}"
    if result.t_phys_s is not None:
        title += f" (t = {result.t_phys_s:.2e} s)"
    ax.set_title(title)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def _get_field(
    *dicts_or_results,
    default=None
) -> np.ndarray | None:
    """
    Helper to retrieve a field from multiple possible locations.

    Tries each key in each dict/result until found.
    Converts torch.Tensor to numpy if needed.
    """
    keys = [k for k in dicts_or_results if isinstance(k, str)]
    sources = [d for d in dicts_or_results if isinstance(d, dict)]

    for key in keys:
        for source in sources:
            if key in source:
                val = source[key]
                if isinstance(val, torch.Tensor):
                    return val.cpu().numpy()
                elif isinstance(val, np.ndarray):
                    return val
                else:
                    return None

    return default


def plot_berry_connection_profiles(
    results: list[DiagnosticResult],
    grid: GridSpec,
    coords: np.ndarray | None = None,
    coord_label: str = "x",
    coord_unit: str = "sim",
    time_labels: list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot Berry connection profiles for multiple timesteps in separate subplots.

    Parameters
    ----------
    results : list[DiagnosticResult]
        List of Berry diagnostic results (one per timestep)
    grid : GridSpec
        Grid specification
    coords : np.ndarray | None
        Coordinate array
    coord_label : str
        Label for coordinate axis
    coord_unit : str
        Unit for coordinate axis
    time_labels : list[str] | None
        Labels for each timestep (if None, use result.t_phys_s)
    figsize : tuple | None
        Figure size (if None, auto-computed based on number of subplots)
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    if coords is None:
        # Use edge coordinates (connection is defined on edges)
        coords = (np.arange(grid.shape[0] - 1) + 0.5) * grid.spacing_sim

    n_times = len(results)

    # Auto-compute figsize if not provided
    # Use vertical stacking (n rows, 1 column) matching standard brane plots
    if figsize is None:
        figsize = (14, max(4, 3 * n_times))

    fig, axes = plt.subplots(n_times, 1, figsize=figsize)
    if n_times == 1:
        axes = [axes]

    # Collect all data to determine common y-axis range
    all_A_axis = []
    for result in results:
        A_axis = _get_field(result.data, "A_axis", "A_x")
        if A_axis is not None:
            all_A_axis.append(A_axis)

    if all_A_axis:
        y_min = min(np.min(A) for A in all_A_axis)
        y_max = max(np.max(A) for A in all_A_axis)
        y_margin = 0.1 * (y_max - y_min) if y_max > y_min else 0.1
        y_lim = (y_min - y_margin, y_max + y_margin)
    else:
        y_lim = None

    for i, (ax, result) in enumerate(zip(axes, results)):
        A_axis = _get_field(result.data, "A_axis", "A_x")
        if A_axis is None:
            continue

        if time_labels:
            label = time_labels[i]
        elif result.t_phys_s is not None:
            label = f"t = {result.t_phys_s:.2e} s"
        else:
            label = f"Snapshot {i}"

        ax.plot(coords, A_axis, 'g-', linewidth=2)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linestyle='-', alpha=0.2)

        if y_lim is not None:
            ax.set_ylim(y_lim)

        ax.set_ylabel(f"A (rad/{coord_unit})", fontsize=11)

        # Time label
        ax.text(0.02, 0.95, label,
                transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    axes[-1].set_xlabel(f"{coord_label} ({coord_unit})", fontsize=12)

    fig.suptitle("Berry Connection Profiles", fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


def plot_berry_phase_profiles(
    results: list[DiagnosticResult],
    grid: GridSpec,
    coords: np.ndarray | None = None,
    coord_label: str = "x",
    coord_unit: str = "sim",
    time_labels: list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    save_path: Path | str | None = None,
) -> Figure:
    """
    Plot Berry phase profiles for multiple timesteps in separate subplots.

    Parameters
    ----------
    results : list[DiagnosticResult]
        List of Berry diagnostic results
    grid : GridSpec
        Grid specification
    coords : np.ndarray | None
        Coordinate array
    coord_label : str
        Label for coordinate axis
    coord_unit : str
        Unit for coordinate axis
    time_labels : list[str] | None
        Labels for each timestep
    figsize : tuple | None
        Figure size (if None, auto-computed based on number of subplots)
    save_path : Path | str | None
        If provided, save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    if coords is None:
        coords = np.arange(grid.shape[0]) * grid.spacing_sim

    n_times = len(results)

    # Auto-compute figsize if not provided
    # Use vertical stacking (n rows, 1 column) matching standard brane plots
    if figsize is None:
        figsize = (14, max(4, 3 * n_times))

    fig, axes = plt.subplots(n_times, 1, figsize=figsize)
    if n_times == 1:
        axes = [axes]

    # Collect all data to determine common y-axis range
    all_gamma = []
    for result in results:
        gamma = _get_field(result.data, "gamma_wrapped", "gamma_unwrapped")
        if gamma is not None:
            all_gamma.append(gamma)

    if all_gamma:
        y_min = min(np.min(g) for g in all_gamma)
        y_max = max(np.max(g) for g in all_gamma)
        y_margin = 0.1 * (y_max - y_min) if y_max > y_min else 0.1
        y_lim = (y_min - y_margin, y_max + y_margin)
    else:
        y_lim = None

    for i, (ax, result) in enumerate(zip(axes, results)):
        gamma = _get_field(result.data, "gamma_wrapped", "gamma_unwrapped")
        if gamma is None:
            continue

        if time_labels:
            label = time_labels[i]
        elif result.t_phys_s is not None:
            label = f"t = {result.t_phys_s:.2e} s"
        else:
            label = f"Snapshot {i}"

        ax.plot(coords, gamma, 'b-', linewidth=2)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linestyle='-', alpha=0.2)

        if y_lim is not None:
            ax.set_ylim(y_lim)

        ax.set_ylabel("γ (rad)", fontsize=11)

        # Time label
        ax.text(0.02, 0.95, label,
                transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    axes[-1].set_xlabel(f"{coord_label} ({coord_unit})", fontsize=12)

    fig.suptitle("Berry Phase Profiles", fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig