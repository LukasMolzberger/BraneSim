"""
Visualization tools for symplectic band structure and Berry phase diagnostics.

This module provides plotting functions for:
- Band structures ω(k)
- Berry phase profiles
- Wilson loop eigenphases
- Polarization vector fields

Note: Unlike core diagnostics, visualization MAY be dimension-specific.
      Different plotting approaches are used for 1D, 2D, and 3D k-spaces.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Optional, List, Tuple

from .symplectic_types import (
    SymplecticBandResult,
    SymplecticWilsonResult,
)


def plot_band_structure(
    result: SymplecticBandResult,
    band_indices: Optional[List[int]] = None,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
    show_degeneracies: bool = True,
) -> Tuple[Figure, Axes]:
    """
    Plot band structure ω(k) along k-path.

    Args:
        result: SymplecticBandResult from solver
        band_indices: Which bands to plot (None = all)
        ax: Existing axes to plot on (None = create new figure)
        title: Custom title (None = auto-generate)
        show_degeneracies: If True, mark degeneracies with shading

    Returns:
        (fig, ax): Matplotlib figure and axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    # Determine which bands to plot
    if band_indices is None:
        band_indices = list(range(result.n_modes))

    # Compute path parameter (distance along k-path)
    k_points = result.kpath.k_points.cpu().numpy()
    if k_points.shape[1] == 1:
        # 1D: use k directly
        k_param = k_points[:, 0]
        xlabel = r"Wave vector $k$ [rad/m]"
    else:
        # Multi-D: use cumulative distance
        dk = np.diff(k_points, axis=0)
        ds = np.linalg.norm(dk, axis=1)
        k_param = np.concatenate([[0], np.cumsum(ds)])
        xlabel = r"Path parameter $s$"

    # Plot each band
    omega = result.omega.cpu().numpy()
    for band_idx in band_indices:
        ax.plot(k_param, omega[:, band_idx], label=f"Band {band_idx}", linewidth=2)

    # Mark degeneracies if requested
    if show_degeneracies and "degeneracy_clusters" in result.meta:
        clusters_all = result.meta["degeneracy_clusters"]
        for k_idx, clusters in enumerate(clusters_all):
            if len(clusters) > 0:
                for cluster in clusters:
                    # Shade degenerate region
                    omega_cluster = omega[k_idx, cluster]
                    k_val = k_param[k_idx]
                    ax.axvspan(k_val - k_param[1]*0.1, k_val + k_param[1]*0.1,
                              alpha=0.1, color='red')

    # Labels and formatting
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(r"Frequency $\omega$ [rad/s]", fontsize=12)

    if title is None:
        title = f"Band Structure: {result.kpath.label}"
    ax.set_title(title, fontsize=14)

    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, ax


def plot_berry_phase_profile(
    result: SymplecticBandResult,
    band_index: int,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
) -> Tuple[Figure, Axes]:
    """
    Plot Berry phase profile along k-path for a single band (U(1) case).

    For a closed path, this should integrate to give the total Berry phase.

    Args:
        result: SymplecticBandResult from solver
        band_index: Index of band to analyze
        ax: Existing axes (None = create new)
        title: Custom title (None = auto-generate)

    Returns:
        (fig, ax): Matplotlib figure and axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    # Import Berry connection function
    from .symplectic_berry import compute_berry_connection_along_path

    # Compute Berry connection A_k
    connection = compute_berry_connection_along_path(result, band_index)
    connection_np = connection.cpu().numpy()

    # Compute path parameter
    k_points = result.kpath.k_points.cpu().numpy()
    if k_points.shape[1] == 1:
        k_param = k_points[:, 0]
        xlabel = r"Wave vector $k$ [rad/m]"
    else:
        dk = np.diff(k_points, axis=0)
        ds = np.linalg.norm(dk, axis=1)
        k_param = np.concatenate([[0], np.cumsum(ds)])
        xlabel = r"Path parameter $s$"

    # Plot Berry connection
    ax.plot(k_param, connection_np, linewidth=2, label="Berry connection $A_k$")

    # If closed path, show accumulated phase
    if result.kpath.closed:
        # Integrate using trapezoidal rule
        accumulated_phase = np.cumsum(connection_np) * (k_param[1] - k_param[0])
        ax_twin = ax.twinx()
        ax_twin.plot(k_param, accumulated_phase, 'r--', alpha=0.7,
                    label="Accumulated phase")
        ax_twin.set_ylabel("Accumulated phase [rad]", fontsize=12, color='r')
        ax_twin.tick_params(axis='y', labelcolor='r')

        total_phase = accumulated_phase[-1]
        ax.text(0.05, 0.95, f"Total Berry phase: {total_phase:.4f} rad\n"
                             f"                  = {total_phase/np.pi:.4f} π",
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Labels
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(r"Berry connection $A_k$", fontsize=12)

    if title is None:
        title = f"Berry Phase Profile: Band {band_index}"
    ax.set_title(title, fontsize=14)

    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, ax


def plot_wilson_eigenphases(
    wilson_result: SymplecticWilsonResult,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
    show_reference_lines: bool = True,
) -> Tuple[Figure, Axes]:
    """
    Plot Wilson loop eigenphases on the unit circle.

    For U(1): single point on circle (Berry phase)
    For U(N): N points on circle (non-Abelian holonomy)

    Args:
        wilson_result: SymplecticWilsonResult from Berry computation
        ax: Existing axes (None = create new)
        title: Custom title (None = auto-generate)
        show_reference_lines: If True, show reference angles (0, π/2, π, etc.)

    Returns:
        (fig, ax): Matplotlib figure and axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    else:
        fig = ax.figure

    # Plot unit circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(theta, np.ones_like(theta), 'k-', alpha=0.3, linewidth=1)

    # Plot eigenphases
    eigenphases = wilson_result.eigenphases
    N = len(eigenphases)

    # Convert to [0, 2π] for polar plot
    phases_positive = np.mod(eigenphases, 2*np.pi)

    for i, phase in enumerate(phases_positive):
        ax.plot([phase], [1.0], 'o', markersize=12, label=f"λ_{i+1}")

    # Reference lines
    if show_reference_lines:
        for ref_angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
            ax.plot([ref_angle, ref_angle], [0, 1], 'k--', alpha=0.2, linewidth=0.5)

    # Labels
    if title is None:
        if N == 1:
            title = f"U(1) Berry Phase: {eigenphases[0]:.4f} rad ({eigenphases[0]/np.pi:.3f} π)"
        else:
            title = f"U({N}) Wilson Loop Eigenphases"
    ax.set_title(title, fontsize=14, pad=20)

    # Add trace and spinorial info
    info_text = f"tr(W) = {wilson_result.trace:.4f}\n"
    if N > 1:
        info_text += f"Spinorial: {wilson_result.meta.get('is_spinorial', False)}"

    ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    if N <= 5:
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

    return fig, ax


def plot_polarization_vectors(
    result: SymplecticBandResult,
    k_index: int,
    band_indices: Optional[List[int]] = None,
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
    component_labels: Optional[List[str]] = None,
) -> Tuple[Figure, Axes]:
    """
    Plot polarization vectors (q-frames) for selected bands at a given k-point.

    Shows the embedding-space components of each polarization vector as a bar chart.

    Args:
        result: SymplecticBandResult from solver
        k_index: Index of k-point to visualize
        band_indices: Which bands to show (None = all)
        ax: Existing axes (None = create new)
        title: Custom title (None = auto-generate)
        component_labels: Labels for embedding components (None = X^0, X^1, ...)

    Returns:
        (fig, ax): Matplotlib figure and axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    if band_indices is None:
        band_indices = list(range(result.n_modes))

    if component_labels is None:
        component_labels = [f"X^{i}" for i in range(result.embedding_dim)]

    # Extract frames at k_index
    frames = result.frames_q[k_index].cpu().numpy()  # [embedding_dim, n_modes]

    # Compute magnitudes (real and imaginary parts)
    x = np.arange(result.embedding_dim)
    width = 0.8 / len(band_indices)

    for i, band_idx in enumerate(band_indices):
        polarization = frames[:, band_idx]

        # Plot real and imaginary parts
        offset = (i - len(band_indices)/2) * width
        ax.bar(x + offset, np.abs(polarization), width,
               label=f"Band {band_idx}", alpha=0.7)

    # Labels
    ax.set_xlabel("Embedding Component", fontsize=12)
    ax.set_ylabel("Polarization Magnitude", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(component_labels)

    if title is None:
        k_val = result.kpath.k_points[k_index].cpu().numpy()
        title = f"Polarization Vectors at k = {k_val}"
    ax.set_title(title, fontsize=14)

    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    return fig, ax


def plot_band_structure_2d(
    result: SymplecticBandResult,
    band_index: int,
    grid_shape: Tuple[int, int],
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
) -> Tuple[Figure, Axes]:
    """
    Plot 2D band structure ω(kx, ky) as a heatmap.

    This requires the k-path to be a 2D grid in k-space.

    Args:
        result: SymplecticBandResult from solver
        band_index: Which band to plot
        grid_shape: Shape of k-space grid (n_kx, n_ky)
        ax: Existing axes (None = create new)
        title: Custom title (None = auto-generate)

    Returns:
        (fig, ax): Matplotlib figure and axes

    Note:
        This function assumes result.kpath.k_points can be reshaped into a 2D grid.
    """
    if result.kpath.d != 2:
        raise ValueError("plot_band_structure_2d requires 2D k-path")

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    else:
        fig = ax.figure

    # Extract frequency data
    omega = result.omega[:, band_index].cpu().numpy()

    # Reshape to grid
    try:
        omega_grid = omega.reshape(grid_shape)
    except ValueError:
        raise ValueError(
            f"Cannot reshape omega with {len(omega)} points into grid {grid_shape}"
        )

    # Extract k-points
    k_points = result.kpath.k_points.cpu().numpy()
    kx = k_points[:, 0].reshape(grid_shape)
    ky = k_points[:, 1].reshape(grid_shape)

    # Plot heatmap
    c = ax.pcolormesh(kx, ky, omega_grid, shading='auto', cmap='viridis')
    cbar = fig.colorbar(c, ax=ax, label=r"Frequency $\omega$ [rad/s]")

    # Labels
    ax.set_xlabel(r"$k_x$ [rad/m]", fontsize=12)
    ax.set_ylabel(r"$k_y$ [rad/m]", fontsize=12)

    if title is None:
        title = f"Band {band_index}: ω(kx, ky)"
    ax.set_title(title, fontsize=14)

    ax.set_aspect('equal')
    plt.tight_layout()

    return fig, ax