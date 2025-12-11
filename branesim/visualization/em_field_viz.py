"""
EM Field Visualization

Generalized visualization tools for electromagnetic fields along centerlines.
Handles both straight waveguides and curved paths (e.g., toroidal geometries).
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Optional, Tuple, List


def visualize_em_fields_along_centerline(
    centerline: np.ndarray,
    E_field: np.ndarray,
    B_field: np.ndarray,
    output_path: str,
    title: str = "EM Fields Along Centerline",
    views: Optional[List[dict]] = None,
    arrow_scale: float = 1.0,
    subsample_step: Optional[int] = None,
    figsize: Tuple[int, int] = (15, 5),
    dpi: int = 300,
    background_draw_func: Optional[callable] = None,
    separate_files: bool = False
) -> str:
    """
    Visualize E and B field vectors along a centerline in 3D.

    Args:
        centerline: (N, 3) array of positions along the centerline
        E_field: (N, 3) array of electric field vectors at each centerline point
        B_field: (N, 3) array of magnetic field vectors at each centerline point
        output_path: Path to save the output figure
        title: Figure title
        views: List of view dictionaries with 'elev', 'azim', 'title' keys
               If None, uses default 3-view layout. Can include 'filename' key for separate files
        arrow_scale: Scale factor for arrow lengths (relative to centerline extent)
        subsample_step: Step size for subsampling arrows (if None, auto-computed)
        figsize: Figure size in inches
        dpi: DPI for saved figure
        background_draw_func: Optional function(ax) to draw custom background elements
        separate_files: If True, save each view as a separate file (requires 'filename' in views)

    Returns:
        Path to saved figure (or first file if separate_files=True)
    """
    num_points = centerline.shape[0]

    # Default views if not provided
    if views is None:
        views = [
            dict(elev=25, azim=-60, title="Oblique view"),
            dict(elev=10, azim=30, title="Side view"),
            dict(elev=80, azim=-90, title="Top view"),
        ]

    # Auto-compute subsample step if not provided
    if subsample_step is None:
        subsample_step = max(1, num_points // 40)

    # Compute arrow length based on centerline extent
    extent_x = np.ptp(centerline[:, 0])
    extent_y = np.ptp(centerline[:, 1])
    extent_z = np.ptp(centerline[:, 2])
    max_extent = max(extent_x, extent_y, extent_z)

    # Arrow length as fraction of domain extent
    arrow_len = arrow_scale * max_extent * 0.05

    # Determine if we're making separate files or one combined figure
    if separate_files:
        saved_paths = []

    # Create figure(s) with multiple views
    for view_idx, view in enumerate(views):
        if separate_files:
            fig = plt.figure(figsize=(8, 8) if separate_files else figsize)
            ax = fig.add_subplot(111, projection="3d")
            subplot_idx = 1
        else:
            if view_idx == 0:
                fig = plt.figure(figsize=figsize)
                fig.suptitle(title, fontsize=16, fontweight='bold')
            ax = fig.add_subplot(1, len(views), view_idx + 1, projection="3d")
            subplot_idx = view_idx + 1

        # Draw custom background elements (e.g., torus hull)
        if background_draw_func is not None:
            background_draw_func(ax)

        # Draw centerline
        ax.plot(
            centerline[:, 0],
            centerline[:, 1],
            centerline[:, 2],
            linewidth=2.0,
            color="k",
            label="Centerline"
        )

        # Draw E and B field arrows at subsampled points
        for idx in range(0, num_points, subsample_step):
            p = centerline[idx]

            # Electric field arrow (blue)
            e = E_field[idx]
            e_norm = np.linalg.norm(e)
            if e_norm > 1e-12:
                e_normalized = e / e_norm
                q_e = p + arrow_len * e_normalized
                ax.plot(
                    [p[0], q_e[0]],
                    [p[1], q_e[1]],
                    [p[2], q_e[2]],
                    linewidth=1.7,
                    color="C0",  # blue
                    alpha=0.8
                )

            # Magnetic field arrow (orange)
            b = B_field[idx]
            b_norm = np.linalg.norm(b)
            if b_norm > 1e-12:
                b_normalized = b / b_norm
                q_b = p + arrow_len * b_normalized
                ax.plot(
                    [p[0], q_b[0]],
                    [p[1], q_b[1]],
                    [p[2], q_b[2]],
                    linewidth=1.7,
                    color="C1",  # orange
                    alpha=0.8
                )

        # Set equal aspect ratio (unless background function handles it)
        if background_draw_func is None:
            _set_equal_aspect_3d(ax, centerline[:, 0], centerline[:, 1], centerline[:, 2])

        # Set view angle
        ax.view_init(elev=view["elev"], azim=view["azim"])
        ax.set_title(view["title"], fontsize=14 if separate_files else 12,
                    fontweight='bold' if separate_files else 'normal')
        ax.set_xlabel("X¹ (brane)" if separate_files else "X [m]",
                     fontsize=12 if separate_files else 10)
        ax.set_ylabel("X² (brane)" if separate_files else "Y [m]",
                     fontsize=12 if separate_files else 10)
        ax.set_zlabel("X³ (brane)" if separate_files else "Z [m]",
                     fontsize=12 if separate_files else 10)

        # Add legend
        if subplot_idx == 1 or separate_files:
            # Create proxy artists for legend
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='k', linewidth=2.0, label='Centerline'),
                Line2D([0], [0], color='C0', linewidth=1.7, label='E-field'),
                Line2D([0], [0], color='C1', linewidth=1.7, label='B-field'),
            ]
            ax.legend(handles=legend_elements, loc='upper right',
                     fontsize=10 if separate_files else 9)

        # Save separate files or continue building combined figure
        if separate_files:
            plt.tight_layout()
            # If output_path is a directory, join with filename from view
            import os
            if os.path.isdir(output_path):
                file_path = os.path.join(output_path, view.get('filename', f"view_{view_idx}.png"))
            else:
                # output_path is a file path, use it as base for numbered files
                file_path = view.get('filename', f"{output_path}_{view_idx}.png")
            plt.savefig(file_path, dpi=dpi, bbox_inches='tight')
            plt.close(fig)
            saved_paths.append(file_path)

    # Save combined figure if not separate files
    if not separate_files:
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        return output_path
    else:
        return saved_paths[0] if saved_paths else output_path


def visualize_em_field_components_2d(
    centerline: np.ndarray,
    E_field: np.ndarray,
    B_field: np.ndarray,
    output_path: str,
    propagation_axis: int = 0,
    title: str = "EM Field Components Along Propagation",
    figsize: Tuple[int, int] = (12, 8),
    dpi: int = 300
) -> str:
    """
    Visualize E and B field components as 2D line plots along the propagation direction.

    Args:
        centerline: (N, 3) array of positions along the centerline
        E_field: (N, 3) array of electric field vectors
        B_field: (N, 3) array of magnetic field vectors
        output_path: Path to save the output figure
        propagation_axis: Index of propagation axis (0=x, 1=y, 2=z)
        title: Figure title
        figsize: Figure size in inches
        dpi: DPI for saved figure

    Returns:
        Path to saved figure
    """
    # Extract propagation coordinate
    s = centerline[:, propagation_axis]

    # Create figure with subplots
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    fig.suptitle(title, fontsize=14, fontweight='bold')

    axis_labels = ['X', 'Y', 'Z']
    colors = ['C0', 'C1', 'C2']

    # Plot E-field components
    for i in range(3):
        axes[0].plot(s, E_field[:, i], label=f'E_{axis_labels[i]}',
                    color=colors[i], linewidth=1.5, alpha=0.8)
    axes[0].set_ylabel('E-field [V/m]', fontsize=12)
    axes[0].legend(loc='upper right', fontsize=10)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Electric Field Components', fontsize=12)

    # Plot B-field components
    for i in range(3):
        axes[1].plot(s, B_field[:, i], label=f'B_{axis_labels[i]}',
                    color=colors[i], linewidth=1.5, alpha=0.8)
    axes[1].set_ylabel('B-field [T]', fontsize=12)
    axes[1].set_xlabel(f'{axis_labels[propagation_axis]} [m]', fontsize=12)
    axes[1].legend(loc='upper right', fontsize=10)
    axes[1].grid(True, alpha=0.3)
    axes[1].set_title('Magnetic Field Components', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    return output_path


def visualize_em_field_volume_3d(
    positions: np.ndarray,
    E_field: np.ndarray,
    B_field: np.ndarray,
    output_dir: str,
    grid_shape: Optional[Tuple[int, int, int]] = None,
    subsample_factor: int = 5,
    arrow_scale: float = 1.0,
    views: Optional[List[dict]] = None,
    dpi: int = 150,
    min_alpha: float = 0.2,
    max_alpha: float = 0.9
) -> Tuple[List[str], List[str]]:
    """
    Visualize E and B field vectors throughout a 3D volume using arrows with heads.

    Creates separate visualizations for E and B fields to avoid clutter.
    Arrow length and opacity both encode field strength.

    Args:
        positions: (N, 3) array of 3D positions where fields are sampled
        E_field: (N, 3) array of electric field vectors [V/m]
        B_field: (N, 3) array of magnetic field vectors [T]
        output_dir: Directory to save output files
        grid_shape: (nx, ny, nz) grid dimensions for proper 3D subsampling
        subsample_factor: Take every Nth point in each dimension (if grid_shape provided),
                         or every Nth point linearly (if grid_shape not provided)
        arrow_scale: Scale factor for arrow lengths
        views: List of view dicts with 'elev', 'azim', 'title', 'filename_e', 'filename_b'
        dpi: DPI for saved figures
        min_alpha: Minimum opacity for weakest fields
        max_alpha: Maximum opacity for strongest fields

    Returns:
        Tuple of (E_field_paths, B_field_paths) - lists of saved file paths
    """
    import os

    # Subsample to reduce clutter
    if grid_shape is not None:
        # Proper 3D subsampling: take every subsample_factor-th point in each dimension
        nx, ny, nz = grid_shape
        positions_3d = positions.reshape(nx, ny, nz, 3)
        E_field_3d = E_field.reshape(nx, ny, nz, 3)
        B_field_3d = B_field.reshape(nx, ny, nz, 3)

        positions_sub = positions_3d[::subsample_factor, ::subsample_factor, ::subsample_factor].reshape(-1, 3)
        E_field_sub = E_field_3d[::subsample_factor, ::subsample_factor, ::subsample_factor].reshape(-1, 3)
        B_field_sub = B_field_3d[::subsample_factor, ::subsample_factor, ::subsample_factor].reshape(-1, 3)

        num_arrows = positions_sub.shape[0]
        print(f"  3D subsampling: {nx}×{ny}×{nz} → {nx//subsample_factor}×{ny//subsample_factor}×{nz//subsample_factor} = {num_arrows} arrows per field")
    else:
        # Linear subsampling (fallback)
        positions_sub = positions[::subsample_factor]
        E_field_sub = E_field[::subsample_factor]
        B_field_sub = B_field[::subsample_factor]
        print(f"  Linear subsampling: {len(positions)} → {len(positions_sub)} arrows per field")

    # Compute field magnitudes for normalization
    E_mag = np.linalg.norm(E_field_sub, axis=1)
    B_mag = np.linalg.norm(B_field_sub, axis=1)

    E_max = E_mag.max() if E_mag.max() > 0 else 1.0
    B_max = B_mag.max() if B_mag.max() > 0 else 1.0

    # Normalize magnitudes to [0, 1] for opacity mapping
    E_mag_norm = E_mag / E_max
    B_mag_norm = B_mag / B_max

    # Map normalized magnitudes to alpha values
    E_alpha = min_alpha + (max_alpha - min_alpha) * E_mag_norm
    B_alpha = min_alpha + (max_alpha - min_alpha) * B_mag_norm

    # Compute arrow scale based on domain extent
    extent_x = np.ptp(positions[:, 0])
    extent_y = np.ptp(positions[:, 1])
    extent_z = np.ptp(positions[:, 2])
    max_extent = max(extent_x, extent_y, extent_z)
    arrow_length_scale = arrow_scale * max_extent * 0.1

    # Default views if not provided
    if views is None:
        views = [
            dict(
                elev=25, azim=-60,
                title_e="E-Field - Oblique View",
                title_b="B-Field - Oblique View",
                filename_e="e_field_volume_oblique.png",
                filename_b="b_field_volume_oblique.png"
            ),
            dict(
                elev=10, azim=30,
                title_e="E-Field - Side View",
                title_b="B-Field - Side View",
                filename_e="e_field_volume_side.png",
                filename_b="b_field_volume_side.png"
            ),
            dict(
                elev=80, azim=-90,
                title_e="E-Field - Top View",
                title_b="B-Field - Top View",
                filename_e="e_field_volume_top.png",
                filename_b="b_field_volume_top.png"
            ),
        ]

    e_paths = []
    b_paths = []

    # Draw centerline for reference (along x-axis at center of y,z)
    x_min, x_max = positions[:, 0].min(), positions[:, 0].max()
    y_center = (positions[:, 1].min() + positions[:, 1].max()) / 2.0
    z_center = (positions[:, 2].min() + positions[:, 2].max()) / 2.0
    centerline_x = np.linspace(x_min, x_max, 100)
    centerline_y = np.full(100, y_center)
    centerline_z = np.full(100, z_center)

    # Create separate figures for E and B fields at each viewing angle
    for view in views:
        # ===== E-FIELD PLOT =====
        fig_e = plt.figure(figsize=(10, 10))
        ax_e = fig_e.add_subplot(111, projection='3d')

        # Draw centerline
        ax_e.plot(centerline_x, centerline_y, centerline_z,
                 color='k', linewidth=2.0, linestyle='--',
                 label='Centerline', alpha=0.7)

        # Plot E-field arrows with varying opacity and length
        for i in range(len(positions_sub)):
            if E_mag[i] > 1e-12:  # Skip near-zero vectors
                pos = positions_sub[i]
                # Arrow length proportional to field magnitude (normalized)
                vec = E_field_sub[i] * (arrow_length_scale / E_max)

                ax_e.quiver(
                    pos[0], pos[1], pos[2],
                    vec[0], vec[1], vec[2],
                    color='C0',  # Blue for E-field
                    alpha=E_alpha[i],
                    arrow_length_ratio=0.3,
                    linewidth=1.5
                )

        # Set equal aspect and view
        _set_equal_aspect_3d(ax_e, positions[:, 0], positions[:, 1], positions[:, 2])
        ax_e.view_init(elev=view['elev'], azim=view['azim'])
        ax_e.set_title(view['title_e'], fontsize=14, fontweight='bold')
        ax_e.set_xlabel('X [m]', fontsize=12)
        ax_e.set_ylabel('Y [m]', fontsize=12)
        ax_e.set_zlabel('Z [m]', fontsize=12)

        # Add colorbar-like legend for field strength
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='C0', linewidth=2, alpha=max_alpha, label='Strong E-field'),
            Line2D([0], [0], color='C0', linewidth=2, alpha=min_alpha, label='Weak E-field'),
        ]
        ax_e.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout()
        e_path = os.path.join(output_dir, view['filename_e'])
        plt.savefig(e_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig_e)
        e_paths.append(e_path)

        # ===== B-FIELD PLOT =====
        fig_b = plt.figure(figsize=(10, 10))
        ax_b = fig_b.add_subplot(111, projection='3d')

        # Draw centerline
        ax_b.plot(centerline_x, centerline_y, centerline_z,
                 color='k', linewidth=2.0, linestyle='--',
                 label='Centerline', alpha=0.7)

        # Plot B-field arrows with varying opacity and length
        for i in range(len(positions_sub)):
            if B_mag[i] > 1e-12:  # Skip near-zero vectors
                pos = positions_sub[i]
                # Arrow length proportional to field magnitude (normalized)
                vec = B_field_sub[i] * (arrow_length_scale / B_max)

                ax_b.quiver(
                    pos[0], pos[1], pos[2],
                    vec[0], vec[1], vec[2],
                    color='C1',  # Orange for B-field
                    alpha=B_alpha[i],
                    arrow_length_ratio=0.3,
                    linewidth=1.5
                )

        # Set equal aspect and view
        _set_equal_aspect_3d(ax_b, positions[:, 0], positions[:, 1], positions[:, 2])
        ax_b.view_init(elev=view['elev'], azim=view['azim'])
        ax_b.set_title(view['title_b'], fontsize=14, fontweight='bold')
        ax_b.set_xlabel('X [m]', fontsize=12)
        ax_b.set_ylabel('Y [m]', fontsize=12)
        ax_b.set_zlabel('Z [m]', fontsize=12)

        # Add legend
        legend_elements = [
            Line2D([0], [0], color='C1', linewidth=2, alpha=max_alpha, label='Strong B-field'),
            Line2D([0], [0], color='C1', linewidth=2, alpha=min_alpha, label='Weak B-field'),
        ]
        ax_b.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout()
        b_path = os.path.join(output_dir, view['filename_b'])
        plt.savefig(b_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig_b)
        b_paths.append(b_path)

    return e_paths, b_paths


def _set_equal_aspect_3d(ax, xs, ys, zs):
    """Set equal aspect ratio for 3D axes based on data ranges."""
    x_min, x_max = np.min(xs), np.max(xs)
    y_min, y_max = np.min(ys), np.max(ys)
    z_min, z_max = np.min(zs), np.max(zs)

    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    ax.set_xlim(x_mid - max_range / 2.0, x_mid + max_range / 2.0)
    ax.set_ylim(y_mid - max_range / 2.0, y_mid + max_range / 2.0)
    ax.set_zlim(z_mid - max_range / 2.0, z_mid + max_range / 2.0)