import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt
import numpy as np

from branesim.geometry.tubular_electron_geometry import (
    TorusKnotParameters,
    sample_torus_knot_centerline,
    compute_frenet_frames,
)
from branesim.geometry.tubular_photon_mode import (
    PhotonModeParameters,
    compute_circular_polarization_EB,
    sample_gaussian_envelope,
)
from branesim.visualization.em_field_viz import visualize_em_fields_along_centerline
from branesim.utils import TestRunManager


# Set this to False if you temporarily want to see ONLY the E/B arrows.
DRAW_CLOUD = False


def _plot_torus_hull(ax, R: float, r: float, num_phi: int = 64, num_theta: int = 32):
    """Draw a semi-transparent torus hull for geometric reference."""
    phi = np.linspace(0.0, 2.0 * np.pi, num_phi)
    theta = np.linspace(0.0, 2.0 * np.pi, num_theta)
    Phi, Theta = np.meshgrid(phi, theta)

    X = (R + r * np.cos(Theta)) * np.cos(Phi)
    Y = (R + r * np.cos(Theta)) * np.sin(Phi)
    Z = r * np.sin(Theta)

    ax.plot_surface(
        X,
        Y,
        Z,
        rstride=2,
        cstride=2,
        linewidth=0.2,
        alpha=0.06,
        edgecolor="k",
    )


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


def visualize_circular_photon_toroidal():
    """
    Visual test: circularly polarized photon mode along a toroidal path.

    - Uses the torus-knot centerline and Frenet frame from the geometry module.
    - Builds a circularly polarized E-field in the local (n, b) plane, with
      total internal phase advance 4π along the loop (spinor-like orientation).
    - Introduces a higher-frequency longitudinal amplitude with 'wave_cycles'
      cos oscillations along the loop (more crests).
    - Constructs a B-field via B ∝ t × E.
    - Saves visualizations to organized test run directory.

    The Möbius / double-loop geometry is present only as the path; the strip
    itself is not drawn here.
    """
    print("=" * 70)
    print("Circularly Polarized Photon - Toroidal Path Visualization")
    print("=" * 70)

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="visualize_circular_photon_toroidal")
    print(run_manager.get_summary())
    # --- Outer geometry (unchanged) -----------------------------------
    torus_params = TorusKnotParameters(
        major_radius=1.0,
        minor_radius=0.36,  # slightly thicker torus
        core_windings=2,
        tube_windings=1,
    )

    num_samples = 900
    centerline = sample_torus_knot_centerline(torus_params, num_samples=num_samples)
    t, n, b = compute_frenet_frames(centerline)

    # --- Photon mode parameters ---------------------------------------
    photon_params = PhotonModeParameters(
        peak_amplitude=1.0,
        sigma_n=0.12,
        sigma_b=0.12,
        extent_sigma=3.5,
        num_radial_samples=8,
        num_angular_samples=32,
        total_phase=4.0 * np.pi,
        phase_offset=0.0,
        wave_cycles=8,   # higher frequency: 8 crests along the loop
        B_over_E=1.0,
    )

    # Electric and magnetic field vectors along the centerline
    e_vectors, b_vectors, phase_internal, a_long = compute_circular_polarization_EB(
        t,
        n,
        b,
        photon_params,
    )

    print(f"\nGenerating EM field visualizations...")
    print(f"  Centerline points: {centerline.shape[0]}")
    print(f"  Total phase advance: {photon_params.total_phase / np.pi:.1f}π")
    print(f"  Wave cycles: {photon_params.wave_cycles}")

    # Define views
    views = [
        dict(elev=25, azim=-60, title="Photon mode - oblique view", filename="tubular_photon_oblique.png"),
        dict(elev=10, azim=30, title="Photon mode - nearly in-plane", filename="tubular_photon_side.png"),
        dict(elev=80, azim=-90, title="Photon mode - top view", filename="tubular_photon_top.png"),
    ]

    # Subsample for drawing E and B vectors
    num = centerline.shape[0]
    step = max(1, num // 40)
    arrow_len = 1.8 * min(photon_params.sigma_n, photon_params.sigma_b)

    # Longitudinal modulation for intensity (for potential cloud drawing)
    longitudinal_modulation = np.abs(a_long)
    field_points, amplitudes = sample_gaussian_envelope(
        centerline,
        n,
        b,
        photon_params,
        longitudinal_modulation=longitudinal_modulation,
    )

    # Normalize amplitudes for color/opacity mapping
    colors = None
    if amplitudes.size > 0 and DRAW_CLOUD:
        amp_min = float(np.min(amplitudes))
        amp_max = float(np.max(amplitudes))
        denom = max(amp_max - amp_min, 1e-12)
        amp_norm = (amplitudes - amp_min) / denom
        cmap = plt.get_cmap()
        colors = cmap(amp_norm)
        gamma = 1.5
        colors[:, 3] = amp_norm ** gamma

    # Compute extent for aspect ratio (add margin for torus hull)
    # Use torus outer radius for extent calculation
    torus_outer_radius = torus_params.major_radius + torus_params.minor_radius + 0.04
    all_x = centerline[:, 0]
    all_y = centerline[:, 1]
    all_z = centerline[:, 2]

    # Extend by torus radius to ensure hull fits
    extent_margin = torus_outer_radius * 1.1

    # Create legend elements once
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='k', linewidth=2.0, label='Centerline'),
        Line2D([0], [0], color='C0', linewidth=1.7, label='E-field'),
        Line2D([0], [0], color='C1', linewidth=1.7, label='B-field'),
    ]

    # Generate separate figure for each view
    for view in views:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")

        # Outer torus hull for context
        _plot_torus_hull(
            ax,
            R=torus_params.major_radius,
            r=torus_params.minor_radius + 0.04,
        )

        # Photon cloud (optional)
        if DRAW_CLOUD and colors is not None:
            ax.scatter(
                field_points[:, 0],
                field_points[:, 1],
                field_points[:, 2],
                s=4,
                c=colors,
                marker="o",
            )

        # Centerline (path of the photon)
        ax.plot(
            centerline[:, 0],
            centerline[:, 1],
            centerline[:, 2],
            linewidth=2.0,
            color="k",
            label="Centerline"
        )

        # E and B field arrows
        for idx in range(0, num, step):
            p = centerline[idx]

            # Electric field arrow (blue)
            e = e_vectors[idx]
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

            # Magnetic field arrow (orange), offset slightly along t to avoid overlap
            b_vec = b_vectors[idx]
            b_norm = np.linalg.norm(b_vec)
            if b_norm > 1e-12:
                b_normalized = b_vec / b_norm
                t_vec = t[idx]
                p_b = p + 0.4 * arrow_len * t_vec
                q_b = p_b + arrow_len * b_normalized
                ax.plot(
                    [p_b[0], q_b[0]],
                    [p_b[1], q_b[1]],
                    [p_b[2], q_b[2]],
                    linewidth=1.7,
                    color="C1",  # orange
                    alpha=0.8
                )

        # Set limits with margin to ensure torus hull fits
        ax.set_xlim(-extent_margin, extent_margin)
        ax.set_ylim(-extent_margin, extent_margin)
        ax.set_zlim(-extent_margin, extent_margin)

        # Set equal aspect ratio
        ax.set_box_aspect([1, 1, 1])

        ax.view_init(elev=view["elev"], azim=view["azim"])
        ax.set_title(view["title"], fontsize=14, fontweight='bold')
        ax.set_xlabel("X¹ (brane)", fontsize=12)
        ax.set_ylabel("X² (brane)", fontsize=12)
        ax.set_zlabel("X³ (brane)", fontsize=12)

        # Add legend
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout()

        # Save to plots directory
        output_path = run_manager.get_plot_path(view["filename"])
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  ✓ Saved: {view['filename']}")
        plt.close(fig)

    # Save configuration
    config = {
        "experiment": "Circularly Polarized Photon - Toroidal Path",
        "torus_major_radius": torus_params.major_radius,
        "torus_minor_radius": torus_params.minor_radius,
        "core_windings": torus_params.core_windings,
        "tube_windings": torus_params.tube_windings,
        "photon_sigma_n": photon_params.sigma_n,
        "photon_sigma_b": photon_params.sigma_b,
        "total_phase": f"{photon_params.total_phase / np.pi:.1f}π",
        "wave_cycles": photon_params.wave_cycles,
        "num_samples": num_samples,
    }
    run_manager.save_config(config)

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    visualize_circular_photon_toroidal()