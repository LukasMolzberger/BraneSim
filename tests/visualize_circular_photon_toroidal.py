import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt
import numpy as np

from branesim.initialization.geometry import (
    TorusKnotParameters,
    sample_torus_knot_centerline,
    compute_frenet_frames,
)
from branesim.initialization.geometry import (
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

    # Define views with filenames for separate file output
    views = [
        dict(elev=25, azim=-60, title="Photon mode - oblique view", filename="tubular_photon_oblique.png"),
        dict(elev=10, azim=30, title="Photon mode - nearly in-plane", filename="tubular_photon_side.png"),
        dict(elev=80, azim=-90, title="Photon mode - top view", filename="tubular_photon_top.png"),
    ]

    # Calculate arrow scale to match previous visualization
    # Previous: arrow_len = 1.8 * 0.12 = 0.216
    # Generic function: arrow_len = arrow_scale * max_extent * 0.05
    # Torus extent ≈ 2*(major_radius + minor_radius) ≈ 2.72
    # So arrow_scale ≈ 0.216 / (2.72 * 0.05) ≈ 1.59
    arrow_scale = 1.59

    # Create background draw function to add torus hull
    def draw_torus_background(ax):
        """Background drawing function for torus hull."""
        _plot_torus_hull(
            ax,
            R=torus_params.major_radius,
            r=torus_params.minor_radius + 0.04,
        )

        # Set axis limits to ensure torus hull fits
        # Compute extent based on torus outer radius with margin
        torus_outer_radius = torus_params.major_radius + torus_params.minor_radius + 0.04
        extent_margin = torus_outer_radius * 1.1

        ax.set_xlim(-extent_margin, extent_margin)
        ax.set_ylim(-extent_margin, extent_margin)
        ax.set_zlim(-extent_margin, extent_margin)

        # Set equal aspect ratio
        ax.set_box_aspect([1, 1, 1])

        # Optionally draw photon cloud
        if DRAW_CLOUD:
            # Longitudinal modulation for intensity
            longitudinal_modulation = np.abs(a_long)
            field_points, amplitudes = sample_gaussian_envelope(
                centerline,
                n,
                b,
                photon_params,
                longitudinal_modulation=longitudinal_modulation,
            )

            # Normalize amplitudes for color/opacity mapping
            if amplitudes.size > 0:
                amp_min = float(np.min(amplitudes))
                amp_max = float(np.max(amplitudes))
                denom = max(amp_max - amp_min, 1e-12)
                amp_norm = (amplitudes - amp_min) / denom
                cmap = plt.get_cmap()
                colors = cmap(amp_norm)
                gamma = 1.5
                colors[:, 3] = amp_norm ** gamma

                ax.scatter(
                    field_points[:, 0],
                    field_points[:, 1],
                    field_points[:, 2],
                    s=4,
                    c=colors,
                    marker="o",
                )

    # Use generic EM field visualization with custom background
    visualize_em_fields_along_centerline(
        centerline=centerline,
        E_field=e_vectors,
        B_field=b_vectors,
        output_path=run_manager.plots_dir,  # Directory path for separate files
        title="EM Fields: Circularly Polarized Photon (Toroidal Path)",
        views=views,
        arrow_scale=arrow_scale,
        subsample_step=None,  # Auto-compute (matches previous: max(1, num // 40))
        figsize=(8, 8),  # Used for separate files
        dpi=150,
        background_draw_func=draw_torus_background,
        separate_files=True
    )

    # Print saved files
    for view in views:
        print(f"  ✓ Saved: {view['filename']}")

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