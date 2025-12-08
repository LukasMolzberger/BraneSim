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
        alpha=0.06,   # slightly thicker but still subtle
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


def visualize_tubular_photon_mode():
    """
    Visual test: circularly polarized photon mode along the toroidal path.

    - Uses the torus-knot centerline and Frenet frame from the geometry module.
    - Builds a circularly polarized E-field in the local (n, b) plane, with
      total internal phase advance 4π along the loop (spinor-like orientation).
    - Introduces a higher-frequency longitudinal amplitude with 'wave_cycles'
      cos oscillations along the loop (more crests).
    - Constructs a B-field via B ∝ t × E.
    - Samples a 2D Gaussian envelope in the transverse directions and maps
      amplitude to transparency: weak field = transparent, strong = opaque.

    The Möbius / double-loop geometry is present only as the path; the strip
    itself is not drawn here.
    """
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

    # Longitudinal modulation for intensity: use |A_long|
    longitudinal_modulation = np.abs(a_long)

    # Gaussian envelope in transverse (n, b) directions, modulated along the path
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
    else:
        amp_norm = amplitudes

    # Map amplitude to RGBA colors with amplitude-dependent alpha
    cmap = plt.get_cmap()  # fixes MatplotlibDeprecationWarning
    colors = cmap(amp_norm)  # (N, 4)
    gamma = 1.5             # controls how fast opacity rises with amplitude
    colors[:, 3] = amp_norm ** gamma

    all_x = np.concatenate([field_points[:, 0], centerline[:, 0]])
    all_y = np.concatenate([field_points[:, 1], centerline[:, 1]])
    all_z = np.concatenate([field_points[:, 2], centerline[:, 2]])

    fig = plt.figure(figsize=(15, 5))

    views = [
        dict(elev=25, azim=-60, title="Photon mode - oblique view"),
        dict(elev=10, azim=30, title="Photon mode - nearly in-plane"),
        dict(elev=80, azim=-90, title="Photon mode - top view"),
    ]

    # Subsample for drawing E and B vectors
    num = centerline.shape[0]
    step = max(1, num // 40)  # slightly denser arrows
    arrow_len = 1.8 * min(photon_params.sigma_n, photon_params.sigma_b)

    for i, view in enumerate(views, start=1):
        ax = fig.add_subplot(1, 3, i, projection="3d")

        # Outer torus hull for context
        _plot_torus_hull(
            ax,
            R=torus_params.major_radius,
            r=torus_params.minor_radius + 0.04,
        )

        # Photon "cloud" with amplitude-dependent transparency
        ax.scatter(
            field_points[:, 0],
            field_points[:, 1],
            field_points[:, 2],
            s=5,
            c=colors,
            marker="o",
        )

        # Centerline (path of the photon)
        ax.plot(
            centerline[:, 0],
            centerline[:, 1],
            centerline[:, 2],
            linewidth=1.7,
            color="k",
        )

        # Polarization vectors (E and B) at a subset of points
        for idx in range(0, num, step):
            p = centerline[idx]

            # Electric field arrow (blue)
            e = e_vectors[idx]
            q_e = p + arrow_len * e
            ax.plot(
                [p[0], q_e[0]],
                [p[1], q_e[1]],
                [p[2], q_e[2]],
                linewidth=1.7,
                color="C0",  # blue
            )

            # Magnetic field arrow (orange), offset slightly along t to avoid overlap
            b_vec = b_vectors[idx]
            t_vec = t[idx]
            p_b = p + 0.4 * arrow_len * t_vec
            q_b = p_b + arrow_len * b_vec
            ax.plot(
                [p_b[0], q_b[0]],
                [p_b[1], q_b[1]],
                [p_b[2], q_b[2]],
                linewidth=1.7,
                color="C1",  # orange
            )

        _set_equal_aspect_3d(ax, all_x, all_y, all_z)

        ax.view_init(elev=view["elev"], azim=view["azim"])
        ax.set_title(view["title"])
        ax.set_xlabel("X¹ (brane)")
        ax.set_ylabel("X² (brane)")
        ax.set_zlabel("X³ (brane)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    visualize_tubular_photon_mode()