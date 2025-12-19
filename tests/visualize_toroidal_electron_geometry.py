import matplotlib.pyplot as plt
import numpy as np

from branesim.initialization.geometry import (
    TorusKnotParameters,
    sample_torus_knot_centerline,
    compute_frenet_frames,
    construct_twisted_strip,
)


def _plot_torus_hull(ax, R: float, r: float, num_phi: int = 64, num_theta: int = 32):
    """Draw a semi-transparent torus hull for visual reference."""
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
        alpha=0.15,
        edgecolor="k",
    )


def _set_equal_aspect_3d(ax, X, Y, Z):
    """Set equal aspect ratio for 3D axes based on data ranges."""
    x_min, x_max = np.min(X), np.max(X)
    y_min, y_max = np.min(Y), np.max(Y)
    z_min, z_max = np.min(Z), np.max(Z)

    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)
    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    ax.set_xlim(x_mid - max_range / 2.0, x_mid + max_range / 2.0)
    ax.set_ylim(y_mid - max_range / 2.0, y_mid + max_range / 2.0)
    ax.set_zlim(z_mid - max_range / 2.0, z_mid + max_range / 2.0)


def visualize_toroidal_electron_geometry():
    """
    Visual test: double-loop strip inside a torus.

    This script constructs
      1. A (2, 1) torus-knot centerline C(z) inside a torus of radii (R, r0).
      2. The associated Frenet–Serret frame (t, n, b).
      3. A narrow strip given by r(z, w) = C(z) + w * n(z).

    It then renders three different 3D views to verify that the strip
    winds twice around the torus core as intended.
    """
    params = TorusKnotParameters(
        major_radius=1.0,
        minor_radius=0.3,
        core_windings=2,
        tube_windings=1,
    )

    centerline = sample_torus_knot_centerline(params, num_samples=1200)
    t, n, b = compute_frenet_frames(centerline)

    strip_X, strip_Y, strip_Z = construct_twisted_strip(
        centerline,
        n,
        b,
        strip_half_width=0.15,
        num_width_samples=32,
    )

    # For aspect-ratio calculation, gather all coordinates.
    all_X = np.concatenate([strip_X.flatten(), centerline[:, 0]])
    all_Y = np.concatenate([strip_Y.flatten(), centerline[:, 1]])
    all_Z = np.concatenate([strip_Z.flatten(), centerline[:, 2]])

    fig = plt.figure(figsize=(15, 5))

    views = [
        dict(elev=25, azim=-60, title="Oblique view"),
        dict(elev=10, azim=30, title="Nearly in-plane view"),
        dict(elev=80, azim=-90, title="Top view"),
    ]

    for i, view in enumerate(views, start=1):
        ax = fig.add_subplot(1, 3, i, projection="3d")

        # Torus hull (slightly thicker than the strip radius for clarity)
        _plot_torus_hull(ax, R=params.major_radius, r=params.minor_radius + 0.05)

        # Centerline
        ax.plot(
            centerline[:, 0],
            centerline[:, 1],
            centerline[:, 2],
            linewidth=2.0,
        )

        # Strip surface
        ax.plot_surface(
            strip_X,
            strip_Y,
            strip_Z,
            rstride=4,
            cstride=1,
            linewidth=0.0,
            alpha=0.8,
        )

        _set_equal_aspect_3d(ax, all_X, all_Y, all_Z)

        ax.view_init(elev=view["elev"], azim=view["azim"])
        ax.set_title(view["title"])
        ax.set_xlabel("X¹ (brane)")
        ax.set_ylabel("X² (brane)")
        ax.set_zlabel("X³ (brane)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    visualize_toroidal_electron_geometry()