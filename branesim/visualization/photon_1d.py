"""
Visualization tools for 1D photon simulations.

Provides reusable plotting functions for amplitude, lateral distortion, velocity,
and lateralization analysis in 1D photon experiments.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from branesim.experiments.common.photon_1d_runner import Photon1DRunData


def _plot_1d_snapshots_generic(
    run,
    data_by_t: dict[float, np.ndarray],
    ylabel: str,
    title: str,
    filename: str,
    ylim_factor: float = 1.5,
    add_zero_line: bool = False,
    add_boundary_markers: bool = True,
) -> None:
    """
    Generic function for plotting 1D snapshots at multiple times.

    Parameters
    ----------
    run : Photon1DRunData
        Run data containing coordinates and run_manager
    data_by_t : dict[float, np.ndarray]
        Dictionary mapping time (fs) to data array [N] in desired units
    ylabel : str
        Y-axis label with units
    title : str
        Plot title
    filename : str
        Output filename
    ylim_factor : float, optional
        Y-limit scaling factor relative to max absolute value
    add_zero_line : bool, optional
        Whether to add a dashed line at y=0
    add_boundary_markers : bool, optional
        Whether to add red markers at boundaries
    """
    x_nm = run.x_coords_phys_m * 1e9  # Convert m → nm
    times_fs = sorted(data_by_t.keys())
    n = len(times_fs)

    fig, axes = plt.subplots(n, 1, figsize=(14, max(4, 3 * n)))
    if n == 1:
        axes = [axes]

    fig.suptitle(title, fontsize=16, fontweight='bold')

    # Compute global y-limits for consistent scaling
    max_val = max([np.abs(data_by_t[t]).max() for t in times_fs])
    ylim = ylim_factor * max_val

    for ax, t_fs in zip(axes, times_fs):
        data = data_by_t[t_fs]

        # Main plot
        ax.plot(x_nm, data, linewidth=2)

        # Optional zero line
        if add_zero_line:
            ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.5)

        # Optional boundary markers
        if add_boundary_markers:
            ax.plot([x_nm[0], x_nm[-1]], [0, 0], 'ro',
                    markersize=8, label='Fixed boundaries' if ax == axes[0] else '')

        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlim(x_nm[0], x_nm[-1])
        ax.set_ylim(-ylim, ylim)
        ax.grid(True, alpha=0.3)

        # Time label
        ax.text(0.02, 0.95, f't = {t_fs:.3f} fs',
                transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        if ax == axes[0] and add_boundary_markers:
            ax.legend(loc='upper right', fontsize=10)

    axes[-1].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run.run_manager.get_plot_path(filename), dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_photon_1d_amplitude_propagation(run: "Photon1DRunData") -> None:
    """
    Plot amplitude ξ (X⁴ displacement) snapshots in nanometers.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with snapshots and coordinates
    """
    # Convert sim → phys → nm
    data_by_t_nm = {}
    for t_phys_s in run.snapshot_times_phys_s:
        t_fs = t_phys_s * 1e15
        xi_sim = run.snapshots_xi[t_phys_s]
        xi_phys_m = run.mapper.to_phys_length(xi_sim)
        data_by_t_nm[t_fs] = xi_phys_m * 1e9  # m → nm

    _plot_1d_snapshots_generic(
        run,
        data_by_t_nm,
        ylabel='ξ [nm]',
        title=f'1D Photon - Amplitude Propagation (c = {run.constants.c:.3e} m/s)',
        filename='photon_1d_example_propagation.png',
        ylim_factor=1.5,
        add_zero_line=False,
        add_boundary_markers=True,
    )


def plot_photon_1d_lateral_distortion(run: "Photon1DRunData") -> None:
    """
    Plot lateral distortion Δx (left/right movement) in picometers.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with snapshots and coordinates
    """
    # Convert sim → phys → pm
    data_by_t_pm = {}
    for t_phys_s in run.snapshot_times_phys_s:
        t_fs = t_phys_s * 1e15
        delta_x_sim = run.snapshots_delta_x[t_phys_s]
        delta_x_phys_m = run.mapper.to_phys_length(delta_x_sim)
        data_by_t_pm[t_fs] = delta_x_phys_m * 1e12  # m → pm

    _plot_1d_snapshots_generic(
        run,
        data_by_t_pm,
        ylabel='Δx [pm]',
        title='1D Photon - Lateral Distortion (Left/Right Movement)',
        filename='photon_1d_example_lateral_distortion.png',
        ylim_factor=1.5,
        add_zero_line=True,
        add_boundary_markers=True,
    )


def plot_photon_1d_amplitude_velocity(run: "Photon1DRunData") -> None:
    """
    Plot amplitude velocity v_ξ (∂ξ/∂t) in m/s.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with snapshots and coordinates
    """
    # Convert sim → phys (m/s)
    data_by_t_ms = {}
    for t_phys_s in run.snapshot_times_phys_s:
        t_fs = t_phys_s * 1e15
        v_xi_sim = run.snapshots_v_xi[t_phys_s]
        v_xi_phys_ms = run.mapper.to_phys_velocity(v_xi_sim)
        data_by_t_ms[t_fs] = v_xi_phys_ms

    _plot_1d_snapshots_generic(
        run,
        data_by_t_ms,
        ylabel='v_ξ [m/s]',
        title='1D Photon - Amplitude Velocity (v_ξ = ∂ξ/∂t)',
        filename='photon_1d_example_amplitude_velocity.png',
        ylim_factor=1.5,
        add_zero_line=True,
        add_boundary_markers=True,
    )


def plot_photon_1d_lateral_velocity(run: "Photon1DRunData") -> None:
    """
    Plot lateral velocity v_x (∂x/∂t) in m/s.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with snapshots and coordinates
    """
    # Convert sim → phys (m/s)
    data_by_t_ms = {}
    for t_phys_s in run.snapshot_times_phys_s:
        t_fs = t_phys_s * 1e15
        v_x_sim = run.snapshots_v_x[t_phys_s]
        v_x_phys_ms = run.mapper.to_phys_velocity(v_x_sim)
        data_by_t_ms[t_fs] = v_x_phys_ms

    _plot_1d_snapshots_generic(
        run,
        data_by_t_ms,
        ylabel='v_x [m/s]',
        title='1D Photon - Lateral Velocity (v_x = ∂x/∂t)',
        filename='photon_1d_example_lateral_velocity.png',
        ylim_factor=1.5,
        add_zero_line=True,
        add_boundary_markers=True,
    )


def plot_photon_1d_lateralization_snapshots(run: "Photon1DRunData") -> None:
    """
    Plot lateralization ratio R_lat (dimensionless, 0 to 1) at snapshot times.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with snapshots and coordinates
    """
    x_nm = run.x_coords_phys_m * 1e9
    # Use original times in seconds to avoid floating-point precision issues
    times_phys_s_sorted = sorted(run.snapshot_times_phys_s)
    n = len(times_phys_s_sorted)

    fig, axes = plt.subplots(n, 1, figsize=(14, max(4, 3 * n)))
    if n == 1:
        axes = [axes]

    fig.suptitle('1D Photon - Lateralization Ratio (R_lat = E_lat / E_total)',
                 fontsize=16, fontweight='bold')

    for ax, t_phys_s in zip(axes, times_phys_s_sorted):
        t_fs = t_phys_s * 1e15
        R_lat = run.snapshots_R_lat[t_phys_s]

        ax.plot(x_nm, R_lat, 'orange', linewidth=2, label='R_lat')
        ax.axhline(y=0.5, color='k', linestyle='--', linewidth=0.5, alpha=0.5,
                   label='R_lat=0.5' if ax == axes[0] else '')
        ax.set_ylabel('R_lat', fontsize=11)
        ax.set_xlim(x_nm[0], x_nm[-1])
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)

        ax.text(0.02, 0.95, f't = {t_fs:.3f} fs',
                transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

        if ax == axes[0]:
            ax.legend(loc='upper right', fontsize=10)

    axes[-1].set_xlabel('Position [nm]', fontsize=12)

    plt.tight_layout()
    plt.savefig(run.run_manager.get_plot_path('photon_1d_example_lateralization.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_photon_1d_lateralization_global(run: "Photon1DRunData") -> None:
    """
    Plot global lateralization ratio vs time.

    Parameters
    ----------
    run : Photon1DRunData
        Run data with tracking history
    """
    times_fs = np.array(run.times_phys_track_s) * 1e15
    R_lat_global = np.array(run.R_lat_global_track)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.plot(times_fs, R_lat_global, 'orange', linewidth=2, label='Global R_lat')
    ax.axhline(y=0.5, color='k', linestyle='--', linewidth=1, alpha=0.5, label='R_lat=0.5')
    ax.set_xlabel('Time [fs]', fontsize=12)
    ax.set_ylabel('Global R_lat', fontsize=12)
    ax.set_title('Global Lateralization Ratio vs Time', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(run.run_manager.get_plot_path('photon_1d_example_lateralization_global.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_photon_1d_tracking_analysis(run: "Photon1DRunData") -> None:
    """
    Plot wave center tracking and energy conservation.

    Creates a 2-panel figure showing:
    1. Wave center position vs time
    2. Normalized energy vs time

    Parameters
    ----------
    run : Photon1DRunData
        Run data with tracking history
    """
    times_fs = np.array(run.times_phys_track_s) * 1e15
    centers_phys_m = run.mapper.to_phys_length(np.array(run.centers_sim_track))
    centers_nm = centers_phys_m * 1e9
    energies = np.array(run.energies_track_J)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Position vs time (in nm and fs)
    axes[0].plot(times_fs, centers_nm, 'b-', linewidth=2, label='Wave center')
    axes[0].set_xlabel('Time [fs]', fontsize=12)
    axes[0].set_ylabel('Wave Center [nm]', fontsize=12)
    axes[0].set_title('Wave Propagation at Speed of Light', fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=10)

    # Energy conservation
    initial_energy = energies[0]
    axes[1].plot(times_fs, energies / initial_energy, 'g-', linewidth=2)
    axes[1].axhline(y=1.0, color='r', linestyle='--', linewidth=1, alpha=0.5)
    axes[1].set_xlabel('Time [fs]', fontsize=12)
    axes[1].set_ylabel('E(t) / E(0)', fontsize=12)
    axes[1].set_title('Energy Conservation', fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(run.run_manager.get_plot_path('photon_1d_example_analysis.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_all_photon_1d_standard(run: "Photon1DRunData") -> None:
    """
    Convenience function to generate all standard 1D photon plots.

    Generates:
    - Amplitude propagation (ξ vs x at different times)
    - Lateral distortion (Δx vs x)
    - Amplitude velocity (v_ξ vs x)
    - Lateral velocity (v_x vs x)
    - Lateralization snapshots (R_lat vs x)
    - Global lateralization (R_lat vs t)
    - Tracking analysis (center and energy vs t)

    Parameters
    ----------
    run : Photon1DRunData
        Run data containing all necessary information

    Example
    -------
    >>> from branesim.experiments.common.photon_1d_runner import run_photon_1d, Photon1DConfig
    >>> from branesim.visualization.photon_1d import plot_all_photon_1d_standard
    >>>
    >>> cfg = Photon1DConfig()
    >>> run = run_photon_1d(cfg)
    >>> plot_all_photon_1d_standard(run)
    """
    print("\nCreating plots...")
    plot_photon_1d_amplitude_propagation(run)
    print("  ✓ Saved: photon_1d_example_propagation.png")

    plot_photon_1d_lateral_distortion(run)
    print("  ✓ Saved: photon_1d_example_lateral_distortion.png")

    plot_photon_1d_amplitude_velocity(run)
    print("  ✓ Saved: photon_1d_example_amplitude_velocity.png")

    plot_photon_1d_lateral_velocity(run)
    print("  ✓ Saved: photon_1d_example_lateral_velocity.png")

    plot_photon_1d_lateralization_snapshots(run)
    print("  ✓ Saved: photon_1d_example_lateralization.png")

    plot_photon_1d_lateralization_global(run)
    print("  ✓ Saved: photon_1d_example_lateralization_global.png")

    plot_photon_1d_tracking_analysis(run)
    print("  ✓ Saved: photon_1d_example_analysis.png")