"""
Visualization tools for Berry phase in 1D systems.

Provides plotting functions for Berry phase profiles and Berry connection along
a 1D spatial coordinate.
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt


def plot_berry_phase_profiles(
    run_manager,
    x_nm: np.ndarray,
    times_fs: list[float],
    gamma_wrapped_by_t: dict[float, np.ndarray],
    title: str = "Berry phase profile γ(x) (wrapped)",
    filename: str = "berry_phase_profile.png",
):
    """
    Plot Berry phase profiles γ(x) at multiple time snapshots.

    Creates a multi-panel figure showing the wrapped Berry phase [-π, π]
    as a function of position at different times.

    Parameters
    ----------
    run_manager : TestRunManager
        Run manager for saving plots
    x_nm : np.ndarray
        Position coordinates in nanometers, shape [N]
    times_fs : list[float]
        List of snapshot times in femtoseconds
    gamma_wrapped_by_t : dict[float, np.ndarray]
        Dictionary mapping time (fs) to Berry phase profile (wrapped to [-π, π]), shape [N]
    title : str, optional
        Figure title
    filename : str, optional
        Output filename (saved to run_manager's plot directory)

    Examples
    --------
    >>> gamma_by_t = {}
    >>> for step in snapshot_steps:
    ...     # ... compute Berry phase ...
    ...     gamma_by_t[t_fs] = result['gamma_wrapped'].cpu().numpy()
    >>>
    >>> plot_berry_phase_profiles(
    ...     run_manager, x_nm, times_fs, gamma_by_t,
    ...     filename="berry_phase_profiles.png"
    ... )
    """
    n = len(times_fs)
    fig, axes = plt.subplots(n, 1, figsize=(14, max(4, 3 * n)))
    if n == 1:
        axes = [axes]

    fig.suptitle(title, fontsize=16, fontweight="bold")

    for ax, t_fs in zip(axes, times_fs):
        y = gamma_wrapped_by_t[t_fs]
        ax.plot(x_nm, y, linewidth=2)
        ax.set_ylabel("γ [rad]")
        ax.grid(True, alpha=0.3)
        ax.text(
            0.02,
            0.92,
            f"t = {t_fs:.3f} fs",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

    axes[-1].set_xlabel("Position [nm]")

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path(filename), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_berry_connection_profiles(
    run_manager,
    x_nm_edges: np.ndarray,
    times_fs: list[float],
    A_x_by_t: dict[float, np.ndarray],
    title: str = "Berry connection A_x(x)",
    filename: str = "berry_connection_profile.png",
):
    """
    Plot Berry connection A_x(x) at multiple time snapshots.

    Creates a multi-panel figure showing the Berry connection (gauge field)
    as a function of position at different times.

    The Berry connection is defined on edges between lattice points, so
    x_nm_edges should contain the midpoint positions.

    Parameters
    ----------
    run_manager : TestRunManager
        Run manager for saving plots
    x_nm_edges : np.ndarray
        Edge position coordinates in nanometers, shape [N-1]
        Typically computed as 0.5 * (x_nm[:-1] + x_nm[1:])
    times_fs : list[float]
        List of snapshot times in femtoseconds
    A_x_by_t : dict[float, np.ndarray]
        Dictionary mapping time (fs) to Berry connection profile, shape [N-1]
        Units: [rad / sim-length] (typically sim-length = 1 grid unit)
    title : str, optional
        Figure title
    filename : str, optional
        Output filename (saved to run_manager's plot directory)

    Notes
    -----
    The Berry connection A_x is related to the phase gradient and represents
    the "gauge field" associated with the adiabatic evolution of the state.

    In the continuum limit: A_x = Im⟨ψ|∂_x ψ⟩

    Examples
    --------
    >>> A_x_by_t = {}
    >>> for step in snapshot_steps:
    ...     # ... compute Berry phase ...
    ...     A_x_by_t[t_fs] = result['A_x'].cpu().numpy()
    >>>
    >>> x_edges_nm = 0.5 * (x_nm[:-1] + x_nm[1:])
    >>> plot_berry_connection_profiles(
    ...     run_manager, x_edges_nm, times_fs, A_x_by_t,
    ...     filename="berry_connection_profiles.png"
    ... )
    """
    n = len(times_fs)
    fig, axes = plt.subplots(n, 1, figsize=(14, max(4, 3 * n)))
    if n == 1:
        axes = [axes]

    fig.suptitle(title, fontsize=16, fontweight="bold")

    for ax, t_fs in zip(axes, times_fs):
        y = A_x_by_t[t_fs]
        ax.plot(x_nm_edges, y, linewidth=2)
        ax.set_ylabel("A_x [rad / sim-length]")
        ax.grid(True, alpha=0.3)
        ax.text(
            0.02,
            0.92,
            f"t = {t_fs:.3f} fs",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

    axes[-1].set_xlabel("Position [nm]")

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path(filename), dpi=150, bbox_inches="tight")
    plt.close(fig)