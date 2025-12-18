"""
CSV export utilities for Berry phase data.

Provides functions to export Berry phase profiles and Berry connection data
to CSV format for analysis in external tools.
"""

from __future__ import annotations
import csv
import numpy as np


def export_berry_phase_csv(
    path: str,
    x_m: np.ndarray,
    gamma_rad: np.ndarray,
    A_x: np.ndarray,
    amp_m: np.ndarray | None = None,
):
    """
    Export per-point Berry phase profile to CSV.

    Writes two CSV files:
    1. Main file (at `path`): per-point Berry phase, position, and optional amplitude
    2. Edge file (path with "_edges" suffix): per-edge Berry connection

    Parameters
    ----------
    path : str
        Output path for main CSV file (e.g., "berry_phase_t_0.000fs.csv")
    x_m : np.ndarray
        Position coordinates in meters, shape [N]
    gamma_rad : np.ndarray
        Berry phase in radians, shape [N]
    A_x : np.ndarray
        Berry connection in [rad / sim-length], shape [N-1]
        Defined on edges between lattice points
    amp_m : np.ndarray, optional
        Amplitude in meters, shape [N]
        If provided, included as an additional column

    Files Created
    -------------
    {path}:
        CSV with columns: point_idx, x_position [m], gamma [rad], amplitude [m] (optional)

    {path with _edges suffix}:
        CSV with columns: edge_idx, x_edge_position [m], A_x [rad / sim-length]

    Notes
    -----
    The Berry connection A_x is defined on edges (between points i and i+1),
    so the edge position is computed as the midpoint: x_edge = (x[i] + x[i+1])/2

    Examples
    --------
    >>> from branesim.io import export_berry_phase_csv
    >>>
    >>> # After computing Berry phase...
    >>> export_berry_phase_csv(
    ...     run_manager.get_data_path("berry_phase_t_0.000fs.csv"),
    ...     x_coords_phys,  # [m]
    ...     gamma_wrapped,  # [rad]
    ...     A_x,  # [rad / sim-length]
    ...     amp_m=amplitude_phys,  # [m], optional
    ... )
    """
    # Write main file (per-point data)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)

        # Header
        header = ["point_idx", "x_position [m]", "gamma [rad]"]
        if amp_m is not None:
            header.append("amplitude [m]")
        w.writerow(header)

        # Data rows
        N = len(x_m)
        for i in range(N):
            row = [i, float(x_m[i]), float(gamma_rad[i])]
            if amp_m is not None:
                row.append(float(amp_m[i]))
            w.writerow(row)

    # Write edge file (per-edge Berry connection)
    edge_path = path.replace(".csv", "_edges.csv")
    with open(edge_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["edge_idx", "x_edge_position [m]", "A_x [rad / sim-length]"])

        for i in range(len(A_x)):
            # Edge position is midpoint between neighboring points
            x_edge = 0.5 * (x_m[i] + x_m[i + 1])
            w.writerow([i, float(x_edge), float(A_x[i])])