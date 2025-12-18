"""
CSV export for 1D photon simulation snapshots.

Exports detailed per-point data including positions, velocities, accelerations,
forces, energies, and lateralization ratios in SI units.
"""

from __future__ import annotations
import csv
import numpy as np
import torch


def export_photon_1d_snapshot_csv(
    filename: str,
    state,
    grid,
    initial_positions: torch.Tensor,
    spring_constant: float,
    lateralization,
    physics,
    h_sim: float,
    mapper,
    rest_length_sim: float = 0.0,
    lateral_dim_idx: int = 0,
    amplitude_dim_idx: int = 3,
) -> None:
    """
    Export detailed CSV snapshot with all brane point data in SI units.

    All quantities are converted from simulation units to physical SI units before export.

    Parameters
    ----------
    filename : str
        Output CSV file path
    state : BraneState
        Current simulation state
    grid : BraneGrid
        Grid topology
    initial_positions : torch.Tensor
        Reference positions for computing displacements, shape [N, 4]
    spring_constant : float
        Spring constant in simulation units
    lateralization : LateralizationMeasurement
        Lateralization measurement object
    physics : SpringForceComputer
        Physics force computer
    h_sim : float
        Grid spacing in simulation units
    mapper : DimensionalMapper
        Unit conversion mapper
    rest_length_sim : float, optional
        Rest length in simulation units (default: 0.0)
    lateral_dim_idx : int, optional
        Index of lateral dimension (default: 0 for x)
    amplitude_dim_idx : int, optional
        Index of amplitude dimension (default: 3 for ξ = X⁴)

    Output Columns
    --------------
    - point_idx: Point index
    - x_position [m]: Lateral position
    - xi_position [m]: Amplitude displacement
    - x_velocity [m/s]: Lateral velocity
    - xi_velocity [m/s]: Amplitude velocity
    - x_acceleration [m/s²]: Lateral acceleration
    - xi_acceleration [m/s²]: Amplitude acceleration
    - delta_x [m]: Lateral displacement from initial
    - delta_xi [m]: Amplitude displacement from initial
    - h_spacing [m]: Grid spacing
    - L_to_next [m]: Current spring length to next neighbor
    - delta_L_to_next [m]: Spring extension beyond reference
    - F_left_x [N], F_left_xi [N]: Force from left neighbor
    - F_right_x [N], F_right_xi [N]: Force from right neighbor
    - E_amp_kin [J], E_amp_pot [J]: Amplitude kinetic and potential energy
    - E_lat_kin [J], E_lat_pot [J]: Lateral kinetic and potential energy
    - R_lat: Lateralization ratio (dimensionless)

    Notes
    -----
    All output is in SI units (meters, seconds, Newtons, Joules).
    The mapper is used to convert from simulation units to physical units.
    """
    nx = state.positions.shape[0]
    neighbors = grid.neighbors

    # Get lateralization measurement
    R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)

    # Convert to numpy (still in sim units)
    positions_sim = state.positions.cpu().numpy()
    velocities_sim = state.velocities.cpu().numpy()
    accelerations_sim = state.accelerations.cpu().numpy()
    initial_pos_sim = initial_positions.cpu().numpy()
    R_lat = R_lat_local.cpu().numpy()
    E_amp_kin = diagnostics['E_amp_kin'].cpu().numpy()  # Already in J
    E_amp_pot = diagnostics['E_amp_pot'].cpu().numpy()  # Already in J
    E_lat_kin = diagnostics['E_lat_kin'].cpu().numpy()  # Already in J
    E_lat_pot = diagnostics['E_lat_pot'].cpu().numpy()  # Already in J

    # Compute spring forces for each point (in sim units)
    k_sim = spring_constant

    # Store forces from left and right neighbors (in sim units)
    F_left_sim = np.zeros((nx, 4))
    F_right_sim = np.zeros((nx, 4))

    for i in range(nx):
        for j_idx in neighbors[i]:
            j = j_idx.item() if isinstance(j_idx, torch.Tensor) else int(j_idx)
            if j < 0:
                continue

            # Vector from i to j (sim units)
            d = positions_sim[j] - positions_sim[i]
            length = np.sqrt(np.sum(d**2))

            if length < 1e-20:
                continue

            # Spring force magnitude: F = k * (L - L0) (sim units)
            extension = length - rest_length_sim
            force_mag = k_sim * extension

            # Force direction (unit vector)
            direction = d / length

            # Force vector on point i from point j (sim units)
            force = force_mag * direction

            # Determine if this is left or right neighbor
            if j == i - 1:
                F_left_sim[i] = force
            elif j == i + 1:
                F_right_sim[i] = force

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header with units
        writer.writerow([
            'point_idx',
            'x_position [m]', 'xi_position [m]',
            'x_velocity [m/s]', 'xi_velocity [m/s]',
            'x_acceleration [m/s^2]', 'xi_acceleration [m/s^2]',
            'delta_x [m]', 'delta_xi [m]',
            'h_spacing [m]',
            'L_to_next [m]', 'delta_L_to_next [m]',
            'F_left_x [N]', 'F_left_xi [N]', 'F_right_x [N]', 'F_right_xi [N]',
            'E_amp_kin [J]', 'E_amp_pot [J]', 'E_lat_kin [J]', 'E_lat_pot [J]', 'R_lat'
        ])

        # Write data for each point (convert to physical units using mapper)
        for i in range(nx):
            # Positions [sim → m]
            x_pos = mapper.to_phys_length(positions_sim[i, lateral_dim_idx])
            xi_pos = mapper.to_phys_length(positions_sim[i, amplitude_dim_idx])

            # Velocities [sim → m/s]
            x_vel = mapper.to_phys_velocity(velocities_sim[i, lateral_dim_idx])
            xi_vel = mapper.to_phys_velocity(velocities_sim[i, amplitude_dim_idx])

            # Accelerations [sim → m/s²]
            x_acc = mapper.to_phys_acceleration(accelerations_sim[i, lateral_dim_idx])
            xi_acc = mapper.to_phys_acceleration(accelerations_sim[i, amplitude_dim_idx])

            # Displacements from initial [sim → m]
            delta_x = mapper.to_phys_length(positions_sim[i, lateral_dim_idx] - initial_pos_sim[i, lateral_dim_idx])
            xi = mapper.to_phys_length(positions_sim[i, amplitude_dim_idx] - initial_pos_sim[i, amplitude_dim_idx])

            # Spring properties (to next neighbor if exists)
            L_to_next_phys = 0.0
            delta_L_to_next_phys = 0.0

            if i < nx - 1:  # Not the last point
                # Find next neighbor
                for j_idx in neighbors[i]:
                    j = j_idx.item() if isinstance(j_idx, torch.Tensor) else int(j_idx)
                    if j == i + 1:  # Next neighbor
                        # Current spring vector (sim units)
                        dX_sim = positions_sim[j] - positions_sim[i]
                        L_to_next_sim = np.sqrt(np.sum(dX_sim**2))

                        # Reference spring vector (sim units)
                        dX0_sim = initial_pos_sim[j] - initial_pos_sim[i]
                        L_ref_sim = np.sqrt(np.sum(dX0_sim**2))

                        # Extension beyond reference (sim units)
                        delta_L_to_next_sim = L_to_next_sim - L_ref_sim

                        # Convert to physical [m]
                        L_to_next_phys = mapper.to_phys_length(L_to_next_sim)
                        delta_L_to_next_phys = mapper.to_phys_length(delta_L_to_next_sim)
                        break

            # Forces from neighbors [sim → N]
            f_left_x = mapper.to_phys_force(F_left_sim[i, lateral_dim_idx])
            f_left_xi = mapper.to_phys_force(F_left_sim[i, amplitude_dim_idx])
            f_right_x = mapper.to_phys_force(F_right_sim[i, lateral_dim_idx])
            f_right_xi = mapper.to_phys_force(F_right_sim[i, amplitude_dim_idx])

            # Energy and lateralization (already in J and dimensionless)
            e_amp_kin = E_amp_kin[i]
            e_amp_pot = E_amp_pot[i]
            e_lat_kin = E_lat_kin[i]
            e_lat_pot = E_lat_pot[i]
            r_lat = R_lat[i]

            # Grid spacing in physical units [m]
            h_phys = mapper.to_phys_length(h_sim)

            # Write row (all in SI units)
            writer.writerow([
                i,
                x_pos, xi_pos,
                x_vel, xi_vel,
                x_acc, xi_acc,
                delta_x, xi,
                h_phys,
                L_to_next_phys, delta_L_to_next_phys,
                f_left_x, f_left_xi, f_right_x, f_right_xi,
                e_amp_kin, e_amp_pot, e_lat_kin, e_lat_pot, r_lat
            ])