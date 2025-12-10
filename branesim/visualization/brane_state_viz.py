"""
3D Brane State Visualization

Generalized visualization tools for 3D brane configurations.
Handles position and velocity fields across different slices (XY, XZ, YZ).
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from typing import Optional, Tuple, List
import torch

from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid
from branesim.physics.dimensional_mapping import DimensionalMapper


def extract_slice_xy(field_3d: np.ndarray, grid_shape: Tuple[int, int, int], z_index: Optional[int] = None) -> np.ndarray:
    """Extract 2D XY slice from 3D field at constant z."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if z_index is None:
        z_index = nz // 2
    return field[:, :, z_index]


def extract_slice_xz(field_3d: np.ndarray, grid_shape: Tuple[int, int, int], y_index: Optional[int] = None) -> np.ndarray:
    """Extract 2D XZ slice from 3D field at constant y."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if y_index is None:
        y_index = ny // 2
    return field[:, y_index, :]


def extract_slice_yz(field_3d: np.ndarray, grid_shape: Tuple[int, int, int], x_index: Optional[int] = None) -> np.ndarray:
    """Extract 2D YZ slice from 3D field at constant x."""
    nx, ny, nz = grid_shape
    field = field_3d.reshape(nx, ny, nz)
    if x_index is None:
        x_index = nx // 2
    return field[x_index, :, :]


class BraneStateVisualizer:
    """
    Visualizer for 3D brane states.

    Generates comprehensive visualizations of position and velocity fields
    across different slice orientations.
    """

    def __init__(
        self,
        state: BraneState,
        grid: BraneGrid,
        mapper: DimensionalMapper,
        initial_positions: Optional[torch.Tensor] = None
    ):
        """
        Initialize visualizer.

        Args:
            state: BraneState object containing positions and velocities
            grid: BraneGrid object with grid information
            mapper: DimensionalMapper for unit conversions
            initial_positions: Initial flat positions (for computing displacements)
        """
        self.state = state
        self.grid = grid
        self.mapper = mapper
        self.initial_positions = initial_positions

        # Extract grid shape
        self.nx, self.ny, self.nz = grid.grid_shape

        # Prepare coordinate arrays (in nm)
        h_phys = mapper.to_phys_length(grid.spacing)
        self.x_coords = np.arange(self.nx) * h_phys * 1e9
        self.y_coords = np.arange(self.ny) * h_phys * 1e9
        self.z_coords = np.arange(self.nz) * h_phys * 1e9

        # Extract and convert data
        self._prepare_data()

    def _prepare_data(self):
        """Extract and convert brane state data to physical units."""
        # Get positions
        if self.initial_positions is not None:
            # Compute displacements from flat configuration
            self.pos_x = (self.state.positions[:, 0] - self.initial_positions[:, 0]).cpu().numpy()
            self.pos_y = (self.state.positions[:, 1] - self.initial_positions[:, 1]).cpu().numpy()
            self.pos_z = (self.state.positions[:, 2] - self.initial_positions[:, 2]).cpu().numpy()
        else:
            # Use absolute positions
            self.pos_x = self.state.positions[:, 0].cpu().numpy()
            self.pos_y = self.state.positions[:, 1].cpu().numpy()
            self.pos_z = self.state.positions[:, 2].cpu().numpy()

        self.pos_ampl = self.state.positions[:, 3].cpu().numpy()

        # Get velocities
        self.vel_x = self.state.velocities[:, 0].cpu().numpy()
        self.vel_y = self.state.velocities[:, 1].cpu().numpy()
        self.vel_z = self.state.velocities[:, 2].cpu().numpy()
        self.vel_ampl = self.state.velocities[:, 3].cpu().numpy()

        # Convert to physical units
        # Positions in picometers
        self.pos_x_phys = self.mapper.to_phys_length(self.pos_x) * 1e12
        self.pos_y_phys = self.mapper.to_phys_length(self.pos_y) * 1e12
        self.pos_z_phys = self.mapper.to_phys_length(self.pos_z) * 1e12
        self.pos_ampl_phys = self.mapper.to_phys_length(self.pos_ampl) * 1e12

        # Velocities as fraction of c
        from branesim.config.physical_constants import PhysicalConstants
        constants = PhysicalConstants()
        self.vel_x_phys = self.mapper.to_phys_velocity(self.vel_x) / constants.c
        self.vel_y_phys = self.mapper.to_phys_velocity(self.vel_y) / constants.c
        self.vel_z_phys = self.mapper.to_phys_velocity(self.vel_z) / constants.c
        self.vel_ampl_phys = self.mapper.to_phys_velocity(self.vel_ampl) / constants.c

    def plot_all_components(
        self,
        output_dir: str = ".",
        filename_prefix: str = "brane_state",
        dpi: int = 150,
        export_csv: bool = True,
        csv_output_dir: str = None
    ) -> List[str]:
        """
        Generate all 24 component plots (3 slices × 4 components × 2 quantities).

        Args:
            output_dir: Directory to save plots
            filename_prefix: Prefix for filenames
            dpi: DPI for saved images
            export_csv: Whether to export CSV files alongside plots
            csv_output_dir: Directory for CSV files (if None, uses output_dir)

        Returns:
            List of saved filenames (both PNG and CSV if enabled)
        """
        import os
        os.makedirs(output_dir, exist_ok=True)

        # Use separate directory for CSV files if specified
        if csv_output_dir is None:
            csv_output_dir = output_dir
        else:
            os.makedirs(csv_output_dir, exist_ok=True)

        component_names = ['x', 'y', 'z', 'amplitude']
        slice_names = ['xy', 'xz', 'yz']

        saved_files = []

        # Position plots
        pos_data = [self.pos_x_phys, self.pos_y_phys, self.pos_z_phys, self.pos_ampl_phys]
        for comp_name, data in zip(component_names, pos_data):
            for slice_name in slice_names:
                filename = f"{filename_prefix}_position_{comp_name}_{slice_name}.png"
                filepath = os.path.join(output_dir, filename)
                self._plot_single_component(
                    data, comp_name.upper(), slice_name.upper(),
                    "Position", "pm", filepath, dpi
                )
                saved_files.append(filename)

                # Export CSV
                if export_csv:
                    csv_filename = f"{filename_prefix}_position_{comp_name}_{slice_name}.csv"
                    csv_filepath = os.path.join(csv_output_dir, csv_filename)
                    self._export_csv(data, slice_name.upper(), csv_filepath)
                    saved_files.append(csv_filename)

        # Velocity plots
        vel_data = [self.vel_x_phys, self.vel_y_phys, self.vel_z_phys, self.vel_ampl_phys]
        for comp_name, data in zip(component_names, vel_data):
            for slice_name in slice_names:
                filename = f"{filename_prefix}_velocity_{comp_name}_{slice_name}.png"
                filepath = os.path.join(output_dir, filename)
                self._plot_single_component(
                    data, comp_name.upper(), slice_name.upper(),
                    "Velocity", "c", filepath, dpi
                )
                saved_files.append(filename)

                # Export CSV
                if export_csv:
                    csv_filename = f"{filename_prefix}_velocity_{comp_name}_{slice_name}.csv"
                    csv_filepath = os.path.join(csv_output_dir, csv_filename)
                    self._export_csv(data, slice_name.upper(), csv_filepath)
                    saved_files.append(csv_filename)

        return saved_files

    def _plot_single_component(
        self,
        data: np.ndarray,
        component_name: str,
        slice_name: str,
        quantity_name: str,
        unit: str,
        filepath: str,
        dpi: int
    ):
        """Plot a single component slice and save to file."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Extract slice
        if slice_name == 'XY':
            slice_data = extract_slice_xy(data, (self.nx, self.ny, self.nz))
            extent = [self.x_coords[0], self.x_coords[-1], self.y_coords[0], self.y_coords[-1]]
            xlabel, ylabel = 'x [nm]', 'y [nm]'
        elif slice_name == 'XZ':
            slice_data = extract_slice_xz(data, (self.nx, self.ny, self.nz))
            extent = [self.x_coords[0], self.x_coords[-1], self.z_coords[0], self.z_coords[-1]]
            xlabel, ylabel = 'x [nm]', 'z [nm]'
        else:  # YZ
            slice_data = extract_slice_yz(data, (self.nx, self.ny, self.nz))
            extent = [self.y_coords[0], self.y_coords[-1], self.z_coords[0], self.z_coords[-1]]
            xlabel, ylabel = 'y [nm]', 'z [nm]'

        # Determine color scale
        vmax = np.abs(slice_data).max()
        if vmax == 0:
            vmax = 1e-10

        im = ax.imshow(
            slice_data.T,
            origin='lower',
            extent=extent,
            cmap='RdBu_r',
            vmin=-vmax * 1.1,
            vmax=vmax * 1.1,
            aspect='equal',
            interpolation='nearest'
        )

        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(f'{quantity_name} - {component_name} ({slice_name} slice)',
                    fontsize=14, fontweight='bold')

        # Colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(unit, fontsize=11)

        plt.tight_layout()
        plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

    def _export_csv(
        self,
        data: np.ndarray,
        slice_name: str,
        filepath: str
    ):
        """Export a single component slice as CSV."""
        import csv

        # Extract slice
        if slice_name == 'XY':
            slice_data = extract_slice_xy(data, (self.nx, self.ny, self.nz))
            row_coords = self.x_coords
            col_coords = self.y_coords
            row_label = 'x[nm]'
            col_label = 'y[nm]'
        elif slice_name == 'XZ':
            slice_data = extract_slice_xz(data, (self.nx, self.ny, self.nz))
            row_coords = self.x_coords
            col_coords = self.z_coords
            row_label = 'x[nm]'
            col_label = 'z[nm]'
        else:  # YZ
            slice_data = extract_slice_yz(data, (self.nx, self.ny, self.nz))
            row_coords = self.y_coords
            col_coords = self.z_coords
            row_label = 'y[nm]'
            col_label = 'z[nm]'

        # Write CSV file
        with open(filepath, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Header row: first column is row coordinate label, then all column coordinates
            header = [f'{row_label}\\{col_label}'] + [f'{coord:.6f}' for coord in col_coords]
            writer.writerow(header)

            # Data rows: first column is row coordinate, then all data values
            for i, row_coord in enumerate(row_coords):
                row = [f'{row_coord:.6f}'] + [f'{slice_data[i, j]:.6e}' for j in range(len(col_coords))]
                writer.writerow(row)

    def print_statistics(self):
        """Print statistics for all components."""
        print(f"\n{'=' * 70}")
        print("Brane State Statistics:")
        print(f"{'=' * 70}")

        print(f"\nPosition Displacements [pm]:")
        print(f"  X: min={self.pos_x_phys.min():.3f}, max={self.pos_x_phys.max():.3f}, "
              f"rms={np.sqrt(np.mean(self.pos_x_phys**2)):.3f}")
        print(f"  Y: min={self.pos_y_phys.min():.3f}, max={self.pos_y_phys.max():.3f}, "
              f"rms={np.sqrt(np.mean(self.pos_y_phys**2)):.3f}")
        print(f"  Z: min={self.pos_z_phys.min():.3f}, max={self.pos_z_phys.max():.3f}, "
              f"rms={np.sqrt(np.mean(self.pos_z_phys**2)):.3f}")
        print(f"  Amplitude: min={self.pos_ampl_phys.min():.3f}, max={self.pos_ampl_phys.max():.3f}, "
              f"rms={np.sqrt(np.mean(self.pos_ampl_phys**2)):.3f}")

        print(f"\nVelocities [fraction of c]:")
        print(f"  X: min={self.vel_x_phys.min():.6f}, max={self.vel_x_phys.max():.6f}, "
              f"rms={np.sqrt(np.mean(self.vel_x_phys**2)):.6f}")
        print(f"  Y: min={self.vel_y_phys.min():.6f}, max={self.vel_y_phys.max():.6f}, "
              f"rms={np.sqrt(np.mean(self.vel_y_phys**2)):.6f}")
        print(f"  Z: min={self.vel_z_phys.min():.6f}, max={self.vel_z_phys.max():.6f}, "
              f"rms={np.sqrt(np.mean(self.vel_z_phys**2)):.6f}")
        print(f"  Amplitude: min={self.vel_ampl_phys.min():.6f}, max={self.vel_ampl_phys.max():.6f}, "
              f"rms={np.sqrt(np.mean(self.vel_ampl_phys**2)):.6f}")


def visualize_brane_state(
    state: BraneState,
    grid: BraneGrid,
    mapper: DimensionalMapper,
    output_dir: str = ".",
    filename_prefix: str = "brane_state",
    initial_positions: Optional[torch.Tensor] = None,
    print_stats: bool = True,
    dpi: int = 150,
    export_csv: bool = True,
    csv_output_dir: str = None
) -> List[str]:
    """
    Convenience function to visualize a brane state.

    Args:
        state: BraneState to visualize
        grid: BraneGrid
        mapper: DimensionalMapper
        output_dir: Directory for plot files
        filename_prefix: Prefix for filenames
        initial_positions: Initial positions (for displacement calculation)
        print_stats: Whether to print statistics
        dpi: DPI for saved images
        export_csv: Whether to export CSV files alongside plots
        csv_output_dir: Directory for CSV files (if None, uses output_dir)

    Returns:
        List of saved filenames (both PNG and CSV if enabled)
    """
    viz = BraneStateVisualizer(state, grid, mapper, initial_positions)

    if print_stats:
        viz.print_statistics()

    saved_files = viz.plot_all_components(output_dir, filename_prefix, dpi, export_csv, csv_output_dir)

    num_plots = len([f for f in saved_files if f.endswith('.png')])
    num_csvs = len([f for f in saved_files if f.endswith('.csv')])

    if csv_output_dir and csv_output_dir != output_dir:
        print(f"\n✓ Saved {num_plots} plots to {output_dir}")
        print(f"✓ Saved {num_csvs} CSV files to {csv_output_dir}")
    else:
        print(f"\n✓ Saved {num_plots} plots and {num_csvs} CSV files to {output_dir}")

    return saved_files