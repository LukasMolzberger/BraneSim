"""
MatplotlibRenderer: Visualization for brane simulations.

This module provides Matplotlib-based rendering for 1D, 2D, and 3D simulations.
"""

import matplotlib.pyplot as plt
import numpy as np
import torch

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid


class MatplotlibRenderer:
    """
    Matplotlib-based visualization for brane simulations.

    Supports:
        - 1D: Line plots of field components
        - 2D: Colormesh/imshow of field slices
        - 3D: Volume rendering or slice views
    """

    def __init__(
        self,
        dimension: Dimensionality,
        component_idx: int = 3,
        figsize: tuple = (10, 6)
    ):
        """
        Initialize renderer.

        Args:
            dimension: Dimensionality enum
            component_idx: Which embedding dimension to visualize (0-3, default 3)
            figsize: Figure size in inches
        """
        self.dimension = dimension
        self.component_idx = component_idx
        self.figsize = figsize

        # Create figure
        self.fig, self.ax = self._create_figure()

        # Storage for colorbar (2D)
        self.colorbar = None

    def _create_figure(self):
        """Create appropriate figure based on dimension."""
        if self.dimension == Dimensionality.ONE_D:
            fig, ax = plt.subplots(figsize=self.figsize)
            ax.set_xlabel('Position [m]')
            ax.set_ylabel(f'$\\xi^{{{self.component_idx}}}$ [m]')
            ax.grid(True, alpha=0.3)

        elif self.dimension == Dimensionality.TWO_D:
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.set_aspect('equal')
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')

        else:  # THREE_D
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure(figsize=self.figsize)
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel(f'$\\xi^{{{self.component_idx}}}$ [m]')

        return fig, ax

    def render_field(
        self,
        state: BraneState,
        grid: BraneGrid,
        title: str = "",
        show: bool = True
    ):
        """
        Render field component.

        Args:
            state: BraneState with positions
            grid: BraneGrid for spatial layout
            title: Plot title
            show: Whether to display the plot
        """
        # Extract field component
        field = state.get_field_component(self.component_idx).cpu().numpy()

        if self.dimension == Dimensionality.ONE_D:
            self._render_1d(field, grid)
        elif self.dimension == Dimensionality.TWO_D:
            self._render_2d(field, grid)
        else:
            self._render_3d(field, grid)

        self.ax.set_title(title)

        if show:
            plt.draw()
            plt.pause(0.001)

    def _render_1d(self, field: np.ndarray, grid: BraneGrid):
        """1D line plot."""
        x = np.arange(grid.grid_shape[0]) * grid.spacing

        self.ax.clear()
        self.ax.plot(x, field, 'b-', linewidth=2, label=f'$\\xi^{{{self.component_idx}}}$')
        self.ax.set_xlabel('Position [m]')
        self.ax.set_ylabel(f'$\\xi^{{{self.component_idx}}}$ [m]')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend()

    def _render_2d(self, field: np.ndarray, grid: BraneGrid):
        """2D heatmap."""
        # Reshape flat field to 2D grid
        field_2d = field.reshape(grid.grid_shape)

        self.ax.clear()
        im = self.ax.imshow(
            field_2d.T,  # Transpose for correct orientation
            extent=[0, grid.grid_shape[0] * grid.spacing,
                    0, grid.grid_shape[1] * grid.spacing],
            origin='lower',
            cmap='RdBu_r',
            interpolation='bilinear',
            aspect='equal'
        )

        if self.colorbar is None:
            self.colorbar = plt.colorbar(im, ax=self.ax)
            self.colorbar.set_label(f'$\\xi^{{{self.component_idx}}}$ [m]')
        else:
            im.set_clim(field_2d.min(), field_2d.max())

        self.ax.set_xlabel('X [m]')
        self.ax.set_ylabel('Y [m]')

    def _render_3d(self, field: np.ndarray, grid: BraneGrid):
        """3D slice visualization."""
        # Show central slice for 3D
        nz = grid.grid_shape[2]
        mid_z = nz // 2

        # Extract slice
        nx, ny = grid.grid_shape[0], grid.grid_shape[1]
        field_3d = field.reshape(grid.grid_shape)
        slice_2d = field_3d[:, :, mid_z]

        self.ax.clear()
        x = np.arange(nx) * grid.spacing
        y = np.arange(ny) * grid.spacing
        X, Y = np.meshgrid(x, y)

        surf = self.ax.plot_surface(
            X, Y, slice_2d.T,
            cmap='RdBu_r',
            linewidth=0,
            antialiased=True
        )

        self.ax.set_xlabel('X [m]')
        self.ax.set_ylabel('Y [m]')
        self.ax.set_zlabel(f'$\\xi^{{{self.component_idx}}}$ [m]')

    def render_intensity(
        self,
        state: BraneState,
        grid: BraneGrid,
        title: str = ""
    ):
        """
        Render energy intensity |ξ|².

        Args:
            state: BraneState
            grid: BraneGrid
            title: Plot title
        """
        intensity = torch.sum(state.positions**2, dim=1).cpu().numpy()

        if self.dimension == Dimensionality.ONE_D:
            x = np.arange(grid.grid_shape[0]) * grid.spacing
            self.ax.clear()
            self.ax.plot(x, intensity, 'r-', linewidth=2, label='$|\\xi|^2$')
            self.ax.set_xlabel('Position [m]')
            self.ax.set_ylabel('Intensity [m²]')
            self.ax.grid(True, alpha=0.3)
            self.ax.legend()

        elif self.dimension == Dimensionality.TWO_D:
            intensity_2d = intensity.reshape(grid.grid_shape)
            self.ax.clear()
            im = self.ax.imshow(
                intensity_2d.T,
                extent=[0, grid.grid_shape[0] * grid.spacing,
                        0, grid.grid_shape[1] * grid.spacing],
                origin='lower',
                cmap='hot',
                interpolation='bilinear'
            )

        self.ax.set_title(title)
        plt.draw()
        plt.pause(0.001)

    def save_frame(self, filename: str, dpi: int = 150):
        """
        Save current frame to file.

        Args:
            filename: Output filename
            dpi: Resolution in dots per inch
        """
        self.fig.savefig(filename, dpi=dpi, bbox_inches='tight')

    def close(self):
        """Close the figure."""
        plt.close(self.fig)

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"MatplotlibRenderer({self.dimension.name}, "
            f"component={self.component_idx})"
        )
