"""
Baseline State Computation for Lateral Distortion Diagnostics

This module provides utilities for computing and storing baseline brane configurations.
The baseline represents the equilibrium state (flat or pre-tensioned) against which
lateral distortions are measured.
"""

import torch
import numpy as np
from typing import Tuple


def compute_flat_baseline_positions(
    grid_shape: Tuple[int, int, int],
    h: float,
    center: Tuple[float, float, float] = None,
    device: str = 'cpu',
    dtype: torch.dtype = torch.float32,
) -> torch.Tensor:
    """
    Compute baseline positions for a flat, uniform grid.

    This creates the reference configuration for measuring lateral distortions.
    The grid is centered at `center` if provided, otherwise at the geometric center.

    Args:
        grid_shape: (nx, ny, nz) number of points in each direction
        h: Grid spacing [m]
        center: (x, y, z) center coordinates [m]. If None, uses geometric center.
        device: Torch device ('cpu' or 'cuda')
        dtype: Torch dtype

    Returns:
        baseline_positions: [nx*ny*nz, 3] tensor of (x, y, z) positions
    """
    nx, ny, nz = grid_shape

    # Create 1D arrays for each dimension
    x = torch.arange(nx, device=device, dtype=dtype) * h
    y = torch.arange(ny, device=device, dtype=dtype) * h
    z = torch.arange(nz, device=device, dtype=dtype) * h

    # Center the grid if requested
    if center is not None:
        cx, cy, cz = center
        x = x - (nx - 1) * h / 2 + cx
        y = y - (ny - 1) * h / 2 + cy
        z = z - (nz - 1) * h / 2 + cz

    # Create 3D meshgrid
    X, Y, Z = torch.meshgrid(x, y, z, indexing='ij')

    # Flatten and stack into [N, 3] array
    baseline_positions = torch.stack([
        X.reshape(-1),
        Y.reshape(-1),
        Z.reshape(-1)
    ], dim=1)

    return baseline_positions


def compute_lateral_distortion(
    positions: torch.Tensor,
    baseline_positions: torch.Tensor,
) -> torch.Tensor:
    """
    Compute lateral distortion magnitude relative to baseline.

    Distortion is the Euclidean distance between current and baseline positions
    in the lateral (X^0, X^1, X^2) coordinates only.

    Args:
        positions: [N, 3] or [N, 4] current positions. If 4D, uses first 3 components.
        baseline_positions: [N, 3] baseline reference positions

    Returns:
        distortion: [N] magnitude of lateral displacement |r - r_baseline|
    """
    # Extract lateral positions (first 3 components)
    if positions.shape[1] == 4:
        pos_lateral = positions[:, :3]
    else:
        pos_lateral = positions

    # Compute displacement vector
    displacement = pos_lateral - baseline_positions

    # Compute magnitude
    distortion = torch.linalg.norm(displacement, dim=1)

    return distortion


def compute_lateral_distortion_grid(
    positions: torch.Tensor,
    baseline_positions: torch.Tensor,
    grid_shape: Tuple[int, int, int],
) -> np.ndarray:
    """
    Compute lateral distortion and reshape to grid for visualization.

    Args:
        positions: [N, 3] or [N, 4] current positions
        baseline_positions: [N, 3] baseline reference positions
        grid_shape: (nx, ny, nz) shape for reshaping

    Returns:
        distortion_grid: [nx, ny, nz] array of distortion magnitudes
    """
    distortion = compute_lateral_distortion(positions, baseline_positions)

    # Convert to numpy and reshape
    distortion_np = distortion.cpu().numpy()
    distortion_grid = distortion_np.reshape(grid_shape)

    return distortion_grid


def initialize_baseline_state(config: dict) -> dict:
    """
    Initialize baseline state from configuration.

    This is the main entry point for setting up baseline positions
    for distortion measurements.

    Args:
        config: Configuration dictionary with keys:
            - grid_shape: (nx, ny, nz)
            - h: grid spacing [m]
            - center: (optional) (x, y, z) center coordinates
            - device: (optional) 'cpu' or 'cuda'
            - dtype: (optional) torch dtype

    Returns:
        Dictionary with baseline information:
            - positions: [N, 3] baseline positions
            - grid_shape: (nx, ny, nz)
            - h: grid spacing
    """
    grid_shape = config['grid_shape']
    h = config['h']
    center = config.get('center', None)
    device = config.get('device', 'cpu')
    dtype = config.get('dtype', torch.float32)

    baseline_positions = compute_flat_baseline_positions(
        grid_shape=grid_shape,
        h=h,
        center=center,
        device=device,
        dtype=dtype,
    )

    return {
        'positions': baseline_positions,
        'grid_shape': grid_shape,
        'h': h,
    }


# ==============================================================================
# Debug and Validation Functions
# ==============================================================================

def validate_baseline_initialization(
    state_positions: torch.Tensor,
    baseline_positions: torch.Tensor,
    tolerance: float = 1e-10,
) -> Tuple[bool, float]:
    """
    Validate that initial state positions match baseline (for states without lateral initialization).

    Args:
        state_positions: [N, 3] or [N, 4] state positions
        baseline_positions: [N, 3] baseline positions
        tolerance: Maximum allowed deviation [m]

    Returns:
        (is_valid, max_deviation): Whether positions match within tolerance, and max deviation
    """
    if state_positions.shape[1] == 4:
        pos_lateral = state_positions[:, :3]
    else:
        pos_lateral = state_positions

    deviation = torch.abs(pos_lateral - baseline_positions).max().item()
    is_valid = deviation < tolerance

    return is_valid, deviation