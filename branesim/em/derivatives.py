"""
Finite-difference operators for EM field solver.

This module implements spatial derivative operators (gradient, divergence, curl,
laplacian) using second-order central differences with periodic and Dirichlet
boundary conditions.
"""

import torch
from branesim.core.state import Dimensionality


def _ndim_from_grid(grid) -> int:
    """Extract spatial dimension from BraneGrid."""
    return grid.dimension.value


def central_diff(f: torch.Tensor, axis: int, h: float, bc: str) -> torch.Tensor:
    """
    First derivative ∂f/∂x_axis with central differences.

    Uses second-order central differences: (f[i+1] - f[i-1]) / (2h)
    For Dirichlet BC, uses one-sided differences at boundaries.

    Args:
        f: Field tensor [..., C] where last dim is channels (untouched)
        axis: Spatial axis to differentiate along (0, 1, or 2)
        h: Grid spacing in meters
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        Derivative tensor of same shape as f
    """
    if bc == "periodic":
        return (torch.roll(f, shifts=-1, dims=axis) - torch.roll(f, shifts=1, dims=axis)) / (2.0 * h)

    # Dirichlet: central interior, one-sided boundaries
    d = torch.zeros_like(f)
    ndim = f.dim()

    # Interior slice 1:-1 on axis
    center = [slice(None)] * ndim
    plus = [slice(None)] * ndim
    minus = [slice(None)] * ndim
    center[axis] = slice(1, -1)
    plus[axis] = slice(2, None)
    minus[axis] = slice(None, -2)

    d[tuple(center)] = (f[tuple(plus)] - f[tuple(minus)]) / (2.0 * h)

    # Left boundary (forward difference)
    left0 = [slice(None)] * ndim
    left1 = [slice(None)] * ndim
    left0[axis] = 0
    left1[axis] = 1
    d[tuple(left0)] = (f[tuple(left1)] - f[tuple(left0)]) / h

    # Right boundary (backward difference)
    rightm1 = [slice(None)] * ndim
    rightm2 = [slice(None)] * ndim
    rightm1[axis] = -1
    rightm2[axis] = -2
    d[tuple(rightm1)] = (f[tuple(rightm1)] - f[tuple(rightm2)]) / h

    return d


def laplacian(f: torch.Tensor, grid, bc: str = "periodic") -> torch.Tensor:
    """
    Compute Laplacian ∇²f for scalar or multi-channel fields.

    Uses second-order finite differences: (f[i+1] - 2f[i] + f[i-1]) / h²
    summed over all spatial dimensions.

    Args:
        f: Field tensor [*grid_shape, C] or [*grid_shape]
           Last dim treated as channels if present
        grid: BraneGrid instance with spacing and dimension info
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        Laplacian tensor of same shape as f
    """
    h = grid.spacing
    ndim = _ndim_from_grid(grid)

    if bc == "periodic":
        lap = torch.zeros_like(f)
        for ax in range(ndim):
            lap = lap + (torch.roll(f, -1, ax) - 2.0 * f + torch.roll(f, 1, ax)) / (h * h)
        return lap

    # Dirichlet: compute interior only; boundaries remain zero (good if boundaries are fixed)
    lap = torch.zeros_like(f)
    for ax in range(ndim):
        # Second derivative along axis ax on interior indices
        ndim_f = f.dim()
        center = [slice(None)] * ndim_f
        plus = [slice(None)] * ndim_f
        minus = [slice(None)] * ndim_f
        center[ax] = slice(1, -1)
        plus[ax] = slice(2, None)
        minus[ax] = slice(None, -2)
        lap[tuple(center)] += (f[tuple(plus)] - 2.0 * f[tuple(center)] + f[tuple(minus)]) / (h * h)
    return lap


def gradient_scalar(phi: torch.Tensor, grid, bc: str = "periodic") -> torch.Tensor:
    """
    Compute gradient of scalar field.

    Returns grad(phi) = (∂φ/∂x, ∂φ/∂y, ∂φ/∂z) as [..., 3].
    Missing dimensions (if ndim < 3) are filled with zeros.

    Args:
        phi: Scalar field tensor [*grid_shape]
        grid: BraneGrid instance with spacing and dimension info
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        Gradient tensor [*grid_shape, 3]
    """
    h = grid.spacing
    ndim = _ndim_from_grid(grid)
    grads = []
    for ax in range(ndim):
        grads.append(central_diff(phi, ax, h, bc))
    # Pad to 3 dimensions
    while len(grads) < 3:
        grads.append(torch.zeros_like(phi))
    return torch.stack(grads[:3], dim=-1)


def divergence(A: torch.Tensor, grid, bc: str = "periodic") -> torch.Tensor:
    """
    Compute divergence of vector field.

    Computes div(A) = ∂Ax/∂x + ∂Ay/∂y + ∂Az/∂z.
    Uses only available spatial dimensions (ignores higher components if ndim < 3).

    Args:
        A: Vector field tensor [*grid_shape, 3]
        grid: BraneGrid instance with spacing and dimension info
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        Divergence scalar field [*grid_shape]
    """
    h = grid.spacing
    ndim = _ndim_from_grid(grid)
    div = torch.zeros_like(A[..., 0])
    for ax in range(ndim):
        div = div + central_diff(A[..., ax], ax, h, bc)
    return div


def curl(A: torch.Tensor, grid, bc: str = "periodic") -> torch.Tensor:
    """
    Compute curl of vector field.

    Computes curl(A) = (∂Az/∂y - ∂Ay/∂z, ∂Ax/∂z - ∂Az/∂x, ∂Ay/∂x - ∂Ax/∂y).
    Missing derivatives (if ndim < 3) are treated as zero.

    Args:
        A: Vector field tensor [*grid_shape, 3]
        grid: BraneGrid instance with spacing and dimension info
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        Curl vector field [*grid_shape, 3]
    """
    h = grid.spacing
    ndim = _ndim_from_grid(grid)

    def d(comp, ax):
        """Helper to compute derivative, returning zero if axis unavailable."""
        if ax >= ndim:
            return torch.zeros_like(A[..., 0])
        return central_diff(comp, ax, h, bc)

    Ax, Ay, Az = A[..., 0], A[..., 1], A[..., 2]
    dAx_dx = d(Ax, 0)
    dAx_dy = d(Ax, 1)
    dAx_dz = d(Ax, 2)
    dAy_dx = d(Ay, 0)
    dAy_dy = d(Ay, 1)
    dAy_dz = d(Ay, 2)
    dAz_dx = d(Az, 0)
    dAz_dy = d(Az, 1)
    dAz_dz = d(Az, 2)

    # curl = (∂Az/∂y - ∂Ay/∂z, ∂Ax/∂z - ∂Az/∂x, ∂Ay/∂x - ∂Ax/∂y)
    cx = dAz_dy - dAy_dz
    cy = dAx_dz - dAz_dx
    cz = dAy_dx - dAx_dy
    return torch.stack([cx, cy, cz], dim=-1)