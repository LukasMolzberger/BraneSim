"""
Dimension-agnostic tensor operations for diagnostics.

This module provides utilities that enable diagnostics to work seamlessly
across 1D, 2D, and 3D grids without dimension-specific branching.

Key operations:
- Reshaping flat arrays to grid structure
- Shifting tensors along arbitrary axes
- Reducing over transverse dimensions
"""

from __future__ import annotations
from typing import Literal
import torch
from .types import GridSpec


def reshape_flat_to_grid(
    x_flat: torch.Tensor,
    grid: GridSpec,
    component_dim: int | None = None
) -> torch.Tensor:
    """
    Reshape flat array to grid structure.

    Parameters
    ----------
    x_flat : torch.Tensor
        Flat tensor, shape [N] or [N, C]
    grid : GridSpec
        Grid specification
    component_dim : int | None
        If not None, indicates x_flat has shape [N, C] with C components

    Returns
    -------
    torch.Tensor
        Reshaped tensor: [*shape] or [*shape, C]

    Examples
    --------
    >>> # Scalar field
    >>> x_flat = torch.randn(100)
    >>> grid = GridSpec(shape=(10, 10), spacing_sim=1.0)
    >>> x_grid = reshape_flat_to_grid(x_flat, grid)
    >>> x_grid.shape
    torch.Size([10, 10])
    >>>
    >>> # Vector field
    >>> x_flat = torch.randn(100, 3)
    >>> x_grid = reshape_flat_to_grid(x_flat, grid, component_dim=3)
    >>> x_grid.shape
    torch.Size([10, 10, 3])
    """
    if x_flat.numel() != grid.num_points and x_flat.shape[0] != grid.num_points:
        raise ValueError(
            f"Tensor size {x_flat.shape} incompatible with grid shape {grid.shape} "
            f"(expected {grid.num_points} points)"
        )

    if component_dim is not None:
        # Vector field [N, C] -> [*shape, C]
        if x_flat.ndim != 2:
            raise ValueError(f"Expected 2D tensor for vector field, got shape {x_flat.shape}")
        return x_flat.reshape(*grid.shape, component_dim)
    else:
        # Scalar field [N] -> [*shape]
        if x_flat.ndim == 2 and x_flat.shape[1] == 1:
            x_flat = x_flat.squeeze(-1)
        return x_flat.reshape(grid.shape)


def shift_along_axis(
    tensor: torch.Tensor,
    axis: int,
    step: int = 1
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Shift tensor along specified axis for neighbor operations.

    Returns overlapping views for computing neighbor differences,
    overlaps, or gradients. Uses safe slicing (no torch.roll) to
    avoid wrapping at boundaries.

    Parameters
    ----------
    tensor : torch.Tensor
        Input tensor, shape [*shape] or [*shape, C]
    axis : int
        Axis along which to shift (0 to D-1)
    step : int, default=1
        Shift step size (positive for forward shift)

    Returns
    -------
    tuple[torch.Tensor, torch.Tensor]
        (left, right) where:
        - left: tensor[..., :-step, ...] (excludes last `step` elements along axis)
        - right: tensor[..., step:, ...] (excludes first `step` elements along axis)

    Examples
    --------
    >>> # 1D case: compute differences
    >>> x = torch.tensor([1.0, 2.0, 3.0, 4.0, 5.0])
    >>> left, right = shift_along_axis(x, axis=0, step=1)
    >>> left
    tensor([1., 2., 3., 4.])
    >>> right
    tensor([2., 3., 4., 5.])
    >>> diff = right - left
    >>> diff
    tensor([1., 1., 1., 1.])
    >>>
    >>> # 2D case: shift along axis=1 (columns)
    >>> x = torch.arange(12).reshape(3, 4).float()
    >>> left, right = shift_along_axis(x, axis=1, step=1)
    >>> left.shape, right.shape
    (torch.Size([3, 3]), torch.Size([3, 3]))
    >>>
    >>> # Vector field: shift along axis=0, preserve component dimension
    >>> x = torch.randn(10, 10, 3)  # [nx, ny, 3]
    >>> left, right = shift_along_axis(x, axis=0, step=1)
    >>> left.shape, right.shape
    (torch.Size([9, 10, 3]), torch.Size([9, 10, 3]))
    """
    if axis < 0 or axis >= tensor.ndim:
        # Handle component dimension case: axis refers to spatial dims only
        # e.g., tensor.shape = [nx, ny, C], axis=0 or 1 (not 2)
        if tensor.ndim > 1 and axis < tensor.ndim - 1:
            pass  # Valid spatial axis
        else:
            raise ValueError(f"Invalid axis {axis} for tensor shape {tensor.shape}")

    # Build slice objects for left and right views
    left_slice = [slice(None)] * tensor.ndim
    right_slice = [slice(None)] * tensor.ndim

    left_slice[axis] = slice(None, -step if step > 0 else None)
    right_slice[axis] = slice(step, None)

    left = tensor[tuple(left_slice)]
    right = tensor[tuple(right_slice)]

    return left, right


def reduce_transverse(
    tensor: torch.Tensor,
    axis: int,
    reduction: Literal["none", "mean", "median", "max"] = "none"
) -> torch.Tensor:
    """
    Reduce tensor over all dimensions except the specified axis.

    This is useful for extracting 1D profiles from higher-dimensional
    data (e.g., averaging Berry phase over transverse coordinates).

    Parameters
    ----------
    tensor : torch.Tensor
        Input tensor, shape [*shape] or [*shape, C]
    axis : int
        Axis to preserve (0 to D-1)
    reduction : {"none", "mean", "median", "max"}
        Reduction operation:
        - "none": return full tensor (no reduction)
        - "mean": average over transverse dimensions
        - "median": median over transverse dimensions
        - "max": maximum over transverse dimensions

    Returns
    -------
    torch.Tensor
        Reduced tensor. If reduction="none", returns input unchanged.
        Otherwise, returns 1D profile along specified axis.

    Examples
    --------
    >>> # 2D tensor: average over axis=1 to get profile along axis=0
    >>> x = torch.randn(10, 20)
    >>> profile = reduce_transverse(x, axis=0, reduction="mean")
    >>> profile.shape
    torch.Size([10])
    >>>
    >>> # 3D tensor: average over axes 1,2 to get profile along axis=0
    >>> x = torch.randn(32, 32, 32)
    >>> profile = reduce_transverse(x, axis=0, reduction="mean")
    >>> profile.shape
    torch.Size([32])
    >>>
    >>> # No reduction
    >>> x = torch.randn(10, 20)
    >>> result = reduce_transverse(x, axis=0, reduction="none")
    >>> result.shape
    torch.Size([10, 20])
    """
    if reduction == "none":
        return tensor

    # Determine dimensions to reduce over (all except axis)
    reduce_dims = [i for i in range(tensor.ndim) if i != axis]

    if not reduce_dims:
        # Already 1D along the specified axis
        return tensor

    if reduction == "mean":
        return tensor.mean(dim=reduce_dims)
    elif reduction == "median":
        # torch.median requires single dimension, so apply iteratively
        result = tensor
        for dim in sorted(reduce_dims, reverse=True):
            result = result.median(dim=dim).values
        return result
    elif reduction == "max":
        result = tensor
        for dim in sorted(reduce_dims, reverse=True):
            result = result.max(dim=dim).values
        return result
    else:
        raise ValueError(f"Unknown reduction: {reduction}")