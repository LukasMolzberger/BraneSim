"""
Dimension-agnostic data types for diagnostics.

This module defines the core data structures used throughout the diagnostic
system, enabling dimension-independent computation (1D/2D/3D).

Key abstractions:
- Axis: Enumeration for spatial axes
- GridSpec: Specification of computational grid (shape, spacing, coordinates)
- Snapshot: Container for field data at a specific time
- DiagnosticResult: Structured output from diagnostic computations
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Any
import torch
import numpy as np


class Axis(IntEnum):
    """Spatial axis enumeration for dimension-agnostic code."""
    X = 0
    Y = 1
    Z = 2


@dataclass(frozen=True)
class GridSpec:
    """
    Specification of a computational grid (dimension-agnostic).

    Attributes
    ----------
    shape : tuple[int, ...]
        Grid shape (length D for D-dimensional grid)
    spacing_sim : float
        Grid spacing in simulation units (assumed uniform)
    coords_sim : tuple[torch.Tensor, ...] | None
        Coordinate arrays in simulation units [optional, can be lazily computed]
    coords_phys : tuple[torch.Tensor, ...] | None
        Coordinate arrays in physical units [optional]

    Examples
    --------
    >>> # 1D grid with 100 points
    >>> grid = GridSpec(shape=(100,), spacing_sim=1.0)
    >>>
    >>> # 2D grid with 50x50 points
    >>> grid = GridSpec(shape=(50, 50), spacing_sim=0.5)
    >>>
    >>> # 3D grid with explicit coordinates
    >>> x = torch.linspace(0, 10, 64)
    >>> y = torch.linspace(0, 10, 64)
    >>> z = torch.linspace(0, 10, 64)
    >>> grid = GridSpec(shape=(64, 64, 64), spacing_sim=10.0/63, coords_sim=(x, y, z))
    """
    shape: tuple[int, ...]
    spacing_sim: float
    coords_sim: tuple[torch.Tensor, ...] | None = None
    coords_phys: tuple[torch.Tensor, ...] | None = None

    @property
    def D(self) -> int:
        """Number of spatial dimensions."""
        return len(self.shape)

    @property
    def num_points(self) -> int:
        """Total number of grid points."""
        result = 1
        for s in self.shape:
            result *= s
        return result

    def is_compatible(self, tensor: torch.Tensor) -> bool:
        """
        Check if a tensor is compatible with this grid.

        Compatible means:
        - tensor.shape == self.shape (scalar field), or
        - tensor.shape == (*self.shape, C) (vector field with C components)
        - tensor.numel() == self.num_points (flat representation)
        """
        if tensor.shape == self.shape:
            return True
        if len(tensor.shape) == self.D + 1 and tensor.shape[:self.D] == self.shape:
            return True
        if tensor.ndim == 1 and tensor.numel() == self.num_points:
            return True
        if tensor.ndim == 2 and tensor.shape[0] == self.num_points:
            return True
        return False


@dataclass
class Snapshot:
    """
    Container for field data at a specific time.

    A Snapshot represents the substrate state (or derived fields) at a single
    moment in time. It provides a uniform interface for diagnostics to access
    field data without depending on solver internals.

    Attributes
    ----------
    t_sim : float | None
        Time in simulation units
    t_phys_s : float | None
        Time in physical units (seconds)
    fields : dict[str, torch.Tensor]
        Dictionary of field data. Common keys:
        - "xi" or "q": position/displacement field
        - "v_xi" or "v": velocity field
        - "pos": full 4D positions [N, 4]
        - "vel": full 4D velocities [N, 4]
    meta : dict[str, Any]
        Metadata (omega, wavelength, device, dtype, etc.)

    Examples
    --------
    >>> # Create snapshot from 1D simulation data
    >>> snapshot = Snapshot(
    ...     t_sim=10.0,
    ...     t_phys_s=1e-15,
    ...     fields={
    ...         "xi": state.positions[:, 3],
    ...         "v_xi": state.velocities[:, 3],
    ...     },
    ...     meta={"omega_sim": 6.28, "device": "cpu"}
    ... )
    """
    t_sim: float | None = None
    t_phys_s: float | None = None
    fields: dict[str, torch.Tensor] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)


@dataclass
class DiagnosticResult:
    """
    Structured output from a diagnostic computation.

    DiagnosticResult separates:
    - Primary computed data (Berry phase, spectrum, energy, etc.)
    - Quality metrics (masks, overlaps, validity flags)
    - Metadata (units, configuration, axis information)

    This structure enables:
    - Systematic export (all fields to NPZ/CSV)
    - Quality-aware visualization (plot only valid regions)
    - Reproducibility (configuration recorded in metadata)

    Attributes
    ----------
    name : str
        Diagnostic name (e.g., "berry_connection_x", "spectrum_1d")
    t_phys_s : float | None
        Time in physical units (seconds)
    t_sim : float | None
        Time in simulation units
    data : dict[str, torch.Tensor | np.ndarray]
        Primary diagnostic outputs. Keys depend on diagnostic type.
        Examples:
        - Berry: "dphi", "A_axis", "gamma_wrapped", "gamma_unwrapped", "curvature"
        - Spectrum: "k", "power", "power_mean", "power_std"
        - Energy: "kinetic", "potential", "total"
    quality : dict[str, torch.Tensor | np.ndarray]
        Quality/validity metrics. Examples:
        - "mask_point": boolean mask for valid points
        - "valid_edge": boolean mask for valid edges/overlaps
        - "overlap_abs": magnitude of wavefunction overlaps
        - "confidence": numerical confidence/weight per point
    meta : dict[str, Any]
        Metadata about the diagnostic. Examples:
        - "units": {"A_axis": "rad/m", "k": "1/m"}
        - "axis": 0 (which axis was analyzed)
        - "config": the configuration object used
        - "method": "berry_plaquette" or "berry_connection"

    Examples
    --------
    >>> # Berry phase result
    >>> result = DiagnosticResult(
    ...     name="berry_connection_x",
    ...     t_sim=10.0,
    ...     data={
    ...         "dphi": dphi_array,
    ...         "A_x": connection_array,
    ...         "gamma_wrapped": phase_array,
    ...     },
    ...     quality={
    ...         "mask_point": amplitude_mask,
    ...         "valid_edge": edge_mask,
    ...         "overlap_abs": overlap_magnitudes,
    ...     },
    ...     meta={
    ...         "axis": 0,
    ...         "spacing_sim": 1.0,
    ...         "config": berry_cfg,
    ...     }
    ... )
    """
    name: str
    t_phys_s: float | None = None
    t_sim: float | None = None
    data: dict[str, torch.Tensor | np.ndarray] = field(default_factory=dict)
    quality: dict[str, torch.Tensor | np.ndarray] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)