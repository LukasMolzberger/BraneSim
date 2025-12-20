"""
Symplectic band structure and Berry phase diagnostics.

This package implements first-order symplectic methods for computing:
- Band structures ω(k) from linearized brane dynamics
- Berry phases and Wilson loops from eigenmode frames
- Dimension-agnostic diagnostics for 1D, 2D, and 3D branes

Key modules:
    symplectic_types: Data structures and configuration
    symplectic_builder: Operator construction with boundary conditions
    symplectic_solver: Eigenvalue solver along k-paths
    symplectic_berry: Berry phase and Wilson loop computation
    symplectic_viz: Visualization tools

Design principles:
    - Everything is prefixed "symplectic_" to distinguish from future second-order methods
    - Boundary conditions (PERIODIC/CLAMPED) affect operator building, not Berry diagnostics
    - All core computations are dimension-agnostic
    - Visualization may be dimension-specific

Example usage:
    >>> import torch
    >>> from branesim.diagnostics.bands import (
    ...     SymplecticBandConfig, BoundaryCondition, KPath,
    ...     solve_symplectic_bands_on_kpath,
    ...     compute_symplectic_wilson_loop,
    ...     plot_band_structure,
    ... )
    >>>
    >>> # Configure 1D lattice
    >>> cfg = SymplecticBandConfig(
    ...     d=1, embedding_dim=4, grid_shape=(64,),
    ...     spacing=1e-12, mass=1e-30, spring_k=1e-6, rest_length=1e-12,
    ...     neighbor_offsets=[(-1,), (1,)],
    ...     boundary=BoundaryCondition.PERIODIC,
    ... )
    >>>
    >>> # Define k-path: Γ → X (closed loop)
    >>> h = cfg.spacing
    >>> k_points = torch.cat([
    ...     torch.linspace(0, torch.pi/h, 50),
    ...     torch.linspace(torch.pi/h, 0, 50)
    ... ]).unsqueeze(1)
    >>> kpath = KPath(k_points=k_points, closed=True, label="Γ→X→Γ")
    >>>
    >>> # Solve bands
    >>> result = solve_symplectic_bands_on_kpath(cfg, kpath, n_modes=4)
    >>>
    >>> # Compute Wilson loop for first two bands (degenerate pair)
    >>> wilson = compute_symplectic_wilson_loop(result, band_indices=[0, 1])
    >>> print(f"Berry phases: {wilson.eigenphases}")
    >>>
    >>> # Visualize
    >>> fig, ax = plot_band_structure(result)
"""

# Core data types
from .symplectic_types import (
    BoundaryCondition,
    SymplecticBandConfig,
    KPath,
    SymplecticBandResult,
    SymplecticWilsonResult,
)

# Operator builder
from .symplectic_builder import (
    build_symplectic_operator_at_k,
    get_nearest_neighbor_offsets,
    validate_operator_symplecticity,
)

# Band solver
from .symplectic_solver import (
    solve_symplectic_bands_on_kpath,
    get_band_velocities,
)

# Berry phase / Wilson loop
from .symplectic_berry import (
    compute_symplectic_wilson_loop,
    compute_berry_connection_along_path,
    verify_gauge_covariance,
)

# Visualization
from .symplectic_viz import (
    plot_band_structure,
    plot_berry_phase_profile,
    plot_wilson_eigenphases,
    plot_polarization_vectors,
    plot_band_structure_2d,
)

__all__ = [
    # Types
    "BoundaryCondition",
    "SymplecticBandConfig",
    "KPath",
    "SymplecticBandResult",
    "SymplecticWilsonResult",
    # Builder
    "build_symplectic_operator_at_k",
    "get_nearest_neighbor_offsets",
    "validate_operator_symplecticity",
    # Solver
    "solve_symplectic_bands_on_kpath",
    "get_band_velocities",
    # Berry
    "compute_symplectic_wilson_loop",
    "compute_berry_connection_along_path",
    "verify_gauge_covariance",
    # Viz
    "plot_band_structure",
    "plot_berry_phase_profile",
    "plot_wilson_eigenphases",
    "plot_polarization_vectors",
    "plot_band_structure_2d",
]