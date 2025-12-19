"""
Dimension-agnostic diagnostic tools for BraneSim.

This module provides:
- Data types: GridSpec, Snapshot, DiagnosticResult
- Complex band state construction
- Berry phase diagnostics (connection, phase, curvature)
- Spectrum analysis
- Energy diagnostics
- Holonomy and degeneracy verification
"""

# Core data types
from branesim.diagnostics.types import (
    Axis,
    GridSpec,
    Snapshot,
    DiagnosticResult,
)

# Analytic signal (ω-free complex state)
from branesim.diagnostics.analytic_signal import (
    analytic_signal_along_axis,
    pointwise_normalize_from_grid,
    pointwise_normalize_scalar,
    pointwise_normalize_vector,
)

# Berry phase diagnostics
from branesim.diagnostics.berry import (
    BerryConfig,
    BerryPhase1DConfig,
    berry_connection_along_axis,
    berry_phase_integrated_along_axis,
    berry_plaquette_curvature,
    berry_phase_profile_along_x,
)

# Tensor operations
from branesim.diagnostics.tensor_ops import (
    reshape_flat_to_grid,
    shift_along_axis,
    reduce_transverse,
)

# Degeneracy verification
from branesim.diagnostics.degeneracy import (
    local_eigen_scan,
    subspace_rank_svd,
    verify_narrowband_preparation,
)

# Holonomy computation
from branesim.diagnostics.holonomy import (
    orthonormalize_frame,
    wilson_loop_holonomy,
    wilson_invariants,
    compute_u1_holonomy,
)

__all__ = [
    # Data types
    "Axis",
    "GridSpec",
    "Snapshot",
    "DiagnosticResult",
    # Analytic signal (ω-free complex state)
    "analytic_signal_along_axis",
    "pointwise_normalize_from_grid",
    "pointwise_normalize_scalar",
    "pointwise_normalize_vector",
    # Berry phase
    "BerryConfig",
    "BerryPhase1DConfig",
    "berry_connection_along_axis",
    "berry_phase_integrated_along_axis",
    "berry_plaquette_curvature",
    "berry_phase_profile_along_x",
    # Tensor ops
    "reshape_flat_to_grid",
    "shift_along_axis",
    "reduce_transverse",
    # Degeneracy verification
    "local_eigen_scan",
    "subspace_rank_svd",
    "verify_narrowband_preparation",
    # Holonomy computation
    "orthonormalize_frame",
    "wilson_loop_holonomy",
    "wilson_invariants",
    "compute_u1_holonomy",
]