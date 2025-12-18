"""
Degeneracy and holonomy diagnostics (no filtering).

This module provides diagnostic tools that verify degeneracy and compute
holonomy WITHOUT relying on spectral band-pass filtering. These are used
to validate preparation-first narrowband initialization.
"""

from branesim.diagnostics.degeneracy import (
    local_eigen_scan,
    subspace_rank_svd,
    verify_narrowband_preparation,
)

from branesim.diagnostics.holonomy import (
    orthonormalize_frame,
    wilson_loop_holonomy,
    wilson_invariants,
    compute_u1_holonomy,
)

__all__ = [
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