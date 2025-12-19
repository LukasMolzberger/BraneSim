"""
Berry phase analysis tools for BraneSim.

DEPRECATED: This module re-exports from branesim.diagnostics for backward compatibility.
New code should import directly from branesim.diagnostics.
"""

# Re-export from new diagnostics location for backward compatibility
from branesim.diagnostics.complex_band_state import (
    complex_band_state_from_quadrature,
    pointwise_normalize,
)
from branesim.diagnostics.berry import (
    BerryPhase1DConfig,
    berry_phase_profile_along_x,
    BerryConfig,
)

__all__ = [
    "complex_band_state_from_quadrature",
    "pointwise_normalize",
    "BerryPhase1DConfig",
    "BerryConfig",
    "berry_phase_profile_along_x",
]