"""Physics module for brane simulations."""

# Import parameters (no dependencies)
from branesim.physics.parameters import (
    compton_calibrated_brane_lattice_params,
    print_calibration_summary
)

# Import forces (requires torch - optional)
try:
    from branesim.physics.forces import SpringForceComputer
    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    SpringForceComputer = None

__all__ = [
    'compton_calibrated_brane_lattice_params',
    'print_calibration_summary',
]

if _HAS_TORCH:
    __all__.append('SpringForceComputer')