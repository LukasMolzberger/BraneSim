"""branesim.berry

Berry-phase diagnostics for the brane simulator.

Design goals
------------
- **Analytical part is dimensionality-agnostic**: it operates on flattened
  lattices with shape (N, E), where N is the number of brane points and E is
  the chosen embedding-component count (e.g. 2 for 1D-in-2D, 3 for 2D-in-3D,
  4 for 3D-in-4D).
- **Visualization is dimensionality-dependent**, similar to the existing
  displacement-field videos.

Hardcoded Compton carrier
-------------------------
This package uses a *hardcoded* carrier angular frequency derived from the
Compton wavelength stored in ``branesim.config.physical_constants``:

    omega0 = 2*pi*c / lambda_C

This follows the convention already used in the photon experiments.

Public API
----------
- :func:`branesim.berry.analysis.compute_berry_time_series`
- :func:`branesim.berry.analysis.compton_omega0`
- :func:`branesim.berry.viz.create_berry_videos`

"""

from .analysis import BerryTimeSeries, compton_omega0, compute_berry_time_series
from .viz import create_berry_videos

__all__ = [
    "BerryTimeSeries",
    "compton_omega0",
    "compute_berry_time_series",
    "create_berry_videos",
]
