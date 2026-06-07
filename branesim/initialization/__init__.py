"""Initialization layer — pure initial-condition generators.

No solver, no diagnostics, no physics beyond setting node positions.
All functions return a full spacelike-slice configuration
    R0 : ndarray, shape (n_nodes, m_ambient), float64
  = ref + displacement.

The only nonlinearity in this project is the geometric Euclidean link-norm
inside the solver.  Seed functions here are pure arithmetic on positions;
they must not introduce additional energy terms, damping, or clamping.
"""

from branesim.initialization.seeds import (
    hedgehog,
    skyrme_twisted_hedgehog,
    axis_triplet,
)
from branesim.initialization.vortex_worldtube import (
    VortexParams,
    inject_vortex_worldtube,
    measure_winding_closure,
    project_carrier_re,
    CARRIER_RE_COMPONENTS,
    CARRIER_RE_WEIGHTS,
    CARRIER_IM,
)

__all__ = [
    "hedgehog",
    "skyrme_twisted_hedgehog",
    "axis_triplet",
    "VortexParams",
    "inject_vortex_worldtube",
    "measure_winding_closure",
    "project_carrier_re",
    "CARRIER_RE_COMPONENTS",
    "CARRIER_RE_WEIGHTS",
    "CARRIER_IM",
]
