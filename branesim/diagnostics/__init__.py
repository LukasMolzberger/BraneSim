"""branesim.diagnostics — read-only diagnostic measurements.

Diagnostic modules compute derived quantities from the world-volume state
(node positions, reference positions) and return plain dicts or arrays.
They must never modify solver state and must never be imported from
simulation/solver code (principles §2, §7.2).
"""

from branesim.diagnostics.confinement import (
    confinement_metrics_per_slice,
    confinement_summary,
    confinement_from_worldvolume,
)

__all__ = [
    "confinement_metrics_per_slice",
    "confinement_summary",
    "confinement_from_worldvolume",
]
