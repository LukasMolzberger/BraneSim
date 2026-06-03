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
from branesim.diagnostics.alpha_separability import (
    projection_operators,
    g_factor,
    closed_form_observables,
    alpha_curve,
    numerical_trace_traceless,
    build_dynamical_block_3d,
    verify_track_a,
    group_velocity_ratio_p1,
    survey_alpha_grid,
)

__all__ = [
    # confinement
    "confinement_metrics_per_slice",
    "confinement_summary",
    "confinement_from_worldvolume",
    # alpha separability (Track A + P1)
    "projection_operators",
    "g_factor",
    "closed_form_observables",
    "alpha_curve",
    "numerical_trace_traceless",
    "build_dynamical_block_3d",
    "verify_track_a",
    "group_velocity_ratio_p1",
    "survey_alpha_grid",
]
