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
from branesim.diagnostics.berry_holonomy import (
    heff_eigenframe,
    plaquette_holonomy_p2,
    verify_p2_all_alpha,
    rotate_and_transport,
    rotation_matrix_spin1,
    rotation_matrix_spin_half,
    spin1_frame_fn,
    spin1_state_fn,
    spin_half_frame_fn,
    spin_half_state_fn,
    verify_p3_so3_holonomy,
    decompose_wz_u1_sun,
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
    # berry / WZ holonomy (Track B, P2 + P3)
    "heff_eigenframe",
    "plaquette_holonomy_p2",
    "verify_p2_all_alpha",
    "rotate_and_transport",
    "rotation_matrix_spin1",
    "rotation_matrix_spin_half",
    "spin1_frame_fn",
    "spin1_state_fn",
    "spin_half_frame_fn",
    "spin_half_state_fn",
    "verify_p3_so3_holonomy",
    "decompose_wz_u1_sun",
]
