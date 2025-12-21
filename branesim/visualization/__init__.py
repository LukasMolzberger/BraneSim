"""Visualization tools for BraneSim."""

from .brane_state_viz import (
    BraneStateVisualizer,
    visualize_brane_state,
    extract_slice_xy,
    extract_slice_xz,
    extract_slice_yz,
)
from .brane_1d_viz import (
    plot_brane_1d_amplitude_propagation,
    plot_brane_1d_lateral_distortion,
    plot_brane_1d_amplitude_velocity,
    plot_brane_1d_lateral_velocity,
    plot_brane_1d_tracking_analysis,
    plot_all_brane_1d_standard,
)

# Displacement field visualization (optional - can import directly from displacement_field_viz)
try:
    from .displacement_field_viz import (
        displacement_from_positions,
        displacement_frames_from_positions_frames,
        create_displacement_arrows_video_3d_in_4d,
        create_displacement_diralpha_slices_videos_3d_in_4d,
        create_displacement_arrows_video_2d_in_3d,
        create_displacement_diralpha_video_2d_in_3d,
        create_displacement_components_video_1d_in_2d,
        create_displacement_magnitude_angle_video_1d_in_2d,
    )
    _DISPLACEMENT_VIZ_AVAILABLE = True
except ImportError:
    _DISPLACEMENT_VIZ_AVAILABLE = False

__all__ = [
    # Brane state visualization
    'BraneStateVisualizer',
    'visualize_brane_state',
    'extract_slice_xy',
    'extract_slice_xz',
    'extract_slice_yz',
    # Generic 1D brane visualization (preferred)
    'plot_brane_1d_amplitude_propagation',
    'plot_brane_1d_lateral_distortion',
    'plot_brane_1d_amplitude_velocity',
    'plot_brane_1d_lateral_velocity',
    'plot_brane_1d_tracking_analysis',
    'plot_all_brane_1d_standard',
]

# Add displacement field viz if available
if _DISPLACEMENT_VIZ_AVAILABLE:
    __all__.extend([
        'displacement_from_positions',
        'displacement_frames_from_positions_frames',
        'create_displacement_arrows_video_3d_in_4d',
        'create_displacement_diralpha_slices_videos_3d_in_4d',
        'create_displacement_arrows_video_2d_in_3d',
        'create_displacement_diralpha_video_2d_in_3d',
        'create_displacement_components_video_1d_in_2d',
        'create_displacement_magnitude_angle_video_1d_in_2d',
    ])