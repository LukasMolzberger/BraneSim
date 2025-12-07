"""
Debug Script: Visualize Electron Cross-Section and Phase Pattern

This script generates two debug plots to verify the electron initialization:
1. Cross-section envelope f(x,y) - should show clean double-lobe dumbbell
2. Phase pattern along centerline - should be smooth with m complete cycles

Run this BEFORE running the full simulation to verify the initialization is correct.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.electron_initialization import (
    calibrate_electron_init_params,
    plot_cross_section_envelope_debug,
    plot_phase_along_centerline_debug,
)

def main():
    """Generate debug plots for electron initialization."""
    print("\n" + "="*70)
    print("ELECTRON INITIALIZATION DEBUG PLOTS")
    print("="*70)
    print("\nThis script verifies:")
    print("  1. Cross-section has correct double-lobe dumbbell shape")
    print("  2. Phase wraps smoothly around torus (no branch cut jumps)")
    print("  3. m-winding pattern is correctly implemented\n")

    # Physical constants
    constants = PhysicalConstants()

    # Grid parameters (from typical experiment setup)
    h_phys = 1e-15  # 1 fm grid spacing
    center = (0.0, 0.0, 0.0)  # Centered at origin

    print("Creating electron initialization parameters...")
    params = calibrate_electron_init_params(
        constants=constants,
        grid_spacing=h_phys,
        center=center,
        amplitude_scale=1e-13,  # Initial guess
    )

    print("\n" + "="*70)
    print("GENERATING DEBUG PLOTS")
    print("="*70)

    # Plot 1: Cross-section envelope
    print("\n1. Cross-section envelope f(x, y)")
    print("   Expected: Two Gaussian lobes separated along y-axis (binormal)")
    print("   Location: y = ±ρ₀")
    plot_cross_section_envelope_debug(
        params,
        n_points=200,
        save_path='debug_cross_section.png'
    )

    # Plot 2: Phase along centerline
    print("\n2. Phase pattern along torus centerline")
    print(f"   Expected: {params.winding_number} complete oscillations")
    print("   Should be: Smooth everywhere, no discontinuities")
    plot_phase_along_centerline_debug(
        params,
        n_points=500,
        save_path='debug_phase_centerline.png'
    )

    print("\n" + "="*70)
    print("DEBUG PLOTS COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print("  - debug_cross_section.png    : Cross-section envelope f(x,y)")
    print("  - debug_phase_centerline.png : Phase pattern along centerline")
    print("\nNext steps:")
    print("  1. Review debug_cross_section.png:")
    print("     ✓ Should see two clear lobes at y = ±ρ₀")
    print("     ✓ Lobes should be centered at x = 0")
    print("     ✓ Should form 'dumbbell' or 'figure-8' shape")
    print("\n  2. Review debug_phase_centerline.png:")
    print(f"     ✓ Should see {params.winding_number} complete cycles")
    print("     ✓ Phase should vary smoothly (no jumps > 0.1 rad)")
    print("     ✓ Amplitude field should be sinusoidal")
    print("\n  3. If plots look correct:")
    print("     → Run full electron experiment")
    print("     → Check that XZ/YZ slices show double-lobe structure")
    print("     → Verify lateral distortion is ~zero at t=0\n")


if __name__ == "__main__":
    main()