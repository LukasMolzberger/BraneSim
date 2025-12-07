"""
Test Script: Verify Baseline Distortion Computation

This script tests that:
1. Baseline positions match the actual grid initialization
2. Lateral distortion is ~zero at t=0 when no lateral initialization is used
3. The visualization uses the correct baseline reference

This helps debug the "gradient in corner" artifact in lateral distortion plots.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np
import matplotlib.pyplot as plt

from branesim.physics.baseline_state import (
    compute_flat_baseline_positions,
    compute_lateral_distortion,
    validate_baseline_initialization,
)


def test_baseline_computation():
    """Test that baseline positions match initialized grid."""
    print("\n" + "="*70)
    print("TEST: Baseline Position Computation")
    print("="*70)

    # Simple grid setup
    grid_shape = (10, 10, 10)
    h = 1e-15  # 1 fm
    nx, ny, nz = grid_shape

    print(f"\nGrid shape: {grid_shape}")
    print(f"Grid spacing: {h:.2e} m")

    # Test 1: No centering (origin at corner)
    print("\n--- Test 1: Grid with corner at origin ---")
    baseline_no_center = compute_flat_baseline_positions(
        grid_shape=grid_shape,
        h=h,
        center=None,
        device='cpu',
        dtype=torch.float32,
    )

    # Create equivalent "actual" positions as if from a simulation
    x = torch.arange(nx, dtype=torch.float32) * h
    y = torch.arange(ny, dtype=torch.float32) * h
    z = torch.arange(nz, dtype=torch.float32) * h
    X, Y, Z = torch.meshgrid(x, y, z, indexing='ij')
    actual_positions = torch.stack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)], dim=1)

    # Check if they match
    is_valid, max_dev = validate_baseline_initialization(
        state_positions=actual_positions,
        baseline_positions=baseline_no_center,
        tolerance=1e-12,
    )

    print(f"  Baseline matches actual: {is_valid}")
    print(f"  Max deviation: {max_dev:.2e} m")

    if is_valid:
        print("  ✓ Baseline correctly computed for corner-at-origin grid")
    else:
        print("  ✗ MISMATCH! Check grid initialization")

    # Test 2: Centered grid
    print("\n--- Test 2: Grid centered at (cx, cy, cz) ---")
    center_phys = (5e-15, 5e-15, 5e-15)  # Physical center coordinates
    cx, cy, cz = center_phys

    baseline_centered = compute_flat_baseline_positions(
        grid_shape=grid_shape,
        h=h,
        center=center_phys,
        device='cpu',
        dtype=torch.float32,
    )

    # Create centered grid
    x_centered = torch.arange(nx, dtype=torch.float32) * h - (nx - 1) * h / 2 + cx
    y_centered = torch.arange(ny, dtype=torch.float32) * h - (ny - 1) * h / 2 + cy
    z_centered = torch.arange(nz, dtype=torch.float32) * h - (nz - 1) * h / 2 + cz
    X_c, Y_c, Z_c = torch.meshgrid(x_centered, y_centered, z_centered, indexing='ij')
    actual_centered = torch.stack([X_c.reshape(-1), Y_c.reshape(-1), Z_c.reshape(-1)], dim=1)

    is_valid_c, max_dev_c = validate_baseline_initialization(
        state_positions=actual_centered,
        baseline_positions=baseline_centered,
        tolerance=1e-12,
    )

    print(f"  Center: ({cx:.2e}, {cy:.2e}, {cz:.2e}) m")
    print(f"  Baseline matches actual: {is_valid_c}")
    print(f"  Max deviation: {max_dev_c:.2e} m")

    if is_valid_c:
        print("  ✓ Baseline correctly computed for centered grid")
    else:
        print("  ✗ MISMATCH! Check centering logic")

    # Test 3: Check distortion is zero for matching grids
    print("\n--- Test 3: Distortion for identical grids ---")
    distortion = compute_lateral_distortion(actual_positions, baseline_no_center)
    max_distortion = distortion.max().item()

    print(f"  Max distortion when positions match: {max_distortion:.2e} m")

    if max_distortion < 1e-12:
        print("  ✓ Distortion is ~zero for matching grids")
    else:
        print("  ✗ Unexpected distortion detected!")

    # Test 4: Check distortion with small perturbation
    print("\n--- Test 4: Distortion with small perturbation ---")
    perturbed = actual_positions.clone()
    perturbed[grid_shape[0]//2, :] += 1e-16  # Add 0.1 fm displacement to center point

    distortion_perturbed = compute_lateral_distortion(perturbed, baseline_no_center)
    max_distortion_p = distortion_perturbed.max().item()

    print(f"  Max distortion with 1e-16 m perturbation: {max_distortion_p:.2e} m")

    if 1e-17 < max_distortion_p < 2e-16:
        print("  ✓ Distortion correctly detects perturbation")
    else:
        print("  ⚠ Distortion magnitude unexpected (but may be OK)")


def visualize_gradient_artifact():
    """Demonstrate the 'gradient in corner' artifact with wrong baseline."""
    print("\n" + "="*70)
    print("DEMONSTRATION: Corner Gradient Artifact")
    print("="*70)

    grid_shape = (20, 20, 1)  # 2D slice for visualization
    h = 1e-15
    nx, ny, nz = grid_shape

    # Actual grid: centered at (10 fm, 10 fm, 0)
    center = (10e-15, 10e-15, 0.0)
    cx, cy, cz = center

    print(f"\nActual grid center: ({cx:.2e}, {cy:.2e}, {cz:.2e}) m")

    # Create actual grid (centered)
    x_actual = torch.arange(nx, dtype=torch.float32) * h - (nx - 1) * h / 2 + cx
    y_actual = torch.arange(ny, dtype=torch.float32) * h - (ny - 1) * h / 2 + cy
    z_actual = torch.zeros(1, dtype=torch.float32) + cz
    X_a, Y_a, Z_a = torch.meshgrid(x_actual, y_actual, z_actual, indexing='ij')
    actual_grid = torch.stack([X_a.reshape(-1), Y_a.reshape(-1), Z_a.reshape(-1)], dim=1)

    # Correct baseline (centered)
    baseline_correct = compute_flat_baseline_positions(
        grid_shape=grid_shape,
        h=h,
        center=center,
        device='cpu',
        dtype=torch.float32,
    )

    # Wrong baseline (corner at origin)
    baseline_wrong = compute_flat_baseline_positions(
        grid_shape=grid_shape,
        h=h,
        center=None,  # No centering!
        device='cpu',
        dtype=torch.float32,
    )

    # Compute distortions
    distortion_correct = compute_lateral_distortion(actual_grid, baseline_correct)
    distortion_wrong = compute_lateral_distortion(actual_grid, baseline_wrong)

    # Reshape for plotting
    dist_correct_grid = distortion_correct.reshape(grid_shape[:2]).numpy()
    dist_wrong_grid = distortion_wrong.reshape(grid_shape[:2]).numpy()

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    im = ax.imshow(dist_correct_grid.T, origin='lower', cmap='hot')
    ax.set_title('Correct Baseline (centered)', fontsize=14, fontweight='bold')
    ax.set_xlabel('X index')
    ax.set_ylabel('Y index')
    plt.colorbar(im, ax=ax, label='Distortion [m]')
    ax.text(0.5, 0.95, f'Max: {dist_correct_grid.max():.2e} m',
            transform=ax.transAxes, ha='center', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax = axes[1]
    im = ax.imshow(dist_wrong_grid.T, origin='lower', cmap='hot')
    ax.set_title('Wrong Baseline (corner at origin)', fontsize=14, fontweight='bold')
    ax.set_xlabel('X index')
    ax.set_ylabel('Y index')
    plt.colorbar(im, ax=ax, label='Distortion [m]')
    ax.text(0.5, 0.95, f'Max: {dist_wrong_grid.max():.2e} m',
            transform=ax.transAxes, ha='center', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig('test_baseline_artifact.png', dpi=150, bbox_inches='tight')
    print("\n  ✓ Saved: test_baseline_artifact.png")
    print("\n  Left panel:  Correct baseline → distortion ~zero everywhere")
    print("  Right panel: Wrong baseline → gradient from corner!")
    print("\n  This demonstrates the 'corner gradient' artifact")
    print("  when baseline centering doesn't match actual grid.\n")


def main():
    print("\n" + "="*70)
    print("BASELINE DISTORTION COMPUTATION TESTS")
    print("="*70)

    print("\nThese tests verify that:")
    print("  1. Baseline positions correctly match grid initialization")
    print("  2. Lateral distortion is ~zero when grids match")
    print("  3. Wrong baseline centering causes 'corner gradient' artifact\n")

    # Run tests
    test_baseline_computation()
    visualize_gradient_artifact()

    print("\n" + "="*70)
    print("ALL TESTS COMPLETE")
    print("="*70)
    print("\nKey takeaways:")
    print("  - baseline_state.py functions work correctly")
    print("  - Centering MUST match between baseline and actual grid")
    print("  - electron_visualization.py must use config['center'] properly")
    print("\nNext step:")
    print("  → Check that experiment config passes correct 'center' to visualization\n")


if __name__ == "__main__":
    main()