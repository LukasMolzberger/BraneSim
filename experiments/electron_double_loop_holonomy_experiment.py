"""
Electron double-loop Wilson loop holonomy experiment.

This experiment demonstrates:
1. Initialization of electron-like double-loop tube geometry
2. Sampling of degenerate 2D polarization frames along the loop
3. Computation of non-Abelian Wilson loop W ∈ U(2)
4. Verification of spinorial signature: W ≈ -I (one circuit), W² ≈ I (two circuits)

Paper reference:
    Section "Preparation-first initialization of narrowband carriers"
    Paragraph "Electron packet (double-loop tube, spinorial transport)"
    Section "Holonomy measurement", paragraph "Degenerate transport (U(N) / Wilczek--Zee)"
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.config.physical_constants import PhysicalConstants
from branesim.initialization.carrier_packets import make_electron_double_loop_packet
from branesim.diagnostics.holonomy import (
    wilson_loop_holonomy,
    verify_gauge_invariance,
)
from branesim.diagnostics.degeneracy import subspace_rank_svd


def run_electron_holonomy_experiment(
    output_dir: str = "outputs/electron_holonomy",
    n_grid: int = 40,  # Grid points per dimension
    torus_major_radius_factor: float = 5.0,  # R / lambda_C
    device: str = "cpu",
):
    """
    Run electron double-loop Wilson loop experiment.

    NOTE: This is a proof-of-concept demonstration of the Wilson loop methodology.
    Full 3D toroidal dynamics would require substantial computational resources
    and is deferred to future numerical work.

    Args:
        output_dir: Directory for outputs
        n_grid: Grid points per dimension (3D grid will be n_grid³)
        torus_major_radius_factor: Major radius in units of λ_C
        device: torch device ('cpu' or 'cuda')
    """
    print("\n" + "="*80)
    print("ELECTRON DOUBLE-LOOP WILSON LOOP HOLONOMY EXPERIMENT")
    print("="*80)
    print("\nExperiment goals:")
    print("  1. Initialize electron-like double-loop tube (half-angle transport)")
    print("  2. Sample degenerate 2D frames along the loop")
    print("  3. Compute Wilson loop W ∈ U(2)")
    print("  4. Verify spinorial signature: W ≈ -I (one circuit)")
    print("="*80 + "\n")

    print("NOTE: This is a proof-of-concept demonstration.")
    print("Full 3D toroidal dynamics deferred to future work.\n")

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Physical constants
    c = PhysicalConstants.c
    lambda_C = PhysicalConstants.lambda_C
    omega_C = PhysicalConstants.omega_C
    m_e = PhysicalConstants.m_e

    k0 = 2 * np.pi / lambda_C
    omega0 = omega_C

    # Grid setup (3D for toroidal geometry)
    h = lambda_C / 5.0  # Grid spacing
    R_major = torus_major_radius_factor * lambda_C
    r_minor = 2 * lambda_C  # Tube radius

    device_torch = torch.device(device)
    grid = BraneGrid(
        grid_shape=(n_grid, n_grid, n_grid),
        dimension=Dimensionality.THREE_D,
        spacing=h,
        device=device_torch,
    )

    # Initialize state
    n_points = n_grid ** 3
    state = BraneState(
        num_points=n_points,
        dimensionality=Dimensionality.THREE_D,
        device=device_torch,
    )

    print(f"\nGrid setup:")
    print(f"  Grid shape:             ({n_grid}, {n_grid}, {n_grid})")
    print(f"  Total points:           {n_points}")
    print(f"  Grid spacing h:         {h:.6e} m ({h/lambda_C:.2f} λ_C)")
    print(f"  Torus major radius R:   {R_major:.6e} m ({R_major/lambda_C:.1f} λ_C)")
    print(f"  Torus minor radius r0:  {r_minor:.6e} m ({r_minor/lambda_C:.1f} λ_C)")

    # Initialize electron double-loop packet
    amplitude_scale = lambda_C / np.sqrt(np.pi)
    torus_center = torch.tensor([
        n_grid * h / 2,
        n_grid * h / 2,
        n_grid * h / 2,
    ], device=device_torch)

    print("\nInitializing electron double-loop packet...")
    make_electron_double_loop_packet(
        state=state,
        grid=grid,
        k0=k0,
        omega0=omega0,
        amplitude=amplitude_scale,
        torus_major_radius=R_major,
        torus_minor_radius=r_minor,
        twist_winding=1,  # One full rotation per loop
        longitudinal_winding=2,  # Mode winds twice per loop
        dof_pair=(2, 3),  # X³ and X⁴ for internal C² representation
        torus_axis=2,  # Aligned with z-axis
        torus_center=torus_center,
    )

    # Sample frames along the loop
    # For a torus aligned with z-axis, loop coordinate is the azimuthal angle φ
    n_samples = 32  # Sample points along loop
    phi_samples = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)

    print(f"\nSampling {n_samples} frames along the loop...")

    frames = []

    for idx, phi in enumerate(phi_samples):
        # Find lattice points near the centerline at this azimuthal angle
        # Centerline: C(φ) = (R cos φ, R sin φ, 0) + center
        x_center = R_major * np.cos(phi) + torus_center[0].item()
        y_center = R_major * np.sin(phi) + torus_center[1].item()
        z_center = torus_center[2].item()

        # Get spatial coordinates of all points
        coords = grid.get_spatial_coordinates()  # [N, 3]

        # Find points within tube radius of this centerline point
        distances = torch.sqrt(
            (coords[:, 0] - x_center) ** 2 +
            (coords[:, 1] - y_center) ** 2 +
            (coords[:, 2] - z_center) ** 2
        )

        # Select points within the tube
        tube_mask = distances < (1.5 * r_minor)
        tube_indices = torch.where(tube_mask)[0]

        if len(tube_indices) < 10:
            print(f"  WARNING: Sample {idx}: only {len(tube_indices)} points found in tube")
            continue

        # Extract 2D frame: [u1, u2] from DOFs (2, 3) at these points
        dof_3 = state.positions[tube_indices, 2]  # [n_tube]
        dof_4 = state.positions[tube_indices, 3]  # [n_tube]

        # Build frame matrix [n_tube, 2]
        frame_data = torch.stack([dof_3, dof_4], dim=1)

        # SVD to get orthonormal 2D subspace
        # frame_data = U Σ V^T, take first 2 columns of U
        try:
            U, S, Vt = torch.linalg.svd(frame_data, full_matrices=False)
            frame_ortho = U[:, :2]  # [n_tube, 2]
        except Exception as e:
            print(f"  WARNING: Sample {idx}: SVD failed ({e})")
            continue

        frames.append(frame_ortho)

        if (idx + 1) % 8 == 0:
            print(f"  Sampled {idx+1}/{n_samples} frames")

    n_frames_actual = len(frames)
    print(f"  Total frames collected: {n_frames_actual}")

    if n_frames_actual < 8:
        print("\n⚠ WARNING: Too few frames collected. Results may not be reliable.")
        print("Consider increasing grid resolution or torus radius.\n")

    # Compute Wilson loop
    print("\n" + "="*80)
    print("WILSON LOOP COMPUTATION")
    print("="*80)

    wilson_result = wilson_loop_holonomy(
        frames=frames,
        reorthonormalize=True,
        gauge_check=True,
    )

    # Analyze results
    print("\n" + "="*80)
    print("SPINORIAL SIGNATURE ANALYSIS")
    print("="*80)

    W = wilson_result.W
    W_squared = W @ W

    # Distance of W² to +I
    I = np.eye(2)
    dist_W2_to_I = np.linalg.norm(W_squared - I, ord='fro') / np.linalg.norm(I, ord='fro')

    print(f"\nWilson loop W:")
    print(f"  W = \n{W}")
    print(f"\nW² (two circuits):")
    print(f"  W² = \n{W_squared}")
    print(f"\nTarget signatures:")
    print(f"  One circuit:  W ≈ -I    →  ||W + I|| / ||I|| = {wilson_result.distance_to_minus_I:.6f}")
    print(f"  Two circuits: W² ≈ I    →  ||W² - I|| / ||I|| = {dist_W2_to_I:.6f}")

    # Success criteria
    success_one_circuit = wilson_result.distance_to_minus_I < 0.2
    success_two_circuits = dist_W2_to_I < 0.2

    print(f"\nVerification:")
    if success_one_circuit:
        print(f"  ✓ One circuit signature (W ≈ -I): PASS")
    else:
        print(f"  ✗ One circuit signature (W ≈ -I): FAIL")

    if success_two_circuits:
        print(f"  ✓ Two circuit signature (W² ≈ I): PASS")
    else:
        print(f"  ✗ Two circuit signature (W² ≈ I): FAIL")

    if success_one_circuit and success_two_circuits:
        print(f"\n✓ SUCCESS: Spinorial (4π) periodicity confirmed!")
    else:
        print(f"\n⚠ PARTIAL: Signatures not fully met (may need higher resolution)")

    # Gauge invariance check
    print("\n" + "="*80)
    print("GAUGE INVARIANCE VERIFICATION")
    print("="*80)

    invariant_check = verify_gauge_invariance(
        frames=frames,
        W_original=W,
        n_random_tests=3,
    )

    # Save results
    results = {
        "grid_shape": (n_grid, n_grid, n_grid),
        "spacing": h,
        "torus_major_radius": R_major,
        "torus_minor_radius": r_minor,
        "n_samples": n_samples,
        "n_frames_collected": n_frames_actual,
        "wilson_loop": {
            "W": W.tolist(),
            "W_squared": W_squared.tolist(),
            "trace": float(wilson_result.trace.real),
            "eigenvalues": wilson_result.eigenvalues.tolist(),
            "eigenphases": wilson_result.eigenphases.tolist(),
            "distance_to_minus_I": wilson_result.distance_to_minus_I,
            "distance_to_plus_I": wilson_result.distance_to_plus_I,
            "W2_distance_to_I": dist_W2_to_I,
            "is_spinorial": wilson_result.is_spinorial,
        },
        "success": {
            "one_circuit": success_one_circuit,
            "two_circuits": success_two_circuits,
            "gauge_invariance": invariant_check,
        },
    }

    np.save(f"{output_dir}/results.npy", results, allow_pickle=True)
    print(f"\nResults saved to {output_dir}/results.npy")

    # Plot eigenphases of W
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    eigenphases = wilson_result.eigenphases

    # Plot on unit circle
    theta = np.linspace(0, 2 * np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k--', alpha=0.3, label="Unit circle")

    # Plot eigenvalues
    for i, phase in enumerate(eigenphases):
        x = np.cos(phase)
        y = np.sin(phase)
        ax.plot([0, x], [0, y], 'o-', markersize=10, linewidth=2, label=f"λ_{i+1}")
        ax.text(1.1*x, 1.1*y, f"e^{{i {phase:.3f}}}", fontsize=10)

    # Mark target: -1 (π phase)
    ax.plot([-1], [0], 'rx', markersize=15, markeredgewidth=3, label="Target (-1)")

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.set_xlabel("Real")
    ax.set_ylabel("Imaginary")
    ax.set_title("Wilson Loop Eigenvalues (one circuit)\nTarget: both at -1 for W ≈ -I")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/wilson_eigenvalues.png", dpi=150)
    print(f"Plot saved to {output_dir}/wilson_eigenvalues.png")

    print("\n" + "="*80)
    print("EXPERIMENT COMPLETE")
    print("="*80)

    return results


if __name__ == "__main__":
    results = run_electron_holonomy_experiment()