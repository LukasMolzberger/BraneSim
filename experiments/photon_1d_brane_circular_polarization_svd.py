"""
Photon circular polarization experiment with degeneracy verification.

This experiment demonstrates the preparation-first initialization approach:
1. Initialize a circularly polarized photon packet (narrowband by construction)
2. Verify 2D polarization subspace via SVD (no FFT filtering)
3. Propagate with periodic boundary conditions
4. Optional FFT cross-check (not required for main results)

Paper reference:
    Section "Preparation-first initialization of narrowband carriers"
    Paragraph "Photon packet (two-dimensional polarization subspace)"
"""

import torch
import numpy as np
import matplotlib.pyplot as plt

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.core.dimensions import MassModel
from branesim.physics.forces import SpringForceComputer
from branesim.config.physical_constants import PhysicalConstants
from branesim.initialization.carrier_packets import make_photon_circular_packet
from branesim.analytics.degeneracy import (
    subspace_rank_svd,
    verify_narrowband_preparation,
)
from branesim.utils import TestRunManager


def run_photon_polarization_experiment(
    nx: int = 200,
    n_cycles: int = 5,
    amplitude_scale: float = 1e-14,  # meters (~ lambda_C / sqrt(pi))
    device: str = "cpu",
):
    """
    Run photon circular polarization experiment.

    Args:
        nx: Number of grid points
        n_cycles: Number of Compton cycles to simulate
        amplitude_scale: Amplitude A (default: ~lambda_C/√π)
        device: torch device ('cpu' or 'cuda')
    """
    print("\n" + "="*80)
    print("PHOTON CIRCULAR POLARIZATION EXPERIMENT")
    print("="*80)
    print("\nExperiment goals:")
    print("  1. Initialize circularly polarized packet (preparation-first)")
    print("  2. Verify 2D polarization subspace via SVD (no FFT)")
    print("  3. Propagate with periodic boundary conditions")
    print("  4. Optional FFT cross-check")
    print("="*80 + "\n")

    # Initialize test run manager
    run_manager = TestRunManager(experiment_name="photon_circular_polarization_svd")
    print(run_manager.get_summary())

    # Physical constants (Compton scale)
    constants = PhysicalConstants()
    c = constants.c
    lambda_C = constants.lambda_C
    m_e = constants.m_e

    # Compute Compton frequency: ω_C = m_e c² / ℏ
    omega_C = m_e * c**2 / constants.hbar

    k0 = 2 * np.pi / lambda_C
    omega0 = omega_C

    # Grid setup with PERIODIC boundary conditions
    h = lambda_C / 10.0  # Grid spacing (10 points per wavelength)
    L = nx * h  # Domain length

    device_torch = torch.device(device)
    grid = BraneGrid(
        grid_shape=(nx,),
        dimension=Dimensionality.ONE_D,
        spacing=h,
        device=device_torch,
        periodic_axes=(True,),  # PERIODIC in x
    )

    # Initialize state
    state = BraneState(
        grid_shape=(nx,),
        dimension=Dimensionality.ONE_D,
        device=device_torch,
    )

    # Calibrate mass density and spring stiffness
    # For 1D chain: need linear density ρ₁ (kg/m)
    # Compton cell mass: m_e ≈ ρ₁ * λ_C  (for 1D)
    # Wave speed: c² = T/ρ₁
    rho_1d = m_e / lambda_C  # Linear density (kg/m) for 1D chain
    T = rho_1d * c ** 2       # Tension (N)

    # Initialize circular polarization packet
    # Packet width: 3 wavelengths (narrowband)
    sigma = 3 * lambda_C
    center = L / 2.0

    print("\nInitializing circularly polarized photon packet...")
    print(f"  Carrier wavelength λ_C: {lambda_C:.6e} m")
    print(f"  Carrier frequency ω_C:  {omega_C:.6e} rad/s")
    print(f"  Envelope width σ:       {sigma:.6e} m ({sigma/lambda_C:.1f} λ_C)")
    print(f"  Amplitude A:            {amplitude_scale:.6e} m")
    print(f"  Grid spacing h:         {h:.6e} m ({h/lambda_C:.2f} λ_C)")
    print(f"  Domain length L:        {L:.6e} m ({L/lambda_C:.1f} λ_C)")
    print(f"  Periodic BCs:           YES")

    make_photon_circular_packet(
        state=state,
        grid=grid,
        k0=k0,
        omega0=omega0,
        amplitude=amplitude_scale,
        center=center,
        sigma=sigma,
        dof_pair=(2, 3),  # X³ and X⁴ for polarization
        propagation_axis=0,
        periodic_axis=True,
    )

    # Verify narrowband preparation (optional FFT cross-check)
    print("\n" + "-"*80)
    narrowband_metrics = verify_narrowband_preparation(
        state=state,
        k0_expected=k0,
        omega0_expected=omega0,
        dofs=[2, 3],
        grid=grid,
    )
    print("-"*80)

    # Initialize physics (spring force computer)
    physics = SpringForceComputer(
        spring_constant=T / h,  # Spring constant from tension
        rest_length=h,  # Rest length equals spacing
    )

    # Create mass model with proper 1D linear density
    mass_model = MassModel.from_density(
        density=rho_1d,      # Linear density in kg/m
        intrinsic_dim=1,     # 1D chain
        spacing=h,
    )

    print(f"\nMass model:")
    print(f"  {mass_model}")

    # Time stepping
    dt = 0.1 * h / c  # CFL condition
    n_steps_per_cycle = int(2 * np.pi / (omega0 * dt))
    n_steps = n_cycles * n_steps_per_cycle

    print(f"\nTime integration settings:")
    print(f"  dt:                     {dt:.6e} s")
    print(f"  CFL number:             {c * dt / h:.4f}")
    print(f"  Steps per cycle:        {n_steps_per_cycle}")
    print(f"  Total steps:            {n_steps}")
    print(f"  Simulation time:        {n_steps * dt:.6e} s ({n_cycles} cycles)")

    # Solver
    solver = VelocityVerletSolver(dt, mass_model, physics, grid)

    # Time series recording for SVD
    record_interval = max(1, n_steps_per_cycle // 20)  # 20 samples per cycle
    time_series_dof3 = []
    time_series_dof4 = []
    times = []

    print("\nRunning simulation...")
    for step in range(n_steps):
        # Record for SVD
        if step % record_interval == 0:
            time_series_dof3.append(state.positions[:, 2].clone().cpu())
            time_series_dof4.append(state.positions[:, 3].clone().cpu())
            times.append(step * dt)

        # Time step
        solver.step(state)

        # Progress
        if (step + 1) % (n_steps // 10) == 0:
            percent = 100 * (step + 1) / n_steps
            print(f"  Progress: {percent:.1f}%")

    print("  Simulation complete.")

    # SVD analysis for degeneracy verification
    print("\n" + "="*80)
    print("DEGENERACY VERIFICATION (SVD of time series)")
    print("="*80)

    # Stack time series: [n_time, n_spatial, 2]
    ts_dof3 = torch.stack(time_series_dof3, dim=0)  # [n_time, nx]
    ts_dof4 = torch.stack(time_series_dof4, dim=0)  # [n_time, nx]
    time_series = torch.stack([ts_dof3, ts_dof4], dim=2)  # [n_time, nx, 2]

    # Take central region (avoid edges if any artifacts)
    nx_center = nx // 2
    center_slice = slice(nx_center - 20, nx_center + 20)
    ts_center = time_series[:, center_slice, :]  # [n_time, 40, 2]

    svd_result = subspace_rank_svd(
        time_series=ts_center,
        energy_threshold=0.99,
        gap_threshold=0.1,
    )

    print("\nSVD Verification Results:")
    print(f"  Dominant dimension:     {svd_result.dominant_dimension}")
    print(f"  Is 2D degenerate:       {'YES' if svd_result.is_degenerate_2d else 'NO'}")
    print(f"  Singular values:        {svd_result.singular_values[:5]}")

    # Check if we have enough energy ratios
    if len(svd_result.energy_ratios) >= 2:
        print(f"  Energy in top 2 modes:  {svd_result.energy_ratios[1]*100:.2f}%")
    elif len(svd_result.energy_ratios) == 1:
        print(f"  Energy in top 1 mode:   {svd_result.energy_ratios[0]*100:.2f}%")
    else:
        print(f"  WARNING: No energy ratios computed")

    if svd_result.is_degenerate_2d:
        print("\n✓ SUCCESS: 2D polarization subspace confirmed (no filtering used)")
    else:
        print("\n⚠ WARNING: Expected 2D subspace not clearly identified")

    # Save results
    results = {
        "grid_shape": (nx,),
        "spacing": h,
        "lambda_C": lambda_C,
        "omega_C": omega_C,
        "amplitude": amplitude_scale,
        "sigma": sigma,
        "n_cycles": n_cycles,
        "narrowband_metrics": narrowband_metrics,
        "svd_result": {
            "singular_values": svd_result.singular_values.tolist(),
            "dominant_dimension": svd_result.dominant_dimension,
            "is_degenerate_2d": svd_result.is_degenerate_2d,
            "gap_ratio": svd_result.gap_ratio,
        },
        "times": times,
    }

    results_path = run_manager.get_data_path('results.npy')
    np.save(results_path, results, allow_pickle=True)
    print(f"\nResults saved to {results_path}")

    # Plot final state
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    x_coords = np.arange(nx) * h / lambda_C  # In units of λ_C

    axes[0].plot(x_coords, state.positions[:, 2].cpu().numpy(), label="DOF 3 (cos)")
    axes[0].set_ylabel("X³ (m)")
    axes[0].set_title("Photon Circular Polarization - Final State")
    axes[0].legend()
    axes[0].grid(True)

    axes[1].plot(x_coords, state.positions[:, 3].cpu().numpy(), label="DOF 4 (sin)", color="orange")
    axes[1].set_xlabel("Position (λ_C)")
    axes[1].set_ylabel("X⁴ (m)")
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path('final_state.png'), dpi=150)
    print(f"Plot saved to final_state.png")

    # Save configuration
    run_manager.save_config({
        "experiment": "Photon Circular Polarization SVD",
        "nx": nx,
        "n_cycles": n_cycles,
        "amplitude_scale": amplitude_scale,
        "device": device,
    })

    print("\n" + "="*80)
    print("EXPERIMENT COMPLETE")
    print("="*80)
    print(f"All outputs saved to: {run_manager.run_dir}")

    return results


if __name__ == "__main__":
    results = run_photon_polarization_experiment()