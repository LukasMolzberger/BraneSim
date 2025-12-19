"""
1D Photon with Berry Phase Analysis

Extends the standard photon_1d_experiment.py with Berry phase computation.
Computes the discrete Berry phase profile γ(x) along the 1D chain at snapshot times.

Mathematical Background:
    For normalized states |u_i⟩ at lattice points i, the discrete Berry phase
    increment between neighbors is:
        Δφ_i = arg⟨u_i|u_{i+1}⟩

    The cumulative Berry phase profile is:
        γ_0 = 0,  γ_{i+1} = γ_i + Δφ_i

    The Berry connection (gauge field) is:
        A_x[i] = Δφ_i / h

    where h is the grid spacing.

Complex State Construction:
    The simulation fields are real (ξ = X^4 displacement and ξ̇ velocity).
    To extract phase information, we construct the analytic signal:
        ψ = ξ + i·ξ̇/ω

    where ω = 2π·c/λ is the carrier frequency of the photon.

This experiment is now a thin wrapper that:
1. Runs the common photon_1d simulation
2. Reuses standard photon plots
3. Adds Berry phase analysis and plots
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from branesim.utils.photon_1d_runner import Photon1DConfig, run_photon_1d
from branesim.visualization import (
    plot_all_brane_1d_standard,
    plot_berry_phase_profiles,
    plot_berry_connection_profiles,
)
from branesim.analysis import (
    complex_band_state_from_quadrature,
    pointwise_normalize,
    BerryPhase1DConfig,
    berry_phase_profile_along_x,
)
from branesim.utils.io import export_berry_phase_csv


def main():
    """Run 1D photon simulation with Berry phase analysis."""
    # Configure simulation (same as standard photon experiment)
    cfg = Photon1DConfig(
        num_steps=20000,
        num_snapshots=7,
        add_snapshot_step1=False,  # Berry phase doesn't need step 1
        points_per_wavelength=20,
        num_wavelengths=100,
        amplitude_factor_h=10.0,
        center_fraction=1.0 / 3.0,
        export_csv_snapshots=True,
        experiment_name="photon_1d_berry_phase_experiment",
        cfl_factor=0.1,
    )

    # Run simulation
    run = run_photon_1d(cfg)

    # Generate all standard photon plots (reuse visualization)
    plot_all_brane_1d_standard(run)

    # ========================================================================
    # Berry Phase Analysis (additional)
    # ========================================================================

    print(f"\nComputing Berry phase profiles...")

    # Berry phase configuration
    berry_cfg = BerryPhase1DConfig(
        spacing=float(run.h_sim),
        amplitude_threshold=1e-6,
        eps=1e-12,
        unwrap=True,
        force_cpu_on_mps=True,
    )

    print(f"\nBerry Phase Configuration:")
    print(f"  Carrier frequency ω_sim = {run.omega_sim:.6e} rad / sim-time")
    print(f"  Grid spacing h_sim = {run.h_sim:.6e}")
    print(f"  Amplitude threshold = {berry_cfg.amplitude_threshold:.6e}")

    # Storage for Berry phase data
    gamma_by_tfs = {}
    Ax_by_tfs = {}

    # Coordinates for plotting
    x_nm = run.x_coords_phys_m * 1e9  # m → nm
    x_edges_nm = 0.5 * (x_nm[:-1] + x_nm[1:])  # Edge midpoints

    # Compute Berry phase at each snapshot time
    for t_phys_s in run.snapshot_times_phys_s:
        t_fs = t_phys_s * 1e15

        # Get snapshots (sim units, numpy arrays)
        xi_sim = run.snapshots_xi[t_phys_s]
        v_xi_sim = run.snapshots_v_xi[t_phys_s]

        # Convert to torch tensors
        import torch
        xi_t = torch.from_numpy(xi_sim).to(run.state.device, run.state.dtype)
        v_xi_t = torch.from_numpy(v_xi_sim).to(run.state.device, run.state.dtype)

        # Build complex band state: ψ = ξ + i·ξ̇/ω
        psi = complex_band_state_from_quadrature(
            xi_t, v_xi_t, run.omega_sim, eps=berry_cfg.eps
        )

        # Normalize pointwise: |ψ̂⟩ = |ψ⟩ / |ψ|
        psi_hat, amp = pointwise_normalize(psi, eps=berry_cfg.eps)

        # Compute Berry phase profile along x
        result = berry_phase_profile_along_x(psi_hat, amp, berry_cfg)

        # Store results (move to CPU for numpy conversion)
        gamma_by_tfs[t_fs] = result["gamma_wrapped"].detach().cpu().numpy()
        Ax_by_tfs[t_fs] = result["A_x"].detach().cpu().numpy()

        # Export to CSV
        csv_path = run.run_manager.get_data_path(f"berry_phase_t_{t_fs:.3f}fs.csv")
        export_berry_phase_csv(
            csv_path,
            run.x_coords_phys_m,  # Position [m]
            result["gamma_wrapped"].detach().cpu().numpy(),  # Berry phase [rad]
            result["A_x"].detach().cpu().numpy(),  # Berry connection [rad/sim-length]
            amp_m=run.mapper.to_phys_length(amp).detach().cpu().numpy(),  # Amplitude [m]
        )

    print(f"  ✓ Computed Berry phase at {len(gamma_by_tfs)} snapshot times")
    print(f"  ✓ Exported CSV files (2 files per snapshot: points + edges)")

    # ========================================================================
    # Berry Phase Visualization
    # ========================================================================

    print(f"\nCreating Berry phase plots...")

    times_fs = sorted(gamma_by_tfs.keys())

    # Plot 1: Berry phase profiles γ(x)
    plot_berry_phase_profiles(
        run.run_manager,
        x_nm,
        times_fs,
        gamma_by_tfs,
        title="1D Photon - Berry Phase Profile γ(x) (wrapped to [-π, π])",
        filename="photon_1d_berry_phase_profiles.png",
    )
    print(f"  ✓ Saved: photon_1d_berry_phase_profiles.png")

    # Plot 2: Berry connection profiles A_x(x)
    plot_berry_connection_profiles(
        run.run_manager,
        x_edges_nm,
        times_fs,
        Ax_by_tfs,
        title="1D Photon - Berry Connection A_x(x)",
        filename="photon_1d_berry_connection_profiles.png",
    )
    print(f"  ✓ Saved: photon_1d_berry_connection_profiles.png")

    # ========================================================================
    # Summary
    # ========================================================================

    print(f"\n{'=' * 70}")
    print("Simulation complete!")
    print(f"{'=' * 70}")
    print(f"\nPhysical Interpretation:")
    print(f"  Domain size: {run.domain_length_phys*1e9:.3f} nm")
    print(f"  Wavelength: {run.wavelength_phys*1e9:.3f} nm")
    print(f"  Simulation time: {run.snapshot_times_phys_s[-1]*1e15:.3f} femtoseconds")
    print(f"  Berry phase snapshots: {len(times_fs)} times")

    print(f"\nBerry Phase Results:")
    print(f"  Carrier frequency ω = {run.omega_sim:.6e} rad/sim-time")
    print(f"  Amplitude threshold = {berry_cfg.amplitude_threshold:.6e}")
    print(f"  Grid spacing h = {run.h_sim:.6e} sim-length")

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run.run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()