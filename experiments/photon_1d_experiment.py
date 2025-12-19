"""
1D Photon with Realistic Physical Scales

Uses actual speed of light c = 299,792,458 m/s and physical length scales
based on the Compton wavelength.

This experiment is now a thin wrapper around the common photon_1d_runner
to eliminate code duplication with photon_1d_berry_phase_experiment.py
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from branesim.utils.photon_1d_runner import Photon1DConfig, run_photon_1d
from branesim.visualization import plot_all_photon_1d_standard


def main():
    """Run 1D photon simulation experiment."""
    # Configure simulation
    cfg = Photon1DConfig(
        num_steps=20000,
        num_snapshots=7,
        add_snapshot_step1=True,
        points_per_wavelength=20,
        num_wavelengths=100,
        amplitude_factor_h=10.0,
        center_fraction=1.0 / 3.0,
        export_csv_snapshots=True,
        experiment_name="photon_1d_experiment",
        cfl_factor=0.1,
    )

    # Run simulation
    run = run_photon_1d(cfg)

    # Generate all standard plots
    plot_all_photon_1d_standard(run)

    print(f"\n{'=' * 70}")
    print(f"All outputs saved to: {run.run_manager.run_dir}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()