"""Test: visualize electron initialization on a 3D brane in 4D (no simulation).

This script seeds an electron-like rotating spherical-harmonic mode and exports
static plots + a camera-orbit video of the initial state.

Run:
    python tests/test_electron_3d_initial_state.py
"""

import sys
import os

# Allow running from repo root
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.utils import TestRunManager
from branesim.visualization.brane_state_viz import visualize_brane_state
from branesim.visualization.brane_3d_viz import (
    subsample_3d_field,
    create_3d_animation,
    camera_orbit,
)

from branesim.electron.electron_initialization import ElectronModeSpec, initialize_electron_mode_3d


def main() -> None:
    run_manager = TestRunManager(experiment_name="electron_3d_initial_state")
    print(run_manager.get_summary())

    constants = PhysicalConstants()

    # Use the same scaling convention as the photon experiments
    points_per_lambdaC = 20
    h_phys = constants.lambda_C / points_per_lambdaC

    # Reference mass: keep consistent with other experiments
    m_point = 2.861821e-27  # kg

    # Wave speed (exactly equals c by construction)
    c_wave = constants.c

    mapper = DimensionalMapper(h_phys=h_phys, c_light=constants.c, mass_reference=m_point)

    # Simulation units
    h_sim = mapper.to_sim_length(h_phys)  # == 1.0
    c_sim = mapper.to_sim_velocity(c_wave)

    # Domain
    nx = ny = nz = 100

    # Device selection
    if torch.cuda.is_available():
        device = torch.device("cuda")
        dtype = torch.float64
        print(f"\n✓ Using NVIDIA GPU: {torch.cuda.get_device_name(0)}")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
        dtype = torch.float32
        print("\n✓ Using Apple Silicon GPU (MPS)")
        print("  Using float32 (MPS doesn't support float64)")
    else:
        device = torch.device("cpu")
        dtype = torch.float64
        print("\n⚠ Using CPU (no GPU detected)")

    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    # Fixed boundaries (same as photon experiments)
    state.set_fixed_boundaries()

    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h_sim, device)

    # Electron mode spec (in sim units)
    # Centered in the domain
    center = (
        (nx - 1) * h_sim / 2.0,
        (ny - 1) * h_sim / 2.0,
        (nz - 1) * h_sim / 2.0,
    )

    # Choose radius ~ reduced Compton wavelength (in sim units)
    radius_sim = float(mapper.to_sim_length(constants.lambda_C))

    spec = ElectronModeSpec(
        l=1, m=1, n=1,
        radius=radius_sim,
        amplitude=0.5,
        center=center,
        wave_speed=c_sim,

        polarization="xy",  # rotating spatial polarization in X^1–X^2
        containment_component=3,  # X^4 trap
        containment_depth=0.25,
        containment_sigma=0.5 * radius_sim,
        smooth_edge=2.0,
    )

    debug = initialize_electron_mode_3d(state=state, grid=grid, spec=spec, return_debug=True)
    print("\nElectron initialization:")
    for k, v in (debug or {}).items():
        print(f"  {k}: {v}")

    # ------------------------------------------------------------------
    # Export static slice plots / CSVs
    # ------------------------------------------------------------------
    visualize_brane_state(
        state=state,
        grid=grid,
        mapper=mapper,
        output_dir=str(run_manager.plots_dir),
        filename_prefix="initial",
        initial_positions=initial_positions,
        print_stats=True,
        dpi=150,
        csv_output_dir=str(run_manager.data_dir),
    )

    # ------------------------------------------------------------------
    # 3D point cloud orbit video (static field, moving camera)
    # ------------------------------------------------------------------
    field_phys = mapper.to_phys_length(state.positions[:, spec.field_component]).detach().cpu().numpy()
    coords, values = subsample_3d_field(
        field_3d=field_phys,
        grid_shape=(nx, ny, nz),
        h_phys=h_phys,
        subsample_factor=2,
    )

    # Reuse the same frame values, but orbit the camera
    num_frames = 180
    frames_data = [(coords, values) for _ in range(num_frames)]
    times = list(np.linspace(0.0, 1.0, num_frames))
    out_path = run_manager.get_plot_path("electron_initial_orbit.mp4")

    print(f"\nRendering 3D orbit video: {out_path}")
    create_3d_animation(
        frames_data=frames_data,
        times=times,
        output_path=out_path,
        title_template="Electron initial state (static) - t = {:.3f} fs",
        fps=30,
        dpi=120,
        camera_motion=camera_orbit,
        gamma=0.5,  # Lower gamma for better visibility (less aggressive opacity mapping)
        alpha_scale=0.9,  # Higher overall opacity
        point_size=8.0,  # Larger points for better visibility
    )

    print("\n✓ Done")


if __name__ == "__main__":
    main()
