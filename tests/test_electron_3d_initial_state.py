"""Test: visualize electron initialization on a 3D brane in 4D (no simulation).

This script seeds an electron-like rotating spherical-harmonic mode using
**spatial polarization** (so the excitation is not only in X^4) and exports
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

    mapper = DimensionalMapper(h_phys=h_phys, c_light=constants.c, mass_reference=m_point)

    # Simulation units
    h_sim = mapper.to_sim_length(h_phys)  # == 1.0

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
    center = (
        (nx - 1) * h_sim / 2.0,
        (ny - 1) * h_sim / 2.0,
        (nz - 1) * h_sim / 2.0,
    )

    # Choose radius ~ one reduced Compton wavelength (in sim units)
    radius_sim = float(mapper.to_sim_length(constants.lambda_C))

    spec = ElectronModeSpec(
        l=1,
        m=1,
        n=1,
        radius=radius_sim,
        amplitude=0.5,
        center=center,
        wave_speed=1.0,
        # key change: excite a spatial polarization plane (not only X^4)
        polarization="spatial_x4",
        smooth_edge=2.0,
    )

    debug = initialize_electron_mode_3d(state=state, grid=grid, spec=spec, return_debug=True)
    print("\nElectron initialization:")
    for k, v in (debug or {}).items():
        if isinstance(v, np.ndarray):
            print(f"  {k}: {v}")
        else:
            print(f"  {k}: {v}")

    # Sanity: show that multiple embedding components are excited at t=0
    disp = (state.positions - initial_positions)
    max_per_comp = torch.max(torch.abs(disp), dim=0).values.detach().cpu().numpy()
    print("\nMax |displacement| per embedding component:")
    for i, mv in enumerate(max_per_comp.tolist()):
        print(f"  comp {i}: {mv:.6e} (sim)")

    # If we have p1/p2, compute the local quadratures a,b and amplitude
    p1 = torch.tensor(debug["p1"], device=device, dtype=state.dtype)
    p2 = torch.tensor(debug["p2"], device=device, dtype=state.dtype)
    a = torch.sum(disp * p1[None, :], dim=1)
    b = torch.sum(disp * p2[None, :], dim=1)
    amp = torch.sqrt(a * a + b * b)

    print(f"\nQuadrature stats (sim units):")
    print(f"  max|a| = {float(torch.max(torch.abs(a))):.6e}")
    print(f"  max|b| = {float(torch.max(torch.abs(b))):.6e}")
    print(f"  max amp = {float(torch.max(amp)):.6e}")

    # ------------------------------------------------------------------
    # Export static plots / CSVs
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
    # We show the *electron quadrature a* (signed) and also X^4 containment.
    # ------------------------------------------------------------------
    num_frames = 180
    times = list(np.linspace(0.0, 1.0, num_frames))

    # 1) Electron quadrature a (signed)
    a_pm = mapper.to_phys_length(a).detach().cpu().numpy() * 1e12  # m -> pm
    coords, values = subsample_3d_field(
        field_3d=a_pm,
        grid_shape=(nx, ny, nz),
        h_phys=h_phys * 1e12,  # m -> pm
        subsample_factor=2,
    )
    frames_data = [(coords, values) for _ in range(num_frames)]
    out_path = run_manager.get_plot_path("electron_initial_orbit_quadrature_a.mp4")
    print(f"\nRendering 3D orbit video (quadrature a): {out_path}")
    create_3d_animation(
        frames_data=frames_data,
        times=times,
        output_path=out_path,
        cmap_name='RdBu_r',
        point_size=6.0,
        gamma=1.0,
        alpha_scale=1.0,
        xlabel='x [pm]',
        ylabel='y [pm]',
        zlabel='z [pm]',
        title_template="Electron initial state (static, quadrature a) - t = {:.3f} as",
        fps=30,
        dpi=120,
        camera_motion=camera_orbit,
    )

    # 2) X^4 containment displacement
    x4_disp = disp[:, 3]
    x4_pm = mapper.to_phys_length(x4_disp).detach().cpu().numpy() * 1e12  # m -> pm
    coords2, values2 = subsample_3d_field(
        field_3d=x4_pm,
        grid_shape=(nx, ny, nz),
        h_phys=h_phys * 1e12,  # m -> pm
        subsample_factor=2,
    )
    frames_data2 = [(coords2, values2) for _ in range(num_frames)]
    out_path2 = run_manager.get_plot_path("electron_initial_orbit_x4_containment.mp4")
    print(f"\nRendering 3D orbit video (X^4 containment): {out_path2}")
    create_3d_animation(
        frames_data=frames_data2,
        times=times,
        output_path=out_path2,
        title_template="Electron initial state (static, X^4) - t = {:.3f} as",
        fps=30,
        dpi=120,
        camera_motion=camera_orbit,
    )

    print("\n✓ Done")


if __name__ == "__main__":
    main()
