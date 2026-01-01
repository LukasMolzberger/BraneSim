"""tests/test_directional_hold_force_physical.py

Directional hold-force measurement (X^4 vs. lateral) with **physical outputs**.

- The brane model uses **SI units internally** (meters, seconds, Newtons).
- Displacements are reported in **nanometers [nm]** for readability.
- Forces are reported in **Newtons [N]**.

Protocol (per direction)
------------------------
1) Start from a perfectly flat brane (fresh state).
2) Quasi-statically ramp the *center node* displacement along a chosen
   embedding-space direction e_dir.
3) After each increment, relax for a fixed number of steps while keeping the
   center node clamped, then measure the internal force at the node.
4) The required holding force is:

        F_hold = -(F_internal(center) · e_dir)

Visualization
-------------
- Two static plots: holding-force curve (X^4) and (lateral)
- Two MP4 videos using the project's **3D point cloud** visualization:
    * coords: deformed (x,y,z) positions in nm
    * color: displacement along e_dir in nm
    * title: shows the current clamp displacement Δ in nm

This test is intentionally opt-in.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import torch
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# TestRunManager
# -----------------------------------------------------------------------------
from branesim.utils.test_run_manager import TestRunManager


# -----------------------------------------------------------------------------
# Project imports
# -----------------------------------------------------------------------------
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.physics.forces import SpringForceComputer
from branesim.core.dimensions import MassModel
from branesim.core.solver import VelocityVerletSolver

from branesim.config.physical_constants import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper

# Point cloud visualization (used for videos)
from branesim.visualization.brane_3d_viz import create_3d_animation, camera_orbit


NM = 1e9  # meters -> nm


@dataclass
class PhysicalConfig:
    """Configuration.

    Notes
    -----
    The simulation is performed directly in SI units (m, s, N) using the project's
    physical parameter calibration. Outputs are reported in nm and N.
    """

    # Domain (lattice points)
    grid_shape: Tuple[int, int, int] = (33, 33, 33)

    # Physical scale choice: lattice spacing h_phys = lambda_C / points_per_lambdaC
    points_per_lambdaC: int = 20

    # Mass per lattice node (physical kg). Used to set the physical density rho = m_point / h^3.
    m_point: float = 2.861821e-27

    # Quasi-static displacement ramp
    d_max_in_spacing: float = 5.0
    ramp_increments: int = 180
    relax_steps_per_increment: int = 160

    # Numerics (CFL time step + simple damping)
    cfl_factor: float = 0.1
    damping_per_step: float = 0.02

    # Video (point cloud)
    fps: int = 20
    frame_stride: int = 2          # keep every Nth displacement step as a frame
    subsample_factor_3d: int = 2   # take every kth point along each intrinsic axis

    # Display units for plots/videos (default: femtometers)
    display_unit_name: str = "fm"
    display_unit_scale: float = 1e15  # meters -> display units

    # Visualization: always include a dense local patch around the center
    patch_radius: int = 8           # in intrinsic grid steps (Chebyshev radius)
    include_coarse_background: bool = True
    center_on_rest: bool = True

    # Point cloud appearance
    point_size: float = 5.5
    gamma: float = 1.1
    alpha_scale: float = 0.9
    min_alpha: float = 0.15


def _center_index_from_coords(grid_coords: torch.Tensor, grid_shape: Tuple[int, ...]) -> int:
    center = torch.tensor([s // 2 for s in grid_shape], device=grid_coords.device, dtype=grid_coords.dtype)
    matches = (grid_coords == center).all(dim=1)
    idx = torch.nonzero(matches, as_tuple=False).squeeze(1)
    if idx.numel() != 1:
        raise RuntimeError(f"Expected exactly one center index, got {idx.numel()}")
    return int(idx.item())


def _subsample_indices_3d(grid_coords: torch.Tensor, subsample_factor: int) -> torch.Tensor:
    """Return indices for points where each intrinsic coordinate is a multiple of subsample_factor."""
    if subsample_factor <= 1:
        return torch.arange(grid_coords.shape[0], device=grid_coords.device)

    mask = (
        (grid_coords[:, 0] % subsample_factor == 0)
        & (grid_coords[:, 1] % subsample_factor == 0)
        & (grid_coords[:, 2] % subsample_factor == 0)
    )
    return torch.nonzero(mask, as_tuple=False).squeeze(1)

def _viz_indices_3d(
    grid_coords: torch.Tensor,
    subsample_factor: int,
    center_idx: int,
    patch_radius: int,
    include_coarse_background: bool = True,
) -> torch.Tensor:
    """Indices for visualization: union of a dense local patch + optional coarse background."""
    center = grid_coords[center_idx]  # [3]
    diffs = (grid_coords - center).abs()
    local_mask = (diffs.max(dim=1).values <= patch_radius)

    if include_coarse_background and subsample_factor > 1:
        coarse = (
            (grid_coords[:, 0] % subsample_factor == 0)
            & (grid_coords[:, 1] % subsample_factor == 0)
            & (grid_coords[:, 2] % subsample_factor == 0)
        )
        mask = local_mask | coarse
    else:
        mask = local_mask

    return torch.nonzero(mask, as_tuple=False).squeeze(1)



def _clamp_center(state: BraneState, center_idx: int, target: torch.Tensor) -> None:
    """Hard clamp the center node to `target` and zero its kinematics."""
    state.positions[center_idx] = target
    state.velocities[center_idx] = 0.0
    state.accelerations[center_idx] = 0.0
    state.new_accelerations[center_idx] = 0.0


def _build_mapper_and_sim_params(cfg: PhysicalConfig) -> Dict[str, Any]:
    """Build **physical (SI)** parameters for the run.

    NOTE: The current branesim core uses SI units internally:
        - positions in meters
        - spring constants in N/m
        - rest lengths in meters
        - time step dt in seconds

    Therefore we compute everything directly in SI and do **not** apply any
    dimensionless mapping here.
    """
    constants = PhysicalConstants()

    # Physical lattice spacing
    h_phys = constants.lambda_C / cfg.points_per_lambdaC

    # Dimensional mapper for stable sim units
    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=cfg.m_point,
    )

    # Density corresponding to m_point on a 3D lattice cell (physical)
    rho_phys = cfg.m_point / (h_phys ** 3)

    # Rest length (physical)
    L0_phys = constants.compute_rest_length(h_phys)

    # Target tension T = rho * c^2 in physical units, map to spring constant k [N/m]
    T_target = constants.compute_target_tension(rho_phys)
    frac = constants.pre_stretch_alpha
    k_phys = float((T_target * h_phys) / max(1e-30, (1.0 - frac)))

    # Sim parameters
    h_sim = float(mapper.to_sim_length(h_phys))  # = 1.0
    rest_length_sim = float(mapper.to_sim_length(L0_phys))
    k_sim = float(mapper.to_sim_spring_constant(k_phys))
    m_sim = float(mapper.to_sim_mass(cfg.m_point))
    rho_sim = float(m_sim / (h_sim ** 3))

    # CFL time step in physical + sim units
    dt_phys = float(cfg.cfl_factor * h_phys / constants.c)
    dt_sim = float(mapper.to_sim_time(dt_phys))

    return {
        "constants": constants,
        "mapper": mapper,
        "h_phys": float(h_phys),
        "h_sim": float(h_sim),
        "rho_phys": float(rho_phys),
        "rho_sim": float(rho_sim),
        "rest_length_phys": float(L0_phys),
        "rest_length_sim": float(rest_length_sim),
        "spring_k_phys": float(k_phys),
        "spring_k_sim": float(k_sim),
        "dt_phys": float(dt_phys),
        "dt_sim": float(dt_sim),
    }


def _make_system(cfg: PhysicalConfig, derived: Dict[str, Any]):
    dim = Dimensionality.THREE_D
    if torch.cuda.is_available():
        device = torch.device('cuda')
        dtype = torch.float32
    elif torch.backends.mps.is_available():
        device = torch.device('mps')
        dtype = torch.float32
    else:
        device = torch.device('cpu')
        dtype = torch.float64

    h_sim = derived["h_sim"]

    state = BraneState(cfg.grid_shape, dim, device=device, dtype=dtype)
    state.initialize_flat_configuration(h_sim)  # spacing [m]
    state.set_fixed_boundaries()  # fix all outer faces

    grid = BraneGrid(cfg.grid_shape, dim, h_sim, device=device)

    physics = SpringForceComputer(
        spring_constant=derived["spring_k_sim"],
        rest_length=derived["rest_length_sim"],
    )

    mass_model = MassModel.from_density(
        density=derived["rho_sim"],
        intrinsic_dim=dim.value,
        spacing=h_sim,
    )

    solver = VelocityVerletSolver(dt=derived["dt_sim"], mass_model=mass_model, physics=physics, grid=grid)
    solver.initialize_accelerations(state)

    return state, grid, physics, solver


def run_measurement(direction: torch.Tensor, name: str, cfg: PhysicalConfig, run: TestRunManager) -> None:
    derived = _build_mapper_and_sim_params(cfg)

    state, grid, physics, solver = _make_system(cfg, derived)

    # Direction must be unit.
    direction = direction.to(state.device, dtype=state.positions.dtype)
    direction = direction / torch.norm(direction)

    center_idx = _center_index_from_coords(state.grid_coords, cfg.grid_shape)
    subs_idx = _viz_indices_3d(state.grid_coords, cfg.subsample_factor_3d, center_idx, cfg.patch_radius, cfg.include_coarse_background)

    mapper = derived["mapper"]

    X0 = state.rest_positions.clone()  # [N,4] (sim)

    # Displacement ramp in sim units
    d_max_sim = cfg.d_max_in_spacing * float(derived["h_sim"])
    disps_sim = np.linspace(0.0, d_max_sim, cfg.ramp_increments, dtype=np.float64)
    hold_sim = np.zeros_like(disps_sim)

    # Point cloud frames
    frames_data_display: List[Tuple[np.ndarray, np.ndarray]] = []
    frame_disp_display: List[float] = []

    # For measuring displacement field along direction
    direction_cpu = direction.detach().cpu().numpy().astype(np.float64)
    X0_cpu = X0.detach().cpu().numpy().astype(np.float64)
    center_rest_cpu = X0_cpu[center_idx]

    for i, d_sim in enumerate(disps_sim):
        target = X0[center_idx] + float(d_sim) * direction

        # Relax to equilibrium under this constraint
        for _ in range(cfg.relax_steps_per_increment):
            _clamp_center(state, center_idx, target)
            solver.step(state)

            # quasi-static damping
            state.velocities *= (1.0 - cfg.damping_per_step)

            # enforce fixed boundaries + exact clamp
            state.apply_fixed_boundaries()
            _clamp_center(state, center_idx, target)

        # Measure internal force and infer required holding force
        F_sim = physics.compute_forces(state, grid)
        hold_sim[i] = -float(torch.dot(F_sim[center_idx], direction).item())

        # Collect point cloud frames
        if (i % cfg.frame_stride) == 0:
            X_sim = state.positions.detach().cpu().numpy().astype(np.float64)

            # displacement field along direction (sim length units)
            delta_dir_sim = (X_sim - X0_cpu) @ direction_cpu  # [N]
            delta_dir_phys = mapper.to_phys_length(delta_dir_sim)

            # coords: deformed xyz positions (sim) -> nm
            coords_m = mapper.to_phys_length(X_sim[subs_idx.cpu().numpy()])

            # Choose a 3D projection so the *moved* direction is visually obvious:
            # - X^4 run: show (X, Y, X^4) so the clamp motion becomes a geometric axis.
            # - Lateral run (we use X direction): show (Y, Z, X) so the clamp motion is "vertical".
            if name == 'x4':
                axes = np.array([0, 1, 3])
            elif name == 'lat':
                axes = np.array([1, 2, 0])
            else:
                axes = np.array([0, 1, 2])

            coords_m = coords_m[:, axes]
            if cfg.center_on_rest:
                coords_m = coords_m - mapper.to_phys_length(center_rest_cpu[axes])

            coords_display = coords_m * cfg.display_unit_scale

            # values: displacement along direction (sim) -> nm
            values_display = delta_dir_phys[subs_idx.cpu().numpy()] * cfg.display_unit_scale

            frames_data_display.append((coords_display, values_display))

            # title parameter: clamp displacement in nm
            d_display = float(mapper.to_phys_length(d_sim) * cfg.display_unit_scale)
            frame_disp_display.append(d_display)

    # Convert final curves to display units + N
    disps_phys_m = mapper.to_phys_length(disps_sim)
    disps_display = disps_phys_m * cfg.display_unit_scale
    hold_N = mapper.to_phys_force(hold_sim)

    X_final = state.positions.detach().cpu().numpy().astype(np.float64)
    delta_dir_final = mapper.to_phys_length((X_final - X0_cpu) @ direction_cpu)
    max_abs_disp_sim = float(np.max(np.abs(delta_dir_final)))
    center_disp_sim = float(delta_dir_final[center_idx])

    print(
        f"[directional_hold_force_physical] {name} max |Δ| = {max_abs_disp_sim * cfg.display_unit_scale:.3f} "
        f"{cfg.display_unit_name}, center Δ = {center_disp_sim * cfg.display_unit_scale:.3f} {cfg.display_unit_name}"
    )

    # Save plot
    fig = plt.figure(figsize=(6.2, 4.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(disps_display, hold_N)
    ax.set_title(f"Holding force curve ({name})")
    ax.set_xlabel(f"Displacement [{cfg.display_unit_name}]")
    ax.set_ylabel("Holding force [N]")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    curve_path = Path(run.get_plot_path(f"force_curve_{name}.png"))
    fig.savefig(curve_path, dpi=220)
    plt.close(fig)

    # Save raw data
    data_path = Path(run.get_data_path(f"force_curve_{name}.npz"))
    np.savez(
        data_path,
        disps_sim=disps_sim,
        hold_sim=hold_sim,
        disps_phys_m=disps_phys_m,
        disps_display=disps_display,
        hold_N=hold_N,
        h_phys_nm=derived["h_phys"] * NM,
        h_sim=derived["h_sim"],
        dt_sim=derived["dt_sim"],
        dt_phys=derived["dt_phys"],
        m_point=cfg.m_point,
        rho_phys=derived["rho_phys"],
        rho_sim=derived["rho_sim"],
        spring_k_phys=derived["spring_k_phys"],
        spring_k_sim=derived["spring_k_sim"],
        rest_length_phys_nm=derived["rest_length_phys"] * NM,
        rest_length_sim=derived["rest_length_sim"]
    )

    # Create 3D point cloud animation
    def camera_motion_func(frame_idx: int, num_frames: int):
        return camera_orbit(
            frame_idx,
            num_frames,
            elev_start=18,
            elev_end=32,
            azim_start=-55,
            azim_revolutions=0.60,
        )

    out_mp4 = run.get_plot_path(f"hold_force_point_cloud_{name}.mp4")
    if name == 'x4':
        xlabel, ylabel, zlabel = f"x [{cfg.display_unit_name}]", f"y [{cfg.display_unit_name}]", f"x^4 [{cfg.display_unit_name}]"
    elif name == 'lat':
        xlabel, ylabel, zlabel = f"y [{cfg.display_unit_name}]", f"z [{cfg.display_unit_name}]", f"x [{cfg.display_unit_name}]"
    else:
        xlabel, ylabel, zlabel = f"x [{cfg.display_unit_name}]", f"y [{cfg.display_unit_name}]", f"z [{cfg.display_unit_name}]"

    create_3d_animation(
        frames_data=frames_data_display,
        times=frame_disp_display,
        output_path=out_mp4,
        cmap_name='RdBu_r',
        point_size=cfg.point_size,
        gamma=cfg.gamma,
        alpha_scale=cfg.alpha_scale,
        min_alpha=cfg.min_alpha,
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel=zlabel,
        title_template=f"{name} pull (Δ = {{:.3f}} {cfg.display_unit_name})",
        time_scale=1.0,
        fps=cfg.fps,
        dpi=110,
        camera_motion=camera_motion_func,
        figsize=(10, 8),
    )

    print(f"[directional_hold_force_physical] saved {curve_path} and {out_mp4}")


def test_directional_hold_force_physical():
    """Pytest entry point (opt-in)."""
    if os.environ.get("RUN_DIRECTIONAL_HOLD_FORCE", "0") != "1":
        import pytest
        pytest.skip("Set RUN_DIRECTIONAL_HOLD_FORCE=1 to run (generates plots/videos).")

    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force_physical")
    cfg = PhysicalConfig()

    derived = _build_mapper_and_sim_params(cfg)

    run.save_config({
        **cfg.__dict__,
        "h_phys_nm": derived["h_phys"] * NM,
        "h_sim": derived["h_sim"],
        "dt_sim": derived["dt_sim"],
        "dt_phys": derived["dt_phys"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_nm": derived["rest_length_phys"] * NM, 
        "rest_length_sim": derived["rest_length_sim"],
        "protocol": "quasi-static displacement clamp of center node; report Δ in display units; point-cloud MP4",
    })
    print(run.get_summary())

    # X^4 direction
    e_x4 = torch.tensor([0.0, 0.0, 0.0, 1.0])

    # Lateral direction (X axis)
    e_lat = torch.tensor([1.0, 0.0, 0.0, 0.0])

    run_measurement(e_x4, name="x4", cfg=cfg, run=run)
    run_measurement(e_lat, name="lat", cfg=cfg, run=run)


def _main():
    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force_physical")
    cfg = PhysicalConfig()

    derived = _build_mapper_and_sim_params(cfg)

    run.save_config({
        **cfg.__dict__,
        "h_phys_nm": derived["h_phys"] * NM,
        "h_sim": derived["h_sim"],
        "dt_sim": derived["dt_sim"],
        "dt_phys": derived["dt_phys"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_nm": derived["rest_length_phys"] * NM, 
        "rest_length_sim": derived["rest_length_sim"],
        "protocol": "quasi-static displacement clamp of center node; report Δ in display units; point-cloud MP4",
    })
    print(run.get_summary())

    e_x4 = torch.tensor([0.0, 0.0, 0.0, 1.0])
    e_lat = torch.tensor([1.0, 0.0, 0.0, 0.0])

    print("[directional_hold_force_physical] running X^4 measurement...")
    run_measurement(e_x4, name="x4", cfg=cfg, run=run)

    print("[directional_hold_force_physical] running lateral measurement...")
    run_measurement(e_lat, name="lat", cfg=cfg, run=run)

    print(f"[directional_hold_force_physical] done. Output in {run.run_dir.resolve()}")


if __name__ == "__main__":
    _main()
