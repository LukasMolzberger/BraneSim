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
    d_max_in_spacing: float = 0.5
    ramp_increments: int = 180
    relax_steps_per_increment: int = 80

    # Numerics (simulation time step + simple damping)
    dt_sim: float = 5e-3
    damping_per_step: float = 0.02

    # Video (point cloud)
    fps: int = 20
    frame_stride: int = 2          # keep every Nth displacement step as a frame
    subsample_factor_3d: int = 2   # take every kth point along each intrinsic axis

    # Visualization: always include a dense local patch around the center
    patch_radius: int = 6           # in intrinsic grid steps (Chebyshev radius)
    include_coarse_background: bool = True

    # Point cloud appearance
    point_size: float = 4.0
    gamma: float = 1.3
    alpha_scale: float = 0.8


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

    # Density corresponding to m_point on a 3D lattice cell (physical)
    rho_phys = cfg.m_point / (h_phys ** 3)

    # Rest length (physical)
    L0_phys = constants.compute_rest_length(h_phys)

    # Target tension T = rho * c^2 in physical units, map to spring constant k [N/m]
    T_target = constants.compute_target_tension(rho_phys)
    frac = constants.rest_length_frac
    k_phys = float((T_target * h_phys) / max(1e-30, (1.0 - frac)))

    # Time step (already in seconds; keep name for backward-compat)
    dt_s = float(cfg.dt_sim)

    return {
        "constants": constants,
        "h_phys": float(h_phys),
        "h_sim": float(h_phys),  # kept for compatibility; this is physical spacing [m]
        "rho_phys": float(rho_phys),
        "rho_sim": float(rho_phys),  # kept for compatibility; this is physical density
        "rest_length_phys": float(L0_phys),
        "rest_length_sim": float(L0_phys),
        "spring_k_phys": float(k_phys),
        "spring_k_sim": float(k_phys),
        "dt_s": float(dt_s),
    }


def _make_system(cfg: PhysicalConfig, derived: Dict[str, Any]):
    dim = Dimensionality.THREE_D
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    h_sim = derived["h_sim"]

    state = BraneState(cfg.grid_shape, dim, device=device, dtype=torch.float32)
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

    solver = VelocityVerletSolver(dt=derived["dt_s"], mass_model=mass_model, physics=physics, grid=grid)
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

    X0 = state.rest_positions.clone()  # [N,4] (sim)

    # Displacement ramp in sim units
    d_max_sim = cfg.d_max_in_spacing * float(derived["h_sim"])
    disps_sim = np.linspace(0.0, d_max_sim, cfg.ramp_increments, dtype=np.float64)
    hold_sim = np.zeros_like(disps_sim)

    # Point cloud frames
    frames_data_nm: List[Tuple[np.ndarray, np.ndarray]] = []
    frame_disp_nm: List[float] = []

    # For measuring displacement field along direction
    direction_cpu = direction.detach().cpu().numpy().astype(np.float64)
    X0_cpu = X0.detach().cpu().numpy().astype(np.float64)

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

            # coords: deformed xyz positions (sim) -> nm
            coords_m = X_sim[subs_idx.cpu().numpy()]

            # Choose a 3D projection so the *moved* direction is visually obvious:
            # - X^4 run: show (X, Y, X^4) so the clamp motion becomes a geometric axis.
            # - Lateral run (we use X direction): show (Y, Z, X) so the clamp motion is "vertical".
            if name == 'x4':
                coords_m = coords_m[:, [0, 1, 3]]
            elif name == 'lat':
                coords_m = coords_m[:, [1, 2, 0]]
            else:
                coords_m = coords_m[:, :3]

            coords_nm = coords_m * NM

            # values: displacement along direction (sim) -> nm
            values_nm = delta_dir_sim[subs_idx.cpu().numpy()] * NM

            frames_data_nm.append((coords_nm, values_nm))

            # title parameter: clamp displacement in nm
            d_nm = float(d_sim * NM)
            frame_disp_nm.append(d_nm)

    # Convert final curves to nm + N
    disps_nm = disps_sim * NM
    hold_N = hold_sim  # already Newtons

    # Save plot
    fig = plt.figure(figsize=(6.2, 4.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(disps_nm, hold_N)
    ax.set_title(f"Holding force curve ({name})")
    ax.set_xlabel("Displacement [nm]")
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
        disps_nm=disps_nm,
        hold_N=hold_N,
        h_phys_nm=derived["h_phys"] * NM,
        h_sim=derived["h_sim"],
        dt_sim=cfg.dt_sim,
        dt_s=derived["dt_s"],
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
    create_3d_animation(
        frames_data=frames_data_nm,
        times=frame_disp_nm,
        output_path=out_mp4,
        cmap_name='RdBu_r',
        point_size=cfg.point_size,
        gamma=cfg.gamma,
        alpha_scale=cfg.alpha_scale,
        xlabel='x [nm]',
        ylabel='y [nm]',
        zlabel='z [nm]',
        title_template=f"{name} pull (Δ = {{:.3f}} nm)",
        fps=cfg.fps,
        dpi=110,
        camera_motion=camera_motion_func,
        figsize=(10, 8),
    )


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
        "dt_s": derived["dt_s"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_nm": derived["rest_length_phys"] * NM, 
        "rest_length_sim": derived["rest_length_sim"],
        "protocol": "quasi-static displacement clamp of center node; report Δ in nm; point-cloud MP4",
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
        "dt_s": derived["dt_s"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_nm": derived["rest_length_phys"] * NM, 
        "rest_length_sim": derived["rest_length_sim"],
        "protocol": "quasi-static displacement clamp of center node; report Δ in nm; point-cloud MP4",
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
