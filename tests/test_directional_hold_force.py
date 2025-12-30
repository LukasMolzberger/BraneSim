"""tests/test_directional_hold_force_physical.py

Directional hold-force measurement (X^4 vs. lateral) with **physical (SI) units**.

This test is intentionally **opt-in** (it generates plots/videos and can take a while).

Run as a script:
    python -m tests.test_directional_hold_force_physical

Run via pytest:
    RUN_DIRECTIONAL_HOLD_FORCE=1 pytest -q -s -k directional_hold_force_physical

Outputs are written via TestRunManager into a timestamped directory:
    test-runs/directional_hold_force_physical_YYYYMMDD_HHMMSS/

What is "physical" here?
------------------------
Internally, the solver still runs in dimensionless simulation units for numerical
stability. However, *every reported quantity* (axes labels, saved data, video
plots) is converted to SI using DimensionalMapper.

So you get:
  - displacement in meters [m]
  - holding force in Newton [N]

The scaling is chosen consistently with the rest of the project (same as the
photon/electron experiments):
  - lattice spacing a = h_phys = lambda_C / points_per_lambdaC
  - mass reference = m_point (physical mass per brane node)
  - time scale = h_phys / c

Holding force definition
------------------------
At each displacement increment, after relaxation:
    F_hold = -(F_internal(center) · e_dir)
where e_dir is the unit direction vector in embedding space.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Any

import numpy as np
import torch
import matplotlib.pyplot as plt

# Optional video writing
try:
    import imageio.v3 as iio
    _HAS_IMAGEIO = True
except Exception:
    _HAS_IMAGEIO = False

# -----------------------------------------------------------------------------
# TestRunManager
# -----------------------------------------------------------------------------
try:
    from branesim.utils.test_run_manager import TestRunManager
except Exception as e:
    raise ImportError(
        "Could not import TestRunManager. Ensure branesim/utils/test_run_manager.py is present."
    ) from e

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


@dataclass
class PhysicalConfig:
    """Configuration in a *physically meaningful* form.

    Notes
    -----
    - The grid spacing is defined physically via lambda_C / points_per_lambdaC.
    - The solver dt is configured in simulation units (dt_sim) because that is
      the stable knob; dt_s is derived.
    """

    # Domain (in lattice points)
    grid_shape: Tuple[int, int, int] = (33, 33, 33)

    # Physical scale choice
    points_per_lambdaC: int = 20

    # Mass per lattice node (kg). Keep consistent with your other calibrated runs.
    # If you prefer, you can compute this from a chosen density and h_phys.
    m_point: float = 2.861821e-27

    # Quasi-static ramp (specified relative to spacing, but outputs in meters)
    d_max_in_spacing: float = 0.5
    ramp_increments: int = 180
    relax_steps_per_increment: int = 80

    # Numerical stability
    dt_sim: float = 5e-3
    damping_per_step: float = 0.02

    # Video
    fps: int = 30
    frame_stride: int = 2

    # Line slice for video (vary this intrinsic axis, hold others at center)
    slice_axis: int = 0


def _center_index_from_coords(grid_coords: torch.Tensor, grid_shape: Tuple[int, ...]) -> int:
    center = torch.tensor([s // 2 for s in grid_shape], device=grid_coords.device, dtype=grid_coords.dtype)
    matches = (grid_coords == center).all(dim=1)
    idx = torch.nonzero(matches, as_tuple=False).squeeze(1)
    if idx.numel() != 1:
        raise RuntimeError(f"Expected exactly one center index, got {idx.numel()}")
    return int(idx.item())


def _center_line_indices(grid_coords: torch.Tensor, grid_shape: Tuple[int, ...], axis: int) -> np.ndarray:
    center = [s // 2 for s in grid_shape]
    mask = torch.ones((grid_coords.shape[0],), device=grid_coords.device, dtype=torch.bool)
    for d in range(len(grid_shape)):
        if d == axis:
            continue
        mask &= (grid_coords[:, d] == center[d])
    idx = torch.nonzero(mask, as_tuple=False).squeeze(1)
    idx = idx[torch.argsort(grid_coords[idx, axis])]
    return idx.detach().cpu().numpy()


def _clamp_center(state: BraneState, center_idx: int, target: torch.Tensor) -> None:
    """Hard clamp the center node to `target` and zero its kinematics."""
    state.positions[center_idx] = target
    state.velocities[center_idx] = 0.0
    state.accelerations[center_idx] = 0.0
    state.new_accelerations[center_idx] = 0.0


def _write_video(frames: List[np.ndarray], out_path_mp4: Path, fps: int) -> Path:
    """Write MP4 if possible, else GIF."""
    out_path_mp4.parent.mkdir(parents=True, exist_ok=True)

    if not frames:
        return out_path_mp4

    if _HAS_IMAGEIO:
        try:
            iio.imwrite(out_path_mp4, frames, fps=fps)
            return out_path_mp4
        except Exception:
            out_gif = out_path_mp4.with_suffix('.gif')
            iio.imwrite(out_gif, frames, fps=fps)
            return out_gif

    return out_path_mp4


def _build_mapper_and_sim_params(cfg: PhysicalConfig) -> Dict[str, Any]:
    """Create DimensionalMapper and derive simulation parameters."""
    constants = PhysicalConstants()

    # physical lattice spacing
    h_phys = constants.lambda_C / cfg.points_per_lambdaC

    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=cfg.m_point,
    )

    # In this convention, h_sim == 1.0
    h_sim = float(mapper.to_sim_length(h_phys))

    # Density corresponding to m_point on a 3D lattice cell of size h_phys
    rho_phys = cfg.m_point / (h_phys ** 3)

    # Convert density to simulation units (for MassModel.from_density)
    # density units: mass / length^3
    density_scale = mapper.mass_scale / (mapper.length_scale ** 3)
    rho_sim = rho_phys / density_scale

    # Rest length from calibrated fraction
    L0_phys = constants.compute_rest_length(h_phys)  # rest_length_frac * h_phys
    L0_sim = float(mapper.to_sim_length(L0_phys))

    # Choose spring constant consistent with target tension (same relation used in the constants docstring)
    # T_target = rho_phys * c^2
    T_target = constants.compute_target_tension(rho_phys)
    # T = (k/a) (1 - L0/a)  =>  k = T * a / (1 - L0/a)
    frac = constants.rest_length_frac
    k_phys = (T_target * h_phys) / max(1e-30, (1.0 - frac))
    k_sim = float(mapper.to_sim_spring_constant(k_phys))

    # Time step
    dt_s = float(mapper.to_phys_time(cfg.dt_sim))

    return {
        "constants": constants,
        "mapper": mapper,
        "h_phys": h_phys,
        "h_sim": h_sim,
        "rho_phys": rho_phys,
        "rho_sim": float(rho_sim),
        "rest_length_phys": L0_phys,
        "rest_length_sim": L0_sim,
        "spring_k_phys": float(k_phys),
        "spring_k_sim": k_sim,
        "dt_s": dt_s,
    }


def _make_system(cfg: PhysicalConfig, derived: Dict[str, Any]):
    dim = Dimensionality.THREE_D

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Simulation spacing
    h_sim = derived["h_sim"]

    state = BraneState(cfg.grid_shape, dim, device=device, dtype=torch.float32)
    state.initialize_flat_configuration(h_sim)
    state.set_fixed_boundaries()  # fixes all outer faces

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

    solver = VelocityVerletSolver(dt=cfg.dt_sim, mass_model=mass_model, physics=physics, grid=grid)
    solver.initialize_accelerations(state)

    return state, grid, physics, solver


def run_measurement(direction: torch.Tensor, name: str, cfg: PhysicalConfig, run: TestRunManager) -> None:
    derived = _build_mapper_and_sim_params(cfg)
    mapper: DimensionalMapper = derived["mapper"]

    state, grid, physics, solver = _make_system(cfg, derived)

    # Direction must be unit.
    direction = direction.to(state.device, dtype=state.positions.dtype)
    direction = direction / torch.norm(direction)

    center_idx = _center_index_from_coords(state.grid_coords, cfg.grid_shape)
    line_idx = _center_line_indices(state.grid_coords, cfg.grid_shape, axis=cfg.slice_axis)

    X0 = state.rest_positions.clone()  # [N,4] (sim)

    # Displacement ramp: define in sim units, report in meters.
    d_max_sim = cfg.d_max_in_spacing * float(derived["h_sim"])
    disps_sim = np.linspace(0.0, d_max_sim, cfg.ramp_increments, dtype=np.float64)

    hold_sim = np.zeros_like(disps_sim)

    # Precompute conversions for x-axis in the slice
    X0_cpu = X0.detach().cpu().numpy()
    x_axis_sim = X0_cpu[line_idx, cfg.slice_axis]
    x_axis_m = mapper.to_phys_length(x_axis_sim)

    frames: List[np.ndarray] = []

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

        # Video frame
        if (i % cfg.frame_stride) == 0:
            X = state.positions.detach().cpu().numpy()

            delta_sim = (X[line_idx] - X0_cpu[line_idx]) @ direction.detach().cpu().numpy()
            delta_m = mapper.to_phys_length(delta_sim)

            # physical curve
            disps_m_so_far = mapper.to_phys_length(disps_sim[: i + 1])
            hold_N_so_far = mapper.to_phys_force(hold_sim[: i + 1])

            fig = plt.figure(figsize=(12, 4))

            ax1 = fig.add_subplot(1, 2, 1)
            ax1.plot(x_axis_m, delta_m)
            ax1.set_title(f"Center line slice: Δ along {name}")
            ax1.set_xlabel("Position along slice [m]")
            ax1.set_ylabel(f"Displacement along {name} [m]")
            ax1.grid(True, alpha=0.3)

            ax2 = fig.add_subplot(1, 2, 2)
            ax2.plot(disps_m_so_far, hold_N_so_far)
            ax2.set_title("Holding force vs displacement")
            ax2.set_xlabel("Displacement [m]")
            ax2.set_ylabel("Holding force [N]")
            ax2.grid(True, alpha=0.3)

            fig.tight_layout()
            fig.canvas.draw()
            buf = np.asarray(fig.canvas.buffer_rgba())
            w, h = fig.canvas.get_width_height()
            # Account for high DPI displays by calculating actual dimensions from buffer size
            actual_size = buf.size // 4  # 4 channels (RGBA)
            actual_h = int(np.sqrt(actual_size * h / w))
            actual_w = actual_size // actual_h
            img = buf.reshape(actual_h, actual_w, 4)[:, :, :3]  # Drop alpha channel to get RGB
            frames.append(img)
            plt.close(fig)

    # Convert final curves to SI
    disps_m = mapper.to_phys_length(disps_sim)
    hold_N = mapper.to_phys_force(hold_sim)

    # Save plot
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(disps_m, hold_N)
    ax.set_title(f"Holding force curve ({name})")
    ax.set_xlabel("Displacement [m]")
    ax.set_ylabel("Holding force [N]")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    curve_path = Path(run.get_plot_path(f"force_curve_{name}.png"))
    fig.savefig(curve_path, dpi=200)
    plt.close(fig)

    # Save raw data (both sim + physical)
    data_path = Path(run.get_data_path(f"force_curve_{name}.npz"))
    np.savez(
        data_path,
        disps_sim=disps_sim,
        hold_sim=hold_sim,
        disps_m=disps_m,
        hold_N=hold_N,
        h_phys=derived["h_phys"],
        h_sim=derived["h_sim"],
        dt_sim=cfg.dt_sim,
        dt_s=derived["dt_s"],
        m_point=cfg.m_point,
        rho_phys=derived["rho_phys"],
        rho_sim=derived["rho_sim"],
        spring_k_phys=derived["spring_k_phys"],
        spring_k_sim=derived["spring_k_sim"],
        rest_length_phys=derived["rest_length_phys"],
        rest_length_sim=derived["rest_length_sim"],
        length_scale=mapper.length_scale,
        time_scale=mapper.time_scale,
        mass_scale=mapper.mass_scale,
    )

    # Save video
    video_path = Path(run.get_plot_path(f"hold_force_{name}.mp4"))
    _write_video(frames, video_path, fps=cfg.fps)


def test_directional_hold_force_physical():
    if os.environ.get("RUN_DIRECTIONAL_HOLD_FORCE", "0") != "1":
        import pytest
        pytest.skip("Set RUN_DIRECTIONAL_HOLD_FORCE=1 to run (generates plots/videos).")

    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force_physical")
    cfg = PhysicalConfig()

    derived = _build_mapper_and_sim_params(cfg)
    mapper: DimensionalMapper = derived["mapper"]

    run.save_config({
        **cfg.__dict__,
        "h_phys_m": derived["h_phys"],
        "h_sim": derived["h_sim"],
        "dt_s": derived["dt_s"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_m": derived["rest_length_phys"],
        "rest_length_sim": derived["rest_length_sim"],
        "mapper": repr(mapper),
        "protocol": "quasi-static displacement clamp of center node; measure holding force in SI",
    })
    print(run.get_summary())

    # X^4 direction
    e_x4 = torch.tensor([0.0, 0.0, 0.0, 1.0])

    # Lateral direction: choose X^1
    e_lat = torch.tensor([1.0, 0.0, 0.0, 0.0])

    run_measurement(e_x4, name="x4", cfg=cfg, run=run)
    run_measurement(e_lat, name="lat", cfg=cfg, run=run)


def _main():
    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force_physical")
    cfg = PhysicalConfig()

    derived = _build_mapper_and_sim_params(cfg)
    mapper: DimensionalMapper = derived["mapper"]

    run.save_config({
        **cfg.__dict__,
        "h_phys_m": derived["h_phys"],
        "h_sim": derived["h_sim"],
        "dt_s": derived["dt_s"],
        "rho_phys_kg_per_m3": derived["rho_phys"],
        "rho_sim": derived["rho_sim"],
        "spring_k_phys_N_per_m": derived["spring_k_phys"],
        "spring_k_sim": derived["spring_k_sim"],
        "rest_length_phys_m": derived["rest_length_phys"],
        "rest_length_sim": derived["rest_length_sim"],
        "mapper": repr(mapper),
        "protocol": "quasi-static displacement clamp of center node; measure holding force in SI",
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
