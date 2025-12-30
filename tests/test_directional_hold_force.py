"""tests/test_directional_hold_force.py

Directional hold-force measurement (X^4 vs. lateral).

This test is intentionally **opt-in** because it generates plots/videos and can
run for a while.

Run as a script:
    python -m tests.test_directional_hold_force

Run via pytest:
    RUN_DIRECTIONAL_HOLD_FORCE=1 pytest -q -s -k directional_hold_force

Outputs are written to a timestamped test run directory via TestRunManager:
    test-runs/directional_hold_force_YYYYMMDD_HHMMSS/

Notes
-----
- Uses the *project's* branesim implementation:
    BraneState, BraneGrid, SpringForceComputer, MassModel, VelocityVerletSolver
- Displacement-controlled quasi-static protocol:
    clamp center node position -> relax -> measure internal force -> holding
    force is negative internal force projected onto displacement direction.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

# Import TestRunManager from branesim.utils
try:
    from branesim.utils.test_run_manager import TestRunManager
except Exception as e:
    raise ImportError("Could not import TestRunManager. Ensure branesim/utils/test_run_manager.py is present.") from e

import numpy as np
import torch
import matplotlib.pyplot as plt

# Video writing (optional)
try:
    import imageio.v3 as iio
    _HAS_IMAGEIO = True
except Exception:
    _HAS_IMAGEIO = False

# Project imports (these are your actual implementation types)
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.physics.forces import SpringForceComputer
from branesim.core.dimensions import MassModel
from branesim.core.solver import VelocityVerletSolver


@dataclass
class Config:
    # Grid
    grid_shape: Tuple[int, int, int] = (33, 33, 33)
    spacing: float = 1.0

    # Spring physics
    spring_constant: float = 1.0
    rest_length: float = 1.0  # usually = spacing for an unstressed flat brane

    # Mass / time step
    density: float = 1.0  # in kg/m^d (d = intrinsic dimension)
    dt: float = 5e-3

    # Quasi-static ramp
    d_max_in_spacing: float = 0.5
    ramp_increments: int = 180
    relax_steps_per_increment: int = 80

    # Additional damping to enforce quasi-static adaptation
    damping_per_step: float = 0.02

    # Video
    fps: int = 30
    frame_stride: int = 2

    # Line-slice for video (vary x-axis, hold other coords at center)
    slice_axis: int = 0


def _center_index_from_coords(grid_coords: torch.Tensor, grid_shape: Tuple[int, ...]) -> int:
    center = torch.tensor([s // 2 for s in grid_shape], device=grid_coords.device, dtype=grid_coords.dtype)
    matches = (grid_coords == center).all(dim=1)
    idx = torch.nonzero(matches, as_tuple=False).squeeze(1)
    if idx.numel() != 1:
        raise RuntimeError(f"Expected exactly one center index, got {idx.numel()}")
    return int(idx.item())


def _center_line_indices(grid_coords: torch.Tensor, grid_shape: Tuple[int, ...], axis: int) -> np.ndarray:
    # Keep all axes except `axis` at center.
    center = [s // 2 for s in grid_shape]
    mask = torch.ones((grid_coords.shape[0],), device=grid_coords.device, dtype=torch.bool)
    for d in range(len(grid_shape)):
        if d == axis:
            continue
        mask &= (grid_coords[:, d] == center[d])
    idx = torch.nonzero(mask, as_tuple=False).squeeze(1)
    # Sort by the varying axis coordinate
    sort_key = grid_coords[idx, axis]
    idx = idx[torch.argsort(sort_key)]
    return idx.detach().cpu().numpy()


def _clamp_center(state: BraneState, center_idx: int, target: torch.Tensor) -> None:
    """Hard clamp the center node to `target` and zero its kinematics."""
    state.positions[center_idx] = target
    state.velocities[center_idx] = 0.0
    state.accelerations[center_idx] = 0.0
    state.new_accelerations[center_idx] = 0.0


def _make_system(cfg: Config) -> Tuple[BraneState, BraneGrid, SpringForceComputer, VelocityVerletSolver]:
    dim = Dimensionality.THREE_D
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    state = BraneState(cfg.grid_shape, dim, device=device, dtype=torch.float32)
    state.initialize_flat_configuration(cfg.spacing)
    state.set_fixed_boundaries()  # fixes all outer faces

    grid = BraneGrid(cfg.grid_shape, dim, cfg.spacing, device=device)

    physics = SpringForceComputer(spring_constant=cfg.spring_constant, rest_length=cfg.rest_length)

    mass_model = MassModel.from_density(cfg.density, intrinsic_dim=dim.value, spacing=cfg.spacing)

    solver = VelocityVerletSolver(dt=cfg.dt, mass_model=mass_model, physics=physics, grid=grid)
    solver.initialize_accelerations(state)

    return state, grid, physics, solver


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
            pass

        # Fallback GIF
        out_gif = out_path_mp4.with_suffix('.gif')
        iio.imwrite(out_gif, frames, fps=fps)
        return out_gif

    # No writer available
    return out_path_mp4


def run_measurement(direction: torch.Tensor, name: str, cfg: Config, run: TestRunManager) -> None:
    """Run one quasi-static displacement-controlled measurement."""
    state, grid, physics, solver = _make_system(cfg)

    # Direction must be unit.
    direction = direction.to(state.device, dtype=state.positions.dtype)
    direction = direction / torch.norm(direction)

    center_idx = _center_index_from_coords(state.grid_coords, cfg.grid_shape)
    line_idx = _center_line_indices(state.grid_coords, cfg.grid_shape, axis=cfg.slice_axis)

    X0 = state.rest_positions.clone()  # type: ignore

    d_max = cfg.d_max_in_spacing * cfg.spacing
    disps = np.linspace(0.0, d_max, cfg.ramp_increments, dtype=np.float64)
    hold = np.zeros_like(disps)

    frames: List[np.ndarray] = []

    for i, d in enumerate(disps):
        target = X0[center_idx] + float(d) * direction

        # Relax to equilibrium under this constraint
        for _ in range(cfg.relax_steps_per_increment):
            _clamp_center(state, center_idx, target)
            solver.step(state)
            # Extra quasi-static damping
            state.velocities *= (1.0 - cfg.damping_per_step)
            # Ensure fixed boundaries and center remain exact
            state.apply_fixed_boundaries()
            _clamp_center(state, center_idx, target)

        # Measure internal force and infer required holding force
        F = physics.compute_forces(state, grid)
        hold[i] = -float(torch.dot(F[center_idx], direction).item())

        # Video frame
        if (i % cfg.frame_stride) == 0:
            # Line slice displacement along direction
            X = state.positions.detach().cpu().numpy()
            Xref = X0.detach().cpu().numpy()

            x_axis = Xref[line_idx, cfg.slice_axis]  # material axis coordinate in embedding component
            delta = (X[line_idx] - Xref[line_idx]) @ direction.detach().cpu().numpy()

            fig = plt.figure(figsize=(12, 4))
            ax1 = fig.add_subplot(1, 2, 1)
            ax1.plot(x_axis, delta)
            ax1.set_title(f"Center line slice: Δ along {name}")
            ax1.set_xlabel(f"X^{cfg.slice_axis+1} (rest)")
            ax1.set_ylabel(f"Δ along {name}")
            ax1.grid(True, alpha=0.3)

            ax2 = fig.add_subplot(1, 2, 2)
            ax2.plot(disps[: i + 1], hold[: i + 1])
            ax2.set_title("Holding force vs displacement")
            ax2.set_xlabel("Displacement")
            ax2.set_ylabel("Holding force")
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
    # Save curve + data via TestRunManager
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(disps, hold)
    ax.set_title(f"Holding force curve ({name})")
    ax.set_xlabel("Displacement")
    ax.set_ylabel("Holding force")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    curve_path = Path(run.get_plot_path(f"force_curve_{name}.png"))
    fig.savefig(curve_path, dpi=200)
    # Save raw data
    data_path = Path(run.get_data_path(f"force_curve_{name}.npz"))
    np.savez(data_path, displacement=disps, hold_force=hold)
    plt.close(fig)

    # Save video
    video_path = Path(run.get_plot_path(f"hold_force_{name}.mp4"))
    _write_video(frames, video_path, fps=cfg.fps)


def test_directional_hold_force_measurement():
    # Opt-in test
    if os.environ.get("RUN_DIRECTIONAL_HOLD_FORCE", "0") != "1":
        import pytest
        pytest.skip("Set RUN_DIRECTIONAL_HOLD_FORCE=1 to run (generates plots/videos).")

    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force")
    cfg = Config()
    run.save_config({
        **cfg.__dict__,
        "protocol": "quasi-static displacement clamp of center node; measure holding force",
    })
    print(run.get_summary())

    # X^4 direction = embedding component index 3
    e_x4 = torch.tensor([0.0, 0.0, 0.0, 1.0])

    # Any lateral direction: choose X^1 (index 0)
    e_lat = torch.tensor([1.0, 0.0, 0.0, 0.0])

    run_measurement(e_x4, name="x4", cfg=cfg, run=run)
    run_measurement(e_lat, name="lat", cfg=cfg, run=run)


def _main():
    run = TestRunManager(base_dir="test-runs", experiment_name="directional_hold_force")
    cfg = Config()
    run.save_config({
        **cfg.__dict__,
        "protocol": "quasi-static displacement clamp of center node; measure holding force",
    })
    print(run.get_summary())

    e_x4 = torch.tensor([0.0, 0.0, 0.0, 1.0])
    e_lat = torch.tensor([1.0, 0.0, 0.0, 0.0])

    print("[directional_hold_force] running X^4 measurement...")
    run_measurement(e_x4, name="x4", cfg=cfg, run=run)

    print("[directional_hold_force] running lateral measurement...")
    run_measurement(e_lat, name="lat", cfg=cfg, run=run)

    print(f"[directional_hold_force] done. Output in {run.run_dir.resolve()}")


if __name__ == "__main__":
    _main()
