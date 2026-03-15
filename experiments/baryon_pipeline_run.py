"""Single-entry orchestration for the modular baryon pipeline.

Runs components in sequence from one JSON config:
1) initialization
2) simulation
3) visualizations (one or more)
4) diagnostics
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import torch

from branesim.baryon_pipeline import (
    BaryonSeedConfig,
    DynamicsConfig,
    LatticeConfig,
    initialize_baryon_triplet_state,
    run_diagnostics_component,
    run_simulation_component,
    run_visualization_component,
    save_initial_state_package,
)


def _load_json(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _resolve_output_dir(config: dict[str, Any], override: str | None) -> Path:
    if override:
        out = Path(override)
    elif config.get("output_dir"):
        out = Path(config["output_dir"])
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out = Path("test-runs") / f"baryon_pipeline_{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _torch_dtype(name: str) -> torch.dtype:
    if name == "float32":
        return torch.float32
    if name == "float64":
        return torch.float64
    raise ValueError(f"Unsupported dtype {name!r}, expected 'float32' or 'float64'")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run full baryon pipeline from one JSON config.")
    parser.add_argument("--config", required=True, help="Pipeline JSON config path")
    parser.add_argument("--output-dir", default=None, help="Override output directory from config")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = _load_json(args.config)

    output_dir = _resolve_output_dir(config, args.output_dir)
    plots_dir = output_dir / "plots"
    diag_dir = output_dir / "diagnostics"
    plots_dir.mkdir(exist_ok=True)
    diag_dir.mkdir(exist_ok=True)

    if "lattice" not in config:
        raise ValueError("Missing required 'lattice' section in config")
    if "dynamics" not in config:
        raise ValueError("Missing required 'dynamics' section in config")

    lattice = LatticeConfig.from_dict(config["lattice"])
    seed = BaryonSeedConfig.from_dict(config.get("seed", {}))
    dynamics = DynamicsConfig.from_dict(config["dynamics"])

    runtime = config.get("runtime", {})
    device = runtime.get("device", "auto")
    dtype = _torch_dtype(runtime.get("dtype", "float64"))

    # 1) Initialization
    init_result = initialize_baryon_triplet_state(
        seed=seed,
        lattice=lattice,
        dynamics=dynamics,
        device=device,
        dtype=dtype,
    )
    initial_state_path = output_dir / "initial_state.npz"
    save_initial_state_package(initial_state_path, init_result.package)

    # 2) Simulation
    sim_cfg = config.get("simulation", {})
    trajectory_path = output_dir / sim_cfg.get("trajectory_name", "trajectory.zip")
    sim_summary = run_simulation_component(
        initial_state_path=initial_state_path,
        output_trajectory_path=trajectory_path,
        num_steps=sim_cfg.get("num_steps"),
        checkpoint_interval=sim_cfg.get("checkpoint_interval"),
        device=sim_cfg.get("device", device),
    )

    # 3) Visualization conversions
    viz_summaries: list[dict[str, Any]] = []
    for idx, viz_cfg in enumerate(config.get("visualizations", [])):
        mode = viz_cfg["mode"]
        output_name = viz_cfg.get("output_name")
        if not output_name:
            ext = ".mp4"
            output_name = f"viz_{idx}_{mode}{ext}"

        kwargs = dict(viz_cfg)
        kwargs.pop("mode", None)
        kwargs.pop("output_name", None)

        summary = run_visualization_component(
            trajectory_path=trajectory_path,
            output_path=plots_dir / output_name,
            mode=mode,
            **kwargs,
        )
        viz_summaries.append(summary)

    # 4) Diagnostics
    diag_cfg = config.get("diagnostics", {})
    diag_summary = run_diagnostics_component(
        trajectory_path=trajectory_path,
        output_dir=diag_dir,
        frame_stride=diag_cfg.get("frame_stride", 1),
        max_frames=diag_cfg.get("max_frames"),
        berry_point_stride=diag_cfg.get("berry_point_stride", 4),
        berry_omega0=diag_cfg.get("berry_omega0"),
        omega_ref=diag_cfg.get("omega_ref", 1.0),
        leakage_radius_factor=diag_cfg.get("leakage_radius_factor", 2.0),
        render_berry_videos=diag_cfg.get("render_berry_videos", False),
    )

    summary = {
        "config": str(Path(args.config).resolve()),
        "output_dir": str(output_dir.resolve()),
        "initial_state": str(initial_state_path.resolve()),
        "trajectory": str(trajectory_path.resolve()),
        "init_debug": init_result.debug,
        "simulation": sim_summary,
        "visualizations": viz_summaries,
        "diagnostics": diag_summary,
    }

    summary_path = output_dir / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Baryon pipeline run complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
