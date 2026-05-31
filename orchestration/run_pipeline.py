"""Pipeline orchestration (outside components).

Runs components as separate processes so components communicate only via files.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


def _load_json(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _output_dir(config: dict[str, Any], override: str | None) -> Path:
    if override:
        out = Path(override)
    elif config.get("output_dir"):
        out = Path(config["output_dir"])
    else:
        out = Path("test-runs") / f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _to_args(options: dict[str, Any]) -> list[str]:
    args: list[str] = []
    for key, value in options.items():
        if value is None:
            continue
        flag = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                args.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            if len(value) == 0:
                continue
            args.extend([flag, ",".join(str(v) for v in value)])
            continue
        args.extend([flag, str(value)])
    return args


def _run(cmd: list[str]) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run full 4-component pipeline from one JSON config")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = _load_json(args.config)
    out_dir = _output_dir(config, args.output_dir)
    plots_dir = out_dir / "plots"
    diag_dir = out_dir / "diagnostics"
    plots_dir.mkdir(exist_ok=True)
    diag_dir.mkdir(exist_ok=True)

    init_cfg = dict(config.get("initialization", {}))
    sim_cfg = dict(config.get("simulation", {}))
    diag_cfg = dict(config.get("diagnostics", {}))
    viz_cfgs = list(config.get("visualizations", []))

    if not init_cfg:
        raise ValueError("Missing required 'initialization' section")
    if not sim_cfg:
        raise ValueError("Missing required 'simulation' section")

    init_module = init_cfg.pop("module", "components.initialization.run")
    diag_module = diag_cfg.pop("module", "components.diagnostics.run") if diag_cfg else None

    initial_state = out_dir / init_cfg.pop("output_name", "initial_state.npz")
    trajectory = out_dir / sim_cfg.pop("output_name", "trajectory.zip")

    commands: list[list[str]] = []

    # 1) initialization
    cmd_init = [
        sys.executable,
        "-m",
        init_module,
        "--output",
        str(initial_state),
    ] + _to_args(init_cfg)
    commands.append(cmd_init)
    _run(cmd_init)

    # 2) simulation
    cmd_sim = [
        sys.executable,
        "-m",
        "components.simulation.run",
        "--input",
        str(initial_state),
        "--output",
        str(trajectory),
    ] + _to_args(sim_cfg)
    commands.append(cmd_sim)
    _run(cmd_sim)

    # 3) visualizations
    viz_outputs: list[str] = []
    for i, vc in enumerate(viz_cfgs):
        vc = dict(vc)
        mode = vc.pop("mode")
        output_name = vc.pop("output_name", f"viz_{i}_{mode}.mp4")
        output = plots_dir / output_name

        cmd_viz = [
            sys.executable,
            "-m",
            "components.visualization.run",
            "--input",
            str(trajectory),
            "--output",
            str(output),
            "--mode",
            str(mode),
        ] + _to_args(vc)
        commands.append(cmd_viz)
        _run(cmd_viz)
        viz_outputs.append(str(output))

    # 4) diagnostics
    if diag_module is not None:
        cmd_diag = [
            sys.executable,
            "-m",
            diag_module,
            "--input",
            str(trajectory),
            "--output-dir",
            str(diag_dir),
        ] + _to_args(diag_cfg)
        commands.append(cmd_diag)
        _run(cmd_diag)

    summary = {
        "config": str(Path(args.config).resolve()),
        "output_dir": str(out_dir.resolve()),
        "initial_state": str(initial_state.resolve()),
        "trajectory": str(trajectory.resolve()),
        "visualizations": viz_outputs,
        "diagnostics": str((diag_dir / "diagnostics_summary.json").resolve()),
        "commands": [" ".join(c) for c in commands],
    }
    summary_path = out_dir / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Pipeline complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
