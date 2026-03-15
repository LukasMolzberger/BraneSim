"""Component 2: pure lattice simulation (input file -> compressed trajectory)."""

from __future__ import annotations

import argparse
import json

from branesim.baryon_pipeline import run_simulation_component


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run baryon simulation from a serialized initial state.")
    parser.add_argument("--input", required=True, help="Input .npz initial-state package")
    parser.add_argument("--output", required=True, help="Output trajectory .zip")
    parser.add_argument("--num-steps", type=int, default=None, help="Override num_steps from init package")
    parser.add_argument("--checkpoint-interval", type=int, default=None, help="Override checkpoint interval")
    parser.add_argument("--device", type=str, default="auto", help="auto/cpu/cuda/mps")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_simulation_component(
        initial_state_path=args.input,
        output_trajectory_path=args.output,
        num_steps=args.num_steps,
        checkpoint_interval=args.checkpoint_interval,
        device=args.device,
    )
    print("Simulation complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
