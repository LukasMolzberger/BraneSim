"""Component 4: diagnostics (Berry + QCD-inspired metrics)."""

from __future__ import annotations

import argparse
import json

from branesim.baryon_pipeline import run_diagnostics_component


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run diagnostics from a compressed trajectory.")
    parser.add_argument("--input", required=True, help="Input trajectory .zip")
    parser.add_argument("--output-dir", required=True, help="Directory for diagnostics outputs")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--berry-point-stride", type=int, default=4)
    parser.add_argument("--berry-omega0", type=float, default=None)
    parser.add_argument("--omega-ref", type=float, default=1.0)
    parser.add_argument("--leakage-radius-factor", type=float, default=2.0)
    parser.add_argument("--render-berry-videos", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = run_diagnostics_component(
        trajectory_path=args.input,
        output_dir=args.output_dir,
        frame_stride=args.frame_stride,
        max_frames=args.max_frames,
        berry_point_stride=args.berry_point_stride,
        berry_omega0=args.berry_omega0,
        omega_ref=args.omega_ref,
        leakage_radius_factor=args.leakage_radius_factor,
        render_berry_videos=args.render_berry_videos,
    )
    print("Diagnostics complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
