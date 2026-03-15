"""Component 3: independent visualization conversion from compressed trajectory."""

from __future__ import annotations

import argparse
import json

from branesim.baryon_pipeline import run_visualization_component


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render visual products from trajectory files.")
    parser.add_argument("--input", required=True, help="Input trajectory .zip")
    parser.add_argument("--output", required=True, help="Output movie path")
    parser.add_argument("--mode", required=True, choices=("volume", "slice"))
    parser.add_argument("--component", type=int, default=3)
    parser.add_argument("--absolute-field", action="store_true", help="Use absolute embedding component instead of displacement")
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=120)

    parser.add_argument("--subsample", type=int, default=2, help="Volume mode")
    parser.add_argument("--plane", type=str, default="xy", choices=("xy", "xz", "yz"), help="Slice mode")
    parser.add_argument("--index", type=int, default=None, help="Slice index for selected plane")
    parser.add_argument("--cmap", type=str, default="RdBu_r", help="Slice mode colormap")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    kwargs = {
        "component": args.component,
        "use_displacement": not args.absolute_field,
        "frame_stride": args.frame_stride,
        "fps": args.fps,
        "dpi": args.dpi,
    }

    if args.mode == "volume":
        kwargs["subsample"] = args.subsample
    else:
        kwargs["plane"] = args.plane
        kwargs["index"] = args.index
        kwargs["cmap"] = args.cmap

    summary = run_visualization_component(
        trajectory_path=args.input,
        output_path=args.output,
        mode=args.mode,
        **kwargs,
    )

    print("Visualization complete")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
