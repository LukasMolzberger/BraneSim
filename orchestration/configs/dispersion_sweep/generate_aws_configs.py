#!/usr/bin/env python3
"""Generate all AWS-scale dispersion sweep JSON configs (96^3 x 512, bvp_chiral).

Run this once locally to produce the per-run config files:
    python orchestration/configs/dispersion_sweep/generate_aws_configs.py

Output: one JSON per (direction, k_index, pol_axis) in the same directory.
The AWS sweep runner (run_dispersion_sweep.sh) iterates over them.
"""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).parent

# Canonical physics parameters (dimensionless units)
BASE = {
    "lattice": {
        "grid_shape": [96, 96, 96],
        "spacing": 1.0,
        "periodic_axes": [True, True, True],
        "axial_weight": 1.0,
    },
    "action": {
        "k_s": 1.0, "alpha": 0.2, "rho": 1.0, "dt": 0.1,
        "n_slices": 512, "m_ambient": 4,
        "r_t": 0.0,
    },
    "solver": {"mode": "bvp_chiral", "chirality": "forward"},
}

POL_VECS = {
    0: [1, 0, 0, 0],
    1: [0, 1, 0, 0],
    2: [0, 0, 1, 0],
}

SPECS = []
# [100]: L (p0) and T (p1) for n in 1..4
for n in [1, 2, 3, 4]:
    SPECS.append(("[100]", [n, 0, 0], 0))
    SPECS.append(("[100]", [n, 0, 0], 1))
# [110]: hard (p0, omega_x=omega_y) and soft (p2, z-out) for n in 1..4
for n in [1, 2, 3, 4]:
    SPECS.append(("[110]", [n, n, 0], 0))
    SPECS.append(("[110]", [n, n, 0], 2))
# [111]: triplet degenerate; two pols for triplet check, n in 1..4
for n in [1, 2, 3, 4]:
    SPECS.append(("[111]", [n, n, n], 0))
    SPECS.append(("[111]", [n, n, n], 1))
# [210]: three distinct eigenvalues; all three pols for n=1,2,3
for n in [1, 2, 3]:
    SPECS.append(("[210]", [2 * n, n, 0], 0))
    SPECS.append(("[210]", [2 * n, n, 0], 1))
    SPECS.append(("[210]", [2 * n, n, 0], 2))


def make_name(direction: str, k_idx: list[int], pol_axis: int) -> str:
    dir_tag = direction.strip("[]")
    k_tag = "_".join(str(v) for v in k_idx)
    return f"aws_{dir_tag}_k{k_tag}_p{pol_axis}"


def make_config(direction: str, k_idx: list[int], pol_axis: int) -> dict:
    import copy
    cfg = copy.deepcopy(BASE)
    cfg["_comment"] = (
        f"AWS sweep: {direction} k={k_idx} pol_axis={pol_axis}. "
        f"96^3 x 512, bvp_chiral. Memory per worldvolume: ~14.5 GB."
    )
    cfg["seed"] = {
        "kind": "plane_wave",
        "amplitude": 1e-3,
        "k_index": k_idx,
        "polarization": POL_VECS[pol_axis],
        "rng_seed": 0,
    }
    return cfg


if __name__ == "__main__":
    names = []
    for direction, k_idx, pol_axis in SPECS:
        name = make_name(direction, k_idx, pol_axis)
        cfg = make_config(direction, k_idx, pol_axis)
        path = HERE / f"{name}.json"
        path.write_text(json.dumps(cfg, indent=2))
        names.append(name)
        print(f"Wrote {path.name}")
    print(f"\nTotal: {len(names)} configs")
