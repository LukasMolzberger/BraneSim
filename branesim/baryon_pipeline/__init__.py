"""Modular baryon simulation pipeline.

Components:
1. `initialization` - build compressed baryon initial states
2. `simulation` - evolve cubic lattice from input package
3. `visualization` - render outputs from stored trajectories
4. `diagnostics` - offline Berry + QCD-inspired analysis
"""

from __future__ import annotations

from .io import (
    InitialStatePackage,
    TrajectoryFrame,
    iter_trajectory_frames,
    load_initial_state_package,
    load_trajectory_manifest,
    save_initial_state_package,
)
from .models import BaryonSeedConfig, DynamicsConfig, LatticeConfig


def initialize_baryon_triplet_state(*args, **kwargs):
    from .initialization import initialize_baryon_triplet_state as _impl

    return _impl(*args, **kwargs)


def run_simulation_component(*args, **kwargs):
    from .simulation import run_simulation_component as _impl

    return _impl(*args, **kwargs)


def run_visualization_component(*args, **kwargs):
    from .visualization import run_visualization_component as _impl

    return _impl(*args, **kwargs)


def run_diagnostics_component(*args, **kwargs):
    from .diagnostics import run_diagnostics_component as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "initialize_baryon_triplet_state",
    "run_simulation_component",
    "run_visualization_component",
    "run_diagnostics_component",
    "InitialStatePackage",
    "TrajectoryFrame",
    "save_initial_state_package",
    "load_initial_state_package",
    "load_trajectory_manifest",
    "iter_trajectory_frames",
    "BaryonSeedConfig",
    "DynamicsConfig",
    "LatticeConfig",
]
