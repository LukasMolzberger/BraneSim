"""electron

Top-level electron initialization module.

This folder is intentionally separate from the `branesim` package so that the
electron seeding logic can evolve independently and be imported from
experiments/tests without depending on internal package layout.
"""

from branesim.electron.electron_initialization import ElectronModeSpec, initialize_electron_mode_3d

__all__ = ["ElectronModeSpec", "initialize_electron_mode_3d"]
