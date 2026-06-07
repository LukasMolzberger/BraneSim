"""Shared plot-style helpers for the diagnostics package.

Extracted here to avoid circular imports between run_measurements and
binding_probe (run_measurements imports binding_probe; binding_probe needs
the same helpers).

Consumers
---------
- branesim.diagnostics.run_measurements
- branesim.diagnostics.binding_probe
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb


_STYLE = {
    "figure.dpi": 120,
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "lines.linewidth": 1.6,
}


def _apply_style() -> None:
    matplotlib.rcParams.update(_STYLE)


def _savefig(fig: "plt.Figure", path: Path) -> None:
    # Only call tight_layout when constrained_layout is NOT active (avoids warning).
    if not fig.get_constrained_layout():
        fig.tight_layout()
    fig.savefig(str(path), dpi=120, bbox_inches="tight")
    plt.close(fig)


def _phase_to_rgb(phase: np.ndarray) -> np.ndarray:
    """Map phase in (-pi, pi) -> HSV hue -> RGB. Shape (...) -> (..., 3)."""
    hue = np.mod((phase + np.pi) / (2.0 * np.pi), 1.0)
    sat = np.ones_like(hue)
    val = np.ones_like(hue)
    return hsv_to_rgb(np.stack([hue, sat, val], axis=-1))