"""File-format contracts for the branesim package (ARCHITECTURE.md §6).

FORMAT_VERSION = "branesim-block-v1"

Two file schemas:

  boundary_problem.npz  — initializer output / solver input
  worldvolume.zip       — solver output / diagnostics+viz input

All JSON blobs are stored as length-1 string arrays (allow_pickle=False
everywhere, matching the validated legacy pattern).

The world-volume zip mirrors the legacy trajectory format but:
  - uses the new FORMAT_VERSION string,
  - stores slices (not frames) under "slices/slice_{l:06d}.npz",
  - carries a solver_report block.
"""

from __future__ import annotations

import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import numpy as np

FORMAT_VERSION = "branesim-block-v1"


# ---------------------------------------------------------------------------
# boundary_problem.npz
# ---------------------------------------------------------------------------


def save_boundary_problem(
    path: str | Path,
    *,
    ref_positions: np.ndarray,
    boundary_slices: np.ndarray,
    boundary_indices: np.ndarray,
    lattice: dict[str, Any],
    action: dict[str, Any],
    boundary_mask: dict[str, Any] | None = None,
    seed: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> Path:
    """Write a boundary problem specification to ``boundary_problem.npz``.

    Parameters
    ----------
    path : path-like
    ref_positions : ndarray, shape (n_nodes, m_ambient)
        Held (unstressed) reference lattice.
    boundary_slices : ndarray, shape (Nb, n_nodes, m_ambient)
        Prescribed slice configurations.
    boundary_indices : ndarray of int, shape (Nb,)
        Which temporal slice index l each boundary slice pins.
    lattice : dict
        Keys: grid_shape, spacing, periodic_axes, axial_weight, dim.
    action : dict
        Keys: k_s, alpha, rho, dt, n_slices, temporal_model, r_t.
    boundary_mask : dict, optional
        Which components/nodes are fixed (for partial BCs / chirality).
    seed : dict, optional
        Seed/ansatz metadata.
    metadata : dict, optional
        Provenance information.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out,
        format_version=np.array([FORMAT_VERSION]),
        ref_positions=np.asarray(ref_positions, dtype=np.float64),
        boundary_slices=np.asarray(boundary_slices, dtype=np.float64),
        boundary_indices=np.asarray(boundary_indices, dtype=np.int64),
        boundary_mask_json=np.array([json.dumps(boundary_mask or {})]),
        lattice_json=np.array([json.dumps(lattice)]),
        action_json=np.array([json.dumps(action)]),
        seed_json=np.array([json.dumps(seed or {})]),
        metadata_json=np.array([json.dumps(metadata or {})]),
    )
    return out


def load_boundary_problem(path: str | Path) -> dict[str, Any]:
    """Load a boundary problem from ``boundary_problem.npz``.

    Returns a dict with keys:
      ref_positions, boundary_slices, boundary_indices,
      boundary_mask (dict), lattice (dict), action (dict),
      seed (dict), metadata (dict).

    Raises
    ------
    ValueError
        If the format version does not match FORMAT_VERSION.
    """
    payload = np.load(path, allow_pickle=False)
    version = str(payload["format_version"][0])
    if version != FORMAT_VERSION:
        raise ValueError(
            f"Unsupported format version {version!r}; "
            f"expected {FORMAT_VERSION!r}"
        )
    return {
        "ref_positions": payload["ref_positions"],
        "boundary_slices": payload["boundary_slices"],
        "boundary_indices": payload["boundary_indices"],
        "boundary_mask": json.loads(str(payload["boundary_mask_json"][0])),
        "lattice": json.loads(str(payload["lattice_json"][0])),
        "action": json.loads(str(payload["action_json"][0])),
        "seed": json.loads(str(payload["seed_json"][0])),
        "metadata": json.loads(str(payload["metadata_json"][0])),
    }


# ---------------------------------------------------------------------------
# worldvolume.zip writer / reader
# ---------------------------------------------------------------------------


@dataclass
class SliceMeta:
    index: int
    time: float
    name: str  # path inside the zip


class WorldVolumeWriter:
    """Write a world-volume zip file (ARCHITECTURE.md §6.2).

    Usage::

        with WorldVolumeWriter(path, manifest_extra) as w:
            for l, positions in enumerate(slices):
                w.write_slice(l, l * dt, positions)
    """

    def __init__(
        self,
        path: str | Path,
        manifest_extra: dict[str, Any] | None = None,
    ) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._zf = zipfile.ZipFile(
            self.path, mode="w", compression=zipfile.ZIP_DEFLATED
        )
        self._slices: list[dict[str, Any]] = []
        self._manifest_extra = manifest_extra or {}

    def write_slice(
        self,
        index: int,
        time: float,
        positions: np.ndarray,
    ) -> None:
        """Append one spacelike slice to the world-volume."""
        name = f"slices/slice_{index:06d}.npz"
        buf = io.BytesIO()
        np.savez_compressed(buf, positions=np.asarray(positions, dtype=np.float64))
        self._zf.writestr(name, buf.getvalue())
        self._slices.append({"index": index, "time": float(time), "name": name})

    def write_npy(self, name: str, array: np.ndarray) -> None:
        """Store an auxiliary numpy array (e.g. ref_positions)."""
        buf = io.BytesIO()
        np.save(buf, np.asarray(array), allow_pickle=False)
        self._zf.writestr(name, buf.getvalue())

    def close(self, solver_report: dict[str, Any] | None = None) -> None:
        manifest = {
            "format_version": FORMAT_VERSION,
            "mode": self._manifest_extra.get("mode", "ivp"),
            "slices": self._slices,
            "solver_report": solver_report or {},
        }
        manifest.update(
            {k: v for k, v in self._manifest_extra.items() if k not in manifest}
        )
        self._zf.writestr(
            "manifest.json",
            json.dumps(manifest, indent=2).encode("utf-8"),
        )
        self._zf.close()

    def __enter__(self) -> "WorldVolumeWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Read manifest.json from a worldvolume.zip.

    Raises ValueError if the format version is wrong.
    """
    with zipfile.ZipFile(path, "r") as zf:
        manifest = json.loads(zf.read("manifest.json").decode("utf-8"))
    version = manifest.get("format_version")
    if version != FORMAT_VERSION:
        raise ValueError(
            f"Unsupported format version {version!r}; expected {FORMAT_VERSION!r}"
        )
    return manifest


def iter_slices(
    path: str | Path,
    stride: int = 1,
) -> Iterator[tuple[int, float, np.ndarray]]:
    """Iterate over slices in a worldvolume.zip.

    Yields
    ------
    (index, time, positions)
        positions has shape (n_nodes, m_ambient).
    """
    if stride < 1:
        raise ValueError("stride must be >= 1")
    manifest = load_manifest(path)
    with zipfile.ZipFile(path, "r") as zf:
        for meta in manifest["slices"][::stride]:
            raw = zf.read(meta["name"])
            payload = np.load(io.BytesIO(raw), allow_pickle=False)
            yield int(meta["index"]), float(meta["time"]), payload["positions"]


def load_npy(path: str | Path, name: str) -> np.ndarray:
    """Load an auxiliary .npy array from a worldvolume.zip."""
    with zipfile.ZipFile(path, "r") as zf:
        raw = zf.read(name)
    return np.load(io.BytesIO(raw), allow_pickle=False)
