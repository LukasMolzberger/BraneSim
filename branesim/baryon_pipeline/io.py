"""Compressed I/O formats for baryon initialization and simulation outputs."""

from __future__ import annotations

import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from .models import DynamicsConfig, LatticeConfig, FORMAT_VERSION


@dataclass(frozen=True)
class InitialStatePackage:
    """Loaded initialization package."""

    positions: np.ndarray
    velocities: np.ndarray
    rest_positions: np.ndarray
    lattice: LatticeConfig
    dynamics: DynamicsConfig
    metadata: dict[str, Any]


@dataclass(frozen=True)
class TrajectoryFrame:
    """One checkpoint frame from a compressed trajectory."""

    index: int
    step: int
    time: float
    positions: np.ndarray
    velocities: np.ndarray


class CompressedTrajectoryWriter:
    """Streaming writer for zip-compressed trajectory checkpoints."""

    def __init__(self, path: str | Path, manifest: dict[str, Any]):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._zip = zipfile.ZipFile(self.path, mode="w", compression=zipfile.ZIP_DEFLATED)
        self._manifest = dict(manifest)
        self._manifest.setdefault("format_version", FORMAT_VERSION)
        self._manifest["frames"] = []
        self._closed = False

    def write_numpy(self, name: str, array: np.ndarray) -> None:
        payload = io.BytesIO()
        np.save(payload, np.asarray(array), allow_pickle=False)
        self._zip.writestr(name, payload.getvalue())

    def write_frame(self, step: int, time: float, positions: np.ndarray, velocities: np.ndarray) -> None:
        if self._closed:
            raise RuntimeError("Cannot write frame after close().")

        frame_index = len(self._manifest["frames"])
        frame_name = f"frames/frame_{frame_index:06d}.npz"

        payload = io.BytesIO()
        np.savez_compressed(
            payload,
            positions=np.asarray(positions),
            velocities=np.asarray(velocities),
        )
        self._zip.writestr(frame_name, payload.getvalue())
        self._manifest["frames"].append(
            {
                "index": frame_index,
                "step": int(step),
                "time": float(time),
                "name": frame_name,
            }
        )

    def close(self) -> None:
        if self._closed:
            return
        self._zip.writestr("manifest.json", json.dumps(self._manifest, indent=2).encode("utf-8"))
        self._zip.close()
        self._closed = True

    def __enter__(self) -> "CompressedTrajectoryWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def save_initial_state_package(path: str | Path, package: InitialStatePackage) -> Path:
    """Save initialization output as one compressed `.npz` package."""

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out_path,
        format_version=np.array([FORMAT_VERSION]),
        positions=np.asarray(package.positions),
        velocities=np.asarray(package.velocities),
        rest_positions=np.asarray(package.rest_positions),
        lattice_json=np.array([json.dumps(package.lattice.to_dict())]),
        dynamics_json=np.array([json.dumps(package.dynamics.to_dict())]),
        metadata_json=np.array([json.dumps(package.metadata)]),
    )
    return out_path


def load_initial_state_package(path: str | Path) -> InitialStatePackage:
    """Load a compressed initial-state package."""

    payload = np.load(Path(path), allow_pickle=False)
    format_version = str(payload["format_version"][0])
    if format_version != FORMAT_VERSION:
        raise ValueError(f"Unsupported initial-state format: {format_version!r}")

    lattice = LatticeConfig.from_dict(json.loads(str(payload["lattice_json"][0])))
    dynamics = DynamicsConfig.from_dict(json.loads(str(payload["dynamics_json"][0])))
    metadata = json.loads(str(payload["metadata_json"][0]))

    return InitialStatePackage(
        positions=payload["positions"],
        velocities=payload["velocities"],
        rest_positions=payload["rest_positions"],
        lattice=lattice,
        dynamics=dynamics,
        metadata=metadata,
    )


def load_trajectory_manifest(path: str | Path) -> dict[str, Any]:
    """Read manifest metadata from compressed trajectory."""

    with zipfile.ZipFile(path, mode="r") as zf:
        manifest = json.loads(zf.read("manifest.json").decode("utf-8"))
    return manifest


def load_trajectory_array(path: str | Path, name: str) -> np.ndarray:
    """Load a named `.npy` array written by `CompressedTrajectoryWriter.write_numpy`."""

    with zipfile.ZipFile(path, mode="r") as zf:
        raw = zf.read(name)
    return np.load(io.BytesIO(raw), allow_pickle=False)


def iter_trajectory_frames(path: str | Path, frame_stride: int = 1) -> Iterator[TrajectoryFrame]:
    """Iterate through trajectory checkpoints in order."""

    if frame_stride <= 0:
        raise ValueError("frame_stride must be >= 1")

    with zipfile.ZipFile(path, mode="r") as zf:
        manifest = json.loads(zf.read("manifest.json").decode("utf-8"))
        frames = manifest.get("frames", [])

        for meta in frames[::frame_stride]:
            raw = zf.read(meta["name"])
            payload = np.load(io.BytesIO(raw), allow_pickle=False)
            yield TrajectoryFrame(
                index=int(meta["index"]),
                step=int(meta["step"]),
                time=float(meta["time"]),
                positions=payload["positions"],
                velocities=payload["velocities"],
            )
