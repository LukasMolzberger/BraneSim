"""Shared file I/O for component communication."""

from __future__ import annotations

import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

import numpy as np

FORMAT_VERSION = "baryon-components-v3"


@dataclass(frozen=True)
class LatticeConfig:
    grid_shape: tuple[int, int, int]
    spacing: float
    periodic_axes: tuple[bool, bool, bool]
    fixed_boundaries: bool
    shell_weights: tuple[float, float, float] = (1.0, 1.0, 1.0)

    def to_dict(self) -> dict[str, Any]:
        return {
            "grid_shape": list(self.grid_shape),
            "spacing": float(self.spacing),
            "periodic_axes": list(self.periodic_axes),
            "fixed_boundaries": bool(self.fixed_boundaries),
            "shell_weights": [float(v) for v in self.shell_weights],
        }

    @staticmethod
    def from_dict(data: dict[str, Any]) -> "LatticeConfig":
        return LatticeConfig(
            grid_shape=tuple(int(v) for v in data["grid_shape"]),
            spacing=float(data["spacing"]),
            periodic_axes=tuple(bool(v) for v in data.get("periodic_axes", (False, False, False))),
            fixed_boundaries=bool(data.get("fixed_boundaries", True)),
            shell_weights=tuple(float(v) for v in data.get("shell_weights", (1.0, 1.0, 1.0))),
        )


@dataclass(frozen=True)
class DynamicsConfig:
    spring_constant: float
    rest_length: float
    mass_density: float
    dt: float
    num_steps: int
    checkpoint_interval: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "spring_constant": float(self.spring_constant),
            "rest_length": float(self.rest_length),
            "mass_density": float(self.mass_density),
            "dt": float(self.dt),
            "num_steps": int(self.num_steps),
            "checkpoint_interval": int(self.checkpoint_interval),
        }

    @staticmethod
    def from_dict(data: dict[str, Any]) -> "DynamicsConfig":
        return DynamicsConfig(
            spring_constant=float(data["spring_constant"]),
            rest_length=float(data["rest_length"]),
            mass_density=float(data["mass_density"]),
            dt=float(data["dt"]),
            num_steps=int(data["num_steps"]),
            checkpoint_interval=int(data.get("checkpoint_interval", 1)),
        )


@dataclass(frozen=True)
class TrajectoryFrame:
    index: int
    step: int
    time: float
    positions: np.ndarray
    velocities: np.ndarray


def save_initial_state(
    path: str | Path,
    *,
    positions: np.ndarray,
    velocities: np.ndarray,
    rest_positions: np.ndarray,
    lattice: LatticeConfig,
    dynamics: DynamicsConfig,
    seed: dict[str, Any],
    metadata: dict[str, Any],
) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out,
        format_version=np.array([FORMAT_VERSION]),
        positions=np.asarray(positions),
        velocities=np.asarray(velocities),
        rest_positions=np.asarray(rest_positions),
        lattice_json=np.array([json.dumps(lattice.to_dict())]),
        dynamics_json=np.array([json.dumps(dynamics.to_dict())]),
        seed_json=np.array([json.dumps(seed)]),
        metadata_json=np.array([json.dumps(metadata)]),
    )
    return out


def load_initial_state(path: str | Path) -> dict[str, Any]:
    payload = np.load(path, allow_pickle=False)
    version = str(payload["format_version"][0])
    if version != FORMAT_VERSION:
        raise ValueError(f"Unsupported format {version!r}, expected {FORMAT_VERSION!r}")

    lattice = LatticeConfig.from_dict(json.loads(str(payload["lattice_json"][0])))
    dynamics = DynamicsConfig.from_dict(json.loads(str(payload["dynamics_json"][0])))

    return {
        "positions": payload["positions"],
        "velocities": payload["velocities"],
        "rest_positions": payload["rest_positions"],
        "lattice": lattice,
        "dynamics": dynamics,
        "seed": json.loads(str(payload["seed_json"][0])),
        "metadata": json.loads(str(payload["metadata_json"][0])),
    }


class TrajectoryWriter:
    def __init__(self, path: str | Path, manifest: dict[str, Any]):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._zf = zipfile.ZipFile(self.path, mode="w", compression=zipfile.ZIP_DEFLATED)
        self.manifest = dict(manifest)
        self.manifest["format_version"] = FORMAT_VERSION
        self.manifest["frames"] = []

    def write_npy(self, name: str, array: np.ndarray) -> None:
        buffer = io.BytesIO()
        np.save(buffer, np.asarray(array), allow_pickle=False)
        self._zf.writestr(name, buffer.getvalue())

    def write_frame(self, *, step: int, time: float, positions: np.ndarray, velocities: np.ndarray) -> None:
        idx = len(self.manifest["frames"])
        name = f"frames/frame_{idx:06d}.npz"
        buff = io.BytesIO()
        np.savez_compressed(buff, positions=np.asarray(positions), velocities=np.asarray(velocities))
        self._zf.writestr(name, buff.getvalue())
        self.manifest["frames"].append({"index": idx, "step": int(step), "time": float(time), "name": name})

    def close(self) -> None:
        self._zf.writestr("manifest.json", json.dumps(self.manifest, indent=2).encode("utf-8"))
        self._zf.close()

    def __enter__(self) -> "TrajectoryWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def load_manifest(path: str | Path) -> dict[str, Any]:
    with zipfile.ZipFile(path, "r") as zf:
        manifest = json.loads(zf.read("manifest.json").decode("utf-8"))
    version = manifest.get("format_version")
    if version != FORMAT_VERSION:
        raise ValueError(f"Unsupported format {version!r}, expected {FORMAT_VERSION!r}")
    return manifest


def load_npy(path: str | Path, name: str) -> np.ndarray:
    with zipfile.ZipFile(path, "r") as zf:
        raw = zf.read(name)
    return np.load(io.BytesIO(raw), allow_pickle=False)


def iter_frames(path: str | Path, frame_stride: int = 1) -> Iterator[TrajectoryFrame]:
    if frame_stride <= 0:
        raise ValueError("frame_stride must be >= 1")
    manifest = load_manifest(path)
    with zipfile.ZipFile(path, "r") as zf:
        for meta in manifest["frames"][::frame_stride]:
            raw = zf.read(meta["name"])
            payload = np.load(io.BytesIO(raw), allow_pickle=False)
            yield TrajectoryFrame(
                index=int(meta["index"]),
                step=int(meta["step"]),
                time=float(meta["time"]),
                positions=payload["positions"],
                velocities=payload["velocities"],
            )
