"""Configuration models for the baryon pipeline components."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any


@dataclass(frozen=True)
class LatticeConfig:
    """Static lattice setup shared by initialization and simulation."""

    grid_shape: tuple[int, int, int]
    spacing: float = 1.0
    periodic_axes: tuple[bool, bool, bool] = (False, False, False)
    fixed_boundaries: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "grid_shape": list(self.grid_shape),
            "spacing": float(self.spacing),
            "periodic_axes": list(self.periodic_axes),
            "fixed_boundaries": bool(self.fixed_boundaries),
        }

    @staticmethod
    def from_dict(payload: dict[str, Any]) -> "LatticeConfig":
        if "grid_shape" not in payload:
            raise ValueError("LatticeConfig requires 'grid_shape'")
        return LatticeConfig(
            grid_shape=tuple(int(v) for v in payload["grid_shape"]),
            spacing=float(payload.get("spacing", 1.0)),
            periodic_axes=tuple(bool(v) for v in payload.get("periodic_axes", (False, False, False))),
            fixed_boundaries=bool(payload.get("fixed_boundaries", True)),
        )


@dataclass(frozen=True)
class DynamicsConfig:
    """Integrator and force parameters for lattice evolution."""

    spring_constant: float
    rest_length: float
    mass_density: float
    dt: float
    num_steps: int
    checkpoint_interval: int = 1

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
    def from_dict(payload: dict[str, Any]) -> "DynamicsConfig":
        required = ("spring_constant", "rest_length", "mass_density", "dt", "num_steps")
        missing = [key for key in required if key not in payload]
        if missing:
            raise ValueError(f"DynamicsConfig missing required keys: {', '.join(missing)}")
        return DynamicsConfig(
            spring_constant=float(payload["spring_constant"]),
            rest_length=float(payload["rest_length"]),
            mass_density=float(payload["mass_density"]),
            dt=float(payload["dt"]),
            num_steps=int(payload["num_steps"]),
            checkpoint_interval=int(payload.get("checkpoint_interval", 1)),
        )


@dataclass(frozen=True)
class BaryonSeedConfig:
    """Spherical-harmonic triplet seed parameters."""

    l: int = 1
    m: int = 1
    n: int = 1
    radius: float = 10.0
    amplitude: float = 0.25
    wave_speed: float = 1.0
    smooth_edge: float = 2.0
    center: tuple[float, float, float] | None = None
    axis_amplitudes: tuple[float, float, float] = (1.0, 1.0, 1.0)
    axis_phase_offsets: tuple[float, float, float] = (0.0, 2.0943951023931953, 4.1887902047863905)
    mixing_strength: float = 0.15
    x4_trace_weight: float = 0.35

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        if self.center is not None:
            payload["center"] = list(self.center)
        payload["axis_amplitudes"] = list(self.axis_amplitudes)
        payload["axis_phase_offsets"] = list(self.axis_phase_offsets)
        return payload

    @staticmethod
    def from_dict(payload: dict[str, Any]) -> "BaryonSeedConfig":
        center = payload.get("center")
        return BaryonSeedConfig(
            l=int(payload.get("l", 1)),
            m=int(payload.get("m", 1)),
            n=int(payload.get("n", 1)),
            radius=float(payload.get("radius", 10.0)),
            amplitude=float(payload.get("amplitude", 0.25)),
            wave_speed=float(payload.get("wave_speed", 1.0)),
            smooth_edge=float(payload.get("smooth_edge", 2.0)),
            center=tuple(float(v) for v in center) if center is not None else None,
            axis_amplitudes=tuple(float(v) for v in payload.get("axis_amplitudes", (1.0, 1.0, 1.0))),
            axis_phase_offsets=tuple(float(v) for v in payload.get("axis_phase_offsets", (0.0, 2.0943951023931953, 4.1887902047863905))),
            mixing_strength=float(payload.get("mixing_strength", 0.15)),
            x4_trace_weight=float(payload.get("x4_trace_weight", 0.35)),
        )


FORMAT_VERSION = "baryon-pipeline-v1"
