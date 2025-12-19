"""
Clean, layered brane initialization pipeline.

This package implements a principled initialization system with strict separation:
- Layer 0: Rest geometry (substrate only)
- Layer 1: Specs (EM-language interface, not physics)
- Layer 2: Carrier compilation (displacements + velocities)

Key principles:
1. Never write carrier values into absolute coordinates
2. Momentum is part of the spec (flows top-down)
3. Diagnostics must not require hardcoded ω
4. Every layer produces artifacts + visualizations
"""

from branesim.init.artifacts import (
    RestGeometryArtifact,
    SpecArtifact,
    CarrierArtifact,
    InitPipelineArtifact,
)
from branesim.init.specs import (
    PhotonSpec,
    ElectronSpec,
    VelocityInit,
)
from branesim.init.compiler import (
    initialize_state_from_spec,
    compile_photon,
    compile_electron,
)

__all__ = [
    "RestGeometryArtifact",
    "SpecArtifact",
    "CarrierArtifact",
    "InitPipelineArtifact",
    "PhotonSpec",
    "ElectronSpec",
    "VelocityInit",
    "initialize_state_from_spec",
    "compile_photon",
    "compile_electron",
]