"""
Artifact dataclasses for the layered initialization pipeline.

Each layer produces an artifact object containing:
- All derived data from that layer
- Metadata for debugging
- Information needed by downstream layers
"""

from dataclasses import dataclass, field
from typing import Optional
import torch


@dataclass
class RestGeometryArtifact:
    """
    Layer 0: Rest geometry artifact (substrate definition).

    This represents the unperturbed brane configuration before
    any wave packet or excitation is added.
    """
    intrinsic_dim: int          # 1, 2, or 3
    embedding_dim: int          # always 4
    grid_shape: tuple[int, ...]
    spacing: float
    rest_positions: torch.Tensor   # [N, 4]
    coords: torch.Tensor           # [N, d] intrinsic coordinates
    fixed_mask: Optional[torch.Tensor] = None  # [N] boolean mask


@dataclass
class SpecArtifact:
    """
    Layer 1: Specification artifact (EM-language interface).

    This is the "what do you want" layer - it describes the desired
    wave packet in physical terms (momentum, polarization, etc.)
    without specifying how to build it on the substrate.
    """
    kind: str               # "photon" or "electron"
    spec: object            # PhotonSpec or ElectronSpec
    k_hat: torch.Tensor     # [d] normalized propagation direction
    k_mag: float            # magnitude of k vector
    notes: dict = field(default_factory=dict)


@dataclass
class CarrierArtifact:
    """
    Layer 2: Compiled carrier artifact (substrate kinematics).

    This is the "how to build it" layer - concrete displacements and
    velocities that will be applied to the substrate.
    """
    envelope: torch.Tensor        # [N] amplitude envelope A(x)
    phase: torch.Tensor           # [N] carrier phase φ(x)
    p1: torch.Tensor              # [4] or [N,4] first polarization basis vector
    p2: torch.Tensor              # [4] or [N,4] second polarization basis vector
    psi: torch.Tensor             # complex [N,4] full complex carrier
    u0: torch.Tensor              # [N,4] displacement field
    v0: torch.Tensor              # [N,4] velocity field
    meta: dict = field(default_factory=dict)  # method-specific metadata


@dataclass
class InitPipelineArtifact:
    """
    Top-level artifact bundle containing outputs from all layers.

    This is returned by initialize_state_from_spec and contains
    everything needed for debugging and validation.
    """
    layer0: RestGeometryArtifact
    layer1: SpecArtifact
    layer2: CarrierArtifact