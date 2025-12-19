"""
Specification layer (Layer 1): EM-facing interface.

These are pure spec dataclasses that describe what you want to initialize,
not how to build it. They use "EM language" as an interface layer but do
not implement any physics - that's handled by the compiler layer.
"""

from dataclasses import dataclass
from typing import Literal
import torch

# Velocity initialization methods
VelocityInit = Literal["time_reversal_shift", "directional_derivative", "complex_quadrature"]


@dataclass
class PhotonSpec:
    """
    Specification for a photon-like wave packet.

    This is the "EM language" interface: you specify what you want
    (center, width, momentum, helicity) and the compiler figures out
    how to build it on the substrate.

    Key principle: k_vector sets BOTH wavelength AND momentum direction.
    There is no separate "Layer 3" for momentum - it flows top-down.
    """
    intrinsic_dim: int               # 1, 2, or 3
    center: torch.Tensor             # [d] center position in intrinsic coords
    sigma: float                     # Gaussian envelope width
    amplitude: float                 # Peak amplitude A
    k_vector: torch.Tensor           # [d] wave vector (sets λ and direction)
    helicity: Literal["L", "R"]      # circular polarization handedness
    prefer_shear: bool = True        # P1: prefer in-brane shear over X⁴
    velocity_init: VelocityInit = "time_reversal_shift"
    shift_cells: int = 1             # for time_reversal_shift mode
    periodic_shift: bool = False     # whether to use periodic wrapping

    def __post_init__(self):
        """Validate spec after initialization."""
        if self.intrinsic_dim not in (1, 2, 3):
            raise ValueError(f"intrinsic_dim must be 1, 2, or 3, got {self.intrinsic_dim}")
        if self.sigma <= 0:
            raise ValueError(f"sigma must be positive, got {self.sigma}")
        if self.amplitude <= 0:
            raise ValueError(f"amplitude must be positive, got {self.amplitude}")
        if len(self.center) != self.intrinsic_dim:
            raise ValueError(
                f"center must have {self.intrinsic_dim} components, "
                f"got {len(self.center)}"
            )
        if len(self.k_vector) != self.intrinsic_dim:
            raise ValueError(
                f"k_vector must have {self.intrinsic_dim} components, "
                f"got {len(self.k_vector)}"
            )


@dataclass
class ElectronSpec:
    """
    Specification for an electron-like tubular excitation (3D only).

    This uses tubular geometry with spinorial transport (half-angle rotation
    of the polarization frame around the loop).
    """
    intrinsic_dim: int = 3           # enforced to be 3D
    center: torch.Tensor = None      # [3] center of the torus
    amplitude: float = 1.0           # peak amplitude
    tube_sigma: float = 1.0          # radial width of tube

    # Torus knot geometry parameters
    torus_major_radius: float = 1.0
    torus_minor_radius: float = 0.3
    p: int = 2                       # core windings
    q: int = 1                       # tube windings
    num_samples: int = 800           # discretization of centerline

    # Carrier parameters
    longitudinal_k: float = 1.0      # phase advance along arclength
    helicity: Literal["L", "R"] = "R"
    velocity_init: VelocityInit = "complex_quadrature"

    def __post_init__(self):
        """Validate spec after initialization."""
        if self.intrinsic_dim != 3:
            raise ValueError("ElectronSpec requires intrinsic_dim=3")
        if self.amplitude <= 0:
            raise ValueError(f"amplitude must be positive, got {self.amplitude}")
        if self.tube_sigma <= 0:
            raise ValueError(f"tube_sigma must be positive, got {self.tube_sigma}")
        if self.torus_major_radius <= 0:
            raise ValueError(f"torus_major_radius must be positive")
        if self.torus_minor_radius <= 0:
            raise ValueError(f"torus_minor_radius must be positive")
        if self.center is not None and len(self.center) != 3:
            raise ValueError("center must have 3 components for electron")