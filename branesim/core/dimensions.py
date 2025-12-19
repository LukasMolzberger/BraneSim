"""
Brane dimensions and mass model: Explicit separation of intrinsic and embedding dimensions.

This module implements a systematic approach to dimensional bookkeeping that prevents
common errors related to mass density units and dimensional confusion.

Key principle:
    Intrinsic dimension (d) determines mass density units:
        d=1: ρ in kg/m  (linear density)   → m_node = ρ * h
        d=2: ρ in kg/m² (surface density)  → m_node = ρ * h²
        d=3: ρ in kg/m³ (volumetric density) → m_node = ρ * h³

    Embedding dimension (N) determines displacement vector size:
        N=4: position vectors are 4D (X⁰, X¹, X², X³)
"""

from dataclasses import dataclass
from typing import Optional
from enum import Enum


class IntrinsicDim(Enum):
    """Intrinsic dimensionality of the brane (d)."""
    ONE_D = 1    # 1D chain/string
    TWO_D = 2    # 2D sheet/membrane
    THREE_D = 3  # 3D volume/bulk


@dataclass(frozen=True)
class BraneDims:
    """
    Explicit separation of intrinsic and embedding dimensions.

    Attributes:
        intrinsic_dim: Intrinsic dimensionality of the brane (d = 1, 2, or 3)
        embedding_dim: Dimension of the embedding space (N, typically 4)

    Examples:
        BraneDims(1, 4): 1D chain living in R^4
        BraneDims(2, 4): 2D sheet living in R^4
        BraneDims(3, 4): 3D volume living in R^4
    """
    intrinsic_dim: int
    embedding_dim: int

    def __post_init__(self):
        """Validate dimensions."""
        if self.intrinsic_dim not in (1, 2, 3):
            raise ValueError(f"Intrinsic dimension must be 1, 2, or 3, got {self.intrinsic_dim}")
        if self.embedding_dim < self.intrinsic_dim:
            raise ValueError(
                f"Embedding dimension ({self.embedding_dim}) must be >= "
                f"intrinsic dimension ({self.intrinsic_dim})"
            )
        if self.embedding_dim > 4:
            raise ValueError(
                f"Embedding dimension must be <= 4 (currently support up to 4D embeddings), "
                f"got {self.embedding_dim}"
            )

    @property
    def d(self) -> int:
        """Intrinsic dimension (convenience accessor)."""
        return self.intrinsic_dim

    @property
    def N(self) -> int:
        """Embedding dimension (convenience accessor)."""
        return self.embedding_dim

    @property
    def num_transverse_dofs(self) -> int:
        """Number of transverse (polarization) degrees of freedom."""
        return self.embedding_dim - self.intrinsic_dim

    def __repr__(self) -> str:
        return f"BraneDims(d={self.intrinsic_dim}, N={self.embedding_dim})"


class MassModel:
    """
    Mass model with explicit density units and conversion to per-node mass.

    This class ensures that mass density is interpreted correctly according to
    the intrinsic brane dimensionality. It prevents the common error of passing
    a volumetric density (kg/m³) to a 1D or 2D brane without proper conversion.

    Attributes:
        m_node: Mass per lattice node [kg]
        density: Mass per intrinsic d-volume (proper units for the brane dimension)
        intrinsic_dim: Intrinsic dimension of the brane
        spacing: Lattice spacing [m]
    """

    def __init__(
        self,
        m_node: float,
        density: float,
        intrinsic_dim: int,
        spacing: float,
    ):
        """
        Initialize mass model.

        Args:
            m_node: Mass per lattice node [kg]
            density: Mass per intrinsic d-volume (units depend on intrinsic_dim)
            intrinsic_dim: Intrinsic dimension (1, 2, or 3)
            spacing: Lattice spacing [m]
        """
        self.m_node = m_node
        self.density = density
        self.intrinsic_dim = intrinsic_dim
        self.spacing = spacing

    @classmethod
    def from_density(
        cls,
        density: float,
        intrinsic_dim: int,
        spacing: float,
    ) -> "MassModel":
        """
        Create mass model from density with proper units.

        Args:
            density: Mass per intrinsic d-volume
                     d=1: kg/m   (linear density)
                     d=2: kg/m²  (surface density)
                     d=3: kg/m³  (volumetric density)
            intrinsic_dim: Intrinsic dimension (1, 2, or 3)
            spacing: Lattice spacing [m]

        Returns:
            MassModel with computed per-node mass

        Examples:
            # 1D chain with linear density
            mass = MassModel.from_density(1e-15, intrinsic_dim=1, spacing=1e-12)
            # m_node = ρ₁ * h

            # 2D sheet with surface density
            mass = MassModel.from_density(1e-6, intrinsic_dim=2, spacing=1e-12)
            # m_node = ρ₂ * h²

            # 3D volume with volumetric density
            mass = MassModel.from_density(1000.0, intrinsic_dim=3, spacing=1e-12)
            # m_node = ρ₃ * h³
        """
        if intrinsic_dim == 1:
            m_node = density * spacing
        elif intrinsic_dim == 2:
            m_node = density * spacing ** 2
        elif intrinsic_dim == 3:
            m_node = density * spacing ** 3
        else:
            raise ValueError(f"Intrinsic dimension must be 1, 2, or 3, got {intrinsic_dim}")

        return cls(
            m_node=m_node,
            density=density,
            intrinsic_dim=intrinsic_dim,
            spacing=spacing,
        )

    @classmethod
    def from_volumetric_density(
        cls,
        rho3: float,
        intrinsic_dim: int,
        spacing: float,
        cross_section: Optional[float] = None,
        thickness: Optional[float] = None,
    ) -> "MassModel":
        """
        Create mass model from volumetric density (kg/m³) with explicit reduction.

        This method allows specifying a 3D density even for lower-dimensional branes,
        but requires explicit geometric parameters for the reduction.

        Args:
            rho3: Volumetric mass density [kg/m³]
            intrinsic_dim: Intrinsic dimension (1, 2, or 3)
            spacing: Lattice spacing [m]
            cross_section: Cross-sectional area [m²] for 1D branes (required if intrinsic_dim=1)
            thickness: Thickness [m] for 2D branes (required if intrinsic_dim=2)

        Returns:
            MassModel with proper per-node mass

        Raises:
            ValueError: If cross_section or thickness is not provided when required

        Examples:
            # 1D string with circular cross-section
            mass = MassModel.from_volumetric_density(
                rho3=1000.0,           # kg/m³
                intrinsic_dim=1,
                spacing=1e-12,         # m
                cross_section=np.pi * (1e-13)**2  # π r²
            )
            # ρ₁ = ρ₃ * A,  m_node = ρ₁ * h

            # 2D sheet with thickness
            mass = MassModel.from_volumetric_density(
                rho3=1000.0,           # kg/m³
                intrinsic_dim=2,
                spacing=1e-12,         # m
                thickness=1e-13        # m
            )
            # ρ₂ = ρ₃ * t,  m_node = ρ₂ * h²
        """
        if intrinsic_dim == 3:
            # Direct volumetric density
            density = rho3
        elif intrinsic_dim == 2:
            # Need thickness to convert to surface density
            if thickness is None:
                raise ValueError(
                    "Must provide 'thickness' parameter when using volumetric density "
                    "with intrinsic_dim=2. Surface density = rho3 * thickness"
                )
            density = rho3 * thickness
        elif intrinsic_dim == 1:
            # Need cross-sectional area to convert to linear density
            if cross_section is None:
                raise ValueError(
                    "Must provide 'cross_section' parameter when using volumetric density "
                    "with intrinsic_dim=1. Linear density = rho3 * cross_section"
                )
            density = rho3 * cross_section
        else:
            raise ValueError(f"Intrinsic dimension must be 1, 2, or 3, got {intrinsic_dim}")

        return cls.from_density(density, intrinsic_dim, spacing)

    def get_density_units(self) -> str:
        """Get the proper units for the density based on intrinsic dimension."""
        if self.intrinsic_dim == 1:
            return "kg/m"
        elif self.intrinsic_dim == 2:
            return "kg/m²"
        elif self.intrinsic_dim == 3:
            return "kg/m³"
        else:
            raise ValueError(f"Unexpected intrinsic dimension: {self.intrinsic_dim}")

    def __repr__(self) -> str:
        units = self.get_density_units()
        return (
            f"MassModel(m_node={self.m_node:.6e} kg, "
            f"ρ={self.density:.6e} {units}, "
            f"d={self.intrinsic_dim}, h={self.spacing:.6e} m)"
        )


# Helper function for backward compatibility
def compute_node_mass(
    density: float,
    intrinsic_dim: int,
    spacing: float,
) -> float:
    """
    Compute mass per node from density (with proper units).

    This is a convenience function for the common case where you just need
    the per-node mass without the full MassModel object.

    Args:
        density: Mass per intrinsic d-volume
        intrinsic_dim: Intrinsic dimension (1, 2, or 3)
        spacing: Lattice spacing [m]

    Returns:
        Mass per node [kg]
    """
    return MassModel.from_density(density, intrinsic_dim, spacing).m_node