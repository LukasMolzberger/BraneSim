"""
Data types for first-order symplectic band diagnostics.

This module defines the configuration and result types for the symplectic
Berry phase pipeline. It uses "symplectic" in all names to distinguish from
potential future second-order (dynamical matrix / phonon) implementations.

Key principle:
    The first-order symplectic approach directly linearizes the equations of
    motion in first-order form (q, p) and extracts eigenmodes from the
    resulting symplectic operator A(k). This is distinct from second-order
    approaches that work with D(k) = K - ω²M.

All types here are dimension-agnostic and work for 1D, 2D, and 3D branes.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Tuple, Optional, Dict, Any
import torch


class BoundaryCondition(Enum):
    """
    Boundary condition types for band structure computation.

    PERIODIC: Use Fourier phases exp(i k·r) for periodic bulk bands.
              Required for k-space band structure computation.

    CLAMPED: Dirichlet boundary conditions (fixed boundaries).
             Not compatible with k-space bands (use for real-space eigenmodes).
    """
    PERIODIC = "periodic"
    CLAMPED = "clamped"


@dataclass
class SymplecticBandConfig:
    """
    Configuration for symplectic band structure computation.

    This config specifies the lattice geometry, material parameters, and
    boundary conditions for building the first-order symplectic operator A(k).

    Attributes:
        d: Intrinsic dimension (1, 2, or 3)
        embedding_dim: Embedding space dimension (typically 4)
        grid_shape: Tuple of grid dimensions, length d
        spacing: Lattice spacing h [m]
        mass: Mass per node [kg]
        spring_k: Spring stiffness [N/m]
        rest_length: Rest length of springs [m]
        neighbor_offsets: List of neighbor displacement vectors (in grid units)
                          Each offset is a tuple of length d (e.g., (1,) or (1,0) or (1,0,0))
        boundary: BoundaryCondition (PERIODIC or CLAMPED)
        device: torch.device for computation
        dtype: torch.dtype for floating-point arrays

    Example:
        # 1D periodic chain
        cfg = SymplecticBandConfig(
            d=1,
            embedding_dim=4,
            grid_shape=(64,),
            spacing=1e-12,
            mass=1e-30,
            spring_k=1e-6,
            rest_length=1e-12,
            neighbor_offsets=[(-1,), (1,)],
            boundary=BoundaryCondition.PERIODIC,
            device=torch.device('cpu'),
            dtype=torch.float64,
        )

        # 2D periodic sheet with nearest neighbors
        cfg = SymplecticBandConfig(
            d=2,
            embedding_dim=4,
            grid_shape=(32, 32),
            spacing=1e-12,
            mass=1e-30,
            spring_k=1e-6,
            rest_length=1e-12,
            neighbor_offsets=[(-1,0), (1,0), (0,-1), (0,1)],
            boundary=BoundaryCondition.PERIODIC,
            device=torch.device('cpu'),
            dtype=torch.float64,
        )
    """
    d: int
    embedding_dim: int
    grid_shape: Tuple[int, ...]
    spacing: float
    mass: float
    spring_k: float
    rest_length: float
    neighbor_offsets: List[Tuple[int, ...]]
    boundary: BoundaryCondition
    device: torch.device = torch.device('cpu')
    dtype: torch.dtype = torch.float64

    def __post_init__(self):
        """Validate configuration."""
        if self.d not in (1, 2, 3):
            raise ValueError(f"Intrinsic dimension d must be 1, 2, or 3, got {self.d}")

        if len(self.grid_shape) != self.d:
            raise ValueError(
                f"grid_shape length {len(self.grid_shape)} must match d={self.d}"
            )

        if self.embedding_dim < self.d:
            raise ValueError(
                f"Embedding dimension ({self.embedding_dim}) must be >= d ({self.d})"
            )

        if self.embedding_dim > 4:
            raise ValueError(
                f"Embedding dimension must be <= 4, got {self.embedding_dim}"
            )

        # Validate neighbor offsets
        for offset in self.neighbor_offsets:
            if len(offset) != self.d:
                raise ValueError(
                    f"Neighbor offset {offset} has wrong dimension "
                    f"(expected {self.d}, got {len(offset)})"
                )

        if self.mass <= 0:
            raise ValueError(f"Mass must be positive, got {self.mass}")

        if self.spring_k <= 0:
            raise ValueError(f"Spring constant must be positive, got {self.spring_k}")

        if self.spacing <= 0:
            raise ValueError(f"Spacing must be positive, got {self.spacing}")

        if self.rest_length <= 0:
            raise ValueError(f"Rest length must be positive, got {self.rest_length}")

    def n_cell_nodes(self) -> int:
        """Number of nodes per unit cell (1 for simple lattices)."""
        return 1

    def n_dof_per_cell(self) -> int:
        """Total DOFs per unit cell."""
        return self.n_cell_nodes() * self.embedding_dim


@dataclass
class KPath:
    """
    Path through k-space for band structure computation.

    Attributes:
        k_points: Tensor of k-vectors [n_k, d] in rad/m
        closed: Whether this is a closed loop (for Wilson loop holonomy)
        label: Human-readable description (e.g., "Γ→X→M→Γ")
        special_point_indices: Optional list of indices marking high-symmetry points

    Example:
        # 1D: Γ to X (Brillouin zone edge)
        k_path = KPath(
            k_points=torch.linspace(0, torch.pi/h, 100).unsqueeze(1),
            closed=False,
            label="Γ→X",
        )

        # 2D: Square Brillouin zone loop Γ→X→M→Γ
        n_seg = 50
        kx = torch.cat([
            torch.linspace(0, torch.pi/h, n_seg),           # Γ→X
            torch.full((n_seg,), torch.pi/h),               # X→M
            torch.linspace(torch.pi/h, 0, n_seg),           # M→Γ
        ])
        ky = torch.cat([
            torch.zeros(n_seg),                             # Γ→X
            torch.linspace(0, torch.pi/h, n_seg),           # X→M
            torch.zeros(n_seg),                             # M→Γ
        ])
        k_path = KPath(
            k_points=torch.stack([kx, ky], dim=1),
            closed=True,
            label="Γ→X→M→Γ",
            special_point_indices=[0, n_seg, 2*n_seg, 3*n_seg-1],
        )
    """
    k_points: torch.Tensor  # [n_k, d]
    closed: bool
    label: str = ""
    special_point_indices: Optional[List[int]] = None

    def __post_init__(self):
        """Validate k-path."""
        if self.k_points.ndim != 2:
            raise ValueError(
                f"k_points must be 2D [n_k, d], got shape {self.k_points.shape}"
            )

        if self.k_points.shape[0] < 2:
            raise ValueError(
                f"Need at least 2 k-points, got {self.k_points.shape[0]}"
            )

    @property
    def n_k(self) -> int:
        """Number of k-points."""
        return self.k_points.shape[0]

    @property
    def d(self) -> int:
        """Dimension of k-space."""
        return self.k_points.shape[1]


@dataclass
class SymplecticBandResult:
    """
    Results from symplectic band structure computation.

    Contains eigenfrequencies and polarization frames (q-part of eigenvectors)
    along a k-path.

    Attributes:
        omega: Eigenfrequencies [n_k, n_modes] in rad/s
               Sorted by increasing frequency at each k
        frames_q: Polarization frames [n_k, embedding_dim, n_modes]
                  The q-part (position) of each eigenmode, normalized
        kpath: Original KPath used for computation
        meta: Dictionary with additional info:
              - normalization: Method used to normalize frames
              - degeneracy_clusters: List of degenerate mode groups per k
              - config: SymplecticBandConfig used

    Note:
        frames_q[j, :, α] is the polarization vector for mode α at k-point j.
        This is the "q" part of the full symplectic eigenvector (q, p).
    """
    omega: torch.Tensor  # [n_k, n_modes]
    frames_q: torch.Tensor  # [n_k, embedding_dim, n_modes]
    kpath: KPath
    meta: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate result."""
        n_k = self.omega.shape[0]
        n_modes = self.omega.shape[1]

        if self.frames_q.shape[0] != n_k:
            raise ValueError(
                f"frames_q has {self.frames_q.shape[0]} k-points, "
                f"but omega has {n_k}"
            )

        if self.frames_q.shape[2] != n_modes:
            raise ValueError(
                f"frames_q has {self.frames_q.shape[2]} modes, "
                f"but omega has {n_modes}"
            )

        if self.kpath.n_k != n_k:
            raise ValueError(
                f"kpath has {self.kpath.n_k} points, but result has {n_k}"
            )

    @property
    def n_k(self) -> int:
        """Number of k-points."""
        return self.omega.shape[0]

    @property
    def n_modes(self) -> int:
        """Number of modes tracked."""
        return self.omega.shape[1]

    @property
    def embedding_dim(self) -> int:
        """Embedding dimension."""
        return self.frames_q.shape[1]


@dataclass
class SymplecticWilsonResult:
    """
    Results from symplectic Wilson loop holonomy computation.

    Computed from the q-frames (polarization vectors) of a degenerate subspace
    along a closed k-path.

    Attributes:
        W: Wilson loop matrix [N, N] (numpy array)
           W = ∏_j U_j† U_{j+1} for the selected band subspace
        trace: tr(W), gauge-invariant
        eigenvalues: Eigenvalues of W (complex)
        eigenphases: Phases of eigenvalues in [-π, π]
        band_indices: Indices of bands included in the subspace
        kpath: Original KPath (must be closed)
        stability_min_overlap: Minimum |⟨u_j, u_{j+1}⟩| along path
        distance_to_identity: Normalized distance ||W - I||_F / ||I||_F
        meta: Additional diagnostics

    Note:
        For U(1) case (N=1), eigenphases[0] is the Berry phase.
        For U(N) case (N>1), eigenphases encode non-Abelian holonomy.
    """
    W: Any  # numpy array [N, N]
    trace: complex
    eigenvalues: Any  # numpy array [N]
    eigenphases: Any  # numpy array [N]
    band_indices: List[int]
    kpath: KPath
    stability_min_overlap: float
    distance_to_identity: float
    meta: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate result."""
        if not self.kpath.closed:
            raise ValueError("Wilson loop requires closed k-path")

        N = len(self.band_indices)
        if self.W.shape != (N, N):
            raise ValueError(
                f"W matrix has shape {self.W.shape}, expected ({N}, {N})"
            )

    @property
    def dimension(self) -> int:
        """Subspace dimension N."""
        return len(self.band_indices)

    def is_u1(self) -> bool:
        """Check if this is U(1) (single band)."""
        return self.dimension == 1

    def berry_phase(self) -> Optional[float]:
        """
        Extract Berry phase (only valid for U(1) case).

        Returns:
            Berry phase in radians (in [-π, π]) if N=1, else None
        """
        if self.is_u1():
            return float(self.eigenphases[0])
        return None