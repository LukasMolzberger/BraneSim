"""
Symplectic operator builder for first-order band structure.

This module constructs the first-order symplectic linear operator A(k) that
governs small-amplitude oscillations of the brane substrate. The operator acts
on the phase-space vector (q, p) where:
    - q: position displacements in embedding space (embedding_dim components)
    - p: momenta (embedding_dim components)

The equations of motion in first-order form are:
    dq/dt = p/m
    dp/dt = F(q) ≈ -K·q  (linearized)

Combined into symplectic form:
    d/dt [q] = A [q]    where A = [    0       1/m  ]
         [p]       [p]              [ -K        0    ]

For a periodic lattice, A(k) is a function of the Bloch wave vector k.

Key principle:
    Boundary conditions affect HOW we build A(k), not Berry diagnostics later.
    - PERIODIC: Use Fourier phases exp(i k·r) for neighbor coupling
    - CLAMPED: Not compatible with k-space bands (error if requested)

All code is dimension-agnostic and works for 1D, 2D, and 3D.
"""

import torch
import numpy as np
from typing import Tuple

from .symplectic_types import SymplecticBandConfig, BoundaryCondition


def build_symplectic_operator_at_k(
    cfg: SymplecticBandConfig,
    k: torch.Tensor,
) -> torch.Tensor:
    """
    Build first-order symplectic operator A(k) for a given k-vector.

    For a single-node unit cell with N = embedding_dim degrees of freedom,
    the operator has size [2N, 2N] and block structure:

        A(k) = [    0       M^{-1}  ]
               [ -K(k)       0      ]

    where:
        - M^{-1} is diagonal: 1/m * I_N
        - K(k) is the stiffness matrix with Fourier phases from periodic BCs

    The stiffness matrix for harmonic springs is:
        K(k)_αβ = (spring_k / rest_length) * Σ_r [δ_αβ - r_α r_β / |r|²] * [1 - cos(k·r)]

    Args:
        cfg: SymplecticBandConfig with lattice parameters and boundary conditions
        k: Wave vector [d] in rad/m

    Returns:
        A(k): Symplectic operator [2*embedding_dim, 2*embedding_dim] (complex)

    Raises:
        ValueError: If boundary condition is CLAMPED (not supported for k-space)
        ValueError: If k dimension doesn't match cfg.d

    Example:
        >>> import torch
        >>> from branesim.diagnostics.bands.symplectic_types import (
        ...     SymplecticBandConfig, BoundaryCondition
        ... )
        >>> cfg = SymplecticBandConfig(
        ...     d=1, embedding_dim=4, grid_shape=(64,),
        ...     spacing=1e-12, mass=1e-30, spring_k=1e-6, rest_length=1e-12,
        ...     neighbor_offsets=[(-1,), (1,)],
        ...     boundary=BoundaryCondition.PERIODIC,
        ... )
        >>> k = torch.tensor([1e10])  # k ~ 1e10 rad/m
        >>> A = build_symplectic_operator_at_k(cfg, k)
        >>> A.shape
        torch.Size([8, 8])  # 2 * embedding_dim = 2 * 4 = 8
    """
    # Validate inputs
    if cfg.boundary == BoundaryCondition.CLAMPED:
        raise ValueError(
            "CLAMPED boundary does not support k-space band computation. "
            "Use PERIODIC or implement finite-domain mode tracking separately."
        )

    if k.shape != (cfg.d,):
        raise ValueError(
            f"k vector has shape {k.shape}, expected ({cfg.d},) to match d={cfg.d}"
        )

    N = cfg.embedding_dim
    device = cfg.device
    dtype = cfg.dtype

    # Complex dtype for A(k)
    if dtype == torch.float64:
        cdtype = torch.complex128
    else:
        cdtype = torch.complex64

    # Initialize A(k) as 2N × 2N complex matrix
    A = torch.zeros(2 * N, 2 * N, device=device, dtype=cdtype)

    # Upper-right block: M^{-1} = (1/m) * I_N
    inv_mass = 1.0 / cfg.mass
    A[:N, N:] = inv_mass * torch.eye(N, device=device, dtype=cdtype)

    # Lower-left block: -K(k)
    # For harmonic springs connecting neighbors via offsets r:
    #   K(k)_αβ = (spring_k / rest_length) * Σ_r [δ_αβ - r̂_α r̂_β] * [1 - cos(k·r)]
    # where r̂ = r / |r| is unit vector along spring

    K = torch.zeros(N, N, device=device, dtype=cdtype)

    for offset in cfg.neighbor_offsets:
        # Convert offset (in grid units) to real-space displacement
        r_vec = torch.tensor(offset, device=device, dtype=dtype) * cfg.spacing

        # Compute |r|
        r_mag = torch.linalg.norm(r_vec)

        if r_mag < 1e-14:
            continue  # Skip self (shouldn't happen with proper offsets)

        # Unit vector r̂
        r_hat = r_vec / r_mag

        # Fourier phase factor: exp(i k·r) for periodic BC
        # k·r (dot product)
        k_dot_r = torch.dot(k.to(dtype), r_vec)

        # Spring contribution matrix: [δ_αβ - r̂_α r̂_β] * [1 - cos(k·r)]
        # This is a projection onto the spring direction
        #
        # For embedding space (4D typically), only the first d components
        # of r̂ are non-zero (lateral directions). The transverse directions
        # (amplitude) have no preferred direction from nearest-neighbor springs.

        # Build r̂ in embedding space (pad with zeros for transverse dimensions)
        r_hat_embed = torch.zeros(N, device=device, dtype=dtype)
        r_hat_embed[:cfg.d] = r_hat

        # Outer product r̂ ⊗ r̂
        r_outer = torch.outer(r_hat_embed, r_hat_embed)

        # Projection matrix P = I - r̂⊗r̂ (perpendicular to spring)
        I_N = torch.eye(N, device=device, dtype=dtype)
        P = I_N - r_outer

        # Fourier weight: 1 - cos(k·r) = 2 sin²(k·r/2)
        # For numerical stability, use 1 - cos directly
        fourier_weight = 1.0 - torch.cos(k_dot_r)

        # Spring stiffness factor
        stiffness_factor = cfg.spring_k / cfg.rest_length

        # Add contribution: K += stiffness * P * weight
        K += stiffness_factor * P.to(cdtype) * fourier_weight

    # Lower-left block: -K(k)
    A[N:, :N] = -K

    # Upper-left and lower-right blocks are zero (already initialized)

    return A


def get_nearest_neighbor_offsets(d: int, include_diagonals: bool = False) -> list:
    """
    Generate nearest-neighbor offset list for d-dimensional lattice.

    Args:
        d: Dimension (1, 2, or 3)
        include_diagonals: If True, include diagonal neighbors (8-connectivity for 2D,
                          26-connectivity for 3D). If False, only axis-aligned (4 for 2D, 6 for 3D).

    Returns:
        List of tuples, each a d-dimensional offset vector

    Example:
        >>> get_nearest_neighbor_offsets(1)
        [(-1,), (1,)]
        >>> get_nearest_neighbor_offsets(2, include_diagonals=False)
        [(-1, 0), (1, 0), (0, -1), (0, 1)]
        >>> len(get_nearest_neighbor_offsets(2, include_diagonals=True))
        8
        >>> len(get_nearest_neighbor_offsets(3, include_diagonals=False))
        6
    """
    offsets = []

    if d == 1:
        # 1D: just left and right
        offsets = [(-1,), (1,)]

    elif d == 2:
        if include_diagonals:
            # 8-connectivity: all 3×3 - 1
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di != 0 or dj != 0:
                        offsets.append((di, dj))
        else:
            # 4-connectivity: axis-aligned only
            offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    elif d == 3:
        if include_diagonals:
            # 26-connectivity: all 3×3×3 - 1
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    for dk in [-1, 0, 1]:
                        if di != 0 or dj != 0 or dk != 0:
                            offsets.append((di, dj, dk))
        else:
            # 6-connectivity: axis-aligned only
            offsets = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]

    else:
        raise ValueError(f"Dimension d must be 1, 2, or 3, got {d}")

    return offsets


def validate_operator_symplecticity(A: torch.Tensor, tol: float = 1e-8) -> Tuple[bool, float]:
    """
    Check if operator A is symplectic: A^T J A = J where J = [[0, I], [-I, 0]].

    For a symplectic operator, the eigenvalues come in ±iω pairs. This function
    checks the algebraic condition directly.

    Args:
        A: Operator [2N, 2N]
        tol: Tolerance for ||A^T J A - J||

    Returns:
        (is_symplectic, error): Boolean and norm of deviation

    Note:
        This is a consistency check during development. For performance-critical
        code, this can be disabled.
    """
    N_half = A.shape[0] // 2
    if A.shape != (2 * N_half, 2 * N_half):
        raise ValueError(f"A must be square with even dimension, got {A.shape}")

    device = A.device
    dtype = A.dtype

    # Build symplectic form J = [[0, I], [-I, 0]]
    I = torch.eye(N_half, device=device, dtype=dtype)
    Z = torch.zeros(N_half, N_half, device=device, dtype=dtype)
    J = torch.cat([
        torch.cat([Z, I], dim=1),
        torch.cat([-I, Z], dim=1)
    ], dim=0)

    # Compute A^T J A
    AtJA = A.T.conj() @ J @ A

    # Check ||A^T J A - J||
    error = torch.linalg.norm(AtJA - J).item()

    is_symplectic = (error < tol)

    return is_symplectic, error