"""
Symplectic band structure solver.

This module computes band structures by solving the eigenvalue problem for the
first-order symplectic operator A(k) along a k-path.

For each k, we find eigenpairs (λ, v) where:
    A(k) v = λ v

The eigenvalues come in ±iω pairs (symplectic structure). We extract the
positive-frequency modes (imag(λ) > 0) and their q-frames (position part of
eigenvector).

The solver is dimension-agnostic and works for 1D, 2D, and 3D.

Key features:
- Degeneracy-aware mode ordering for stable band tracking
- Proper normalization of q-frames
- Metadata export for diagnostics
"""

import torch
import numpy as np
from typing import Optional, List, Tuple
from dataclasses import dataclass

from .symplectic_types import (
    SymplecticBandConfig,
    KPath,
    SymplecticBandResult,
)
from .symplectic_builder import build_symplectic_operator_at_k


def solve_symplectic_eigenvalue_problem(
    A: torch.Tensor,
    select_positive_freq: bool = True,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Solve eigenvalue problem for symplectic operator A(k).

    Eigenvalues of a symplectic matrix come in pairs ±λ. For our first-order
    system with real frequencies, eigenvalues are ±iω, so we identify:
        ω = |imag(λ)|

    Args:
        A: Symplectic operator [2N, 2N] (complex)
        select_positive_freq: If True, return only positive-frequency modes
                             (imag(λ) > 0). Otherwise return all.

    Returns:
        (eigenvalues, eigenvectors):
            eigenvalues: [n_modes] complex tensor
            eigenvectors: [2N, n_modes] complex tensor (column vectors)

    Note:
        For positive frequencies, we expect imag(λ) > 0, which corresponds to
        oscillations e^{iωt} with ω > 0.
    """
    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = torch.linalg.eig(A)

    # eigenvalues: [2N] complex
    # eigenvectors: [2N, 2N] complex (column vectors)

    if select_positive_freq:
        # Select modes with imag(λ) > 0
        # Use small tolerance to handle numerical noise
        pos_freq_mask = eigenvalues.imag > 1e-12

        if pos_freq_mask.sum() == 0:
            # Fallback: if no positive frequencies, take all and sort by |imag|
            # This can happen for zero-frequency modes or numerical issues
            abs_imag = torch.abs(eigenvalues.imag)
            sorted_indices = torch.argsort(abs_imag, descending=True)
            n_half = len(eigenvalues) // 2
            selected_indices = sorted_indices[:n_half]
        else:
            selected_indices = torch.where(pos_freq_mask)[0]

        eigenvalues = eigenvalues[selected_indices]
        eigenvectors = eigenvectors[:, selected_indices]

    return eigenvalues, eigenvectors


def extract_q_frames(
    eigenvectors: torch.Tensor,
    embedding_dim: int,
    normalize: bool = True,
) -> torch.Tensor:
    """
    Extract q-frames (position part) from symplectic eigenvectors.

    The eigenvector has structure v = [q, p] where both q and p have length
    embedding_dim. We extract the q part and optionally normalize.

    Args:
        eigenvectors: [2*embedding_dim, n_modes] complex tensor
        embedding_dim: Embedding dimension (N)
        normalize: If True, normalize each q-frame to unit Euclidean norm

    Returns:
        frames_q: [embedding_dim, n_modes] complex tensor
                 Each column is a normalized polarization vector
    """
    # Extract q-part (first half)
    frames_q = eigenvectors[:embedding_dim, :]

    if normalize:
        # Normalize each column
        norms = torch.linalg.norm(frames_q, dim=0, keepdim=True)
        # Avoid division by zero
        norms = torch.clamp(norms, min=1e-14)
        frames_q = frames_q / norms

    return frames_q


def order_modes_by_frequency(
    eigenvalues: torch.Tensor,
    frames_q: torch.Tensor,
    degeneracy_threshold: float = 1e-3,
) -> Tuple[torch.Tensor, torch.Tensor, List[List[int]]]:
    """
    Order modes by increasing frequency with degeneracy detection.

    Args:
        eigenvalues: [n_modes] complex, expect ±iω structure
        frames_q: [embedding_dim, n_modes] polarization frames
        degeneracy_threshold: Relative frequency difference |ω_i - ω_j|/ω_mean < threshold
                             to consider modes degenerate

    Returns:
        (omega_sorted, frames_sorted, degeneracy_clusters):
            omega_sorted: [n_modes] real frequencies ω in increasing order
            frames_sorted: [embedding_dim, n_modes] reordered frames
            degeneracy_clusters: List of index lists for degenerate mode groups
    """
    # Extract frequencies: ω = |imag(λ)|
    omega = torch.abs(eigenvalues.imag)

    # Sort by increasing frequency
    sorted_indices = torch.argsort(omega)
    omega_sorted = omega[sorted_indices]
    frames_sorted = frames_q[:, sorted_indices]

    # Detect degenerate clusters
    degeneracy_clusters = []
    if len(omega_sorted) > 1:
        current_cluster = [0]
        for i in range(1, len(omega_sorted)):
            omega_mean = (omega_sorted[i] + omega_sorted[i-1]) / 2.0
            if omega_mean > 1e-12:
                rel_diff = torch.abs(omega_sorted[i] - omega_sorted[i-1]) / omega_mean
                if rel_diff < degeneracy_threshold:
                    # Degenerate with previous mode
                    current_cluster.append(i)
                else:
                    # End of cluster
                    if len(current_cluster) > 1:
                        degeneracy_clusters.append(current_cluster)
                    current_cluster = [i]
            else:
                # Both frequencies near zero, consider degenerate
                current_cluster.append(i)

        # Don't forget last cluster
        if len(current_cluster) > 1:
            degeneracy_clusters.append(current_cluster)

    return omega_sorted, frames_sorted, degeneracy_clusters


def solve_symplectic_bands_on_kpath(
    cfg: SymplecticBandConfig,
    kpath: KPath,
    n_modes: Optional[int] = None,
    degeneracy_threshold: float = 1e-3,
    verbose: bool = True,
) -> SymplecticBandResult:
    """
    Solve symplectic band structure along a k-path.

    For each k-point, constructs A(k), solves eigenvalue problem, and extracts
    positive-frequency modes with their q-frames (polarization vectors).

    Args:
        cfg: SymplecticBandConfig with lattice and material parameters
        kpath: KPath defining the path through k-space
        n_modes: Number of modes to keep (None = all positive-frequency modes)
        degeneracy_threshold: Threshold for detecting degenerate modes
        verbose: If True, print progress

    Returns:
        SymplecticBandResult with omega[n_k, n_modes] and frames_q[n_k, embedding_dim, n_modes]

    Example:
        >>> import torch
        >>> from branesim.diagnostics.bands.symplectic_types import (
        ...     SymplecticBandConfig, BoundaryCondition, KPath
        ... )
        >>> cfg = SymplecticBandConfig(
        ...     d=1, embedding_dim=4, grid_shape=(64,),
        ...     spacing=1e-12, mass=1e-30, spring_k=1e-6, rest_length=1e-12,
        ...     neighbor_offsets=[(-1,), (1,)],
        ...     boundary=BoundaryCondition.PERIODIC,
        ... )
        >>> k_points = torch.linspace(0, torch.pi/cfg.spacing, 50).unsqueeze(1)
        >>> kpath = KPath(k_points=k_points, closed=False, label="Γ→X")
        >>> result = solve_symplectic_bands_on_kpath(cfg, kpath, n_modes=4)
        >>> result.omega.shape
        torch.Size([50, 4])
        >>> result.frames_q.shape
        torch.Size([50, 4, 4])  # [n_k, embedding_dim, n_modes]
    """
    # Validate inputs
    if kpath.d != cfg.d:
        raise ValueError(
            f"k-path dimension ({kpath.d}) must match config dimension ({cfg.d})"
        )

    n_k = kpath.n_k
    embedding_dim = cfg.embedding_dim

    # Determine number of modes to track
    if n_modes is None:
        # All positive-frequency modes (half of total DOFs)
        n_modes = embedding_dim
    else:
        if n_modes > embedding_dim:
            raise ValueError(
                f"n_modes ({n_modes}) cannot exceed embedding_dim ({embedding_dim})"
            )

    if verbose:
        print(f"\nSolving symplectic bands on k-path:")
        print(f"  Path: {kpath.label}")
        print(f"  Number of k-points: {n_k}")
        print(f"  Embedding dimension: {embedding_dim}")
        print(f"  Modes to track: {n_modes}")
        print(f"  Boundary condition: {cfg.boundary.value}")

    # Allocate storage
    omega_all = torch.zeros(n_k, n_modes, device=cfg.device, dtype=cfg.dtype)
    frames_q_all = torch.zeros(
        n_k, embedding_dim, n_modes,
        device=cfg.device,
        dtype=torch.complex128 if cfg.dtype == torch.float64 else torch.complex64
    )

    # Store degeneracy info per k
    all_degeneracy_clusters = []

    # Loop over k-points
    for k_idx in range(n_k):
        k = kpath.k_points[k_idx]

        # Build A(k)
        A = build_symplectic_operator_at_k(cfg, k)

        # Solve eigenvalue problem
        eigenvalues, eigenvectors = solve_symplectic_eigenvalue_problem(
            A, select_positive_freq=True
        )

        # Extract q-frames
        frames_q = extract_q_frames(eigenvectors, embedding_dim, normalize=True)

        # Order by frequency
        omega_sorted, frames_sorted, degeneracy_clusters = order_modes_by_frequency(
            eigenvalues, frames_q, degeneracy_threshold
        )

        # Keep only first n_modes
        omega_all[k_idx] = omega_sorted[:n_modes].real
        frames_q_all[k_idx] = frames_sorted[:, :n_modes]
        all_degeneracy_clusters.append(degeneracy_clusters)

        if verbose and (k_idx % max(1, n_k // 10) == 0):
            print(f"  Progress: {k_idx+1}/{n_k} k-points")

    if verbose:
        print(f"  Done!")
        print(f"\n  Frequency range: [{omega_all.min():.3e}, {omega_all.max():.3e}] rad/s")

        # Count degeneracies
        n_degenerate_points = sum(1 for clusters in all_degeneracy_clusters if len(clusters) > 0)
        if n_degenerate_points > 0:
            print(f"  Degeneracies detected at {n_degenerate_points}/{n_k} k-points")

    # Build result
    result = SymplecticBandResult(
        omega=omega_all,
        frames_q=frames_q_all,
        kpath=kpath,
        meta={
            "normalization": "euclidean",
            "degeneracy_clusters": all_degeneracy_clusters,
            "degeneracy_threshold": degeneracy_threshold,
            "config": cfg,
        }
    )

    return result


def get_band_velocities(
    result: SymplecticBandResult,
) -> torch.Tensor:
    """
    Compute group velocities v_g = dω/dk for each band.

    Uses finite differences along the k-path. For closed loops, uses periodic
    differences.

    Args:
        result: SymplecticBandResult from solve_symplectic_bands_on_kpath

    Returns:
        velocities: [n_k, n_modes] tensor of group velocity magnitudes |v_g|

    Note:
        For multi-dimensional k-paths, this returns |dω/dk| where k is the
        path parameter (not a directional derivative in k-space).
    """
    n_k = result.n_k
    n_modes = result.n_modes

    omega = result.omega  # [n_k, n_modes]
    k_points = result.kpath.k_points  # [n_k, d]

    # Compute path parameter s = |k - k_0| (distance along path)
    dk = torch.diff(k_points, dim=0)  # [n_k-1, d]
    ds = torch.linalg.norm(dk, dim=1)  # [n_k-1]

    # Frequency differences
    domega = torch.diff(omega, dim=0)  # [n_k-1, n_modes]

    # Group velocities
    vg = domega / ds.unsqueeze(1)  # [n_k-1, n_modes]

    # Pad to match original size (use forward difference at boundaries)
    if result.kpath.closed:
        # Periodic: wrap around
        vg_padded = torch.cat([vg, vg[:1]], dim=0)
    else:
        # Open: duplicate last value
        vg_padded = torch.cat([vg, vg[-1:]], dim=0)

    return vg_padded.abs()