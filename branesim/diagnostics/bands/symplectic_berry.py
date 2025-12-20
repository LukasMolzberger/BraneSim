"""
Berry phase and Wilson loop diagnostics for symplectic bands.

This module computes gauge-invariant holonomy (Berry phase / Wilson loop)
from the q-frames extracted from symplectic band structure calculations.

Key principle:
    Berry transport is computed from EIGENFRAMES, not from time-domain signals.
    The frames come from the band solver, which provides the polarization
    vectors (q-part of eigenmodes) at each k-point along a path.

For U(1) (single band):
    Berry phase γ = arg(∏_j ⟨u_j | u_{j+1}⟩)

For U(N) (degenerate subspace):
    Wilson loop W = ∏_j U_j† U_{j+1}
    Eigenphases of W encode non-Abelian holonomy

All computations are dimension-agnostic and work for 1D, 2D, and 3D k-spaces.
"""

import torch
import numpy as np
from typing import List, Optional

from .symplectic_types import (
    SymplecticBandResult,
    SymplecticWilsonResult,
)

# Import existing holonomy functions (reuse where possible)
from branesim.diagnostics.holonomy import (
    wilson_loop_holonomy,
    orthonormalize_frame,
    WilsonLoopResult,
)


def compute_symplectic_wilson_loop(
    band_result: SymplecticBandResult,
    band_indices: List[int],
    reorthonormalize: bool = True,
    verbose: bool = True,
) -> SymplecticWilsonResult:
    """
    Compute Wilson loop holonomy for a subspace of symplectic bands.

    Takes the q-frames (polarization vectors) for the selected bands along
    a closed k-path and computes the non-Abelian Wilson loop:

        W = ∏_j M_j  where  M_j = U_j† U_{j+1}

    For U(1) case (single band), W is a scalar and its phase is the Berry phase.
    For U(N) case (N > 1), W is an N×N unitary matrix encoding spinorial transport.

    Args:
        band_result: SymplecticBandResult with frames_q from band solver
        band_indices: List of band indices to include in subspace (e.g., [0, 1] for
                     first two bands). Must have at least 1 element.
        reorthonormalize: If True, QR-orthonormalize frames before computing overlaps
                         (recommended for numerical stability)
        verbose: If True, print diagnostic information

    Returns:
        SymplecticWilsonResult with Wilson loop matrix, eigenphases, and diagnostics

    Raises:
        ValueError: If k-path is not closed
        ValueError: If band_indices contains invalid indices

    Example:
        >>> # Single band (U(1) Berry phase)
        >>> result = compute_symplectic_wilson_loop(band_result, band_indices=[2])
        >>> berry_phase = result.berry_phase()
        >>> print(f"Berry phase: {berry_phase:.6f} rad")

        >>> # Degenerate pair (U(2) Wilson loop)
        >>> result = compute_symplectic_wilson_loop(band_result, band_indices=[0, 1])
        >>> print(f"Wilson eigenphases: {result.eigenphases}")
        >>> print(f"Spinorial: {result.W}")
    """
    # Validate inputs
    if not band_result.kpath.closed:
        raise ValueError(
            "Wilson loop requires a closed k-path. "
            "Set kpath.closed=True or use an open path for Berry connection only."
        )

    n_bands = len(band_indices)
    if n_bands == 0:
        raise ValueError("Must provide at least one band index")

    for idx in band_indices:
        if idx < 0 or idx >= band_result.n_modes:
            raise ValueError(
                f"Band index {idx} out of range [0, {band_result.n_modes})"
            )

    # Extract frames for selected bands
    # frames_q has shape [n_k, embedding_dim, n_modes]
    # We want frames for bands in band_indices: [n_k, embedding_dim, n_bands]
    frames_selected = band_result.frames_q[:, :, band_indices]
    # frames_selected.shape = [n_k, embedding_dim, n_bands]

    n_k = frames_selected.shape[0]
    embedding_dim = frames_selected.shape[1]

    if verbose:
        print(f"\n{'='*70}")
        print(f"Computing Wilson loop holonomy")
        print(f"{'='*70}")
        print(f"  k-path: {band_result.kpath.label}")
        print(f"  Number of k-points: {n_k}")
        print(f"  Subspace dimension: {n_bands}")
        print(f"  Band indices: {band_indices}")
        print(f"  Embedding dimension: {embedding_dim}")

    # Convert to list of frames for holonomy function
    # Each frame is [embedding_dim, n_bands] for compatibility with holonomy.py
    frames_list: List[torch.Tensor] = []
    for k_idx in range(n_k):
        frame_k = frames_selected[k_idx, :, :]  # [embedding_dim, n_bands]
        frames_list.append(frame_k)

    # Compute Wilson loop using existing holonomy function
    if verbose:
        print(f"\n  Calling wilson_loop_holonomy...")

    wilson_result: WilsonLoopResult = wilson_loop_holonomy(
        frames=frames_list,
        reorthonormalize=reorthonormalize,
        gauge_check=True,
    )

    # Compute stability metric: minimum overlap along path
    min_overlap = compute_frame_overlap_stability(frames_list)

    # Compute distance to identity
    I = np.eye(n_bands)
    dist_to_identity = np.linalg.norm(wilson_result.W - I, ord='fro') / np.linalg.norm(I, ord='fro')

    if verbose:
        print(f"\n  Wilson loop computed successfully!")
        print(f"  Minimum frame overlap: {min_overlap:.6f}")
        print(f"  Distance to identity: {dist_to_identity:.6f}")

        if n_bands == 1:
            berry_phase_rad = float(wilson_result.eigenphases[0])
            berry_phase_pi = berry_phase_rad / np.pi
            print(f"\n  U(1) Berry phase: {berry_phase_rad:.6f} rad ({berry_phase_pi:.4f} π)")
        else:
            print(f"\n  U({n_bands}) Wilson loop:")
            print(f"    Trace: {wilson_result.trace:.6f}")
            print(f"    Eigenphases: {wilson_result.eigenphases} rad")
            print(f"    Spinorial (W ≈ -I): {wilson_result.is_spinorial}")

    # Build result
    result = SymplecticWilsonResult(
        W=wilson_result.W,
        trace=wilson_result.trace,
        eigenvalues=wilson_result.eigenvalues,
        eigenphases=wilson_result.eigenphases,
        band_indices=band_indices,
        kpath=band_result.kpath,
        stability_min_overlap=min_overlap,
        distance_to_identity=dist_to_identity,
        meta={
            "reorthonormalized": reorthonormalize,
            "dimension": n_bands,
            "is_spinorial": wilson_result.is_spinorial,
            "distance_to_minus_I": wilson_result.distance_to_minus_I,
        }
    )

    return result


def compute_frame_overlap_stability(frames: List[torch.Tensor]) -> float:
    """
    Compute minimum overlap |⟨u_j, u_{j+1}⟩| along the path.

    For a well-conditioned Berry phase computation, frames should vary smoothly,
    meaning adjacent frames should have large overlap.

    Args:
        frames: List of frames [d, N] at each k-point

    Returns:
        min_overlap: Minimum |⟨u_j, u_{j+1}⟩| over all j

    Note:
        For U(1) case (N=1), this is just |⟨u_j | u_{j+1}⟩|.
        For U(N) case, we use the minimum singular value of U_j† U_{j+1}.
    """
    n_points = len(frames)
    overlaps = []

    for j in range(n_points):
        U_j = frames[j]
        U_j1 = frames[(j + 1) % n_points]  # Periodic

        # Overlap matrix: U_j† U_{j+1}
        M = U_j.T.conj() @ U_j1

        if M.shape == (1, 1):
            # U(1) case: scalar overlap
            overlap = torch.abs(M[0, 0]).item()
        else:
            # U(N) case: minimum singular value
            singular_values = torch.linalg.svdvals(M)
            overlap = singular_values.min().item()

        overlaps.append(overlap)

    return float(min(overlaps))


def compute_berry_connection_along_path(
    band_result: SymplecticBandResult,
    band_index: int,
) -> torch.Tensor:
    """
    Compute Berry connection A_k = i⟨u | d/dk u⟩ along k-path for a single band.

    This is the local (gauge-dependent) Berry connection. Its integral gives
    the Berry phase for a closed loop.

    Args:
        band_result: SymplecticBandResult with frames_q
        band_index: Index of band to analyze (must be non-degenerate)

    Returns:
        connection: [n_k] tensor of Berry connection values at each k-point

    Note:
        For U(N) case, use compute_symplectic_wilson_loop instead.
        This function is only for single bands (U(1)).
    """
    if band_index < 0 or band_index >= band_result.n_modes:
        raise ValueError(
            f"Band index {band_index} out of range [0, {band_result.n_modes})"
        )

    # Extract frame for this band
    frames = band_result.frames_q[:, :, band_index]  # [n_k, embedding_dim]

    n_k = frames.shape[0]

    # Compute finite-difference derivative du/dk
    # Use central differences in interior, forward/backward at boundaries
    connection = torch.zeros(n_k, device=frames.device, dtype=frames.real.dtype)

    for j in range(n_k):
        u_j = frames[j]

        # Get adjacent k-points (with periodic or open boundary handling)
        if band_result.kpath.closed:
            # Periodic
            j_next = (j + 1) % n_k
            j_prev = (j - 1) % n_k
            u_next = frames[j_next]
            u_prev = frames[j_prev]

            # Central difference
            du_dk = (u_next - u_prev) / 2.0
        else:
            # Open boundaries
            if j == 0:
                # Forward difference
                u_next = frames[j + 1]
                du_dk = u_next - u_j
            elif j == n_k - 1:
                # Backward difference
                u_prev = frames[j - 1]
                du_dk = u_j - u_prev
            else:
                # Central difference
                u_next = frames[j + 1]
                u_prev = frames[j - 1]
                du_dk = (u_next - u_prev) / 2.0

        # Berry connection: A = i ⟨u | du/dk⟩
        # Note: imaginary part of ⟨u | du/dk⟩ gives the connection
        overlap = torch.vdot(u_j, du_dk)  # ⟨u | du/dk⟩
        connection[j] = overlap.imag

    return connection


def verify_gauge_covariance(
    band_result: SymplecticBandResult,
    band_indices: List[int],
    n_tests: int = 3,
) -> bool:
    """
    Verify that Wilson loop invariants are gauge-covariant.

    Applies random local unitary transformations to frames and checks that
    Wilson loop eigenvalues remain unchanged (up to numerical tolerance).

    Args:
        band_result: SymplecticBandResult with frames
        band_indices: Bands to include in subspace
        n_tests: Number of random gauge transformations to test

    Returns:
        True if all tests pass, False otherwise

    Note:
        This is a development/debugging tool. For production, disable or
        use sparingly as it's expensive (recomputes Wilson loop n_tests times).
    """
    print(f"\n{'='*70}")
    print(f"Gauge covariance verification")
    print(f"{'='*70}")

    # Compute original Wilson loop
    result_orig = compute_symplectic_wilson_loop(
        band_result, band_indices, verbose=False
    )

    eigenphases_orig = sorted(result_orig.eigenphases)
    trace_orig = result_orig.trace

    print(f"  Original trace: {trace_orig:.6f}")
    print(f"  Original eigenphases: {eigenphases_orig}")

    all_pass = True

    for test_idx in range(n_tests):
        # Apply random gauge transformation
        frames_transformed = apply_random_gauge_transformation(
            band_result.frames_q[:, :, band_indices]
        )

        # Create modified band result
        band_result_transformed = SymplecticBandResult(
            omega=band_result.omega,
            frames_q=frames_transformed,
            kpath=band_result.kpath,
            meta=band_result.meta,
        )

        # Compute Wilson loop with transformed frames
        result_transformed = compute_symplectic_wilson_loop(
            band_result_transformed,
            band_indices=list(range(len(band_indices))),  # Indices in transformed array
            verbose=False,
        )

        # Compare invariants
        eigenphases_new = sorted(result_transformed.eigenphases)
        trace_new = result_transformed.trace

        trace_error = abs(trace_new - trace_orig)
        eigenphase_error = np.max(np.abs(np.array(eigenphases_new) - np.array(eigenphases_orig)))

        print(f"\n  Test {test_idx + 1}:")
        print(f"    Δtrace: {trace_error:.2e}")
        print(f"    Δeigenphase_max: {eigenphase_error:.2e}")

        if trace_error > 1e-6 or eigenphase_error > 1e-6:
            print(f"    FAIL: Invariants changed significantly")
            all_pass = False
        else:
            print(f"    PASS")

    if all_pass:
        print(f"\n  ✓ All {n_tests} gauge covariance tests passed")
    else:
        print(f"\n  ✗ Some gauge covariance tests FAILED")

    return all_pass


def apply_random_gauge_transformation(
    frames: torch.Tensor,
) -> torch.Tensor:
    """
    Apply random local gauge transformations to frames.

    For each k-point, applies U_k → U_k @ G_k where G_k ∈ U(N) is a random
    unitary matrix.

    Args:
        frames: [n_k, embedding_dim, n_bands] frames to transform

    Returns:
        transformed_frames: [n_k, embedding_dim, n_bands] with random gauge applied
    """
    n_k, embedding_dim, n_bands = frames.shape
    device = frames.device
    dtype = frames.dtype

    transformed_frames = frames.clone()

    for k_idx in range(n_k):
        # Generate random unitary G = exp(iH) where H is Hermitian
        H = torch.randn(n_bands, n_bands, device=device, dtype=dtype)
        H = (H + H.T.conj()) / 2.0  # Make Hermitian
        G = torch.matrix_exp(1j * H)

        # Apply transformation: U_k → U_k @ G
        transformed_frames[k_idx] = frames[k_idx] @ G

    return transformed_frames