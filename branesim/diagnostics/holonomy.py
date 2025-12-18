"""
Holonomy computation: U(1) and non-Abelian Wilson loops.

This module implements gauge-invariant holonomy measurements for both
non-degenerate (U(1)) and degenerate (U(N) / Wilczek-Zee) transport.

Key principle:
    Holonomy is computed from overlap matrices between frames at successive
    points along a loop. For degenerate subspaces, this yields a non-Abelian
    Wilson loop matrix W ∈ U(N) whose gauge-invariant properties (trace,
    eigenvalues) encode spinorial transport signatures.
"""

import torch
import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class WilsonLoopResult:
    """Results from non-Abelian Wilson loop computation."""
    W: np.ndarray  # Wilson loop matrix [N×N]
    trace: complex  # tr(W), gauge-invariant
    eigenvalues: np.ndarray  # Eigenvalues of W
    eigenphases: np.ndarray  # Phases of eigenvalues (in [-π, π])
    distance_to_minus_I: float  # ||W + I||_F / ||I||_F
    distance_to_plus_I: float  # ||W - I||_F / ||I||_F
    is_spinorial: bool  # True if W ≈ -I (one circuit sign flip)
    dimension: int  # Dimension N of the subspace


@dataclass
class U1HolonomyResult:
    """Results from U(1) holonomy computation."""
    gamma: float  # Holonomy angle Γ ∈ [-π, π]
    total_phase: float  # Accumulated phase (may be >2π for multiple loops)
    path_length: float  # Total path length
    average_curvature: float  # Γ / (enclosed area), if applicable


def orthonormalize_frame(
    frame: torch.Tensor,
    method: str = "qr",
) -> torch.Tensor:
    """
    Orthonormalize a frame of vectors.

    Takes a matrix U = [u1, u2, ..., uN] where columns are vectors in the
    subspace, and returns an orthonormal basis via QR decomposition or
    Gram-Schmidt.

    Args:
        frame: Input frame [d × N] where d is ambient dimension, N is subspace dim
        method: "qr" for QR decomposition or "gs" for Gram-Schmidt

    Returns:
        Orthonormal frame [d × N]

    Note:
        QR is more numerically stable than Gram-Schmidt for nearly-parallel vectors.
    """
    if method == "qr":
        Q, R = torch.linalg.qr(frame)
        # Ensure consistent phase: make diagonal of R positive
        signs = torch.sign(torch.diag(R))
        signs[signs == 0] = 1.0
        Q = Q * signs.unsqueeze(0)
        return Q
    elif method == "gs":
        # Gram-Schmidt (less stable but explicit)
        d, N = frame.shape
        Q = torch.zeros_like(frame)
        for i in range(N):
            v = frame[:, i].clone()
            for j in range(i):
                v = v - torch.dot(Q[:, j], frame[:, i]) * Q[:, j]
            norm = torch.linalg.norm(v)
            if norm < 1e-12:
                raise ValueError(f"Frame vector {i} is linearly dependent")
            Q[:, i] = v / norm
        return Q
    else:
        raise ValueError(f"Unknown orthonormalization method: {method}")


def wilson_loop_holonomy(
    frames: List[torch.Tensor],
    reorthonormalize: bool = True,
    gauge_check: bool = True,
) -> WilsonLoopResult:
    """
    Compute non-Abelian Wilson loop holonomy W ∈ U(N).

    For a degenerate N-dimensional subspace sampled at points along a closed
    loop, computes the ordered product of overlap matrices:

        W = ∏_j M_j  where  M_j = U_j† U_{j+1}

    and U_j is an orthonormal frame at sample point j.

    The Wilson loop is gauge-invariant: tr(W) and the eigenvalues of W are
    independent of the local basis choice at each point.

    Args:
        frames: List of orthonormal frames [d × N] at each sample point
                along the loop. Last frame should be close to first (periodic).
        reorthonormalize: If True, QR-orthonormalize each frame before computing overlaps
        gauge_check: If True, verify gauge invariance by checking ||U_j† U_j - I|| < tol

    Returns:
        WilsonLoopResult with W matrix and gauge-invariant diagnostics

    Paper reference:
        Section "Holonomy measurement", paragraph "Degenerate transport (U(N) / Wilczek--Zee)"
        Equation: W = ∏_j M_j ∈ U(N)

    Note:
        For electron spinorial transport, target signature is:
            - One circuit: W ≈ -I  (state sign flip)
            - Two circuits: W² ≈ I (return to original state)
    """
    n_points = len(frames)
    if n_points < 2:
        raise ValueError("Need at least 2 frames to compute Wilson loop")

    device = frames[0].device
    dtype = frames[0].dtype
    d, N = frames[0].shape

    # Validate and optionally reorthonormalize
    frames_ortho = []
    for i, U in enumerate(frames):
        if U.shape != (d, N):
            raise ValueError(f"Frame {i} has shape {U.shape}, expected ({d}, {N})")

        if reorthonormalize:
            U_ortho = orthonormalize_frame(U, method="qr")
        else:
            U_ortho = U

        # Gauge check: verify orthonormality
        if gauge_check:
            overlap = U_ortho.T.conj() @ U_ortho
            identity = torch.eye(N, device=device, dtype=dtype)
            error = torch.linalg.norm(overlap - identity).item()
            if error > 1e-6:
                print(f"  WARNING: Frame {i} not orthonormal, error = {error:.2e}")

        frames_ortho.append(U_ortho)

    # Compute Wilson loop: W = ∏_j M_j where M_j = U_j† U_{j+1}
    print(f"\nComputing Wilson loop holonomy...")
    print(f"  Number of sample points: {n_points}")
    print(f"  Subspace dimension N: {N}")

    W = torch.eye(N, device=device, dtype=dtype)

    for j in range(n_points):
        U_j = frames_ortho[j]
        U_j1 = frames_ortho[(j + 1) % n_points]  # Periodic: wrap around

        # Overlap matrix M_j = U_j† U_{j+1}
        M_j = U_j.T.conj() @ U_j1

        # Accumulate: W = W * M_j
        W = W @ M_j

    # Extract gauge-invariant quantities
    W_np = W.cpu().numpy()

    # Trace (gauge-invariant)
    trace = np.trace(W_np)

    # Eigenvalues (gauge-invariant)
    eigenvalues = np.linalg.eigvals(W_np)
    eigenphases = np.angle(eigenvalues)  # Phases in [-π, π]

    # Distance to ±I
    I = np.eye(N)
    dist_minus_I = np.linalg.norm(W_np + I, ord='fro') / np.linalg.norm(I, ord='fro')
    dist_plus_I = np.linalg.norm(W_np - I, ord='fro') / np.linalg.norm(I, ord='fro')

    # Check for spinorial signature: W ≈ -I
    is_spinorial = (dist_minus_I < 0.1)  # Threshold can be adjusted

    print(f"\n  Wilson loop results:")
    print(f"    tr(W) = {trace:.6f}")
    print(f"    Eigenvalues: {eigenvalues}")
    print(f"    Eigenphases: {eigenphases} rad")
    print(f"    ||W + I|| / ||I||: {dist_minus_I:.6e}")
    print(f"    ||W - I|| / ||I||: {dist_plus_I:.6e}")
    print(f"    Spinorial (W ≈ -I): {'YES' if is_spinorial else 'NO'}")

    return WilsonLoopResult(
        W=W_np,
        trace=trace,
        eigenvalues=eigenvalues,
        eigenphases=eigenphases,
        distance_to_minus_I=dist_minus_I,
        distance_to_plus_I=dist_plus_I,
        is_spinorial=is_spinorial,
        dimension=N,
    )


def wilson_invariants(W: np.ndarray) -> Dict[str, any]:
    """
    Compute gauge-invariant diagnostics from a Wilson loop matrix.

    Args:
        W: Wilson loop matrix [N × N]

    Returns:
        Dictionary with invariants:
            - trace: tr(W)
            - determinant: det(W) (should be ≈ 1 for unitary W)
            - eigenvalues: eigenvalues of W
            - eigenphases: phases of eigenvalues
            - frobenius_norm: ||W||_F
    """
    trace = np.trace(W)
    det = np.linalg.det(W)
    eigenvalues = np.linalg.eigvals(W)
    eigenphases = np.angle(eigenvalues)
    frob_norm = np.linalg.norm(W, ord='fro')

    return {
        "trace": trace,
        "determinant": det,
        "eigenvalues": eigenvalues,
        "eigenphases": eigenphases,
        "frobenius_norm": frob_norm,
    }


def compute_u1_holonomy(
    phases: torch.Tensor,
    closed_loop: bool = True,
) -> U1HolonomyResult:
    """
    Compute U(1) holonomy from discrete link phases.

    For non-degenerate (N=1) transport, the holonomy is simply the accumulated
    phase Γ = Σ_j arg(U_j) where U_j are the link variables along the path.

    Args:
        phases: Link phases [n_links] in radians, each in [-π, π]
        closed_loop: If True, expect phases to wrap back (Γ mod 2π = 0)

    Returns:
        U1HolonomyResult with total phase and diagnostics

    Paper reference:
        Section "Holonomy measurement", paragraph "Non-degenerate transport (U(1))"
        Equation: Γ[γ] ≈ arg ∏_j U_{μ_j}
    """
    # Sum phases
    total_phase = torch.sum(phases).item()

    # Wrap to [-π, π] for the final holonomy angle
    gamma = np.angle(np.exp(1j * total_phase))

    # Path length (in number of steps)
    path_length = len(phases)

    # Average curvature (holonomy per unit length)
    avg_curvature = gamma / path_length if path_length > 0 else 0.0

    print(f"\nU(1) Holonomy results:")
    print(f"  Number of links: {path_length}")
    print(f"  Total phase: {total_phase:.6f} rad  ({total_phase / np.pi:.3f} π)")
    print(f"  Holonomy Γ: {gamma:.6f} rad  ({gamma / np.pi:.3f} π)")
    print(f"  Average curvature: {avg_curvature:.6e} rad/step")

    if closed_loop:
        closure_error = abs(gamma)
        print(f"  Closure error |Γ|: {closure_error:.6e}")
        if closure_error > 0.1:
            print(f"  WARNING: Loop not closed, |Γ| = {closure_error:.3f} > 0.1")

    return U1HolonomyResult(
        gamma=gamma,
        total_phase=total_phase,
        path_length=path_length,
        average_curvature=avg_curvature,
    )


def verify_gauge_invariance(
    frames: List[torch.Tensor],
    W_original: np.ndarray,
    n_random_tests: int = 5,
) -> bool:
    """
    Verify that Wilson loop invariants are gauge-invariant.

    Applies random local unitary transformations G_j to each frame:
        U_j → U_j @ G_j
    and checks that the Wilson loop invariants (trace, eigenvalues) remain unchanged.

    Args:
        frames: List of orthonormal frames [d × N]
        W_original: Original Wilson loop matrix
        n_random_tests: Number of random gauge transformations to test

    Returns:
        True if invariants are stable under gauge transformations

    Note:
        This is a crucial consistency check to ensure the implementation is correct.
    """
    device = frames[0].device
    dtype = frames[0].dtype
    d, N = frames[0].shape

    inv_original = wilson_invariants(W_original)
    trace_orig = inv_original["trace"]
    eigenvals_orig = sorted(np.abs(inv_original["eigenvalues"]))

    print(f"\nGauge invariance verification:")
    print(f"  Original tr(W) = {trace_orig:.6f}")

    all_pass = True

    for test_idx in range(n_random_tests):
        # Generate random local gauge transformations G_j ∈ U(N)
        frames_transformed = []
        for U_j in frames:
            # Random unitary: G = exp(iH) where H is Hermitian
            H = torch.randn(N, N, device=device, dtype=dtype)
            H = (H + H.T.conj()) / 2.0  # Make Hermitian
            G = torch.matrix_exp(1j * H)

            # Transform: U_j' = U_j @ G
            U_j_transformed = U_j @ G
            frames_transformed.append(U_j_transformed)

        # Compute Wilson loop with transformed frames
        result_transformed = wilson_loop_holonomy(
            frames_transformed,
            reorthonormalize=True,
            gauge_check=False,  # Skip checks for speed
        )

        # Compare invariants
        trace_new = result_transformed.trace
        eigenvals_new = sorted(np.abs(result_transformed.eigenvalues))

        trace_error = abs(trace_new - trace_orig)
        eigenval_error = np.max(np.abs(np.array(eigenvals_new) - np.array(eigenvals_orig)))

        print(f"  Test {test_idx+1}: Δtr = {trace_error:.2e}, Δλ_max = {eigenval_error:.2e}")

        if trace_error > 1e-6 or eigenval_error > 1e-6:
            print(f"    FAIL: Invariants changed under gauge transformation")
            all_pass = False

    if all_pass:
        print(f"  ✓ Gauge invariance verified ({n_random_tests} tests passed)")
    else:
        print(f"  ✗ Gauge invariance FAILED")

    return all_pass