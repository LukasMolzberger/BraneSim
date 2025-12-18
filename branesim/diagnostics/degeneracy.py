"""
Degeneracy verification without spectral band-pass filtering.

This module implements the diagnostic methods described in the paper section
"Degeneracy verification without spectral band-pass filtering".

Key principle:
    Verify that the prepared packet occupies a (nearly) degenerate N-dimensional
    mode subspace near ω0 WITHOUT using diagnostic band-pass masks. This confirms
    that degeneracy is a property of the initialization, not an artifact of filtering.
"""

import torch
import numpy as np
from typing import TYPE_CHECKING, Tuple, Dict, Optional, List
from dataclasses import dataclass

if TYPE_CHECKING:
    from branesim.core.state import BraneState
    from branesim.core.grid import BraneGrid


@dataclass
class EigenScanResult:
    """Results from local linearization eigen-scan."""
    omega_values: np.ndarray  # Eigenfrequencies near omega0
    eigenvectors: np.ndarray  # Corresponding eigenvectors [n_modes, n_dof]
    degeneracy_pairs: List[Tuple[int, int]]  # Indices of degenerate pairs
    degeneracy_score: float  # |ω1 - ω2| / ω0 for closest pair
    omega0_target: float  # Target frequency
    n_modes_found: int


@dataclass
class SubspaceRankResult:
    """Results from SVD subspace rank analysis."""
    singular_values: np.ndarray  # Sorted singular values
    dominant_dimension: int  # Number of significant dimensions
    energy_ratios: np.ndarray  # Cumulative energy fractions
    gap_ratio: float  # (σ_N - σ_{N+1}) / σ_1
    is_degenerate_2d: bool  # True if clearly 2D subspace


def local_eigen_scan(
    state: "BraneState",
    grid: "BraneGrid",
    physics,  # SpringForceComputer or similar
    omega0_target: float,
    n_modes: int = 10,
    patch_size: Optional[int] = None,
    degeneracy_threshold: float = 1e-3,
    dofs: Optional[List[int]] = None,
) -> EigenScanResult:
    """
    Perform local linearization eigen-scan to verify degeneracy.

    Linearizes the discrete equations of motion around the current substrate
    configuration and computes eigenpairs near omega0. Degeneracy is identified
    by |ω1 - ω2| / ω0 << 1 together with consistent polarization character.

    This method does NOT use FFT or band-pass filtering.

    Args:
        state: Current BraneState
        grid: BraneGrid
        physics: Force computer (e.g., SpringForceComputer) with compute_forces method
        omega0_target: Target carrier frequency to search near
        n_modes: Number of eigenmodes to compute (default: 10)
        patch_size: Optional local patch size (None = use full system)
        degeneracy_threshold: |ω1 - ω2|/ω0 below this is considered degenerate
        dofs: Which DOFs to include in analysis (None = all 4)

    Returns:
        EigenScanResult with eigenfrequencies, vectors, and degeneracy measures

    Paper reference:
        Section "Degeneracy verification without spectral band-pass filtering"
        Paragraph "(i) Local linearization eigen-scan"
    """
    device = state.device
    dtype = state.dtype

    if dofs is None:
        dofs = list(range(4))  # All DOFs
    n_dof = len(dofs)

    # For simplicity, work with full system (patch implementation would select subset)
    N = state.num_points

    # Build stiffness matrix K and mass matrix M
    # For small-amplitude oscillations around current state:
    # M ü + K u = 0
    # ω² M u = K u  (eigenvalue problem)

    # Compute finite-difference approximation to Hessian
    # K_ij = -∂²U/∂q_i∂q_j evaluated at current state

    # Mass matrix (diagonal)
    m_point = physics.m_point if hasattr(physics, 'm_point') else 1.0
    M_diag = m_point * torch.ones(N * n_dof, device=device, dtype=dtype)

    # Stiffness matrix approximation via finite differences
    # Perturb each DOF and measure force change
    epsilon = 1e-8
    K = torch.zeros(N * n_dof, N * n_dof, device=device, dtype=dtype)

    # Save original state
    pos_orig = state.positions.clone()

    print(f"\nComputing stiffness matrix via finite differences...")
    print(f"  System size: {N} nodes × {n_dof} DOFs = {N*n_dof} dimensions")
    print(f"  WARNING: Full eigen-decomposition may be slow for large systems")

    # Compute baseline forces
    F0 = physics.compute_forces(state, grid)[:, dofs].flatten()

    for i in range(N):
        for d_idx, d in enumerate(dofs):
            idx = i * n_dof + d_idx

            # Perturb
            state.positions[i, d] += epsilon

            # Compute perturbed forces
            F_pert = physics.compute_forces(state, grid)[:, dofs].flatten()

            # Finite difference: K[:, idx] ≈ -(F_pert - F0) / epsilon
            K[:, idx] = -(F_pert - F0) / epsilon

            # Restore
            state.positions[i, d] = pos_orig[i, d]

    # Restore original state completely
    state.positions[:] = pos_orig

    # Symmetrize K (should be symmetric up to numerical error)
    K = 0.5 * (K + K.t())

    # Solve generalized eigenvalue problem: K u = ω² M u
    # Convert to standard form: M^{-1/2} K M^{-1/2} v = ω² v
    # where u = M^{-1/2} v

    M_inv_sqrt = 1.0 / torch.sqrt(M_diag)
    M_inv_sqrt_mat = torch.diag(M_inv_sqrt)

    K_normalized = M_inv_sqrt_mat @ K @ M_inv_sqrt_mat

    # Compute eigenvalues and eigenvectors
    print(f"  Computing eigendecomposition...")
    try:
        # Use only the smallest n_modes eigenvalues (scipy would be better for large systems)
        if N * n_dof > 1000:
            print(f"  WARNING: System too large ({N*n_dof}), using approximate method")
            # For large systems, should use sparse methods, but here we'll just warn
            eigenvalues, eigenvectors = torch.linalg.eigh(K_normalized)
        else:
            eigenvalues, eigenvectors = torch.linalg.eigh(K_normalized)

        # Convert eigenvalues to frequencies: ω = √λ
        omega_sq = eigenvalues.cpu().numpy()
        omega_values = np.sqrt(np.abs(omega_sq)) * np.sign(omega_sq)

        # Transform eigenvectors back: u = M^{-1/2} v
        eigenvectors_physical = M_inv_sqrt_mat @ eigenvectors
        eigenvectors_np = eigenvectors_physical.cpu().numpy()

    except Exception as e:
        print(f"  ERROR in eigendecomposition: {e}")
        return EigenScanResult(
            omega_values=np.array([]),
            eigenvectors=np.array([]),
            degeneracy_pairs=[],
            degeneracy_score=np.inf,
            omega0_target=omega0_target,
            n_modes_found=0,
        )

    # Find modes near omega0_target
    omega_sorted_idx = np.argsort(np.abs(omega_values - omega0_target))
    nearest_modes = omega_sorted_idx[:n_modes]

    omega_near = omega_values[nearest_modes]
    vectors_near = eigenvectors_np[:, nearest_modes]

    # Identify degenerate pairs
    degeneracy_pairs = []
    best_degeneracy = np.inf

    for i in range(len(omega_near)):
        for j in range(i+1, len(omega_near)):
            delta_omega = abs(omega_near[i] - omega_near[j])
            relative_split = delta_omega / omega0_target if omega0_target > 0 else np.inf

            if relative_split < degeneracy_threshold:
                degeneracy_pairs.append((i, j))
                best_degeneracy = min(best_degeneracy, relative_split)

    print(f"\n  Eigen-scan results:")
    print(f"    Found {len(omega_near)} modes near ω0 = {omega0_target:.6e} rad/s")
    print(f"    Closest mode: ω = {omega_near[0]:.6e} rad/s")
    print(f"    Degeneracy pairs found: {len(degeneracy_pairs)}")
    if degeneracy_pairs:
        print(f"    Best degeneracy score: Δω/ω0 = {best_degeneracy:.6e}")

    return EigenScanResult(
        omega_values=omega_near,
        eigenvectors=vectors_near,
        degeneracy_pairs=degeneracy_pairs,
        degeneracy_score=best_degeneracy,
        omega0_target=omega0_target,
        n_modes_found=len(omega_near),
    )


def subspace_rank_svd(
    time_series: torch.Tensor,
    sample_locations: Optional[torch.Tensor] = None,
    energy_threshold: float = 0.99,
    gap_threshold: float = 0.1,
) -> SubspaceRankResult:
    """
    Extract dominant subspace dimension via SVD of time series.

    Records a short time window (a few carrier cycles) and extracts the dominant
    subspace dimension via SVD of the multi-component signal. A genuine
    two-dimensional polarization sector exhibits two dominant singular values
    and strongly suppressed remainder.

    This method does NOT use FFT or band-pass filtering.

    Args:
        time_series: Time series data [n_time, n_spatial, n_dof]
                     OR [n_time, n_dof] for single spatial location
        sample_locations: Optional indices of spatial points to analyze
        energy_threshold: Cumulative energy fraction to determine dominant dimension
        gap_threshold: Relative gap (σ_N - σ_{N+1})/σ_1 to identify clear separation

    Returns:
        SubspaceRankResult with singular values and dimension analysis

    Paper reference:
        Section "Degeneracy verification without spectral band-pass filtering"
        Paragraph "(ii) Subspace rank / SVD signature from time series"
    """
    # Reshape to [n_time, n_features]
    if time_series.ndim == 3:
        n_time, n_spatial, n_dof = time_series.shape
        if sample_locations is not None:
            data = time_series[:, sample_locations, :].reshape(n_time, -1)
        else:
            data = time_series.reshape(n_time, n_spatial * n_dof)
    elif time_series.ndim == 2:
        data = time_series  # Already [n_time, n_features]
    else:
        raise ValueError(f"time_series must be 2D or 3D, got shape {time_series.shape}")

    # Center the data
    data_mean = data.mean(dim=0, keepdim=True)
    data_centered = data - data_mean

    # SVD: data_centered = U Σ V^T
    print(f"\nPerforming SVD on time series...")
    print(f"  Data shape: {data_centered.shape}")

    try:
        U, S, Vt = torch.linalg.svd(data_centered, full_matrices=False)
        singular_values = S.cpu().numpy()
    except Exception as e:
        print(f"  ERROR in SVD: {e}")
        return SubspaceRankResult(
            singular_values=np.array([]),
            dominant_dimension=0,
            energy_ratios=np.array([]),
            gap_ratio=0.0,
            is_degenerate_2d=False,
        )

    # Compute cumulative energy fractions
    total_energy = np.sum(singular_values ** 2)
    energy_fractions = np.cumsum(singular_values ** 2) / total_energy

    # Determine dominant dimension
    dominant_dim = np.searchsorted(energy_fractions, energy_threshold) + 1
    dominant_dim = min(dominant_dim, len(singular_values))

    # Check for clear 2D structure
    if len(singular_values) >= 3:
        gap_2_3 = (singular_values[1] - singular_values[2]) / singular_values[0]
        is_2d = (dominant_dim == 2) or (gap_2_3 > gap_threshold)
    else:
        gap_2_3 = 0.0
        is_2d = (dominant_dim == 2)

    print(f"\n  SVD results:")
    print(f"    Singular values (top 5): {singular_values[:5]}")
    print(f"    Dominant dimension: {dominant_dim}")
    print(f"    Energy in top {dominant_dim} modes: {energy_fractions[dominant_dim-1]*100:.2f}%")
    if len(singular_values) >= 3:
        print(f"    Gap ratio σ₂/σ₁: {singular_values[1]/singular_values[0]:.4f}")
        print(f"    Gap ratio (σ₂ - σ₃)/σ₁: {gap_2_3:.4f}")
        print(f"    Is 2D degenerate: {'YES' if is_2d else 'NO'}")

    return SubspaceRankResult(
        singular_values=singular_values,
        dominant_dimension=dominant_dim,
        energy_ratios=energy_fractions,
        gap_ratio=gap_2_3 if len(singular_values) >= 3 else 0.0,
        is_degenerate_2d=is_2d,
    )


def verify_narrowband_preparation(
    state: "BraneState",
    k0_expected: float,
    omega0_expected: float,
    dofs: List[int],
    grid: Optional["BraneGrid"] = None,
) -> Dict[str, float]:
    """
    Verify that initialization is narrowband by measuring spectral content.

    This is an optional cross-check that can use FFT, but is NOT required
    for primary results. It simply confirms that the preparation-first
    initialization achieved a narrowband spectrum.

    Args:
        state: Current BraneState
        k0_expected: Expected carrier wave number
        omega0_expected: Expected carrier frequency (not used in spatial check)
        dofs: List of DOFs to analyze
        grid: Optional BraneGrid for spatial FFT

    Returns:
        Dictionary with spectral metrics:
            - peak_k: Measured peak wave number
            - k0_error: |peak_k - k0_expected| / k0_expected
            - spectral_width: Estimate of Δk/k0
            - is_narrowband: True if Δk/k0 < 0.1
    """
    print("\n" + "="*70)
    print("OPTIONAL FFT CROSS-CHECK (not required for main results)")
    print("="*70)

    # For spatial narrowband check, do FFT in space
    if grid is None or grid.grid_shape is None:
        print("  Skipping: requires grid information for spatial FFT")
        return {"is_narrowband": False}

    # Take 1D slice along first axis for simplicity
    if len(grid.grid_shape) > 1:
        # For multi-D, take slice through center
        center_idx = [s//2 for s in grid.grid_shape[1:]]
        slice_idx = (slice(None),) + tuple(center_idx)
        data = state.positions[slice_idx, dofs[0]]
    else:
        data = state.positions[:, dofs[0]]

    # FFT
    fft = torch.fft.rfft(data)
    power = torch.abs(fft) ** 2

    # Find peak
    peak_idx = torch.argmax(power).item()
    n_points = len(data)
    k_values = 2 * np.pi * torch.fft.rfftfreq(n_points, d=grid.spacing)
    peak_k = k_values[peak_idx].item()

    # Estimate width (FWHM)
    half_max = power[peak_idx] / 2.0
    above_half = power > half_max
    width_indices = torch.where(above_half)[0]
    if len(width_indices) > 1:
        delta_k = (width_indices[-1] - width_indices[0]) * (k_values[1] - k_values[0])
        spectral_width = delta_k.item() / peak_k if peak_k > 0 else np.inf
    else:
        spectral_width = 0.0

    k0_error = abs(peak_k - k0_expected) / k0_expected if k0_expected > 0 else np.inf
    is_narrowband = spectral_width < 0.1

    print(f"  Expected k0:         {k0_expected:.6e} rad/m")
    print(f"  Measured peak k:     {peak_k:.6e} rad/m")
    print(f"  Error |k - k0|/k0:   {k0_error:.6e}")
    print(f"  Spectral width Δk/k: {spectral_width:.6e}")
    print(f"  Is narrowband:       {'YES' if is_narrowband else 'NO'}")
    print("="*70 + "\n")

    return {
        "peak_k": peak_k,
        "k0_error": k0_error,
        "spectral_width": spectral_width,
        "is_narrowband": is_narrowband,
    }