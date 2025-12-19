"""
P1 Polarization basis selection: shear-first, no privileged X⁴.

This module implements the P1 rule for selecting polarization basis vectors:
- In 3D: both are tangential shear vectors perpendicular to k_hat
- In 2D: one tangential, second from remaining embedding axes
- In 1D: two embedding axes chosen deterministically

Key principle: NO special treatment of X⁴. The selection is purely geometric
and deterministic, based on the propagation direction and dimensionality.
"""

import torch
from typing import Tuple


def photon_polarization_basis(
    intrinsic_dim: int,
    k_hat: torch.Tensor,      # [d] normalized propagation direction
    prefer_shear: bool = True,
    device=None,
    dtype=None,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Select polarization basis vectors p1, p2 for a photon packet.

    Returns two 4D unit vectors representing the 2D degenerate polarization
    subspace. The selection follows the P1 rule: prefer in-brane shear modes
    over amplitude excitation, with deterministic fallback ordering.

    P1 Rule implementation:
    - 3D: Both p1, p2 are tangential (perpendicular to k_hat) and embedded as (e, 0)
    - 2D: p1 is tangential, p2 chosen from remaining embedding dimensions
    - 1D: Both chosen from embedding dimensions (e.g., axes 1 and 2)

    Args:
        intrinsic_dim: Intrinsic dimensionality (1, 2, or 3)
        k_hat: Normalized propagation direction [d]
        prefer_shear: If True, apply P1 rule (default)
        device: torch.device for output tensors
        dtype: torch.dtype for output tensors

    Returns:
        (p1, p2): Tuple of two [4] tensors, orthonormal polarization basis
    """
    if device is None:
        device = k_hat.device
    if dtype is None:
        dtype = k_hat.dtype

    k_hat = k_hat.to(device=device, dtype=dtype)

    if intrinsic_dim == 3:
        return _polarization_3d(k_hat, prefer_shear, device, dtype)
    elif intrinsic_dim == 2:
        return _polarization_2d(k_hat, prefer_shear, device, dtype)
    elif intrinsic_dim == 1:
        return _polarization_1d(k_hat, prefer_shear, device, dtype)
    else:
        raise ValueError(f"intrinsic_dim must be 1, 2, or 3, got {intrinsic_dim}")


def _polarization_3d(
    k_hat: torch.Tensor,
    prefer_shear: bool,
    device,
    dtype,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    3D P1 polarization: both vectors are tangential shear modes.

    We find two orthonormal vectors perpendicular to k_hat in R³,
    then embed them as (e, 0) in R⁴.
    """
    # Find reference axis least aligned with k_hat
    k_abs = torch.abs(k_hat)
    min_idx = torch.argmin(k_abs)

    # Build reference vector along that axis
    ref = torch.zeros(3, device=device, dtype=dtype)
    ref[min_idx] = 1.0

    # First transverse: e1 = normalize(k_hat × ref)
    e1 = torch.cross(k_hat, ref)
    e1 = e1 / torch.linalg.norm(e1)

    # Second transverse: e2 = k_hat × e1 (automatically normalized)
    e2 = torch.cross(k_hat, e1)
    e2 = e2 / torch.linalg.norm(e2)  # normalize for numerical safety

    # Embed into 4D: (e, 0)
    p1 = torch.zeros(4, device=device, dtype=dtype)
    p1[:3] = e1

    p2 = torch.zeros(4, device=device, dtype=dtype)
    p2[:3] = e2

    return p1, p2


def _polarization_2d(
    k_hat: torch.Tensor,
    prefer_shear: bool,
    device,
    dtype,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    2D P1 polarization: one tangential, one from remaining embedding axes.

    k_hat is in R² (x, y). We choose:
    - p1: tangential in-plane, embedded as (e, 0, 0)
    - p2: deterministically from remaining embedding axes (axis 2, then axis 3)
    """
    # Tangential in 2D: rotate k_hat by 90°
    # k = (kx, ky) → e1 = (-ky, kx) / norm
    e1 = torch.tensor([-k_hat[1], k_hat[0]], device=device, dtype=dtype)
    e1 = e1 / torch.linalg.norm(e1)

    # Embed into 4D
    p1 = torch.zeros(4, device=device, dtype=dtype)
    p1[:2] = e1

    # Second basis vector: choose from remaining embedding axes
    # Deterministic choice: axis 2 (index 2)
    p2 = torch.zeros(4, device=device, dtype=dtype)
    p2[2] = 1.0

    # Gram-Schmidt: ensure p2 ⊥ p1
    p2 = p2 - torch.dot(p2, p1) * p1
    p2 = p2 / torch.linalg.norm(p2)

    return p1, p2


def _polarization_1d(
    k_hat: torch.Tensor,
    prefer_shear: bool,
    device,
    dtype,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    1D P1 polarization: choose two embedding axes deterministically.

    k_hat is in R¹ (just direction, always along x-axis).
    We choose axes 1 and 2 (indices 1 and 2) as the polarization plane.
    This treats all transverse embedding directions equally.
    """
    p1 = torch.zeros(4, device=device, dtype=dtype)
    p1[1] = 1.0  # Y-axis (X¹)

    p2 = torch.zeros(4, device=device, dtype=dtype)
    p2[2] = 1.0  # Z-axis (X²)

    return p1, p2


def validate_polarization_basis(
    p1: torch.Tensor,
    p2: torch.Tensor,
    k_hat: torch.Tensor,
    intrinsic_dim: int,
    atol: float = 1e-6,
) -> dict:
    """
    Validate that polarization basis satisfies P1 requirements.

    Returns:
        dict with validation results and diagnostics
    """
    results = {}

    # Check normalization
    p1_norm = torch.linalg.norm(p1).item()
    p2_norm = torch.linalg.norm(p2).item()
    results["p1_norm"] = p1_norm
    results["p2_norm"] = p2_norm
    results["normalized"] = abs(p1_norm - 1.0) < atol and abs(p2_norm - 1.0) < atol

    # Check orthogonality
    dot_p1_p2 = torch.dot(p1, p2).item()
    results["dot_p1_p2"] = dot_p1_p2
    results["orthogonal"] = abs(dot_p1_p2) < atol

    # Check transversality (for 3D)
    if intrinsic_dim == 3:
        dot_p1_k = torch.dot(p1[:3], k_hat).item()
        dot_p2_k = torch.dot(p2[:3], k_hat).item()
        results["dot_p1_k"] = dot_p1_k
        results["dot_p2_k"] = dot_p2_k
        results["transverse"] = abs(dot_p1_k) < atol and abs(dot_p2_k) < atol
    else:
        results["transverse"] = None  # not applicable

    results["valid"] = results["normalized"] and results["orthogonal"]
    if intrinsic_dim == 3:
        results["valid"] = results["valid"] and results["transverse"]

    return results