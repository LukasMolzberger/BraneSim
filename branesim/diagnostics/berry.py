"""
Dimension-agnostic Berry phase diagnostics.

Computes Berry connection, phase profiles, and gauge-invariant curvature
for 1D, 2D, and 3D systems.

Key concepts:
- Berry connection A_μ: gauge-dependent, computed from phase increments
- Berry phase γ: gauge-dependent cumulative phase (integrated connection)
- Berry curvature F: gauge-invariant, computed from plaquette holonomies

Gauge dependence:
- Connection and phase depend on the choice of phase origin
- Only curvature and closed-loop holonomies are gauge-invariant
- Use amplitude and overlap thresholds to mask unreliable regions
"""

from __future__ import annotations
import math
from dataclasses import dataclass
import numpy as np
import torch

from .tensor_ops import shift_along_axis


@dataclass(frozen=True)
class BerryConfig:
    """
    Configuration for Berry phase diagnostics.

    Attributes
    ----------
    spacing_sim : float
        Grid spacing in simulation units
    amplitude_threshold : float
        Minimum amplitude to include a point in Berry phase calculation.
        Points with |ψ| < threshold are masked out.
        This is a *definition* of the valid region, not a numerical hack.
    overlap_threshold : float
        Minimum overlap magnitude |⟨u_i|u_{i+1}⟩| to trust the phase increment.
        Edges with |overlap| < threshold indicate nodes, mode crossings, or
        regions where the phase becomes ill-defined.
    eps : float
        Small constant for numerical stability
    unwrap : bool
        If True, apply phase unwrapping to cumulative phase (removes 2π jumps)
    force_cpu_on_mps : bool
        If True and device is MPS, move tensors to CPU for complex operations.
        MPS complex support can be unreliable.
    """
    spacing_sim: float
    amplitude_threshold: float = 1e-8
    overlap_threshold: float = 1e-3
    eps: float = 1e-12
    unwrap: bool = True
    force_cpu_on_mps: bool = True


# Backward compatibility wrapper
def BerryPhase1DConfig(
    spacing: float = None,
    spacing_sim: float = None,
    amplitude_threshold: float = 1e-8,
    overlap_threshold: float = 1e-3,
    eps: float = 1e-12,
    unwrap: bool = True,
    force_cpu_on_mps: bool = True,
) -> BerryConfig:
    """
    Backward compatibility wrapper for BerryPhase1DConfig.

    Accepts both 'spacing' (old API) and 'spacing_sim' (new API).
    """
    # Handle backward compatibility: 'spacing' -> 'spacing_sim'
    if spacing is not None and spacing_sim is None:
        spacing_sim = spacing
    elif spacing_sim is None and spacing is None:
        raise ValueError("Either 'spacing' or 'spacing_sim' must be provided")

    return BerryConfig(
        spacing_sim=spacing_sim,
        amplitude_threshold=amplitude_threshold,
        overlap_threshold=overlap_threshold,
        eps=eps,
        unwrap=unwrap,
        force_cpu_on_mps=force_cpu_on_mps,
    )


def berry_connection_along_axis(
    u_hat: torch.Tensor,
    amp: torch.Tensor,
    axis: int,
    cfg: BerryConfig
) -> dict[str, torch.Tensor]:
    """
    Compute Berry connection along a specified axis.

    The discrete Berry connection is defined by the phase increment between
    neighboring points:
        A_μ[i] ≈ (1/h) · arg⟨u_i|u_{i+1}⟩

    This is gauge-dependent (changes under u → exp(iχ)u).

    Parameters
    ----------
    u_hat : torch.Tensor
        Normalized complex state
        Shape: [*grid_shape] (scalar) or [*grid_shape, C] (vector)
    amp : torch.Tensor
        Real amplitude before normalization
        Shape: [*grid_shape]
    axis : int
        Axis along which to compute connection (0 to D-1)
    cfg : BerryConfig
        Configuration parameters

    Returns
    -------
    dict
        Dictionary containing:
        - dphi: phase increment between neighbors, shape reduced by 1 along axis
        - A_axis: Berry connection (dphi / spacing), same shape as dphi
        - overlap_abs: magnitude |⟨u_i|u_{i+1}⟩|, same shape as dphi
        - mask_point: boolean mask for valid points, shape [*grid_shape]
        - valid_edge: boolean mask for valid edges, shape reduced by 1 along axis

    Examples
    --------
    >>> # 1D case
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> cfg = BerryConfig(spacing_sim=1.0)
    >>> result = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    >>> A_x = result['A_axis']  # Berry connection along x
    >>>
    >>> # 2D case: compute connection along x-axis
    >>> psi = complex_band_state_from_quadrature(xi_2d, v_2d, omega)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> result = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    >>> A_x = result['A_axis']  # Shape: [nx-1, ny]
    """
    # Handle MPS device limitations with complex numbers
    if cfg.force_cpu_on_mps and u_hat.device.type == "mps":
        u_hat = u_hat.to("cpu")
        amp = amp.to("cpu")

    # Mask low-amplitude regions
    mask_point = amp > cfg.amplitude_threshold

    # Get neighboring slices along the specified axis
    u_left, u_right = shift_along_axis(u_hat, axis=axis, step=1)
    mask_left, mask_right = shift_along_axis(mask_point, axis=axis, step=1)

    # Edge is valid if both endpoints are valid
    valid_edge = mask_left & mask_right

    # Compute overlap ⟨u_i|u_{i+1}⟩
    if u_hat.ndim == len(amp.shape):
        # Scalar case: simple conjugate multiply
        overlap = torch.conj(u_left) * u_right
    else:
        # Vector case: last dimension is components
        # u_hat.shape = [*grid_shape, C]
        # u_left.shape = [*reduced_shape, C]
        overlap = torch.sum(torch.conj(u_left) * u_right, dim=-1)

    # Compute overlap magnitude
    overlap_abs = torch.abs(overlap)

    # Update valid_edge: require both amplitude and overlap thresholds
    valid_edge = valid_edge & (overlap_abs > cfg.overlap_threshold)

    # Set invalid edges to unity overlap (phase increment = 0)
    overlap = torch.where(valid_edge, overlap, torch.ones_like(overlap))

    # Extract phase increment Δφ = arg⟨u_i|u_{i+1}⟩
    dphi = torch.angle(overlap)

    # Berry connection: A_axis = Δφ / h
    A_axis = dphi / float(cfg.spacing_sim)

    return {
        "dphi": dphi,
        "A_axis": A_axis,
        "overlap_abs": overlap_abs,
        "mask_point": mask_point,
        "valid_edge": valid_edge,
    }


def berry_phase_integrated_along_axis(
    dphi: torch.Tensor,
    axis: int,
    cfg: BerryConfig
) -> dict[str, torch.Tensor]:
    """
    Integrate phase increments to get cumulative Berry phase profiles.

    For each transverse coordinate, integrates dphi along the specified axis
    starting from index 0. In 1D, this gives a single profile γ(x).
    In 2D/3D, this gives a family of profiles (one for each transverse point).

    The Berry phase is gauge-dependent (it changes under local phase rotations).

    Parameters
    ----------
    dphi : torch.Tensor
        Phase increments from berry_connection_along_axis
        Shape: reduced by 1 along the specified axis
    axis : int
        Axis along which to integrate (0 to D-1)
    cfg : BerryConfig
        Configuration (for unwrap setting)

    Returns
    -------
    dict
        Dictionary containing:
        - gamma_unwrapped: cumulative phase (may grow beyond [-π, π])
        - gamma_wrapped: wrapped to [-π, π]

    Notes
    -----
    The output shape matches the original grid shape along the integration axis
    (we prepend a zero for the starting point).

    Examples
    --------
    >>> # 1D case
    >>> result_conn = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    >>> dphi = result_conn['dphi']  # [N-1]
    >>> result_phase = berry_phase_integrated_along_axis(dphi, axis=0, cfg=cfg)
    >>> gamma = result_phase['gamma_unwrapped']  # [N]
    >>>
    >>> # 2D case: integrate along x-axis
    >>> result_conn = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    >>> dphi = result_conn['dphi']  # [nx-1, ny]
    >>> result_phase = berry_phase_integrated_along_axis(dphi, axis=0, cfg=cfg)
    >>> gamma = result_phase['gamma_unwrapped']  # [nx, ny]
    """
    # Determine the shape for the output (add 1 along axis for the zero initial condition)
    output_shape = list(dphi.shape)
    output_shape[axis] += 1

    # Initialize gamma with zeros
    gamma = torch.zeros(output_shape, device=dphi.device, dtype=dphi.dtype)

    # Compute cumulative sum along the specified axis
    if cfg.unwrap:
        # True phase unwrapping
        # For multi-dimensional case, we need to unwrap each 1D slice independently
        gamma_raw = torch.cumsum(dphi, dim=axis)

        # Apply unwrapping along the specified axis
        # Move axis to front, reshape to [axis_len, ...], unwrap each slice
        gamma_raw_moved = torch.moveaxis(gamma_raw, axis, 0)
        original_shape = gamma_raw_moved.shape
        gamma_raw_flat = gamma_raw_moved.reshape(original_shape[0], -1)

        # Unwrap each transverse slice
        gamma_unwrapped_flat = np.unwrap(gamma_raw_flat.cpu().numpy(), axis=0)
        gamma_unwrapped = torch.from_numpy(gamma_unwrapped_flat).to(
            device=dphi.device, dtype=dphi.dtype
        ).reshape(original_shape)

        # Move axis back
        gamma_unwrapped = torch.moveaxis(gamma_unwrapped, 0, axis)

        # Insert the initial zero
        # Build slice to insert at position 1: along axis
        insert_slice = [slice(None)] * len(output_shape)
        insert_slice[axis] = slice(1, None)
        gamma[tuple(insert_slice)] = gamma_unwrapped
    else:
        # Simple cumsum without unwrapping
        gamma_cumsum = torch.cumsum(dphi, dim=axis)
        insert_slice = [slice(None)] * len(output_shape)
        insert_slice[axis] = slice(1, None)
        gamma[tuple(insert_slice)] = gamma_cumsum

    # Wrapped version: fold back into [-π, π]
    two_pi = 2.0 * math.pi
    gamma_wrapped = (gamma + math.pi) % two_pi - math.pi

    return {
        "gamma_unwrapped": gamma,
        "gamma_wrapped": gamma_wrapped,
    }


def berry_plaquette_curvature(
    u_hat: torch.Tensor,
    amp: torch.Tensor,
    axes: tuple[int, int],
    cfg: BerryConfig
) -> dict[str, torch.Tensor]:
    """
    Compute gauge-invariant Berry curvature via plaquette holonomy.

    For a 2D slice defined by axes (a, b), the curvature is computed from
    the U(1) phase around each elementary plaquette:

        U = ⟨u_{00}|u_{10}⟩ ⟨u_{10}|u_{11}⟩ ⟨u_{11}|u_{01}⟩ ⟨u_{01}|u_{00}⟩

    The Berry curvature is:
        F_{ab} = arg(U)

    This quantity is gauge-invariant (invariant under u → exp(iχ)u).

    Parameters
    ----------
    u_hat : torch.Tensor
        Normalized complex state
        Shape: [*grid_shape] (scalar) or [*grid_shape, C] (vector)
    amp : torch.Tensor
        Real amplitude before normalization
        Shape: [*grid_shape]
    axes : tuple[int, int]
        Pair of axes defining the 2D slice (e.g., (0, 1) for xy-plane)
    cfg : BerryConfig
        Configuration parameters

    Returns
    -------
    dict
        Dictionary containing:
        - curvature: Berry curvature F_{ab}, shape reduced by 1 along both axes
        - valid_plaquette: boolean mask for valid plaquettes
        - edge_overlap_abs: dict with overlap magnitudes for each edge direction

    Examples
    --------
    >>> # 2D case: compute curvature in xy-plane
    >>> psi = complex_band_state_from_quadrature(xi_2d, v_2d, omega)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> cfg = BerryConfig(spacing_sim=1.0)
    >>> result = berry_plaquette_curvature(psi_hat, amp, axes=(0, 1), cfg=cfg)
    >>> F_xy = result['curvature']  # Shape: [nx-1, ny-1]
    >>>
    >>> # 3D case: compute curvature in xy-plane at each z
    >>> psi = complex_band_state_from_quadrature(xi_3d, v_3d, omega)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> result = berry_plaquette_curvature(psi_hat, amp, axes=(0, 1), cfg=cfg)
    >>> F_xy = result['curvature']  # Shape: [nx-1, ny-1, nz]
    """
    # Handle MPS device limitations
    if cfg.force_cpu_on_mps and u_hat.device.type == "mps":
        u_hat = u_hat.to("cpu")
        amp = amp.to("cpu")

    axis_a, axis_b = axes

    # Mask low-amplitude regions
    mask_point = amp > cfg.amplitude_threshold

    # Get all four corners of each plaquette
    # Corner notation: (i_a, i_b)
    # 00: (i, j)
    # 10: (i+1, j)
    # 01: (i, j+1)
    # 11: (i+1, j+1)

    # Build slice objects to extract corners
    def get_corner(da: int, db: int) -> torch.Tensor:
        """Get corner of plaquette with offset (da, db)."""
        # Build slice list
        slices = [slice(None)] * u_hat.ndim
        # For spatial dimensions that aren't in axes, keep all
        # For axis_a and axis_b, apply offset and trim
        shape_a = u_hat.shape[axis_a]
        shape_b = u_hat.shape[axis_b]

        if da == 0:
            slices[axis_a] = slice(None, -1)
        else:
            slices[axis_a] = slice(1, None)

        if db == 0:
            slices[axis_b] = slice(None, -1)
        else:
            slices[axis_b] = slice(1, None)

        return u_hat[tuple(slices)]

    def get_corner_mask(da: int, db: int) -> torch.Tensor:
        """Get mask corner."""
        slices = [slice(None)] * mask_point.ndim

        if da == 0:
            slices[axis_a] = slice(None, -1)
        else:
            slices[axis_a] = slice(1, None)

        if db == 0:
            slices[axis_b] = slice(None, -1)
        else:
            slices[axis_b] = slice(1, None)

        return mask_point[tuple(slices)]

    u_00 = get_corner(0, 0)
    u_10 = get_corner(1, 0)
    u_01 = get_corner(0, 1)
    u_11 = get_corner(1, 1)

    mask_00 = get_corner_mask(0, 0)
    mask_10 = get_corner_mask(1, 0)
    mask_01 = get_corner_mask(0, 1)
    mask_11 = get_corner_mask(1, 1)

    # Plaquette is valid if all four corners are valid
    valid_plaquette = mask_00 & mask_10 & mask_01 & mask_11

    # Compute overlaps for each edge
    def compute_overlap(u_left: torch.Tensor, u_right: torch.Tensor) -> torch.Tensor:
        """Compute ⟨u_left|u_right⟩."""
        if u_hat.ndim == len(amp.shape):
            # Scalar case
            return torch.conj(u_left) * u_right
        else:
            # Vector case
            return torch.sum(torch.conj(u_left) * u_right, dim=-1)

    overlap_00_10 = compute_overlap(u_00, u_10)
    overlap_10_11 = compute_overlap(u_10, u_11)
    overlap_11_01 = compute_overlap(u_11, u_01)
    overlap_01_00 = compute_overlap(u_01, u_00)

    # Compute overlap magnitudes
    overlap_abs_00_10 = torch.abs(overlap_00_10)
    overlap_abs_10_11 = torch.abs(overlap_10_11)
    overlap_abs_11_01 = torch.abs(overlap_11_01)
    overlap_abs_01_00 = torch.abs(overlap_01_00)

    # Update valid_plaquette: all edges must have sufficient overlap
    valid_plaquette = valid_plaquette & (
        (overlap_abs_00_10 > cfg.overlap_threshold) &
        (overlap_abs_10_11 > cfg.overlap_threshold) &
        (overlap_abs_11_01 > cfg.overlap_threshold) &
        (overlap_abs_01_00 > cfg.overlap_threshold)
    )

    # Compute plaquette holonomy: U = product of four overlaps around loop
    U = overlap_00_10 * overlap_10_11 * overlap_11_01 * overlap_01_00

    # Set invalid plaquettes to unity (curvature = 0)
    U = torch.where(valid_plaquette, U, torch.ones_like(U))

    # Berry curvature: F = arg(U)
    curvature = torch.angle(U)

    return {
        "curvature": curvature,
        "valid_plaquette": valid_plaquette,
        "edge_overlap_abs": {
            "00_10": overlap_abs_00_10,
            "10_11": overlap_abs_10_11,
            "11_01": overlap_abs_11_01,
            "01_00": overlap_abs_01_00,
        },
    }


def berry_phase_profile_along_x(
    u_hat: torch.Tensor,
    amp: torch.Tensor,
    cfg: BerryConfig
) -> dict[str, torch.Tensor]:
    """
    Compute Berry phase profile along x-axis (1D compatibility wrapper).

    This function provides backward compatibility with the original 1D API.
    It calls berry_connection_along_axis with axis=0 and integrates the result.

    Parameters
    ----------
    u_hat : torch.Tensor
        Normalized complex state, shape [N] or [N, C]
    amp : torch.Tensor
        Real amplitude, shape [N]
    cfg : BerryConfig
        Configuration parameters

    Returns
    -------
    dict
        Dictionary containing:
        - gamma_unwrapped [N]: cumulative Berry phase
        - gamma_wrapped [N]: wrapped to [-π, π]
        - dphi [N-1]: phase increment between neighbors
        - A_x [N-1]: Berry connection along x
        - overlap_abs [N-1]: magnitude |⟨u_i|u_{i+1}⟩|
        - mask [N]: boolean mask for valid points
        - valid_edge [N-1]: boolean mask for valid edges

    Examples
    --------
    >>> from branesim.diagnostics import (
    ...     complex_band_state_from_quadrature,
    ...     pointwise_normalize,
    ...     berry_phase_profile_along_x,
    ...     BerryPhase1DConfig
    ... )
    >>>
    >>> psi = complex_band_state_from_quadrature(xi, xi_dot, omega_sim)
    >>> psi_hat, amp = pointwise_normalize(psi)
    >>> cfg = BerryPhase1DConfig(spacing=h_sim)
    >>> result = berry_phase_profile_along_x(psi_hat, amp, cfg)
    >>> gamma = result['gamma_wrapped']
    """
    # Compute connection along x-axis (axis=0)
    conn_result = berry_connection_along_axis(u_hat, amp, axis=0, cfg=cfg)

    # Integrate to get phase profile
    phase_result = berry_phase_integrated_along_axis(conn_result["dphi"], axis=0, cfg=cfg)

    # Combine results in the expected format
    return {
        "gamma_unwrapped": phase_result["gamma_unwrapped"],
        "gamma_wrapped": phase_result["gamma_wrapped"],
        "dphi": conn_result["dphi"],
        "A_x": conn_result["A_axis"],
        "overlap_abs": conn_result["overlap_abs"],
        "mask": conn_result["mask_point"],
        "valid_edge": conn_result["valid_edge"],
    }