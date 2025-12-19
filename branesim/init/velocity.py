"""
Velocity initialization methods for wave packets.

This module implements three methods for computing initial velocities:
1. time_reversal_shift: Use brane forces to time-reverse a spatial shift
2. directional_derivative: Simple v = -c * (k̂ · ∇u) formula
3. complex_quadrature: v = Re(-iω * ψ) for narrowband carriers
"""

import torch
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from branesim.core.state import BraneState
    from branesim.core.grid import BraneGrid
    from branesim.physics.forces import SpringForceComputer
    from branesim.init.artifacts import RestGeometryArtifact


def velocities_time_reversal_shift(
    geom: "RestGeometryArtifact",
    u0: torch.Tensor,              # [N, 4] initial displacements
    physics: "SpringForceComputer",
    grid: "BraneGrid",
    m_point: float,
    wave_speed: float,
    shift_cells: int,
    k_hat: torch.Tensor,           # [d] propagation direction
    periodic: bool = False,
) -> torch.Tensor:
    """
    Compute velocities by time-reversing a spatial shift along k_hat.

    This uses the ACTUAL brane forces to determine consistent initial
    velocities. Given a displacement field u0, we:
    1. Create a target field u_target by shifting u0 along k_hat direction
    2. Compute forces at the initial state
    3. Solve: u_target ≈ u0 + v0*dt + 0.5*a0*dt² for v0

    This method works for ALL 4 components simultaneously and shifts
    along the dominant direction of k_hat.

    Args:
        geom: Rest geometry artifact
        u0: [N, 4] displacement field
        physics: Force computer
        grid: Grid topology
        m_point: Point mass
        wave_speed: Target wave speed (for computing dt)
        shift_cells: Number of cells to shift
        k_hat: [d] normalized propagation direction
        periodic: Whether to use periodic wrapping

    Returns:
        v0: [N, 4] initial velocities
    """
    from branesim.core.state import BraneState

    device = u0.device
    dtype = u0.dtype
    N = u0.shape[0]
    grid_shape = geom.grid_shape
    ndim = geom.intrinsic_dim

    # Compute time step for shift
    dt = (shift_cells * geom.spacing) / wave_speed

    # Determine dominant axis (largest |k_hat| component)
    k_abs = torch.abs(k_hat)
    dom_axis = torch.argmax(k_abs).item()

    # Create target displacement by shifting along dominant axis
    u_target = torch.zeros_like(u0)

    # Reshape to grid for shifting
    u0_grid = u0.view(*grid_shape, 4)  # [*grid_shape, 4]
    u_target_grid = torch.zeros_like(u0_grid)

    # Shift along dominant axis for all components
    if ndim == 1:
        nx = grid_shape[0]
        if shift_cells < nx:
            if periodic:
                u_target_grid[:, :] = torch.roll(u0_grid, -shift_cells, dims=0)
            else:
                u_target_grid[shift_cells:, :] = u0_grid[:-shift_cells, :]

    elif ndim == 2:
        nx, ny = grid_shape
        if dom_axis == 0 and shift_cells < nx:
            if periodic:
                u_target_grid[:, :, :] = torch.roll(u0_grid, -shift_cells, dims=0)
            else:
                u_target_grid[shift_cells:, :, :] = u0_grid[:-shift_cells, :, :]
        elif dom_axis == 1 and shift_cells < ny:
            if periodic:
                u_target_grid[:, :, :] = torch.roll(u0_grid, -shift_cells, dims=1)
            else:
                u_target_grid[:, shift_cells:, :] = u0_grid[:, :-shift_cells, :]

    elif ndim == 3:
        nx, ny, nz = grid_shape
        if dom_axis == 0 and shift_cells < nx:
            if periodic:
                u_target_grid[:, :, :, :] = torch.roll(u0_grid, -shift_cells, dims=0)
            else:
                u_target_grid[shift_cells:, :, :, :] = u0_grid[:-shift_cells, :, :, :]
        elif dom_axis == 1 and shift_cells < ny:
            if periodic:
                u_target_grid[:, :, :, :] = torch.roll(u0_grid, -shift_cells, dims=1)
            else:
                u_target_grid[:, shift_cells:, :, :] = u0_grid[:, :-shift_cells, :, :]
        elif dom_axis == 2 and shift_cells < nz:
            if periodic:
                u_target_grid[:, :, :, :] = torch.roll(u0_grid, -shift_cells, dims=2)
            else:
                u_target_grid[:, :, shift_cells:, :] = u0_grid[:, :, :-shift_cells, :]

    u_target = u_target_grid.reshape(N, 4)

    # Create temporary state at initial displaced configuration
    from branesim.core.state import Dimensionality
    if ndim == 1:
        dim = Dimensionality.ONE_D
    elif ndim == 2:
        dim = Dimensionality.TWO_D
    else:
        dim = Dimensionality.THREE_D

    temp_state = BraneState(grid_shape, dim, device, dtype)
    temp_state.rest_positions = geom.rest_positions.clone()
    temp_state.positions = geom.rest_positions + u0
    temp_state.velocities = torch.zeros_like(u0)

    # Compute forces at initial state
    forces = physics.compute_forces(temp_state, grid)  # [N, 4]
    a0 = forces / m_point  # [N, 4]

    # Solve for v0: u_target ≈ u0 + v0*dt + 0.5*a0*dt²
    # → v0 = (u_target - u0 - 0.5*a0*dt²) / dt
    half_dt2 = 0.5 * dt * dt
    v0 = (u_target - u0 - half_dt2 * a0) / dt

    return v0


def velocities_directional_derivative(
    geom: "RestGeometryArtifact",
    u0: torch.Tensor,           # [N, 4] displacements
    wave_speed: float,
    k_hat: torch.Tensor,        # [d] propagation direction
) -> torch.Tensor:
    """
    Compute velocities using directional derivative: v = -c * (k̂ · ∇u).

    This is a simple kinematic formula assuming wave propagation
    along k_hat with speed c. Works for any dimensionality.

    Args:
        geom: Rest geometry artifact
        u0: [N, 4] displacement field
        wave_speed: Wave speed c
        k_hat: [d] normalized propagation direction

    Returns:
        v0: [N, 4] initial velocities
    """
    device = u0.device
    dtype = u0.dtype
    grid_shape = geom.grid_shape
    ndim = geom.intrinsic_dim
    h = geom.spacing
    N = u0.shape[0]

    v0 = torch.zeros_like(u0)

    # Reshape to grid
    u_grid = u0.view(*grid_shape, 4)

    # Compute gradient for each component
    for comp in range(4):
        u_comp = u_grid[..., comp]
        grad_components = []

        # Compute ∂u/∂x_i for each spatial dimension
        for axis in range(ndim):
            grad_axis = torch.zeros_like(u_comp)

            # Central differences for interior
            if grid_shape[axis] > 2:
                # Build slice tuples
                center_slice = [slice(None)] * ndim
                center_slice[axis] = slice(1, -1)

                plus_slice = [slice(None)] * ndim
                plus_slice[axis] = slice(2, None)

                minus_slice = [slice(None)] * ndim
                minus_slice[axis] = slice(None, -2)

                grad_axis[tuple(center_slice)] = (
                    u_comp[tuple(plus_slice)] - u_comp[tuple(minus_slice)]
                ) / (2.0 * h)

            # One-sided at boundaries
            left_slice = [slice(None)] * ndim
            left_slice[axis] = 0
            left_plus = [slice(None)] * ndim
            left_plus[axis] = 1
            grad_axis[tuple(left_slice)] = (
                u_comp[tuple(left_plus)] - u_comp[tuple(left_slice)]
            ) / h

            right_slice = [slice(None)] * ndim
            right_slice[axis] = -1
            right_minus = [slice(None)] * ndim
            right_minus[axis] = -2
            grad_axis[tuple(right_slice)] = (
                u_comp[tuple(right_slice)] - u_comp[tuple(right_minus)]
            ) / h

            grad_components.append(grad_axis)

        # Stack gradients: [ndim, *grid_shape]
        grad = torch.stack(grad_components, dim=0)

        # Compute directional derivative: k̂ · ∇u
        k_view = k_hat.view(ndim, *([1]*ndim))
        dir_deriv = (k_view * grad).sum(dim=0)  # [*grid_shape]

        # v = -c * (k̂ · ∇u)
        v0[:, comp] = (-wave_speed * dir_deriv).reshape(N)

    return v0


def velocities_from_complex_quadrature(
    psi: torch.Tensor,    # complex [N, 4]
    omega: float,         # carrier frequency
) -> torch.Tensor:
    """
    Compute velocities from complex carrier using quadrature formula.

    For a narrowband carrier ψ(x) = A(x) exp(iφ(x)) p with frequency ω:
        v(x,0) = Re(-iω * ψ(x))

    This assumes the carrier was prepared with known ω. For photons, ω
    should be derived from k and the dispersion relation. For electrons,
    ω comes from the longitudinal phase advance.

    Args:
        psi: complex [N, 4] carrier field
        omega: Carrier angular frequency

    Returns:
        v0: [N, 4] initial velocities
    """
    # v = Re(-i * ω * ψ)
    #   = Re(ω * (-i) * ψ)
    #   = Re(ω * ψ * e^(-iπ/2))
    #   = ω * Im(ψ)  (since -i * z = Im(z) - i*Re(z), and we take real part)
    #
    # Wait, let me recalculate:
    # -i * ψ = -i * (a + ib) = -i*a - i²*b = -i*a + b = b - i*a
    # Re(-i * ψ) = Re(b - i*a) = b = Im(ψ)

    v0 = omega * psi.imag

    return v0


def estimate_omega_from_spectrum(
    u0: torch.Tensor,         # [N, 4] displacement
    coords: torch.Tensor,     # [N, d] coordinates
    k_hat: torch.Tensor,      # [d] propagation direction
    wave_speed: float,
) -> float:
    """
    Estimate carrier frequency ω from the displacement spectrum.

    This is an OPTIONAL helper for cases where ω is not known a priori.
    It computes the FFT along the propagation direction and finds the
    dominant frequency, then uses ω ≈ c * k.

    Args:
        u0: Displacement field
        coords: Spatial coordinates
        k_hat: Propagation direction
        wave_speed: Wave speed c

    Returns:
        omega: Estimated carrier frequency
    """
    # Project coordinates along propagation direction
    x_proj = torch.matmul(coords, k_hat)  # [N]

    # Use first non-zero component for spectrum
    comp_idx = 0
    for i in range(4):
        if torch.abs(u0[:, i]).max() > 1e-10:
            comp_idx = i
            break

    # Sort by projected coordinate
    sort_idx = torch.argsort(x_proj)
    x_sorted = x_proj[sort_idx]
    u_sorted = u0[sort_idx, comp_idx]

    # FFT (only if evenly spaced, approximate)
    fft = torch.fft.rfft(u_sorted)
    fft_mag = torch.abs(fft)

    # Find dominant frequency
    k_idx = torch.argmax(fft_mag[1:]) + 1  # skip DC
    N = len(u_sorted)
    L = x_sorted[-1] - x_sorted[0]
    k_mag = 2 * np.pi * k_idx / L

    omega = wave_speed * k_mag

    return omega.item() if torch.is_tensor(omega) else omega