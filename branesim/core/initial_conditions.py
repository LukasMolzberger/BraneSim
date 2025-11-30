"""
Initial conditions for brane simulations.

This module provides dimension-independent utilities for setting up initial
conditions, particularly for wave packets that propagate with speed c.

KEY PRINCIPLE:
    For a right-moving wave ξ(x,t) = f(n·x - ct), the initial velocity is:
    v(x) = ∂_t ξ(x,0) = -c (n · ∇ξ(x,0))

This works in 1D, 2D, and 3D - only the gradient computation changes.
"""

import torch
import numpy as np
from typing import TYPE_CHECKING
from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid

if TYPE_CHECKING:
    from branesim.physics.forces import SpringForceComputer


def initialize_right_moving_velocities(
    state: BraneState,
    grid: BraneGrid,
    wave_speed: float,
    direction: torch.Tensor | None = None,
    field_component: int = 3,
) -> None:
    """
    Initialize velocities for a right-moving wave with speed c.

    Given an already initialized field ξ in positions[:, field_component],
    sets velocities so that the pattern travels with speed `wave_speed`
    along `direction`.

    Physics:
        For a wave traveling in direction n with speed c:
            ξ(x,t) = f(n·x - ct)
            ∂_t ξ = -c (n · ∇ξ)

    Implementation uses central finite differences for the gradient,
    with one-sided differences at boundaries.

    Args:
        state: BraneState with positions already initialized
        grid: BraneGrid with grid_shape and spacing
        wave_speed: Propagation speed (should equal √(T/ρ_m) = c)
        direction: Optional direction vector in R^D (default: +x axis)
        field_component: Which embedding component stores the wave (default: 3 for X⁴)

    Note:
        - Works for 1D, 2D, and 3D automatically
        - Respects fixed boundary conditions (sets v=0 for fixed points)
        - Direction is automatically normalized
    """
    device = state.device
    dtype = state.dtype
    grid_shape = grid.grid_shape
    ndim = len(grid_shape)
    N = state.num_points

    # Validate grid consistency
    expected_N = int(np.prod(grid_shape))
    if N != expected_N:
        raise ValueError(
            f"Grid shape {grid_shape} implies {expected_N} points, "
            f"but state has {N} points"
        )

    # Extract scalar field and reshape to grid
    field = state.positions[:, field_component].view(*grid_shape)

    h = grid.spacing

    # Default direction: +x axis
    if direction is None:
        direction = torch.zeros(ndim, device=device, dtype=dtype)
        direction[0] = 1.0  # +x direction
    else:
        direction = direction.to(device=device, dtype=dtype)

    # Normalize direction
    norm = torch.linalg.norm(direction)
    if norm < 1e-10:
        raise ValueError("Direction vector must have non-zero norm")
    direction = direction / norm

    # Compute spatial gradient ∇ξ using finite differences
    grads = []
    for ax in range(ndim):
        # Initialize gradient component
        d = torch.zeros_like(field)

        # Central differences for interior points
        if grid_shape[ax] > 2:
            # Build slice tuples for indexing
            center = [slice(None)] * ndim
            center[ax] = slice(1, -1)

            plus = [slice(None)] * ndim
            plus[ax] = slice(2, None)

            minus = [slice(None)] * ndim
            minus[ax] = slice(None, -2)

            # Central difference: (f(x+h) - f(x-h)) / (2h)
            d[tuple(center)] = (field[tuple(plus)] - field[tuple(minus)]) / (2.0 * h)

        # One-sided difference at left boundary
        left0 = [slice(None)] * ndim
        left0[ax] = 0
        left1 = [slice(None)] * ndim
        left1[ax] = 1
        d[tuple(left0)] = (field[tuple(left1)] - field[tuple(left0)]) / h

        # One-sided difference at right boundary
        rightm1 = [slice(None)] * ndim
        rightm1[ax] = -1
        rightm2 = [slice(None)] * ndim
        rightm2[ax] = -2
        d[tuple(rightm1)] = (field[tuple(rightm1)] - field[tuple(rightm2)]) / h

        grads.append(d)

    # Stack gradient components: [D, *grid_shape]
    grad = torch.stack(grads, dim=0)

    # Compute directional derivative: n · ∇ξ
    # Reshape direction for broadcasting: [D, 1, 1, ...]
    view_shape = [ndim] + [1] * ndim
    dir_view = direction.view(*view_shape)
    directional_derivative = (dir_view * grad).sum(dim=0)  # [*grid_shape]

    # Apply formula: v = -c * (n · ∇ξ)
    v = -wave_speed * directional_derivative.reshape(N)

    # Set velocities in the field component
    state.velocities[:, field_component] = v

    # Respect fixed boundary conditions
    if state.fixed_mask is not None:
        state.velocities[state.fixed_mask, field_component] = 0.0

    print(f"\nInitialized right-moving velocities:")
    print(f"  Wave speed: {wave_speed:.6e} m/s")
    print(f"  Direction: {direction.cpu().numpy()}")
    print(f"  Max |v|: {torch.abs(v).max().item():.6e} m/s")


def initialize_right_moving_velocities_time_reversed(
    state: BraneState,
    grid: BraneGrid,
    physics: "SpringForceComputer",
    m_point: float,
    wave_speed: float,
    field_component: int = 3,
    shift_cells: int = 1,
) -> None:
    """
    Initialize velocities for a right-moving wave by *time-reversing* the
    brane dynamics over Δt = (shift_cells * h) / c.

    This method uses the ACTUAL brane forces (including pretension, geometric
    coupling, and nonlinear saturation) to compute initial velocities, rather
    than assuming a simple scalar wave model. This is more consistent with
    the full 4D geometry and should reduce artificial splitting artifacts.

    Physical principle:
        - Current state has initial shape ξ₀(x) in positions[:, field_component]
        - We define target shape ξ_target(x) which is ξ₀ shifted by
          'shift_cells' grid cells to the right
        - Using the full brane forces F (from SpringForceComputer), we compute
          the acceleration a₀ at the initial state
        - We solve the 2nd-order Taylor expansion:
              ξ_target ≈ ξ₀ + v₀ Δt + 0.5 a₀ Δt²
          for v₀, giving:
              v₀ = (ξ_target - ξ₀ - 0.5 a₀ Δt²) / Δt

    This gives initial velocities consistent with the actual brane model
    (up to O(Δt³) errors), accounting for lateral distortion via the 4D
    geometric coupling.

    Args:
        state: BraneState with initial positions already set (shape only)
        grid: BraneGrid (currently 1D) with spacing h
        physics: SpringForceComputer used in the simulation
        m_point: Point mass (kg) for each brane node
        wave_speed: Target propagation speed (should equal √(T/ρ_m) = c)
        field_component: Index of embedding component containing the wave
                        (default: 3 for X⁴)
        shift_cells: How many grid cells the packet should move in Δt (default: 1)

    Note:
        - Currently implemented for 1D grids
        - Only the amplitude component velocities are explicitly set; lateral
          distortions emerge naturally from the forces during evolution
        - This method respects the substrate-only evolution principle: it uses
          only the microscopic spring forces, with no back-reaction from
          emergent fields
    """
    device = state.device
    dtype = state.dtype
    grid_shape = grid.grid_shape

    # Currently only 1D is implemented
    if len(grid_shape) != 1:
        raise ValueError(
            "initialize_right_moving_velocities_time_reversed currently "
            "supports only 1D grids."
        )

    nx = grid_shape[0]
    N = state.num_points
    if N != nx:
        raise ValueError(f"State and grid size mismatch: N={N}, nx={nx}")

    h = grid.spacing
    dt_shift = (shift_cells * h) / wave_speed  # time to move 'shift_cells' cells at speed c

    # --- 1) Extract current amplitude field ξ₀ -----------------------
    xi0 = state.positions[:, field_component].to(device=device, dtype=dtype)

    # --- 2) Build target field ξ_target = ξ₀ shifted right -----------
    # Shift by 'shift_cells' cells; fill left boundary with 0 (fixed boundary)
    xi_target = torch.zeros_like(xi0)
    if shift_cells < nx:
        xi_target[shift_cells:] = xi0[:-shift_cells]

    # Enforce fixed boundaries if present
    if state.fixed_mask is not None:
        xi_target[state.fixed_mask] = 0.0

    # --- 3) Compute accelerations a₀ from full brane forces ----------
    # Forces: [N, 4] from all springs, including pretension and 4D geometry
    forces = physics.compute_forces(state, grid)  # [N, 4]

    # Acceleration in amplitude component
    a_xi = forces[:, field_component] / m_point  # [N]

    # --- 4) Solve Taylor expansion for v₀ in amplitude component -----
    dt = torch.tensor(dt_shift, device=device, dtype=dtype)
    half_dt2 = 0.5 * dt * dt

    # v₀ = (ξ_target - ξ₀ - 0.5 a₀ Δt²) / Δt
    v_xi = (xi_target - xi0 - half_dt2 * a_xi) / dt

    # --- 5) Write velocities back to state ----------------------------
    # Only set the amplitude component; leave lateral components as they were
    # (typically zero). Lateral distortion will be generated dynamically by
    # the forces.
    state.velocities[:, field_component] = v_xi

    # Respect fixed boundaries: clamp velocities at fixed nodes
    if state.fixed_mask is not None:
        state.velocities[state.fixed_mask, field_component] = 0.0

    # Diagnostics
    print("\nInitialized right-moving velocities (time-reversed):")
    print(f"  dt_shift       = {dt_shift:.6e} s")
    print(f"  grid spacing h = {h:.6e} m")
    print(f"  nominal speed  = {wave_speed:.6e} m/s")
    print(f"  shift_cells    = {shift_cells}")
    print(f"  max |v_xi|     = {torch.abs(v_xi).max().item():.6e} m/s")


def measure_wave_speed(
    state: BraneState,
    grid: BraneGrid,
    field_component: int = 3,
    threshold_fraction: float = 0.5,
) -> tuple[float, float]:
    """
    Measure the center of mass of a wave packet.

    This can be used to track wave propagation and verify it moves at speed c.

    Args:
        state: Current brane state
        grid: Grid with coordinate information
        field_component: Which component contains the wave (default: 3 for X⁴)
        threshold_fraction: Only include points above this fraction of max amplitude

    Returns:
        Tuple of (center_x, center_y) in meters
        For 1D, center_y will be 0.0
    """
    field = state.positions[:, field_component]

    # Find significant points (above threshold)
    max_amp = torch.abs(field).max()
    threshold = threshold_fraction * max_amp
    significant = torch.abs(field) > threshold

    if not significant.any():
        # No significant amplitude
        return 0.0, 0.0

    # Get coordinates of significant points
    coords = grid.get_spatial_coordinates()  # [N, D]
    sig_coords = coords[significant]
    sig_amp = torch.abs(field[significant])

    # Weighted center
    total_weight = sig_amp.sum()
    center_x = (sig_amp * sig_coords[:, 0]).sum() / total_weight

    if coords.shape[1] > 1:
        center_y = (sig_amp * sig_coords[:, 1]).sum() / total_weight
        return center_x.item(), center_y.item()
    else:
        return center_x.item(), 0.0


def verify_wave_propagation(
    state: BraneState,
    grid: BraneGrid,
    wave_speed_expected: float,
    time_elapsed: float,
    initial_center_x: float,
    field_component: int = 3,
) -> None:
    """
    Verify that a wave packet has propagated the expected distance.

    Prints diagnostic information comparing expected vs actual displacement.

    Args:
        state: Current brane state
        grid: Grid with coordinate information
        wave_speed_expected: Expected wave speed c [m/s]
        time_elapsed: Time since initialization [s]
        initial_center_x: Initial x-coordinate of wave center [m]
        field_component: Which component contains the wave (default: 3)
    """
    center_x, _ = measure_wave_speed(state, grid, field_component)

    displacement = center_x - initial_center_x
    expected_displacement = wave_speed_expected * time_elapsed

    if expected_displacement > 0:
        actual_speed = displacement / time_elapsed
        relative_error = abs(actual_speed - wave_speed_expected) / wave_speed_expected

        print(f"\n--- Wave Propagation Verification ---")
        print(f"Time elapsed: {time_elapsed:.6e} s")
        print(f"Initial center: {initial_center_x:.6e} m")
        print(f"Current center: {center_x:.6e} m")
        print(f"Displacement: {displacement:.6e} m")
        print(f"Expected displacement: {expected_displacement:.6e} m")
        print(f"Expected speed: {wave_speed_expected:.6e} m/s")
        print(f"Measured speed: {actual_speed:.6e} m/s")
        print(f"Relative error: {relative_error * 100:.4f}%")

        if relative_error < 0.01:
            print("✓ Wave speed matches c within 1%")
        elif relative_error < 0.05:
            print("⚠ Wave speed matches c within 5%")
        else:
            print(f"✗ Wave speed error {relative_error*100:.1f}% exceeds tolerance")