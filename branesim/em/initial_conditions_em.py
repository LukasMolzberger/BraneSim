"""
Initial conditions for EM field simulations.

This module provides functions to initialize the four-potential A^μ with
various physically meaningful configurations (plane waves, Gaussian pulses, etc.).
"""

import torch
import math


def initialize_plane_wave_Ay(
    state,
    grid,
    c: float,
    amplitude: float,
    wavelength: float,
    phase: float = 0.0,
    direction_axis: int = 0,
):
    """
    Initialize a vacuum plane wave with A along y-direction.

    Creates a plane wave propagating along the specified axis with:
        A_y(x, t) = A₀ sin(kx - ωt + φ)
        A^0 = 0 (scalar potential)

    This automatically satisfies the Lorenz gauge (∇·A = 0 in this case).

    At t=0:
        A_y = A₀ sin(kx + φ)
        ∂A_y/∂t = -ω A₀ cos(kx + φ)

    The electric field will be E_y = -∂A_y/∂t and the magnetic field
    will be B_z = -∂A_y/∂x (for propagation along x).

    Args:
        state: EMState instance to initialize
        grid: BraneGrid instance with spatial coordinates
        c: Speed of light in m/s
        amplitude: Wave amplitude A₀ (in units such that E ~ ω A₀)
        wavelength: Wavelength λ in meters
        phase: Initial phase φ in radians
        direction_axis: Axis along which wave propagates (0=x, 1=y, 2=z)
    """
    k = 2.0 * math.pi / wavelength
    omega = c * k

    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()  # [N, D] in meters
    x = coords[:, direction_axis]

    # Compute initial values at t=0
    Ay = amplitude * torch.sin(k * x + phase).to(state.device, state.dtype)
    Ay_dot = (-omega * amplitude * torch.cos(k * x + phase)).to(state.device, state.dtype)

    # Zero out all components
    state.potential.zero_()
    state.velocity.zero_()

    # Set only y-component (index 2: A^0=0, Ax=1, Ay=2, Az=3)
    state.potential[:, 2] = Ay
    state.velocity[:, 2] = Ay_dot

    state.apply_fixed_boundaries()


def initialize_plane_wave_general(
    state,
    grid,
    c: float,
    amplitude: float,
    wavelength: float,
    polarization: torch.Tensor,
    propagation_dir: torch.Tensor,
    phase: float = 0.0,
):
    """
    Initialize a general plane wave with arbitrary polarization and direction.

    Creates a plane wave: A(x, t) = A₀ ê_pol sin(k·x - ωt + φ)
    where ê_pol is the polarization direction (must be perpendicular to k).

    Args:
        state: EMState instance to initialize
        grid: BraneGrid instance
        c: Speed of light in m/s
        amplitude: Wave amplitude A₀
        wavelength: Wavelength λ in meters
        polarization: Polarization vector [3] (will be normalized)
        propagation_dir: Propagation direction [3] (will be normalized)
        phase: Initial phase φ in radians

    Raises:
        ValueError: If polarization is not perpendicular to propagation direction
    """
    # Normalize directions
    pol = polarization / torch.norm(polarization)
    prop = propagation_dir / torch.norm(propagation_dir)

    # Check orthogonality (gauge condition)
    if abs(torch.dot(pol, prop).item()) > 1e-6:
        raise ValueError("Polarization must be perpendicular to propagation direction")

    k_mag = 2.0 * math.pi / wavelength
    omega = c * k_mag
    k_vec = k_mag * prop  # [3]

    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()  # [N, D]

    # Pad coords to 3D if needed
    if coords.shape[1] < 3:
        coords_3d = torch.zeros((coords.shape[0], 3), device=coords.device, dtype=coords.dtype)
        coords_3d[:, :coords.shape[1]] = coords
        coords = coords_3d

    # Compute k·x
    k_dot_x = torch.sum(coords * k_vec.to(coords.device), dim=1)  # [N]

    # Compute wave profile
    phi = torch.sin(k_dot_x + phase).to(state.device, state.dtype)
    phi_dot = (-omega * torch.cos(k_dot_x + phase)).to(state.device, state.dtype)

    # Zero out all components
    state.potential.zero_()
    state.velocity.zero_()

    # Set vector potential: A = A₀ ê_pol sin(...)
    for i in range(3):
        state.potential[:, i + 1] = amplitude * pol[i] * phi
        state.velocity[:, i + 1] = amplitude * pol[i] * phi_dot

    state.apply_fixed_boundaries()


def initialize_gaussian_pulse(
    state,
    grid,
    center: torch.Tensor,
    width: float,
    amplitude: float,
    polarization: torch.Tensor,
):
    """
    Initialize a localized Gaussian pulse in the vector potential.

    Creates: A(x) = A₀ ê_pol exp(-|x - x₀|²/(2σ²))

    This is not a wave solution but useful for studying pulse propagation
    and dispersion.

    Args:
        state: EMState instance to initialize
        grid: BraneGrid instance
        center: Pulse center position [D] in meters
        width: Gaussian width σ in meters
        amplitude: Peak amplitude A₀
        polarization: Polarization direction [3] (will be normalized)
    """
    pol = polarization / torch.norm(polarization)

    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()  # [N, D]

    # Compute distance from center
    center_padded = center.to(coords.device)[:coords.shape[1]]
    r_squared = torch.sum((coords - center_padded) ** 2, dim=1)

    # Gaussian envelope
    envelope = amplitude * torch.exp(-r_squared / (2.0 * width * width))
    envelope = envelope.to(state.device, state.dtype)

    # Zero out all components
    state.potential.zero_()
    state.velocity.zero_()

    # Set vector potential: A = A₀ ê_pol exp(...)
    for i in range(3):
        state.potential[:, i + 1] = pol[i] * envelope

    # Initial velocity is zero for static pulse
    state.apply_fixed_boundaries()