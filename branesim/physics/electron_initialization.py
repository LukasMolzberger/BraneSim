"""
Electron Initialization for W&vdM Toroidal Double-Loop Model

This module implements the Williamson & van der Mark toroidal electron model
adapted to the brane simulation framework. It provides:

1. Analytical ansatz for the electron structure (tubular coordinates, double-loop)
2. Physical parameter calibration (Compton scale, charge, energy, spin)
3. Optimization-ready parameterization for stability refinement

The electron is modeled as a toroidal soliton with:
- Centerline at Compton scale (R ~ λ_C / 2π)
- Double-loop cross-section (two opposite lobes)
- Internal Compton-frequency oscillation (ω_C = m_e c² / ℏ)
- Self-consistent bound mode via geometric nonlinearity

Implementation follows PROJECT_PRINCIPLES.md:
- Substrate-only evolution (positions and velocities only)
- No back-reaction from emergent fields
- Pure geometric coupling via 4D distance
"""

import torch
import math
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional

from branesim.core.state import BraneState
from branesim.config.simulation_config import PhysicalConstants


@dataclass
class ElectronInitParams:
    """
    Parameters for electron initialization in rest frame.

    These parameters define the geometry and field structure of the
    toroidal electron. Most can be optimized for stability.

    Geometric Parameters:
        center: (x, y, z) coordinates of torus center in physical units [m]
        R: Major radius of torus (centerline radius) [m]
        rho0: Peak radius for cross-section envelope [m]
        sigma_r: Radial width of cross-section envelope [m]
        sigma_theta: Angular width of each tunnel [rad]

    Twist Parameters:
        l_twist: Twist winding number (ℓ ∈ ℤ).
                 ℓ=0 → symmetric double lobe (no twist)
                 ℓ=1 → W&vdM intertwined double loop (one rotation per torus revolution)
        alpha0: Initial twist angle offset [rad]

    Field Parameters:
        A: Amplitude scale for X⁴ Compton mode [m]
        phase_offset: Global phase offset for internal oscillation [rad]

    Physical Parameters:
        compton_omega: Compton frequency ω_C [rad/s]
        wave_speed: Effective wave speed c in brane [m/s]
        spin_axis: (x, y, z) unit vector defining expected spin axis [unitless]

    Numerical Controls:
        tube_max_radius: Maximum transverse distance from centerline [m]
        lateral_spin_scale: Scale factor for lateral spin velocity (0 = no spin)
    """
    # Geometric parameters
    center: Tuple[float, float, float]
    R: float
    rho0: float
    sigma_r: float

    # Field parameters
    A: float

    # Physical parameters
    compton_omega: float
    wave_speed: float

    # Numerical controls
    tube_max_radius: float

    # Fields with defaults (must come last in dataclass)
    sigma_theta: float = 0.5  # Angular width in radians (default ~30°)
    l_twist: int = 1  # Twist winding number (ℓ=1 for W&vdM-like)
    alpha0: float = 0.0  # Initial twist angle offset
    phase_offset: float = 0.0
    lateral_spin_scale: float = 0.0
    spin_axis: Tuple[float, float, float] = (0.0, 0.0, 1.0)
    winding_number: int = 2  # m=2 for longitudinal winding

    def __repr__(self) -> str:
        return (
            f"ElectronInitParams(\n"
            f"  center={self.center},\n"
            f"  R={self.R:.6e} m,\n"
            f"  rho0={self.rho0:.6e} m,\n"
            f"  sigma_r={self.sigma_r:.6e} m,\n"
            f"  A={self.A:.6e} m,\n"
            f"  ω_C={self.compton_omega:.6e} rad/s,\n"
            f"  tube_max_radius={self.tube_max_radius:.6e} m\n"
            f")"
        )


def compute_tubular_coords_for_point(
    pos: torch.Tensor,
    center: torch.Tensor,
    R: float,
) -> Tuple[float, float, float]:
    """
    Map a 3D brane position to tubular coordinates (z, x, y) around a circular
    centerline of radius R in the X1-X2 plane (with X3 as vertical/binormal).

    Tubular coordinates:
        z: Arclength along the circle centerline [0, 2πR)
        x: Transverse coordinate in radial direction (in-plane)
        y: Transverse coordinate along binormal (X3 direction)

    The torus centerline is a circle in the X1-X2 plane centered at `center`.

    Args:
        pos: Position tensor [3] in physical units (X1, X2, X3)
        center: Center of torus [3] in physical units
        R: Major radius of torus [m]

    Returns:
        (z, x, y): Tubular coordinates as Python floats
    """
    # Shift to frame where center is origin
    v = pos - center  # [3]

    # Project into X1-X2 plane
    vx, vy, vz = v[0], v[1], v[2]
    r_xy = torch.sqrt(vx * vx + vy * vy + 1e-30)

    # Angular coordinate around circle
    phi = torch.atan2(vy, vx)  # range (-π, π]
    if phi < 0:
        phi = phi + 2.0 * torch.pi  # [0, 2π)

    # Centerline point on circle at this angle
    cx = R * torch.cos(phi)
    cy = R * torch.sin(phi)
    cz = torch.tensor(0.0, dtype=pos.dtype, device=pos.device)
    C = torch.stack((cx, cy, cz)) + center

    # Frenet frame for circle in X1-X2 plane:
    # Normal n points radially outward in the plane
    # Binormal b points along X3
    n = torch.stack([
        torch.cos(phi),
        torch.sin(phi),
        torch.tensor(0.0, device=pos.device, dtype=pos.dtype)
    ])
    b = torch.tensor([0.0, 0.0, 1.0], dtype=pos.dtype, device=pos.device)

    # Transverse vector from centerline to pos
    delta = pos - C

    # x = component along n (radial in-plane)
    x = torch.dot(delta, n)

    # y = component along binormal (X3)
    y = torch.dot(delta, b)

    # Arclength coordinate z = R * φ
    z = R * phi

    return z.item(), x.item(), y.item()


def compute_tubular_coords_vectorized(
    positions: torch.Tensor,
    center: torch.Tensor,
    R: float,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Vectorized version of tubular coordinate computation.

    Args:
        positions: [N, 3] tensor of positions in physical units
        center: [3] tensor of torus center
        R: Major radius

    Returns:
        (z, x, y): Each is [N] tensor of tubular coordinates
    """
    # Shift to frame where center is origin
    v = positions - center.unsqueeze(0)  # [N, 3]

    vx, vy, vz = v[:, 0], v[:, 1], v[:, 2]
    r_xy = torch.sqrt(vx * vx + vy * vy + 1e-30)

    # Angular coordinate
    phi = torch.atan2(vy, vx)
    phi = torch.where(phi < 0, phi + 2.0 * torch.pi, phi)

    # Centerline points
    cx = R * torch.cos(phi)
    cy = R * torch.sin(phi)
    cz = torch.zeros_like(cx)
    C = torch.stack([cx, cy, cz], dim=1) + center.unsqueeze(0)  # [N, 3]

    # Normal and binormal
    n = torch.stack([torch.cos(phi), torch.sin(phi), torch.zeros_like(phi)], dim=1)  # [N, 3]
    b = torch.tensor([0.0, 0.0, 1.0], device=positions.device, dtype=positions.dtype).expand(positions.shape[0], 3)

    # Transverse displacement
    delta = positions - C  # [N, 3]

    # Transverse coordinates
    x = (delta * n).sum(dim=1)  # [N]
    y = (delta * b).sum(dim=1)  # [N]

    # Arclength coordinate
    z = R * phi  # [N]

    return z, x, y


def twisted_tunnel_envelope(
    x: torch.Tensor,
    y: torch.Tensor,
    z: torch.Tensor,
    R: float,
    rho0: float,
    sigma_r: float,
    sigma_theta: float,
    l_twist: int,
    alpha0: float,
) -> torch.Tensor:
    """
    Unified twisted-tunnel cross-section envelope f(ρ,θ;z) implementing the
    W&vdM model with adjustable twist parameter.

    This implements the unified tubular ansatz:
        f(ρ,θ;z) = G(ρ) · [exp(-(θ-α(z))²/(2σ_θ²)) + exp(-(θ-α(z)-π)²/(2σ_θ²))]

    where:
        - (ρ,θ) are polar coordinates on the cross-section
        - G(ρ) = exp(-(ρ-ρ₀)²/(2σ_r²)) is the radial envelope
        - α(z) = α₀ + ℓ·φ(z) = α₀ + ℓ·(z/R) controls tunnel twist

    The electron ansatz in this work uses ℓ=1 (W&vdM intertwined double loop),
    so that the ridges rotate once per torus revolution. Other integer ℓ are
    allowed in principle (ℓ=0 would correspond to a symmetric double-lobe
    profile) but are used only for low-level tests, not for the physical
    electron model.

    The two ridges are located at angular positions θ = α(z) and θ = α(z)+π.

    Args:
        x: Transverse coordinate (radial from centerline) [N] or [*grid_shape]
        y: Transverse coordinate (binormal, vertical) [N] or [*grid_shape]
        z: Arclength coordinate along torus [N] or [*grid_shape]
        R: Major radius of torus [m]
        rho0: Peak radius for radial envelope [m]
        sigma_r: Radial width of envelope [m]
        sigma_theta: Angular width of each tunnel [rad]
        l_twist: Twist winding number ℓ ∈ ℤ
        alpha0: Initial twist angle offset α₀ [rad]

    Returns:
        Envelope f(ρ,θ;z) with same shape as input
    """
    # Convert Cartesian (x, y) to polar (ρ, θ) on cross-section
    rho = torch.sqrt(x ** 2 + y ** 2 + 1e-30)
    theta = torch.atan2(y, x)  # range (-π, π]

    # Radial envelope G(ρ) = exp(-(ρ-ρ₀)²/(2σ_r²))
    G_rho = torch.exp(-0.5 * ((rho - rho0) / sigma_r) ** 2)

    # Twist angle α(z) = α₀ + ℓ·φ(z) where φ(z) = z/R
    phi = z / R
    alpha = alpha0 + l_twist * phi

    # Angular envelope with two tunnels at θ = α(z) and θ = α(z)+π
    # Need to handle angular wrapping: use cos distance for periodicity
    # Angular deviation from tunnel 1 at θ = α(z)
    delta_theta_1 = theta - alpha
    # Wrap to [-π, π]
    delta_theta_1 = torch.atan2(torch.sin(delta_theta_1), torch.cos(delta_theta_1))
    tunnel_1 = torch.exp(-0.5 * (delta_theta_1 / sigma_theta) ** 2)

    # Angular deviation from tunnel 2 at θ = α(z)+π
    delta_theta_2 = theta - (alpha + torch.pi)
    # Wrap to [-π, π]
    delta_theta_2 = torch.atan2(torch.sin(delta_theta_2), torch.cos(delta_theta_2))
    tunnel_2 = torch.exp(-0.5 * (delta_theta_2 / sigma_theta) ** 2)

    # Combine: f(ρ,θ;z) = G(ρ) · [tunnel_1 + tunnel_2]
    return G_rho * (tunnel_1 + tunnel_2)


def double_loop_envelope(
    x: torch.Tensor,
    y: torch.Tensor,
    rho0: float,
    sigma_r: float
) -> torch.Tensor:
    """
    Legacy symmetric double-lobe envelope (ℓ=0 case of twisted_tunnel_envelope).

    This is kept for backward compatibility with existing test code.
    New code should use twisted_tunnel_envelope with l_twist=0.

    Cross-section envelope f(x, y) with two lobes along the binormal (y) direction.
    Two Gaussian peaks separated along the y-axis at y = ±rho0.

    Args:
        x: Transverse coordinate (radial from centerline) [N] or [*grid_shape]
        y: Transverse coordinate (binormal, vertical) [N] or [*grid_shape]
        rho0: Half-separation between lobes [m]
        sigma_r: Width of each lobe [m]

    Returns:
        Envelope f(x, y) with same shape as input
    """
    # Two Gaussian lobes centered at y = +rho0 and y = -rho0
    # Both lobes have small x extent (close to centerline)
    lobe_plus = torch.exp(-0.5 * ((y - rho0) / sigma_r) ** 2 - 0.5 * (x / sigma_r) ** 2)
    lobe_minus = torch.exp(-0.5 * ((y + rho0) / sigma_r) ** 2 - 0.5 * (x / sigma_r) ** 2)

    return lobe_plus + lobe_minus


def init_electron_amplitude(
    state: BraneState,
    params: ElectronInitParams,
    mask: torch.Tensor,
) -> None:
    """
    Initialize X⁴ (amplitude) and dX⁴/dt for the electron region using the
    unified twisted-tunnel W&vdM ansatz.

    The field is initialized as:
        ξ(z, x, y, t=0) = A · f(ρ,θ;z) · cos(m·k_C·z + φ₀)
        ∂_t ξ|_{t=0} = -A · ω_C · f(ρ,θ;z) · sin(m·k_C·z + φ₀)

    where:
        - f(ρ,θ;z) is the twisted-tunnel cross-section envelope
        - k_C = ω_C / c is the Compton wave number
        - m is the longitudinal winding number
        - z is the arclength coordinate along the torus

    Args:
        state: BraneState to initialize
        params: ElectronInitParams with geometry and field parameters
        mask: Boolean mask [N] indicating electron tube region
    """
    center = torch.tensor(params.center, dtype=state.dtype, device=state.device)
    R = params.R
    rho0 = params.rho0
    sigma_r = params.sigma_r
    sigma_theta = params.sigma_theta
    l_twist = params.l_twist
    alpha0 = params.alpha0
    A = params.A
    omega_C = params.compton_omega
    c_eff = params.wave_speed
    k_C = omega_C / c_eff
    phase_offset = params.phase_offset

    # Get lateral coordinates (X^0, X^1, X^2)
    X_lat = state.positions[:, :3]  # [N, 3]

    # Compute tubular coordinates for all points (vectorized)
    z, x, y = compute_tubular_coords_vectorized(X_lat, center, R)

    # Build twisted-tunnel cross-section envelope f(ρ,θ;z)
    f_xyz = twisted_tunnel_envelope(
        x, y, z, R, rho0, sigma_r, sigma_theta, l_twist, alpha0
    )

    # Phase at t=0 with m-fold winding: φ = -m * k_C * z + φ₀
    # For m=2 (double loop), phase completes TWO cycles around torus
    # Combined with l_twist, this creates the (m,ℓ) torus-knot pattern
    m = params.winding_number
    phase0 = -m * k_C * z + phase_offset

    # Amplitude displacement X⁴
    xi0 = A * f_xyz * torch.cos(phase0)

    # Amplitude velocity dX⁴/dt
    dxi_dt0 = -A * omega_C * f_xyz * torch.sin(phase0)

    # Apply only on masked region
    amp_idx = 3  # X⁴ is the 4th component (index 3)
    state.positions[:, amp_idx] = torch.where(mask, xi0, state.positions[:, amp_idx])
    state.velocities[:, amp_idx] = torch.where(mask, dxi_dt0, state.velocities[:, amp_idx])

    # Determine twist description
    if l_twist == 0:
        twist_desc = "symmetric (no twist)"
    elif l_twist == 1:
        twist_desc = "W&vdM intertwined (one rotation per revolution)"
    else:
        twist_desc = f"{l_twist} rotations per revolution"

    print(f"\n=== Initialized Electron Amplitude Field ===")
    print(f"  Model: Unified twisted-tunnel W&vdM ansatz")
    print(f"  Torus major radius R = {R:.6e} m")
    print(f"  Cross-section peak radius ρ₀ = {rho0:.6e} m")
    print(f"  Radial width σ_r = {sigma_r:.6e} m")
    print(f"  Angular width σ_θ = {sigma_theta:.3f} rad")
    print(f"  Twist winding ℓ = {l_twist} ({twist_desc})")
    print(f"  Twist offset α₀ = {alpha0:.3f} rad")
    print(f"  Longitudinal winding m = {m} (Compton phase cycles)")
    print(f"  Torus-knot pattern: ({m},{l_twist})-type")
    print(f"  Amplitude scale A = {A:.6e} m")
    print(f"  Compton frequency ω_C = {omega_C:.6e} rad/s")
    print(f"  Wave number k_C = {k_C:.6e} rad/m")
    print(f"  Effective wave number m*k_C = {m*k_C:.6e} rad/m")
    print(f"  Max |ξ| = {torch.abs(xi0[mask]).max().item():.6e} m")
    print(f"  Max |v_ξ| = {torch.abs(dxi_dt0[mask]).max().item():.6e} m/s")


def init_electron_lateral_geometry(
    state: BraneState,
    params: ElectronInitParams,
    mask: torch.Tensor,
) -> None:
    """
    Optionally imprint a tubular lateral deformation into X^0, X^1, X^2.

    This pre-shapes the brane to have a toroidal "tunnel" structure modulated
    by the same cross-section envelope. This can help stabilize the electron
    by providing a pre-formed geometric channel for the mode to live in.

    The lateral displacement is:
        Δr = α · f(x, y) · n̂

    where n̂ is the radial normal direction and α ~ sigma_r.

    Args:
        state: BraneState to modify
        params: ElectronInitParams
        mask: Boolean mask [N] for electron tube region
    """
    center = torch.tensor(params.center, dtype=state.dtype, device=state.device)
    R = params.R
    rho0 = params.rho0
    sigma_r = params.sigma_r

    X_lat = state.positions[:, :3]  # [N, 3]

    # Compute tubular coords
    z, x, y = compute_tubular_coords_vectorized(X_lat, center, R)

    # Cross-section envelope
    f_xy = double_loop_envelope(x, y, rho0=rho0, sigma_r=sigma_r)

    # Local radial direction n̂(φ) in X1-X2 plane
    vx = X_lat[:, 0] - center[0]
    vy = X_lat[:, 1] - center[1]
    phi = torch.atan2(vy, vx)

    n_x = torch.cos(phi)
    n_y = torch.sin(phi)
    n_z = torch.zeros_like(n_x)
    n_vec = torch.stack([n_x, n_y, n_z], dim=1)  # [N, 3]

    # Lateral displacement scale (small fraction of sigma_r)
    alpha = 0.5 * sigma_r

    delta_lat = alpha * f_xy.unsqueeze(1) * n_vec  # [N, 3]

    # Apply displacement only in masked region
    state.positions[:, :3] = torch.where(
        mask.unsqueeze(1),
        X_lat + delta_lat,
        X_lat
    )

    print(f"\n=== Applied Lateral Geometry Deformation ===")
    print(f"  Displacement scale α = {alpha:.6e} m")
    print(f"  Max |Δr| = {torch.norm(delta_lat[mask], dim=1).max().item():.6e} m")


def init_electron_lateral_spin_velocity(
    state: BraneState,
    params: ElectronInitParams,
    mask: torch.Tensor,
) -> None:
    """
    Optionally add a lateral spin velocity field inside the tube to give the
    electron internal angular momentum ~ ℏ/2.

    The velocity field is tangent to the torus centerline and localized by
    the same radial envelope. This adds circulation around the torus.

    Args:
        state: BraneState to modify
        params: ElectronInitParams
        mask: Boolean mask [N] for electron tube region
    """
    if params.lateral_spin_scale == 0.0:
        return  # No spin requested

    center = torch.tensor(params.center, dtype=state.dtype, device=state.device)
    R = params.R
    rho0 = params.rho0
    sigma_r = params.sigma_r

    X_lat = state.positions[:, :3]  # [N, 3]

    # Tangent direction along the circle (azimuthal)
    vx = X_lat[:, 0] - center[0]
    vy = X_lat[:, 1] - center[1]
    phi = torch.atan2(vy, vx)

    t_x = -torch.sin(phi)
    t_y = torch.cos(phi)
    t_z = torch.zeros_like(t_x)
    t_vec = torch.stack([t_x, t_y, t_z], dim=1)  # [N, 3]

    # Radial envelope to localize spin near rho0
    z, x, y = compute_tubular_coords_vectorized(X_lat, center, R)
    rho = torch.sqrt(x * x + y * y + 1e-30)
    radial = torch.exp(-0.5 * ((rho - rho0) / sigma_r) ** 2)

    spin_scale = params.lateral_spin_scale
    v_spin = spin_scale * radial.unsqueeze(1) * t_vec  # [N, 3]

    # Add to velocities in masked region
    state.velocities[:, :3] = torch.where(
        mask.unsqueeze(1),
        state.velocities[:, :3] + v_spin,
        state.velocities[:, :3]
    )

    print(f"\n=== Applied Lateral Spin Velocity ===")
    print(f"  Spin scale = {spin_scale:.6e} m/s")
    print(f"  Max |v_spin| = {torch.norm(v_spin[mask], dim=1).max().item():.6e} m/s")


def init_electron_state(
    state: BraneState,
    params: ElectronInitParams,
) -> BraneState:
    """
    Initialize a single electron soliton in the given BraneState using the
    W&vdM toroidal double-loop ansatz.

    This is the main entry point for electron initialization. It:
    1. Computes the tube mask (region where electron lives)
    2. Initializes amplitude field X⁴ with Compton oscillation
    3. Optionally pre-shapes lateral geometry (tunnel)
    4. Optionally adds lateral spin velocity

    The electron is initialized in the rest frame with:
    - Centerline at Compton scale R ~ λ_C / 2π
    - Double-loop cross-section with two opposite lobes
    - Internal Compton-frequency mode ω_C = m_e c² / ℏ

    Args:
        state: BraneState to initialize (must have flat grid positions set)
        params: ElectronInitParams with all geometry and field parameters

    Returns:
        Modified state with electron initialized
    """
    print(f"\n{'='*60}")
    print(f"Initializing Electron (W&vdM Toroidal Double-Loop Model)")
    print(f"{'='*60}")
    print(params)

    center = torch.tensor(params.center, dtype=state.dtype, device=state.device)
    X_lat = state.positions[:, :3]  # [N, 3]

    # Compute mask of points belonging to the tube region:
    # Points within tube_max_radius of the circular centerline
    v = X_lat - center
    vx, vy, vz = v[:, 0], v[:, 1], v[:, 2]
    r_xy = torch.sqrt(vx * vx + vy * vy + 1e-30)

    # Distance from circle centerline in XY plane
    R = params.R
    radial_offset = r_xy - R

    # Transverse distance (radial in plane + vertical)
    transverse = torch.sqrt(radial_offset * radial_offset + vz * vz)
    mask = transverse <= params.tube_max_radius

    num_electron = mask.sum().item()
    num_total = state.num_points
    print(f"\n=== Electron Tube Region ===")
    print(f"  Tube max radius = {params.tube_max_radius:.6e} m")
    print(f"  Points in electron region: {num_electron}/{num_total} ({100*num_electron/num_total:.2f}%)")

    # Initialize amplitude field and velocities
    init_electron_amplitude(state, params, mask)

    # Optional: imprint lateral geometry
    # init_electron_lateral_geometry(state, params, mask)

    # Optional: add lateral spin
    # init_electron_lateral_spin_velocity(state, params, mask)

    # Optional: enforce net momentum ~ 0 by subtracting center-of-mass velocity
    # (Not needed in rest frame if initialization is symmetric)

    print(f"\n{'='*60}")
    print(f"Electron Initialization Complete")
    print(f"{'='*60}\n")

    return state


def default_electron_geometry(
    constants: PhysicalConstants,
    grid_spacing: float
) -> Tuple[float, float, float]:
    """
    Compute default geometric parameters for electron from Compton scale.

    Default choices:
        - Torus circumference L = λ_C (one Compton wavelength)
        - R = λ_C / (2π)
        - rho0 = 0.3 * R (lobes at 30% of major radius)
        - sigma_r = 0.1 * R (narrow lobes)

    Args:
        constants: PhysicalConstants with m_e, c, ℏ
        grid_spacing: Grid spacing h [m] for validation

    Returns:
        (R, rho0, sigma_r): Torus major radius, lobe radius, lobe width [m]
    """
    lambda_C = constants.lambda_C

    R_phys = lambda_C / (2.0 * math.pi)
    rho0 = 0.3 * R_phys
    sigma_r = 0.1 * R_phys

    print(f"\n=== Default Electron Geometry ===")
    print(f"  Compton wavelength λ_C = {lambda_C:.6e} m")
    print(f"  Major radius R = λ_C/(2π) = {R_phys:.6e} m")
    print(f"  Lobe radius rho0 = 0.3R = {rho0:.6e} m")
    print(f"  Lobe width sigma_r = 0.1R = {sigma_r:.6e} m")
    print(f"  Grid spacing h = {grid_spacing:.6e} m")
    print(f"  R/h = {R_phys/grid_spacing:.2f}")

    return R_phys, rho0, sigma_r


def plot_cross_section_envelope_debug(
    params: ElectronInitParams,
    n_points: int = 200,
    save_path: str = 'debug_cross_section.png'
) -> None:
    """
    Debug plot: visualize the cross-section envelope f(x, y).

    This should show a clean dumbbell shape with two lobes separated
    along the y (binormal) axis.

    Args:
        params: ElectronInitParams with geometry
        n_points: Resolution for plotting
        save_path: Path to save the figure
    """
    import matplotlib.pyplot as plt

    # Create grid in transverse (x, y) plane
    extent = 3.0 * max(params.rho0, params.sigma_r)  # Plot ±3σ range
    x = np.linspace(-extent, extent, n_points)
    y = np.linspace(-extent, extent, n_points)
    X, Y = np.meshgrid(x, y)

    # Convert to tensors
    X_t = torch.from_numpy(X).float()
    Y_t = torch.from_numpy(Y).float()

    # Compute envelope
    envelope = double_loop_envelope(X_t, Y_t, params.rho0, params.sigma_r)
    envelope_np = envelope.numpy()

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 2D heatmap
    ax = axes[0]
    im = ax.contourf(X * 1e12, Y * 1e12, envelope_np, levels=20, cmap='viridis')
    ax.set_xlabel('x (radial) [pm]')
    ax.set_ylabel('y (binormal) [pm]')
    ax.set_title('Cross-Section Envelope f(x, y)')
    ax.axhline(y=params.rho0 * 1e12, color='r', linestyle='--', alpha=0.5, label=f'y = +ρ₀')
    ax.axhline(y=-params.rho0 * 1e12, color='r', linestyle='--', alpha=0.5, label=f'y = -ρ₀')
    ax.axvline(x=0, color='w', linestyle='--', alpha=0.3)
    ax.legend()
    plt.colorbar(im, ax=ax)
    ax.set_aspect('equal')

    # 1D cuts
    ax = axes[1]
    # Cut along y axis (x=0)
    y_cut = envelope_np[n_points//2, :]
    ax.plot(y * 1e12, y_cut, 'b-', linewidth=2, label='Cut along y (x=0)')

    # Cut along x axis (y=0)
    x_cut = envelope_np[:, n_points//2]
    ax.plot(x * 1e12, x_cut, 'r-', linewidth=2, label='Cut along x (y=0)')

    ax.set_xlabel('Coordinate [pm]')
    ax.set_ylabel('Envelope amplitude')
    ax.set_title('1D Cuts Through Center')
    ax.axvline(x=params.rho0 * 1e12, color='b', linestyle='--', alpha=0.3)
    ax.axvline(x=-params.rho0 * 1e12, color='b', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Debug plot saved: {save_path}")
    print(f"    Expected: Two lobes at y = ±{params.rho0:.2e} m = ±{params.rho0*1e12:.2f} pm")
    plt.close()


def plot_phase_along_centerline_debug(
    params: ElectronInitParams,
    n_points: int = 500,
    save_path: str = 'debug_phase_centerline.png'
) -> None:
    """
    Debug plot: visualize the phase pattern along the torus centerline.

    This should show a smooth cos(m*φ) pattern with NO discontinuity,
    where m is the winding number.

    Args:
        params: ElectronInitParams with geometry and wave parameters
        n_points: Number of points around the torus
        save_path: Path to save the figure
    """
    import matplotlib.pyplot as plt

    # Angular coordinate around torus
    phi = np.linspace(0, 2 * np.pi, n_points)

    # Arclength coordinate
    z = params.R * phi

    # Wave number
    k_C = params.compton_omega / params.wave_speed

    # Phase with m-fold winding
    m = params.winding_number
    phase = -m * k_C * z + params.phase_offset

    # Field amplitude at centerline (x=0, y=0) at t=0
    # Envelope at centerline depends on cross-section shape
    # For double-lobe at y=±rho0, centerline (x=0, y=0) has intermediate envelope value
    x_center = torch.zeros(n_points)
    y_center = torch.zeros(n_points)
    envelope_center = double_loop_envelope(x_center, y_center, params.rho0, params.sigma_r)
    envelope_np = envelope_center.numpy()

    xi = params.A * envelope_np * np.cos(phase)

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    # Phase vs angle
    ax = axes[0]
    ax.plot(phi, phase, 'b-', linewidth=2)
    ax.set_xlabel('φ [rad]')
    ax.set_ylabel('Phase [rad]')
    ax.set_title(f'Phase Pattern (m={m} winding) - Should be smooth, no jumps')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.grid(True, alpha=0.3)

    # Expected: m cycles of 2π
    expected_phase_range = m * 2 * np.pi
    ax.text(0.02, 0.95, f'Expected phase range: {expected_phase_range/np.pi:.1f}π',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Field amplitude vs angle
    ax = axes[1]
    ax.plot(phi, xi * 1e14, 'r-', linewidth=2)
    ax.set_xlabel('φ [rad]')
    ax.set_ylabel('ξ [×10⁻¹⁴ m]')
    ax.set_title(f'Amplitude Field Along Centerline - Should show {m} complete cycles')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Debug plot saved: {save_path}")
    print(f"    Expected: {m} complete oscillations, smooth everywhere")
    print(f"    Phase range: {phase.min():.2f} to {phase.max():.2f} rad")
    print(f"    Phase jump check: max(|Δφ|) = {np.abs(np.diff(phase)).max():.2e} rad")
    if np.abs(np.diff(phase)).max() > 0.1:
        print(f"    ⚠ WARNING: Large phase jump detected! Check φ wrapping.")
    plt.close()


def calibrate_electron_init_params(
    constants: PhysicalConstants,
    grid_spacing: float,
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    amplitude_scale: float = 1e-13,  # Initial guess for A [m]
) -> ElectronInitParams:
    """
    Create calibrated electron initialization parameters from physical constants.

    This provides a reasonable starting point for the electron ansatz using:
    - Geometry derived from Compton wavelength
    - Amplitude scale based on expected electron size
    - Compton frequency from m_e c² / ℏ

    The amplitude scale is an initial guess; it should be refined by:
    1. Measuring effective charge and adjusting A to match -e
    2. Measuring total energy and adjusting geometry to match m_e c²
    3. Running optimization for stability

    Args:
        constants: PhysicalConstants with m_e, c, ℏ
        grid_spacing: Grid spacing h [m]
        center: Center of torus (x, y, z) [m]
        amplitude_scale: Initial guess for amplitude A [m]

    Returns:
        ElectronInitParams with calibrated values
    """
    # Geometry from Compton scale
    R, rho0, sigma_r = default_electron_geometry(constants, grid_spacing)

    # Tube radius: 3σ envelope
    tube_max_radius = 3.0 * sigma_r

    # Compton frequency and wave speed
    omega_C = constants.m_e * constants.c ** 2 / constants.hbar
    c_eff = constants.c

    # Create params
    params = ElectronInitParams(
        center=center,
        R=R,
        rho0=rho0,
        sigma_r=sigma_r,
        A=amplitude_scale,
        phase_offset=0.0,
        compton_omega=omega_C,
        wave_speed=c_eff,
        tube_max_radius=tube_max_radius,
        lateral_spin_scale=0.0,
    )

    print(f"\n=== Calibrated Electron Parameters ===")
    print(params)

    return params