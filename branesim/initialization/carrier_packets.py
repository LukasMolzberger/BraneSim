"""
Preparation-first narrowband carrier initialization.

This module implements the initialization protocols described in the paper
(Section on "Preparation-first initialization of narrowband carriers").

Key principle:
    The band is defined by PHYSICAL PREPARATION of the substrate initial
    conditions, not by post-hoc filtering. We initialize narrowband carrier
    wave packets with chosen (k0, omega0) and controlled internal polarization
    subspace. Diagnostics only READ OUT the phase and degeneracy structure.
"""

import torch
import numpy as np
from typing import Optional, Tuple
from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


def make_photon_circular_packet(
    state: BraneState,
    grid: BraneGrid,
    k0: float,
    omega0: float,
    amplitude: float,
    center: float | torch.Tensor,
    sigma: float,
    dof_pair: Tuple[int, int] = (2, 3),
    propagation_axis: int = 0,
    periodic_axis: bool = True,
) -> None:
    """
    Initialize a circularly polarized photon-like wave packet.

    Creates a narrowband carrier with 2D polarization subspace, narrowband
    by construction (no diagnostic filtering required).

    Physical setup:
        u(x,0) = Re{A(x) exp(i k0·x) p}
        v(x,0) = Re{-i omega0 A(x) exp(i k0·x) p}

    where:
        - A(x) is a slowly varying Gaussian envelope
        - p ∝ (1, i) gives circular polarization in the two DOFs
        - The real-valued fields are:
            DOF a: A(x) * cos(k0 x)
            DOF b: A(x) * sin(k0 x)

    This creates quadrature components that form a circularly polarized
    carrier without requiring any FFT or band-pass filtering.

    Args:
        state: BraneState to initialize (positions and velocities will be set)
        grid: BraneGrid defining the lattice
        k0: Carrier wave number (2π/λ_C for Compton scale)
        omega0: Carrier angular frequency (typically m_e c^2 / ℏ)
        amplitude: Peak amplitude A (typically λ_C / √π from calibration)
        center: Center position of the envelope (scalar or ndarray for higher dims)
        sigma: Width of Gaussian envelope (should be >> λ_C for narrowband)
        dof_pair: Tuple (a, b) specifying which two DOFs form the polarization plane
                  Default: (2, 3) means X³ and X⁴
        propagation_axis: Spatial axis along which the packet propagates (default: 0 = x)
        periodic_axis: If True, the propagation axis should have periodic boundary conditions

    Notes:
        - Spectrum is narrow by construction: Δk ~ 1/sigma << k0
        - The two DOFs form a genuinely 2D degenerate polarization subspace
        - For verification: use SVD on time series to confirm 2 dominant singular values
        - No FFT masking is used or required

    Paper reference:
        Section "Preparation-first initialization of narrowband carriers"
        Paragraph "Photon packet (two-dimensional polarization subspace)"
    """
    device = state.device
    dtype = state.dtype
    ndim = len(grid.grid_shape)

    # Validate inputs
    dof_a, dof_b = dof_pair
    if dof_a == dof_b:
        raise ValueError("DOF pair must contain two distinct components")
    if not (0 <= dof_a < 4 and 0 <= dof_b < 4):
        raise ValueError("DOF indices must be in [0,3]")
    if propagation_axis >= ndim:
        raise ValueError(f"Propagation axis {propagation_axis} invalid for {ndim}D grid")

    # Get spatial coordinates
    coords = grid.get_spatial_coordinates()  # [N, ndim]
    x_prop = coords[:, propagation_axis]  # Propagation coordinate

    # Build envelope A(x) - Gaussian centered at 'center'
    if isinstance(center, (int, float)):
        center_tensor = torch.tensor([center], device=device, dtype=dtype)
    else:
        center_tensor = torch.as_tensor(center, device=device, dtype=dtype)

    # For multidimensional grids, compute radial distance from center
    if ndim == 1:
        r_sq = (coords[:, 0] - center_tensor[0]) ** 2
    else:
        # Compute distance in all dimensions from center
        center_full = torch.zeros(ndim, device=device, dtype=dtype)
        if len(center_tensor) == 1:
            center_full[propagation_axis] = center_tensor[0]
        else:
            center_full = center_tensor
        r_sq = torch.sum((coords - center_full) ** 2, dim=1)

    envelope = amplitude * torch.exp(-r_sq / (2 * sigma ** 2))

    # Carrier phase: φ(x) = k0 * x_prop
    carrier_phase = k0 * x_prop

    # Build quadrature components for circular polarization
    # DOF a: A(x) * cos(k0 x)  [Real part of complex polarization (1, 0)]
    # DOF b: A(x) * sin(k0 x)  [Imaginary part becomes sin for p ∝ (1, i)]
    state.positions[:, dof_a] = envelope * torch.cos(carrier_phase)
    state.positions[:, dof_b] = envelope * torch.sin(carrier_phase)

    # Velocities: v = Re{-i omega0 A(x) exp(i k0 x) p}
    # -i * exp(i φ) = -i * (cos φ + i sin φ) = sin φ - i cos φ
    # So: v_a = omega0 * A * sin(k0 x)
    #     v_b = -omega0 * A * cos(k0 x)
    state.velocities[:, dof_a] = omega0 * envelope * torch.sin(carrier_phase)
    state.velocities[:, dof_b] = -omega0 * envelope * torch.cos(carrier_phase)

    # Respect fixed boundary conditions if present
    if state.fixed_mask is not None:
        state.positions[state.fixed_mask, dof_a] = 0.0
        state.positions[state.fixed_mask, dof_b] = 0.0
        state.velocities[state.fixed_mask, dof_a] = 0.0
        state.velocities[state.fixed_mask, dof_b] = 0.0

    # Diagnostics
    lambda_c = 2 * np.pi / k0
    omega_c = omega0
    narrowband_ratio = sigma / lambda_c

    print("\n" + "="*70)
    print("PHOTON CIRCULAR PACKET INITIALIZATION (preparation-first)")
    print("="*70)
    print(f"Carrier wavelength λ_C:  {lambda_c:.6e} m")
    print(f"Carrier wave number k0:  {k0:.6e} rad/m")
    print(f"Carrier frequency ω0:    {omega0:.6e} rad/s")
    print(f"Envelope width σ:        {sigma:.6e} m  ({narrowband_ratio:.1f} λ_C)")
    print(f"Peak amplitude A:        {amplitude:.6e} m")
    print(f"Spectral width Δk/k0:    ~{1.0/(k0*sigma):.6e} (narrowband)")
    print(f"Polarization DOFs:       ({dof_a}, {dof_b})")
    print(f"Propagation axis:        {propagation_axis} ({'periodic' if periodic_axis else 'non-periodic'})")
    print(f"Max |position|:          {torch.abs(state.positions[:, dof_a]).max().item():.6e} m")
    print(f"Max |velocity|:          {torch.abs(state.velocities[:, dof_a]).max().item():.6e} m/s")

    if narrowband_ratio < 3:
        print("⚠ WARNING: Envelope width < 3λ_C may not be sufficiently narrowband")

    print("="*70)
    print("Verification targets:")
    print("  1. SVD of time series should show 2 dominant singular values")
    print("  2. Optional FFT cross-check: power concentrated near k0")
    print("  3. No post-hoc band-pass filtering required")
    print("="*70 + "\n")


def make_electron_double_loop_packet(
    state: BraneState,
    grid: BraneGrid,
    k0: float,
    omega0: float,
    amplitude: float,
    torus_major_radius: float,
    torus_minor_radius: float,
    twist_winding: int = 1,
    longitudinal_winding: int = 2,
    dof_pair: Tuple[int, int] = (2, 3),
    torus_axis: int = 2,
    torus_center: Optional[torch.Tensor] = None,
) -> None:
    """
    Initialize an electron-like double-loop tube packet with spinorial transport.

    Creates a narrowband carrier confined to a tubular neighborhood of a closed
    double-loop centerline with half-angle rotation of the internal polarization
    frame (spinorial 4π periodicity).

    Physical setup:
        - Toroidal centerline C(z) with major radius R and minor radius r0
        - Carrier: ξ(t,z,ρ,θ) = A f(ρ,θ;z) cos(ω0 t - m k0 z)
        - Double-loop ridge profile with twist α(z) = ℓ z/R
        - Half-angle polarization: p(z) ∝ (cos(α/2), sin(α/2))
        - One circuit: sign flip (W ≈ -I)
        - Two circuits: return (W² ≈ I)

    The ridge profile creates two narrow amplitude peaks at θ = α(z) and
    θ = α(z) + π that spiral around each other as z increases.

    Args:
        state: BraneState to initialize
        grid: BraneGrid (must be 3D for toroidal geometry)
        k0: Compton wave number (2π/λ_C)
        omega0: Compton frequency (m_e c²/ℏ)
        amplitude: Peak amplitude (typically λ_C/√π)
        torus_major_radius: Major radius R of the torus
        torus_minor_radius: Minor radius r0 (tube thickness)
        twist_winding: Twist winding number ℓ (default: 1, each ridge rotates once per loop)
        longitudinal_winding: Longitudinal winding m (default: 2, mode winds twice per loop)
        dof_pair: Tuple (q1, q2) for the internal C² polarization representation
        torus_axis: Which spatial axis the torus is aligned with (default: 2 = z-axis)
        torus_center: Center position of the torus (default: grid center)

    Notes:
        - Target signature: Wilson loop W ≈ -I after one circuit, W² ≈ I after two
        - The half-angle construction (α/2) encodes spinorial behavior
        - This is the scaffold described in Appendix B of the paper
        - Degeneracy is built-in via the double-loop ridge geometry

    Paper reference:
        Section "Preparation-first initialization of narrowband carriers"
        Paragraph "Electron packet (double-loop tube, spinorial transport)"
        Appendix B: "Tubular world-tube description"
    """
    device = state.device
    dtype = state.dtype
    ndim = len(grid.grid_shape)

    if ndim != 3:
        raise ValueError(f"Electron double-loop requires 3D grid, got {ndim}D")

    dof_q1, dof_q2 = dof_pair
    if dof_q1 == dof_q2:
        raise ValueError("DOF pair must contain two distinct components")

    # Get coordinates
    coords = grid.get_spatial_coordinates()  # [N, 3]

    # Default: center of grid
    if torus_center is None:
        torus_center = torch.tensor([
            grid.grid_shape[0] * grid.spacing / 2,
            grid.grid_shape[1] * grid.spacing / 2,
            grid.grid_shape[2] * grid.spacing / 2,
        ], device=device, dtype=dtype)
    else:
        torus_center = torch.as_tensor(torus_center, device=device, dtype=dtype)

    # Shift coordinates relative to torus center
    coords_centered = coords - torus_center

    # Torus parametrization:
    # We need to map each lattice point to (z, ρ, θ) tubular coordinates
    # where z is arclength along the centerline
    #
    # For a torus aligned with the z-axis (torus_axis=2):
    # Centerline C(z) = R (cos(z/R), sin(z/R), 0) in the xy-plane
    # z ∈ [0, 2πR)
    #
    # For simplicity, we'll work in Cartesian and map to toroidal coordinates

    # Identify which axes form the torus plane
    if torus_axis == 0:  # yz plane
        ax1, ax2, ax_axial = 1, 2, 0
    elif torus_axis == 1:  # xz plane
        ax1, ax2, ax_axial = 0, 2, 1
    else:  # xy plane (default, torus_axis==2)
        ax1, ax2, ax_axial = 0, 1, 2

    x1 = coords_centered[:, ax1]
    x2 = coords_centered[:, ax2]
    x_axial = coords_centered[:, ax_axial]

    # Cylindrical coordinates in the torus plane
    r_plane = torch.sqrt(x1**2 + x2**2)
    phi_plane = torch.atan2(x2, x1)  # Angle around torus axis

    # Arclength coordinate: z = R * phi_plane (mod 2πR)
    R = torus_major_radius
    z = R * phi_plane

    # Radial distance from the centerline (distance from ideal circle at radius R)
    # ρ = distance from centerline point C(z) to the lattice point
    # In the plane: distance from circle is |r_plane - R|
    # Out-of-plane: x_axial
    rho = torch.sqrt((r_plane - R)**2 + x_axial**2)

    # Angular coordinate θ around the tube cross-section
    # θ = 0 points radially outward from torus axis
    # θ increases counterclockwise when viewed from +z along the tube
    theta = torch.atan2(x_axial, r_plane - R)

    # Radial envelope: confine to tube of radius r0
    r0 = torus_minor_radius
    sigma_radial = r0 / 3.0  # Gaussian width
    G_rho = torch.exp(-rho**2 / (2 * sigma_radial**2))

    # Twist angle: α(z) = ℓ z / R
    ell = twist_winding
    alpha = ell * z / R

    # Double-loop ridge profile:
    # Two peaks at θ = α(z) and θ = α(z) + π
    # Angular width σ_θ
    sigma_theta = 0.5  # radians, adjust for narrowness
    ridge1 = torch.exp(-((theta - alpha) % (2*np.pi))**2 / (2 * sigma_theta**2))
    ridge2 = torch.exp(-((theta - (alpha + np.pi)) % (2*np.pi))**2 / (2 * sigma_theta**2))

    # Cross-section profile
    f_profile = G_rho * (ridge1 + ridge2)

    # Longitudinal phase: m k0 z where m = longitudinal_winding
    m = longitudinal_winding
    longitudinal_phase = m * k0 * z

    # Full amplitude field envelope
    envelope = amplitude * f_profile

    # Half-angle polarization for spinorial transport:
    # Internal frame rotates by α/2 (not α)
    # p(z) ∝ (cos(α/2), sin(α/2))
    # After 2π rotation of torus: α → α + 2π, so α/2 → α/2 + π → p → -p (sign flip)
    # After 4π rotation: α/2 → α/2 + 2π → p → +p (return)

    alpha_half = alpha / 2.0

    # Real-valued fields for the two DOFs:
    # q1 ~ A f(ρ,θ;z) cos(ω0 t - m k0 z) cos(α/2)
    # q2 ~ A f(ρ,θ;z) cos(ω0 t - m k0 z) sin(α/2)
    #
    # At t=0:
    cos_phase = torch.cos(longitudinal_phase)
    sin_phase = torch.sin(longitudinal_phase)

    state.positions[:, dof_q1] = envelope * cos_phase * torch.cos(alpha_half)
    state.positions[:, dof_q2] = envelope * cos_phase * torch.sin(alpha_half)

    # Velocities: ∂_t [envelope cos(ω0 t - m k0 z)] at t=0
    # = -ω0 envelope sin(-m k0 z) = ω0 envelope sin(m k0 z)
    state.velocities[:, dof_q1] = omega0 * envelope * sin_phase * torch.cos(alpha_half)
    state.velocities[:, dof_q2] = omega0 * envelope * sin_phase * torch.sin(alpha_half)

    # Respect fixed boundaries
    if state.fixed_mask is not None:
        state.positions[state.fixed_mask, dof_q1] = 0.0
        state.positions[state.fixed_mask, dof_q2] = 0.0
        state.velocities[state.fixed_mask, dof_q1] = 0.0
        state.velocities[state.fixed_mask, dof_q2] = 0.0

    # Diagnostics
    lambda_c = 2 * np.pi / k0
    loop_circumference = 2 * np.pi * R
    num_wavelengths = loop_circumference / lambda_c

    print("\n" + "="*70)
    print("ELECTRON DOUBLE-LOOP PACKET INITIALIZATION (preparation-first)")
    print("="*70)
    print(f"Compton wavelength λ_C:     {lambda_c:.6e} m")
    print(f"Compton wave number k0:     {k0:.6e} rad/m")
    print(f"Compton frequency ω0:       {omega0:.6e} rad/s")
    print(f"Torus major radius R:       {R:.6e} m")
    print(f"Torus minor radius r0:      {r0:.6e} m")
    print(f"Loop circumference:         {loop_circumference:.6e} m  ({num_wavelengths:.2f} λ_C)")
    print(f"Twist winding ℓ:            {ell}")
    print(f"Longitudinal winding m:     {m}")
    print(f"Peak amplitude A:           {amplitude:.6e} m")
    print(f"Polarization DOFs:          ({dof_q1}, {dof_q2})")
    print(f"Max |position|:             {torch.abs(state.positions[:, dof_q1]).max().item():.6e} m")
    print(f"Max |velocity|:             {torch.abs(state.velocities[:, dof_q1]).max().item():.6e} m/s")
    print("="*70)
    print("Spinorial transport structure:")
    print(f"  Half-angle rotation:      α/2 along loop")
    print(f"  After 1 circuit (2π):     W ≈ -I  (sign flip)")
    print(f"  After 2 circuits (4π):    W² ≈ I (return)")
    print("="*70)
    print("Verification targets:")
    print("  1. Wilson loop holonomy W ∈ U(2)")
    print("  2. Eigenvalues of W after one loop: λ ≈ -1")
    print("  3. Eigenvalues of W² after two loops: λ ≈ +1")
    print("  4. Degeneracy: SVD should show 2D polarization subspace")
    print("="*70 + "\n")