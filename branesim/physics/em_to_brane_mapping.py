"""
EM to Brane Mapping: Inverted dictionary from EM fields to brane state.

This module implements the mathematical mapping from electromagnetic fields
(E, B) in physical units to brane state (amplitude and lateral velocities).

The mapping is based on:
1. Energy matching: EM energy density → brane amplitude via oscillation energy
2. Poynting flow: EM energy flux → lateral velocity of brane points

This provides the "inverted dictionary" for initializing brane states from
classical EM field configurations.
"""

import torch
from dataclasses import dataclass
from typing import Tuple


@dataclass
class EMMaterialParams:
    """
    Effective EM material parameters for the brane-EM correspondence.

    These represent the effective permittivity, permeability, and mass density
    that relate EM fields to brane dynamics.
    """
    epsilon_eff: float  # Effective permittivity (≈ ε₀) [F/m]
    mu_eff: float       # Effective permeability (≈ μ₀) [H/m]
    rho_mass: float     # Effective mass density of brane continuum [kg/m³]


def compute_em_energy_and_poynting(
    E: torch.Tensor,
    B: torch.Tensor,
    params: EMMaterialParams,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Compute EM energy density and Poynting vector from field vectors.

    Energy density:
        u = (1/2)(εE² + B²/μ)

    Poynting vector:
        S = (1/μ) E × B

    Args:
        E: (N, 3) electric field in SI units [V/m]
        B: (N, 3) magnetic field in SI units [T]
        params: EMMaterialParams with ε, μ constants

    Returns:
        u: (N,) EM energy density [J/m³]
        S: (N, 3) Poynting vector [W/m²]
    """
    epsilon = params.epsilon_eff
    mu = params.mu_eff

    # Compute squared magnitudes
    E2 = (E ** 2).sum(dim=-1)  # |E|²
    B2 = (B ** 2).sum(dim=-1)  # |B|²

    # Energy density: u = (1/2)(εE² + B²/μ)
    u = 0.5 * (epsilon * E2 + B2 / mu)

    # Poynting vector: S = (1/μ) E × B
    S = torch.cross(E, B, dim=-1) / mu

    return u, S


def compute_brane_amplitude_from_em_energy(
    u_em: torch.Tensor,
    omega_phys: float,
    rho_mass: float,
    max_amplitude: float = None,
) -> torch.Tensor:
    """
    Map EM energy density to brane amplitude magnitude via energy matching.

    The mapping comes from equating time-averaged energies:
        ⟨u_brane⟩ = (1/2) ρ ω² A²
        ⟨u_EM⟩ = u_em

    Solving for amplitude:
        A = √(2 u_em / (ρ ω²))

    Args:
        u_em: (N,) EM energy density [J/m³]
        omega_phys: Physical angular frequency [rad/s]
        rho_mass: Mass density of brane [kg/m³]
        max_amplitude: Optional maximum amplitude [m]; if provided, A is clamped

    Returns:
        A: (N,) Amplitude magnitude [m]
    """
    # Avoid division by zero
    denom = max(rho_mass * (omega_phys ** 2), 1e-60)

    # A = √(2u / (ρω²))
    A = torch.sqrt(2.0 * u_em / denom)

    # Optional clamping to maximum amplitude
    if max_amplitude is not None:
        A = torch.clamp(A, max=max_amplitude)

    return A


def initialize_brane_from_em_fields(
    *,
    state,
    grid,
    mapper,
    m_point_phys: float,
    h_phys: float,
    omega_phys: float,
    E_field_phys: torch.Tensor,
    B_field_phys: torch.Tensor,
    epsilon_eff: float,
    mu_eff: float,
    c_light: float,
    field_component: int = 3,
    max_amplitude_fraction_of_h: float = 0.1,
    velocity_clip_to_c: bool = True,
) -> None:
    """
    Initialize brane positions and velocities from EM fields.

    This is the main entry point for the inverted EM→brane mapping.
    It performs the following steps:

    1. Compute effective mass density from point mass and lattice spacing
    2. Calculate EM energy density and Poynting vector
    3. Map energy density to brane amplitude via energy matching
    4. Map Poynting vector to lateral velocities
    5. Write amplitude to state.positions[:, field_component]
    6. Write lateral velocities to state.velocities[:, 0:3]

    Args:
        state: BraneState (positions and velocities modified in-place)
        grid: BraneGrid (for N and spacing)
        mapper: DimensionalMapper (for SI → sim unit conversion)
        m_point_phys: Point mass [kg]
        h_phys: Lattice spacing [m]
        omega_phys: Photon angular frequency [rad/s]
        E_field_phys: (N, 3) E-field in SI units at t=0 [V/m]
        B_field_phys: (N, 3) B-field in SI units at t=0 [T]
        epsilon_eff: Effective permittivity (≈ ε₀) [F/m]
        mu_eff: Effective permeability (≈ μ₀) [H/m]
        c_light: Speed of light [m/s]
        field_component: Index of amplitude dimension in state.positions (default: 3)
        max_amplitude_fraction_of_h: Clamp amplitude to this fraction of h_phys
        velocity_clip_to_c: If True, clip |v| ≤ c
    """

    # Validate input dimensions
    N = grid.num_points
    assert E_field_phys.shape == (N, 3), f"E field shape {E_field_phys.shape} != ({N}, 3)"
    assert B_field_phys.shape == (N, 3), f"B field shape {B_field_phys.shape} != ({N}, 3)"

    # 1) Compute effective mass density
    #    For 3D: rho_D = m_point / h³
    rho_mass = m_point_phys / (h_phys ** 3)

    params = EMMaterialParams(
        epsilon_eff=epsilon_eff,
        mu_eff=mu_eff,
        rho_mass=rho_mass,
    )

    # 2) Compute EM energy density and Poynting vector
    u_em, S_phys = compute_em_energy_and_poynting(E_field_phys, B_field_phys, params)

    # 3) Map energy density to amplitude
    max_amp_phys = max_amplitude_fraction_of_h * h_phys
    A_phys = compute_brane_amplitude_from_em_energy(
        u_em=u_em,
        omega_phys=omega_phys,
        rho_mass=rho_mass,
        max_amplitude=max_amp_phys,
    )

    # Convert amplitude to simulation units
    A_sim = mapper.to_sim_length(A_phys)

    # 4) Extract spatial phase from EM fields (for wave structure)
    #    For circularly/linearly polarized waves, extract phase from field direction
    #    Use the first non-zero transverse component to get spatial phase
    E_magnitudes = torch.linalg.norm(E_field_phys, dim=-1)
    max_E_mag = E_magnitudes.max()

    if max_E_mag > 1e-10:  # If fields are non-trivial
        # Find dominant transverse components (not aligned with propagation)
        # Assume propagation along x (could be generalized)
        # For circular polarization: E_y = E₀ cos(kx), E_z = E₀ sin(kx)
        # Extract phase: φ = atan2(E_z, E_y) gives kx modulo 2π
        phase = torch.atan2(E_field_phys[:, 2], E_field_phys[:, 1] + 1e-30)
        spatial_modulation = torch.cos(phase)
    else:
        spatial_modulation = torch.ones_like(A_sim)

    # 5) Write amplitude with spatial phase modulation to 4th coordinate
    #    Position = A(r) * cos(φ(x)) captures both envelope and wave structure
    with torch.no_grad():
        state.positions[:, field_component] = A_sim * spatial_modulation

    # 5) Compute lateral velocities from Poynting flow
    #    Energy transport velocity: v_energy = S / u
    #    For an EM wave, this gives the direction and speed of energy flow
    eps = 1e-30  # Regularization to avoid division by zero
    u_expanded = (u_em + eps).unsqueeze(-1)  # (N, 1)
    v_energy_phys = S_phys / u_expanded      # (N, 3) [m/s]

    # Optional: clip velocity magnitude to speed of light
    if velocity_clip_to_c:
        v_norm = torch.linalg.norm(v_energy_phys, dim=-1, keepdim=True)
        factor = torch.clamp(v_norm, max=c_light) / (v_norm + eps)
        v_energy_phys = v_energy_phys * factor

    # Convert to simulation units
    v_energy_sim = mapper.to_sim_velocity(v_energy_phys)

    # Write lateral velocities (x, y, z components)
    with torch.no_grad():
        state.velocities[:, 0:3] = v_energy_sim

    # Print diagnostic info
    print(f"\nEM → Brane Mapping:")
    print(f"  Mass density ρ = {rho_mass:.6e} kg/m³")
    print(f"  Angular frequency ω = {omega_phys:.6e} rad/s")
    print(f"  Max EM energy density: {u_em.max().item():.6e} J/m³")
    print(f"  Max amplitude (phys): {A_phys.max().item():.6e} m ({A_phys.max().item()/h_phys:.3f} × h)")
    print(f"  Max amplitude (sim): {A_sim.max().item():.6e}")
    print(f"  Max Poynting magnitude: {torch.linalg.norm(S_phys, dim=-1).max().item():.6e} W/m²")
    print(f"  Max lateral velocity (phys): {torch.linalg.norm(v_energy_phys, dim=-1).max().item():.6e} m/s")
    print(f"  Max lateral velocity / c: {(torch.linalg.norm(v_energy_phys, dim=-1).max().item() / c_light):.6f}")