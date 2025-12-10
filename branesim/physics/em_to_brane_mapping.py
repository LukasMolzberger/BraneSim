"""
EM to Brane Mapping: Field-theoretic dictionary from EM fields to brane state.

This module implements the mathematical mapping from electromagnetic fields
(E, B) in physical units to brane state (amplitude and lateral velocities).

The mapping is based on:
1. Energy matching: EM energy density → brane amplitude via field energy density
2. Energy-momentum matching: EM Poynting vector → lateral velocity of brane points

This provides the "inverse dictionary" for initializing brane states from
classical EM field configurations using only field-theoretic quantities.
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
    K_eff: float,
    max_amplitude: float = None,
) -> torch.Tensor:
    """
    Map EM energy density to brane amplitude magnitude via energy matching.

    The mapping comes from equating energy densities using a field-theoretic
    calibration constant:
        u_brane = (1/2) K_eff A²
        u_EM = u_em

    Solving for amplitude:
        A = √(2 u_em / K_eff)

    The calibration constant K_eff encodes the brane material parameters
    (ρ_m, T, etc.) and is fixed once by Compton-scale calibration.

    Args:
        u_em: (N,) EM energy density [J/m³]
        K_eff: Effective calibration constant [J/m⁵] = [kg/m³·s²]
        max_amplitude: Optional maximum amplitude [m]; if provided, A is clamped

    Returns:
        A: (N,) Amplitude magnitude [m]
    """
    # Avoid division by zero
    K_safe = max(K_eff, 1e-60)

    # A = √(2u / K_eff)
    A = torch.sqrt(2.0 * u_em / K_safe)

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
    T_phys: float,
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
    Initialize brane positions and velocities from EM fields using field-theoretic mapping.

    This is the main entry point for the inverse EM→brane mapping.
    It performs the following steps:

    1. Compute effective mass density and calibration constant K_eff
    2. Calculate EM energy density and Poynting vector
    3. Map energy density to brane amplitude via energy matching (no phase/trig)
    4. Set normal displacement: A(r) → state.positions[:, field_component]
    5. Set normal velocity to zero: ∂_t ξ_⊥ = 0 (gauge choice)
    6. Map Poynting vector to lateral velocities → state.velocities[:, 0:3]

    Args:
        state: BraneState (positions and velocities modified in-place)
        grid: BraneGrid (for N and spacing)
        mapper: DimensionalMapper (for SI → sim unit conversion)
        m_point_phys: Point mass [kg]
        h_phys: Lattice spacing [m]
        T_phys: Brane tension [N]
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

    # 1) Compute effective mass density and calibration constant
    #    For 3D: rho_D = m_point / h³
    rho_mass = m_point_phys / (h_phys ** 3)

    # K_eff = ρ_m c² (from traveling wave energy density)
    # This encodes the brane material parameters fixed by calibration
    K_eff = rho_mass * (c_light ** 2)

    params = EMMaterialParams(
        epsilon_eff=epsilon_eff,
        mu_eff=mu_eff,
        rho_mass=rho_mass,
    )

    # 2) Compute EM energy density and Poynting vector
    u_em, S_phys = compute_em_energy_and_poynting(E_field_phys, B_field_phys, params)

    # 3) Map energy density to amplitude (field-theoretic, no oscillator phase)
    max_amp_phys = max_amplitude_fraction_of_h * h_phys
    A_phys = compute_brane_amplitude_from_em_energy(
        u_em=u_em,
        K_eff=K_eff,
        max_amplitude=max_amp_phys,
    )

    # Convert amplitude to simulation units
    A_sim = mapper.to_sim_length(A_phys)

    # 4) Write amplitude directly to 4th coordinate (no trig modulation)
    #    Position = A(r) is the normal displacement
    with torch.no_grad():
        state.positions[:, field_component] = A_sim

    # 5) Set normal velocity to zero (gauge choice: all energy initially in "potential")
    #    ∂_t ξ_⊥(x, 0) = 0
    with torch.no_grad():
        state.velocities[:, field_component] = 0.0

    # 6) Compute lateral velocities from Poynting flow (energy-momentum matching)
    #    Energy transport velocity: v_energy = S / u
    #    This is the local energy-flow velocity of the electromagnetic field
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
    print(f"\nEM → Brane Mapping (Field-Theoretic):")
    print(f"  Mass density ρ = {rho_mass:.6e} kg/m³")
    print(f"  Calibration constant K_eff = {K_eff:.6e} kg/(m³·s²)")
    print(f"  Max EM energy density: {u_em.max().item():.6e} J/m³")
    print(f"  Max amplitude (phys): {A_phys.max().item():.6e} m ({A_phys.max().item()/h_phys:.3f} × h)")
    print(f"  Max amplitude (sim): {A_sim.max().item():.6e}")
    print(f"  Normal velocity: 0.0 (gauge choice)")
    print(f"  Max Poynting magnitude: {torch.linalg.norm(S_phys, dim=-1).max().item():.6e} W/m²")
    print(f"  Max lateral velocity (phys): {torch.linalg.norm(v_energy_phys, dim=-1).max().item():.6e} m/s")
    print(f"  Max lateral velocity / c: {(torch.linalg.norm(v_energy_phys, dim=-1).max().item() / c_light):.6f}")