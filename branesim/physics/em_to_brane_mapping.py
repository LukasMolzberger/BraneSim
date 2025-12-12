"""
Electrostatic Brane-EM Mapping (Bidirectional)

This module implements the minimal, self-consistent electrostatic part of the
brane ↔ EM mapping:

FORWARD (brane → EM):
    Φ(x) = κ_EM * X^4(x)     (scalar potential from brane amplitude)
    E(x) = -∇Φ(x)            (electric field from potential gradient)
    ρ(x) = ε₀ ∇·E(x)         (charge density from Gauss's law)

INVERSE (EM → brane):
    Given E(x) or ρ(x), solve ∇²Φ = -∇·E (or -ρ/ε₀) via FFT
    Then X^4(x) = Φ(x) / κ_EM

Both directions use periodic boundary conditions and FFT-based Poisson solver.
"""

import torch
from typing import Tuple


class ElectrostaticMapping:
    """
    Minimal EM-to-brane mapping (electrostatic part only).

    Relates the brane amplitude field X^4(x) to electromagnetic quantities
    through the electrostatic potential:
        Φ(x)  = κ_EM * X^4(x)
        E(x)  = -∇Φ(x)
        ρ(x)  = ε₀ ∇·E(x)

    This implements the forward brane → EM direction using pure field theory.
    """

    def __init__(
        self,
        kappa_EM: float = 1.0,
        epsilon_0: float = 8.854187817e-12,
        dx: float = 1.0,
        device: torch.device = None,
        dtype: torch.dtype = torch.float64
    ):
        """
        Initialize the electrostatic mapping.

        Args:
            kappa_EM: Coupling constant relating X^4 to Φ [dimensionless or V/m]
            epsilon_0: Vacuum permittivity [F/m]
            dx: Spatial grid spacing for finite differences [m]
            device: torch device for computations
            dtype: torch dtype for computations
        """
        self.kappa_EM = kappa_EM
        self.epsilon_0 = epsilon_0
        self.dx = dx
        self.device = device if device is not None else torch.device('cpu')
        self.dtype = dtype

    def compute_potential(self, X4: torch.Tensor) -> torch.Tensor:
        """
        Compute the macroscopic electric potential field.

        Φ(x) = κ_EM * X^4(x)

        Args:
            X4: (nx, ny, nz) brane amplitude field (4th coordinate)

        Returns:
            Φ: (nx, ny, nz) electric potential [V]
        """
        return self.kappa_EM * X4

    def compute_electric_field(self, Phi: torch.Tensor) -> torch.Tensor:
        """
        Compute E = -∇Φ using central finite differences.

        Uses second-order accurate central differences in the interior and
        one-sided differences at boundaries.

        Args:
            Phi: (nx, ny, nz) electric potential [V]

        Returns:
            E: (nx, ny, nz, 3) electric field [V/m]
        """
        # Compute gradients along each axis
        # torch gradient uses central differences in interior, forward/backward at edges
        grad_x = torch.gradient(Phi, spacing=self.dx, dim=0)[0]
        grad_y = torch.gradient(Phi, spacing=self.dx, dim=1)[0]
        grad_z = torch.gradient(Phi, spacing=self.dx, dim=2)[0]

        # E = -∇Φ
        E = -torch.stack([grad_x, grad_y, grad_z], dim=-1)

        return E

    def compute_charge_density(self, E: torch.Tensor) -> torch.Tensor:
        """
        Compute charge density ρ = ε₀ ∇·E using Gauss's law.

        Args:
            E: (nx, ny, nz, 3) electric field [V/m]

        Returns:
            rho: (nx, ny, nz) charge density [C/m³]
        """
        # Extract components
        Ex = E[..., 0]
        Ey = E[..., 1]
        Ez = E[..., 2]

        # Compute divergence: ∇·E = ∂Ex/∂x + ∂Ey/∂y + ∂Ez/∂z
        div_E = (
            torch.gradient(Ex, spacing=self.dx, dim=0)[0] +
            torch.gradient(Ey, spacing=self.dx, dim=1)[0] +
            torch.gradient(Ez, spacing=self.dx, dim=2)[0]
        )

        # Gauss's law: ρ = ε₀ ∇·E
        rho = self.epsilon_0 * div_E

        return rho

    def map_from_brane(
        self,
        X4: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Given the brane amplitude X^4(x), compute the emergent electromagnetic quantities.

        This implements the complete forward electrostatic mapping:
            Φ = κ_EM * X^4
            E = -∇Φ
            ρ = ε₀ ∇·E

        Args:
            X4: (nx, ny, nz) brane amplitude field (4th coordinate)

        Returns:
            Phi: (nx, ny, nz) electric potential [V]
            E: (nx, ny, nz, 3) electric field [V/m]
            rho: (nx, ny, nz) charge density [C/m³]
        """
        # Step 1: Φ = κ_EM * X^4
        Phi = self.compute_potential(X4)

        # Step 2: E = -∇Φ
        E = self.compute_electric_field(Phi)

        # Step 3: ρ = ε₀ ∇·E
        rho = self.compute_charge_density(E)

        return Phi, E, rho

    def compute_electric_field_energy_density(self, E: torch.Tensor) -> torch.Tensor:
        """
        Compute the electric field energy density.

        u_E = (1/2) ε₀ |E|²

        Args:
            E: (nx, ny, nz, 3) electric field [V/m]

        Returns:
            u_E: (nx, ny, nz) energy density [J/m³]
        """
        E_squared = torch.sum(E**2, dim=-1)
        return 0.5 * self.epsilon_0 * E_squared


# ============================================================================
# INVERSE MAPPING (EM → brane) with Periodic Boundary Conditions
# ============================================================================

def _central_diff_periodic(field: torch.Tensor, dx: float, axis: int) -> torch.Tensor:
    """
    Compute central difference along axis with periodic boundary conditions.

    Uses second-order central differences: df/dx ≈ (f[i+1] - f[i-1]) / (2*dx)

    Args:
        field: Input tensor (nx, ny, nz) or (nx, ny, nz, 3)
        dx: Grid spacing
        axis: Axis along which to differentiate (0, 1, or 2)

    Returns:
        Derivative along specified axis (same shape as input)
    """
    # Roll forward and backward with periodic wrapping
    field_forward = torch.roll(field, shifts=-1, dims=axis)
    field_backward = torch.roll(field, shifts=1, dims=axis)

    return (field_forward - field_backward) / (2.0 * dx)


def divergence_periodic(
    vector_field: torch.Tensor,
    dx: float
) -> torch.Tensor:
    """
    Compute divergence ∇·V using periodic boundary conditions.

    Args:
        vector_field: (nx, ny, nz, 3) vector field
        dx: Grid spacing

    Returns:
        divergence: (nx, ny, nz) scalar field
    """
    Vx = vector_field[..., 0]
    Vy = vector_field[..., 1]
    Vz = vector_field[..., 2]

    dVx_dx = _central_diff_periodic(Vx, dx, axis=0)
    dVy_dy = _central_diff_periodic(Vy, dx, axis=1)
    dVz_dz = _central_diff_periodic(Vz, dx, axis=2)

    return dVx_dx + dVy_dy + dVz_dz


def gradient_periodic(
    scalar_field: torch.Tensor,
    dx: float
) -> torch.Tensor:
    """
    Compute gradient ∇f using periodic boundary conditions.

    Args:
        scalar_field: (nx, ny, nz) scalar field
        dx: Grid spacing

    Returns:
        gradient: (nx, ny, nz, 3) vector field
    """
    grad_x = _central_diff_periodic(scalar_field, dx, axis=0)
    grad_y = _central_diff_periodic(scalar_field, dx, axis=1)
    grad_z = _central_diff_periodic(scalar_field, dx, axis=2)

    return torch.stack([grad_x, grad_y, grad_z], dim=-1)


def solve_poisson_periodic_fft(
    rhs: torch.Tensor,
    dx: float,
    device: torch.device = None,
    dtype: torch.dtype = torch.float64
) -> torch.Tensor:
    """
    Solve Poisson equation ∇²φ = rhs with periodic boundary conditions using FFT.

    The solution is obtained by transforming to Fourier space where:
        -k² φ̂ = F[rhs]
        φ̂ = -F[rhs] / k²

    The zero mode (k=0) is set to zero, fixing φ up to an additive constant.

    Args:
        rhs: (nx, ny, nz) right-hand side of Poisson equation
        dx: Grid spacing
        device: torch device for computations
        dtype: torch dtype for computations

    Returns:
        phi: (nx, ny, nz) solution to Poisson equation
    """
    if device is None:
        device = rhs.device
    if dtype is None:
        dtype = rhs.dtype

    nx, ny, nz = rhs.shape

    # FFT of right-hand side
    rhs_hat = torch.fft.fftn(rhs, dim=(0, 1, 2))

    # Build k² in Fourier space
    # Frequencies: k_i = 2π n_i / (N_i * dx) for n_i ∈ [0, N_i-1]
    # with wrapping for negative frequencies
    kx = 2.0 * torch.pi * torch.fft.fftfreq(nx, d=dx, device=device, dtype=dtype)
    ky = 2.0 * torch.pi * torch.fft.fftfreq(ny, d=dx, device=device, dtype=dtype)
    kz = 2.0 * torch.pi * torch.fft.fftfreq(nz, d=dx, device=device, dtype=dtype)

    # Meshgrid for 3D
    KX, KY, KZ = torch.meshgrid(kx, ky, kz, indexing='ij')
    k_squared = KX**2 + KY**2 + KZ**2

    # Avoid division by zero at k=0 (set to 1, will be zeroed out anyway)
    k_squared[0, 0, 0] = 1.0

    # Solve in Fourier space: φ̂ = -F[rhs] / k²
    phi_hat = -rhs_hat / k_squared

    # Zero out the mean (k=0 mode)
    phi_hat[0, 0, 0] = 0.0

    # Inverse FFT to get solution
    phi = torch.fft.ifftn(phi_hat, dim=(0, 1, 2)).real

    return phi


def potential_from_E_periodic(
    E_field: torch.Tensor,
    dx: float,
    device: torch.device = None,
    dtype: torch.dtype = torch.float64
) -> torch.Tensor:
    """
    Reconstruct potential Φ from electric field E using ∇²Φ = -∇·E.

    This uses the electrostatic relation E = -∇Φ, which implies that
    the Laplacian of Φ is minus the divergence of E.

    Args:
        E_field: (nx, ny, nz, 3) electric field [V/m]
        dx: Grid spacing [m]
        device: torch device for computations
        dtype: torch dtype for computations

    Returns:
        Phi: (nx, ny, nz) electric potential [V]
    """
    # Compute divergence of E field
    div_E = divergence_periodic(E_field, dx)

    # Solve ∇²Φ = -∇·E
    Phi = solve_poisson_periodic_fft(-div_E, dx, device=device, dtype=dtype)

    return Phi


def potential_from_rho_periodic(
    rho: torch.Tensor,
    epsilon_0: float,
    dx: float,
    device: torch.device = None,
    dtype: torch.dtype = torch.float64
) -> torch.Tensor:
    """
    Reconstruct potential Φ from charge density ρ using ∇²Φ = -ρ/ε₀.

    This is Poisson's equation for the electrostatic potential.

    Args:
        rho: (nx, ny, nz) charge density [C/m³]
        epsilon_0: Vacuum permittivity [F/m]
        dx: Grid spacing [m]
        device: torch device for computations
        dtype: torch dtype for computations

    Returns:
        Phi: (nx, ny, nz) electric potential [V]
    """
    # Solve ∇²Φ = -ρ/ε₀
    Phi = solve_poisson_periodic_fft(-rho / epsilon_0, dx, device=device, dtype=dtype)

    return Phi


def initialize_brane_from_electrostatics(
    state,
    kappa_EM: float,
    dx: float,
    epsilon_0: float = 8.854187817e-12,
    E_field: torch.Tensor = None,
    rho: torch.Tensor = None,
    Phi: torch.Tensor = None,
    field_component: int = 3,
    device: torch.device = None,
    dtype: torch.dtype = torch.float64
):
    """
    Initialize brane state from electrostatic quantities: E, ρ, or Φ.

    The inverse electrostatic mapping is:
        1. If E is given: Solve ∇²Φ = -∇·E
        2. If ρ is given: Solve ∇²Φ = -ρ/ε₀
        3. If Φ is given: Use directly
        4. Set X^4 = Φ / κ_EM

    Args:
        state: BraneState object to initialize
        kappa_EM: Coupling constant relating X^4 to Φ
        dx: Grid spacing [m]
        epsilon_0: Vacuum permittivity [F/m]
        E_field: (nx, ny, nz, 3) electric field [V/m] (optional)
        rho: (nx, ny, nz) charge density [C/m³] (optional)
        Phi: (nx, ny, nz) electric potential [V] (optional)
        field_component: Which component of state.positions to write to (default: 3 for X^4)
        device: torch device for computations
        dtype: torch dtype for computations

    Raises:
        ValueError: If none of E_field, rho, or Phi are provided
    """
    if device is None:
        device = state.positions.device
    if dtype is None:
        dtype = state.positions.dtype

    # Determine which input was provided and compute Φ
    if Phi is not None:
        # Direct input of potential
        pass
    elif E_field is not None:
        # Reconstruct Φ from E via ∇²Φ = -∇·E
        Phi = potential_from_E_periodic(E_field, dx, device=device, dtype=dtype)
    elif rho is not None:
        # Reconstruct Φ from ρ via ∇²Φ = -ρ/ε₀
        Phi = potential_from_rho_periodic(rho, epsilon_0, dx, device=device, dtype=dtype)
    else:
        raise ValueError("Must provide at least one of: E_field, rho, or Phi")

    # Inverse electrostatic mapping: X^4 = Φ / κ_EM
    X4 = Phi / kappa_EM

    # Write to brane state
    state.positions[:, field_component] = X4.flatten()

    # Zero out velocities in the normal direction (gauge choice)
    state.velocities[:, field_component] = 0.0


def initialize_brane_from_em_fields(
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
):
    """
    Initialize brane state from electromagnetic fields (backward compatibility wrapper).

    This function provides backward compatibility with the old interface by wrapping
    the new simplified electrostatic mapping and adding energy-momentum matching for
    lateral velocities.

    Args:
        state: BraneState object to initialize
        grid: BraneGrid object
        mapper: DimensionalMapper for unit conversions
        m_point_phys: Point mass [kg]
        h_phys: Lattice spacing [m]
        T_phys: Brane tension [N]
        E_field_phys: (nx, ny, nz, 3) electric field [V/m]
        B_field_phys: (nx, ny, nz, 3) magnetic field [T]
        epsilon_eff: Effective permittivity [F/m]
        mu_eff: Effective permeability [H/m]
        c_light: Speed of light [m/s]
        field_component: Which component of state.positions to write to (default: 3 for X^4)
        max_amplitude_fraction_of_h: Maximum amplitude as fraction of h_phys
        velocity_clip_to_c: Whether to clip velocities to speed of light
    """
    device = state.positions.device
    dtype = state.positions.dtype

    # Get grid shape
    nx, ny, nz = grid.grid_shape

    # Reshape E and B fields to 3D if flattened
    if E_field_phys.ndim == 2:
        E_field_phys = E_field_phys.reshape(nx, ny, nz, 3)
    if B_field_phys.ndim == 2:
        B_field_phys = B_field_phys.reshape(nx, ny, nz, 3)

    # Compute electromagnetic energy density
    # u_EM = (1/2)(ε|E|² + |B|²/μ)
    E_squared = torch.sum(E_field_phys**2, dim=-1)
    B_squared = torch.sum(B_field_phys**2, dim=-1)
    u_EM = 0.5 * (epsilon_eff * E_squared + B_squared / mu_eff)

    # Compute Poynting vector S = E × B / μ
    S_x = (E_field_phys[..., 1] * B_field_phys[..., 2] -
           E_field_phys[..., 2] * B_field_phys[..., 1]) / mu_eff
    S_y = (E_field_phys[..., 2] * B_field_phys[..., 0] -
           E_field_phys[..., 0] * B_field_phys[..., 2]) / mu_eff
    S_z = (E_field_phys[..., 0] * B_field_phys[..., 1] -
           E_field_phys[..., 1] * B_field_phys[..., 0]) / mu_eff
    S_mag = torch.sqrt(S_x**2 + S_y**2 + S_z**2)

    # Compute lateral velocities from Poynting vector
    # v_lateral = S / u_EM (where u_EM > 0)
    v_lateral_x = torch.zeros_like(S_x)
    v_lateral_y = torch.zeros_like(S_y)
    v_lateral_z = torch.zeros_like(S_z)

    mask = u_EM > 0
    v_lateral_x[mask] = S_x[mask] / u_EM[mask]
    v_lateral_y[mask] = S_y[mask] / u_EM[mask]
    v_lateral_z[mask] = S_z[mask] / u_EM[mask]

    # Clip velocities to c if requested
    if velocity_clip_to_c:
        v_mag = torch.sqrt(v_lateral_x**2 + v_lateral_y**2 + v_lateral_z**2)
        clip_mask = v_mag > c_light
        if clip_mask.any():
            scale = c_light / v_mag[clip_mask]
            v_lateral_x[clip_mask] *= scale
            v_lateral_y[clip_mask] *= scale
            v_lateral_z[clip_mask] *= scale

    # Convert lateral velocities to simulation units
    v_lateral_x_sim = mapper.to_sim_velocity(v_lateral_x)
    v_lateral_y_sim = mapper.to_sim_velocity(v_lateral_y)
    v_lateral_z_sim = mapper.to_sim_velocity(v_lateral_z)

    # Set lateral velocities in state
    state.velocities[:, 0] = v_lateral_x_sim.flatten()
    state.velocities[:, 1] = v_lateral_y_sim.flatten()
    state.velocities[:, 2] = v_lateral_z_sim.flatten()

    # Compute normal amplitude from energy density
    # From energy matching: u_EM = (1/2) K_eff * A²
    # where K_eff depends on brane material parameters
    # For traveling wave: K_eff ≈ ρ_m * c²
    rho_m = m_point_phys / (h_phys ** 3)
    K_eff = rho_m * c_light**2

    # Amplitude in physical units
    A_phys = torch.sqrt(2.0 * u_EM / K_eff)

    # Apply amplitude clipping
    A_max = max_amplitude_fraction_of_h * h_phys
    A_phys = torch.clamp(A_phys, max=A_max)

    # Convert to simulation units
    A_sim = mapper.to_sim_length(A_phys)

    # Set normal displacement (X^4 component)
    state.positions[:, field_component] = A_sim.flatten()

    # Zero out normal velocity (gauge choice: ∂_t X^4 = 0 at t=0)
    state.velocities[:, field_component] = 0.0


# Example usage function
def example_usage():
    """
    Example demonstrating how to use the ElectrostaticMapping class.

    Shows both forward (brane → EM) and inverse (EM → brane) mappings.
    """
    # Create a synthetic brane amplitude field (64³ grid)
    nx, ny, nz = 64, 64, 64
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dtype = torch.float64

    # Synthetic X^4 field (small random amplitudes)
    X4_original = torch.randn(nx, ny, nz, device=device, dtype=dtype) * 0.01

    # Initialize the mapping
    dx = 1e-15  # 1 fm grid spacing
    kappa_EM = 1.0
    epsilon_0 = 8.854187817e-12

    mapper = ElectrostaticMapping(
        kappa_EM=kappa_EM,
        epsilon_0=epsilon_0,
        dx=dx,
        device=device,
        dtype=dtype
    )

    print("=" * 70)
    print("FORWARD MAPPING (brane → EM)")
    print("=" * 70)

    # Compute emergent EM fields
    Phi, E, rho = mapper.map_from_brane(X4_original)

    # Compute energy density
    u_E = mapper.compute_electric_field_energy_density(E)

    print(f"Potential range: [{Phi.min():.6e}, {Phi.max():.6e}] V")
    print(f"E-field magnitude: [{E.norm(dim=-1).min():.6e}, {E.norm(dim=-1).max():.6e}] V/m")
    print(f"Charge density range: [{rho.min():.6e}, {rho.max():.6e}] C/m³")
    print(f"Energy density: [{u_E.min():.6e}, {u_E.max():.6e}] J/m³")

    print("\n" + "=" * 70)
    print("INVERSE MAPPING (EM → brane)")
    print("=" * 70)

    # Test inverse mapping from E field
    Phi_reconstructed = potential_from_E_periodic(E, dx, device=device, dtype=dtype)
    X4_reconstructed = Phi_reconstructed / kappa_EM

    # Compare original and reconstructed X^4
    diff = X4_reconstructed - X4_original
    rel_error = torch.abs(diff).mean() / torch.abs(X4_original).mean()

    print(f"X^4 reconstruction from E field:")
    print(f"  Mean absolute difference: {torch.abs(diff).mean():.6e}")
    print(f"  Relative error: {rel_error:.6e}")
    print(f"  Max absolute difference: {torch.abs(diff).max():.6e}")

    # Test inverse mapping from charge density
    Phi_from_rho = potential_from_rho_periodic(rho, epsilon_0, dx, device=device, dtype=dtype)
    X4_from_rho = Phi_from_rho / kappa_EM

    diff_rho = X4_from_rho - X4_original
    rel_error_rho = torch.abs(diff_rho).mean() / torch.abs(X4_original).mean()

    print(f"\nX^4 reconstruction from ρ field:")
    print(f"  Mean absolute difference: {torch.abs(diff_rho).mean():.6e}")
    print(f"  Relative error: {rel_error_rho:.6e}")
    print(f"  Max absolute difference: {torch.abs(diff_rho).max():.6e}")

    print("\n" + "=" * 70)
    print("Round-trip mapping successful!")
    print("=" * 70)

    return Phi, E, rho, u_E, X4_reconstructed


if __name__ == "__main__":
    # Run example if executed as script
    example_usage()