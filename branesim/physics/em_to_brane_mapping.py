"""
Electrostatic Brane-EM Mapping with Mixed Boundary Conditions

This module implements the minimal, self-consistent electrostatic part of the
brane ↔ EM mapping with support for mixed boundary conditions (clamped/periodic):

FORWARD (brane → EM):
    Φ(x) = κ_EM * X^4(x)     (scalar potential from brane amplitude)
    E(x) = -∇Φ(x)            (electric field from potential gradient)
    ρ(x) = ε₀ ∇·E(x)         (charge density from Gauss's law)

INVERSE (EM → brane):
    Given E(x), solve via optimization (CG or L-BFGS):
        min_{X^4} (1/2) |κ_EM ∇X^4 + E|² + (λ/2)|X^4|²
    Then X^4(x) = Φ(x) / κ_EM

Boundary Conditions:
    - "periodic": Wrap indices across boundaries (periodic repetition)
    - "dirichlet0": Treat values outside domain as 0 (clamped boundaries)
    - Can be mixed per axis: e.g., ("periodic", "dirichlet0", "dirichlet0")

No FFT is used. Inverse mapping uses iterative solvers (CG/L-BFGS).
"""

import torch
from typing import Tuple, Literal

# Type alias for boundary condition specification
# Order: (x, y, z) corresponding to tensor dims (0, 1, 2)
BC = Tuple[Literal["periodic", "dirichlet0"],
          Literal["periodic", "dirichlet0"],
          Literal["periodic", "dirichlet0"]]


# ============================================================================
# Boundary-Aware Shift Operator
# ============================================================================

def shift_mixed(
    f: torch.Tensor,
    dim: int,
    step: int,
    bc: BC
) -> torch.Tensor:
    """
    Shift tensor along dimension with boundary condition handling.

    Args:
        f: Input tensor (nx, ny, nz) or (nx, ny, nz, 3)
        dim: Dimension to shift (0, 1, or 2)
        step: Shift amount (+1 for forward, -1 for backward)
        bc: Boundary conditions for each axis

    Returns:
        Shifted tensor with boundary conditions applied
    """
    bc_type = bc[dim]

    if bc_type == "periodic":
        # Periodic: wrap around
        # torch.roll with positive shifts moves data "to the right" (later indices)
        # We want: step=+1 to get f[i-1], step=-1 to get f[i+1]
        # torch.roll(f, shifts=1) at index i gives f[i-1]
        # torch.roll(f, shifts=-1) at index i gives f[i+1]
        # So we use shifts=step (not -step)
        return torch.roll(f, shifts=step, dims=dim)

    elif bc_type == "dirichlet0":
        # Dirichlet with zero outside: shift and fill with zeros
        out = torch.zeros_like(f)

        if step == 1:
            # Forward shift: out[1:] = f[:-1], out[0] = 0
            if dim == 0:
                out[1:, ...] = f[:-1, ...]
            elif dim == 1:
                out[:, 1:, ...] = f[:, :-1, ...]
            elif dim == 2:
                out[:, :, 1:, ...] = f[:, :, :-1, ...]

        elif step == -1:
            # Backward shift: out[:-1] = f[1:], out[-1] = 0
            if dim == 0:
                out[:-1, ...] = f[1:, ...]
            elif dim == 1:
                out[:, :-1, ...] = f[:, 1:, ...]
            elif dim == 2:
                out[:, :, :-1, ...] = f[:, :, 1:, ...]
        else:
            raise ValueError(f"Only step=±1 supported, got {step}")

        return out

    else:
        raise ValueError(f"Unknown boundary condition: {bc_type}")


# ============================================================================
# Mixed Finite-Difference Operators
# ============================================================================

def grad_mixed(
    phi: torch.Tensor,
    h: float,
    bc: BC
) -> torch.Tensor:
    """
    Compute gradient ∇φ with mixed boundary conditions.

    Uses second-order central differences:
        ∂φ/∂x ≈ (φ[i+1] - φ[i-1]) / (2h)

    Args:
        phi: Scalar field (nx, ny, nz)
        h: Grid spacing [m]
        bc: Boundary conditions per axis

    Returns:
        Gradient field (nx, ny, nz, 3)
    """
    # Compute partial derivatives using central differences
    grad_x = (shift_mixed(phi, 0, -1, bc) - shift_mixed(phi, 0, 1, bc)) / (2.0 * h)
    grad_y = (shift_mixed(phi, 1, -1, bc) - shift_mixed(phi, 1, 1, bc)) / (2.0 * h)
    grad_z = (shift_mixed(phi, 2, -1, bc) - shift_mixed(phi, 2, 1, bc)) / (2.0 * h)

    return torch.stack([grad_x, grad_y, grad_z], dim=-1)


def div_mixed(
    E: torch.Tensor,
    h: float,
    bc: BC
) -> torch.Tensor:
    """
    Compute divergence ∇·E with mixed boundary conditions.

    Args:
        E: Vector field (nx, ny, nz, 3)
        h: Grid spacing [m]
        bc: Boundary conditions per axis

    Returns:
        Divergence field (nx, ny, nz)
    """
    Ex, Ey, Ez = E[..., 0], E[..., 1], E[..., 2]

    # Compute divergence: ∇·E = ∂Ex/∂x + ∂Ey/∂y + ∂Ez/∂z
    dEx_dx = (shift_mixed(Ex, 0, -1, bc) - shift_mixed(Ex, 0, 1, bc)) / (2.0 * h)
    dEy_dy = (shift_mixed(Ey, 1, -1, bc) - shift_mixed(Ey, 1, 1, bc)) / (2.0 * h)
    dEz_dz = (shift_mixed(Ez, 2, -1, bc) - shift_mixed(Ez, 2, 1, bc)) / (2.0 * h)

    return dEx_dx + dEy_dy + dEz_dz


def laplacian_mixed(
    phi: torch.Tensor,
    h: float,
    bc: BC
) -> torch.Tensor:
    """
    Compute Laplacian ∇²φ with mixed boundary conditions.

    Uses 7-point stencil:
        ∇²φ ≈ (φ[i+1] - 2φ[i] + φ[i-1])/h² + ... (for each dimension)

    Args:
        phi: Scalar field (nx, ny, nz)
        h: Grid spacing [m]
        bc: Boundary conditions per axis

    Returns:
        Laplacian field (nx, ny, nz)
    """
    laplacian = torch.zeros_like(phi)

    for dim in range(3):
        phi_forward = shift_mixed(phi, dim, -1, bc)
        phi_backward = shift_mixed(phi, dim, 1, bc)
        laplacian += (phi_forward - 2.0 * phi + phi_backward) / (h * h)

    return laplacian


# ============================================================================
# Reshape Helpers (Flat ↔ Grid)
# ============================================================================

def reshape_flat_to_grid_scalar(
    f_flat: torch.Tensor,
    nx: int,
    ny: int,
    nz: int
) -> torch.Tensor:
    """
    Reshape flat scalar field to 3D grid.

    Args:
        f_flat: Flat field (N,) where N = nx*ny*nz
        nx, ny, nz: Grid dimensions

    Returns:
        Grid field (nx, ny, nz)
    """
    return f_flat.reshape(nx, ny, nz)


def reshape_grid_to_flat_scalar(f_xyz: torch.Tensor) -> torch.Tensor:
    """
    Reshape 3D grid scalar field to flat.

    Args:
        f_xyz: Grid field (nx, ny, nz)

    Returns:
        Flat field (N,)
    """
    return f_xyz.flatten()


def reshape_flat_to_grid_vec(
    v_flat: torch.Tensor,
    nx: int,
    ny: int,
    nz: int
) -> torch.Tensor:
    """
    Reshape flat vector field to 3D grid.

    Args:
        v_flat: Flat vector field (N, 3) where N = nx*ny*nz
        nx, ny, nz: Grid dimensions

    Returns:
        Grid vector field (nx, ny, nz, 3)
    """
    return v_flat.reshape(nx, ny, nz, 3)


def reshape_grid_to_flat_vec(v_xyz: torch.Tensor) -> torch.Tensor:
    """
    Reshape 3D grid vector field to flat.

    Args:
        v_xyz: Grid vector field (nx, ny, nz, 3)

    Returns:
        Flat vector field (N, 3)
    """
    nx, ny, nz = v_xyz.shape[:3]
    return v_xyz.reshape(nx * ny * nz, 3)


# ============================================================================
# Forward Electrostatic Mapping (Brane → EM)
# ============================================================================

def brane_X4_to_E_field(
    state,
    grid,
    mapper,
    kappa_EM: float = 1.0,
    bc: BC = ("dirichlet0", "dirichlet0", "dirichlet0"),
    return_potential: bool = False,
    return_charge_density: bool = False,
    epsilon_0: float = 8.854187817e-12
) -> Tuple[torch.Tensor, ...]:
    """
    Compute electric field E from brane amplitude X^4.

    Forward electrostatic mapping:
        Φ(x) = κ_EM * X^4(x)
        E(x) = -∇Φ(x) = -κ_EM * ∇X^4(x)
        ρ(x) = ε₀ ∇·E(x)  (optional)

    Args:
        state: BraneState with positions[:, 3] = X^4 in simulation units
        grid: BraneGrid with grid_shape and spacing
        mapper: DimensionalMapper for unit conversion
        kappa_EM: Coupling constant [V/m] relating X^4 to Φ
        bc: Boundary conditions per axis (x, y, z)
        return_potential: If True, also return potential Φ
        return_charge_density: If True, also return charge density ρ
        epsilon_0: Vacuum permittivity [F/m]

    Returns:
        E_phys_flat: Electric field (N, 3) in physical units [V/m]
        If return_potential=True: also returns Phi_phys_flat (N,) [V]
        If return_charge_density=True: also returns rho_phys_flat (N,) [C/m³]
    """
    # Extract X^4 in simulation units
    X4_sim_flat = state.positions[:, 3]

    # Convert to physical meters
    X4_phys_flat = mapper.to_phys_length(X4_sim_flat)

    # Get grid shape and spacing
    nx, ny, nz = grid.grid_shape
    h_phys = grid.spacing  # Already in physical units

    # Reshape to 3D grid
    X4_xyz = reshape_flat_to_grid_scalar(X4_phys_flat, nx, ny, nz)

    # Compute potential: Φ = κ_EM * X^4
    Phi_xyz = kappa_EM * X4_xyz

    # Compute electric field: E = -∇Φ = -κ_EM * ∇X^4
    E_xyz = -kappa_EM * grad_mixed(X4_xyz, h_phys, bc)

    # Flatten for output
    E_phys_flat = reshape_grid_to_flat_vec(E_xyz)

    outputs = [E_phys_flat]

    if return_potential:
        Phi_phys_flat = reshape_grid_to_flat_scalar(Phi_xyz)
        outputs.append(Phi_phys_flat)

    if return_charge_density:
        # ρ = ε₀ ∇·E
        rho_xyz = epsilon_0 * div_mixed(E_xyz, h_phys, bc)
        rho_phys_flat = reshape_grid_to_flat_scalar(rho_xyz)
        outputs.append(rho_phys_flat)

    if len(outputs) == 1:
        return outputs[0]
    else:
        return tuple(outputs)


def brane_X4_to_Phi(
    state,
    grid,
    mapper,
    kappa_EM: float = 1.0
) -> torch.Tensor:
    """
    Compute electric potential Φ from brane amplitude X^4.

    Φ(x) = κ_EM * X^4(x)

    Args:
        state: BraneState with positions[:, 3] = X^4 in simulation units
        grid: BraneGrid with grid_shape
        mapper: DimensionalMapper for unit conversion
        kappa_EM: Coupling constant [V/m]

    Returns:
        Phi_phys_flat: Electric potential (N,) in physical units [V]
    """
    # Extract X^4 in simulation units
    X4_sim_flat = state.positions[:, 3]

    # Convert to physical meters
    X4_phys_flat = mapper.to_phys_length(X4_sim_flat)

    # Compute potential: Φ = κ_EM * X^4
    Phi_phys_flat = kappa_EM * X4_phys_flat

    return Phi_phys_flat


def brane_X4_to_rho(
    state,
    grid,
    mapper,
    kappa_EM: float = 1.0,
    bc: BC = ("dirichlet0", "dirichlet0", "dirichlet0"),
    epsilon_0: float = 8.854187817e-12
) -> torch.Tensor:
    """
    Compute charge density ρ from brane amplitude X^4.

    ρ(x) = ε₀ ∇·E(x) where E = -κ_EM ∇X^4

    Args:
        state: BraneState with positions[:, 3] = X^4 in simulation units
        grid: BraneGrid with grid_shape and spacing
        mapper: DimensionalMapper for unit conversion
        kappa_EM: Coupling constant [V/m]
        bc: Boundary conditions per axis
        epsilon_0: Vacuum permittivity [F/m]

    Returns:
        rho_phys_flat: Charge density (N,) in physical units [C/m³]
    """
    # Get E field
    E_phys_flat = brane_X4_to_E_field(state, grid, mapper, kappa_EM, bc,
                                       return_potential=False,
                                       return_charge_density=False)

    # Get grid parameters
    nx, ny, nz = grid.grid_shape
    h_phys = grid.spacing

    # Reshape to grid
    E_xyz = reshape_flat_to_grid_vec(E_phys_flat, nx, ny, nz)

    # Compute divergence: ρ = ε₀ ∇·E
    rho_xyz = epsilon_0 * div_mixed(E_xyz, h_phys, bc)

    # Flatten for output
    rho_phys_flat = reshape_grid_to_flat_scalar(rho_xyz)

    return rho_phys_flat


# ============================================================================
# Inverse Mapping Stub (EM → Brane)
# ============================================================================

def inverse_E_to_X4_opt(
    E_phys_flat: torch.Tensor,
    grid,
    mapper,
    kappa_EM: float = 1.0,
    bc: BC = ("dirichlet0", "dirichlet0", "dirichlet0"),
    regularization: float = 1e-6,
    max_iterations: int = 1000,
    tolerance: float = 1e-8
) -> torch.Tensor:
    """
    Inverse electrostatic mapping: reconstruct X^4 from E field.

    Solves the optimization problem:
        min_{X^4} (1/2) |κ_EM ∇X^4 + E|² + (λ/2)|X^4|²

    using CG or L-BFGS (no FFT).

    NOTE: The solution is unique up to an additive constant (gauge freedom).
          The regularization term λ|X^4|² fixes this gauge by preferring
          solutions with small mean value.

    Args:
        E_phys_flat: Electric field (N, 3) in physical units [V/m]
        grid: BraneGrid with grid_shape and spacing
        mapper: DimensionalMapper for unit conversion
        kappa_EM: Coupling constant [V/m]
        bc: Boundary conditions per axis
        regularization: Regularization parameter λ for gauge fixing
        max_iterations: Maximum number of optimization iterations
        tolerance: Convergence tolerance

    Returns:
        X4_sim_flat: Brane amplitude (N,) in simulation units

    Raises:
        NotImplementedError: This function is not yet implemented

    Implementation Notes:
        When implementing, use one of:
        1. Conjugate Gradient (CG) on the normal equations:
           (κ_EM² ∇† ∇ + λI) X^4 = -κ_EM ∇† E

        2. L-BFGS minimization of the objective:
           f(X^4) = (1/2) |κ_EM ∇X^4 + E|² + (λ/2)|X^4|²

        Both approaches must respect the boundary conditions via the
        mixed finite-difference operators defined above.

    Boundary Condition Consistency:
        If using periodic BC, the system is consistent only if ∮ E·dl = 0
        around any closed loop (i.e., E must be curl-free and have zero
        circulation). Otherwise, no potential Φ exists.

        If using Dirichlet BC, the problem is always well-posed, and the
        regularization term ensures uniqueness by fixing X^4 = 0 on average.
    """
    raise NotImplementedError(
        "Inverse electrostatic mapping (EM → brane) not yet implemented.\n"
        "\n"
        "This function will solve for X^4 given E using an optimization-based\n"
        "approach (CG or L-BFGS) with the boundary conditions specified in 'bc'.\n"
        "\n"
        "Implementation options:\n"
        "  1. Conjugate Gradient on normal equations:\n"
        "     (κ_EM² ∇† ∇ + λI) X^4 = -κ_EM ∇† E\n"
        "\n"
        "  2. L-BFGS minimization:\n"
        "     min_{X^4} (1/2) |κ_EM ∇X^4 + E|² + (λ/2)|X^4|²\n"
        "\n"
        "Both must use the boundary-aware operators (grad_mixed, div_mixed, etc.)\n"
        "to ensure consistency with the forward mapping.\n"
        "\n"
        "Gauge issue: The solution is unique up to an additive constant.\n"
        "The regularization term λ|X^4|² fixes this by preferring small mean.\n"
    )


# ============================================================================
# Backward Compatibility: Initialization from EM Fields
# ============================================================================

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

    This function provides backward compatibility with the old interface by computing
    lateral velocities from the Poynting vector and normal amplitude from energy density.

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


# ============================================================================
# Legacy Compatibility Layer
# ============================================================================

class ElectrostaticMapping:
    """
    Legacy interface for electrostatic mapping.

    DEPRECATED: Use functional API (brane_X4_to_E_field, etc.) instead.

    This class provides backward compatibility with the old class-based interface.
    New code should use the functional API directly for better flexibility with
    boundary conditions.
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
            kappa_EM: Coupling constant [V/m] relating X^4 to Φ
            epsilon_0: Vacuum permittivity [F/m]
            dx: Spatial grid spacing [m]
            device: torch device for computations
            dtype: torch dtype for computations
        """
        self.kappa_EM = kappa_EM
        self.epsilon_0 = epsilon_0
        self.dx = dx
        self.device = device if device is not None else torch.device('cpu')
        self.dtype = dtype

        # Default to Dirichlet BC (clamped boundaries)
        self.bc: BC = ("dirichlet0", "dirichlet0", "dirichlet0")

    def set_boundary_conditions(self, bc: BC):
        """
        Set boundary conditions for all operations.

        Args:
            bc: Tuple of ("periodic" or "dirichlet0") for (x, y, z) axes
        """
        self.bc = bc

    def compute_potential(self, X4: torch.Tensor) -> torch.Tensor:
        """
        Compute electric potential Φ = κ_EM * X^4.

        Args:
            X4: (nx, ny, nz) brane amplitude field

        Returns:
            Φ: (nx, ny, nz) electric potential [V]
        """
        return self.kappa_EM * X4

    def compute_electric_field(self, Phi: torch.Tensor) -> torch.Tensor:
        """
        Compute E = -∇Φ using mixed boundary conditions.

        Args:
            Phi: (nx, ny, nz) electric potential [V]

        Returns:
            E: (nx, ny, nz, 3) electric field [V/m]
        """
        return -grad_mixed(Phi, self.dx, self.bc)

    def compute_charge_density(self, E: torch.Tensor) -> torch.Tensor:
        """
        Compute charge density ρ = ε₀ ∇·E.

        Args:
            E: (nx, ny, nz, 3) electric field [V/m]

        Returns:
            rho: (nx, ny, nz) charge density [C/m³]
        """
        return self.epsilon_0 * div_mixed(E, self.dx, self.bc)

    def map_from_brane(
        self,
        X4: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward mapping: brane → EM.

        Φ = κ_EM * X^4
        E = -∇Φ
        ρ = ε₀ ∇·E

        Args:
            X4: (nx, ny, nz) brane amplitude field

        Returns:
            Phi: (nx, ny, nz) electric potential [V]
            E: (nx, ny, nz, 3) electric field [V/m]
            rho: (nx, ny, nz) charge density [C/m³]
        """
        Phi = self.compute_potential(X4)
        E = self.compute_electric_field(Phi)
        rho = self.compute_charge_density(E)
        return Phi, E, rho

    def compute_electric_field_energy_density(self, E: torch.Tensor) -> torch.Tensor:
        """
        Compute electric field energy density u_E = (1/2) ε₀ |E|².

        Args:
            E: (nx, ny, nz, 3) electric field [V/m]

        Returns:
            u_E: (nx, ny, nz) energy density [J/m³]
        """
        E_squared = torch.sum(E**2, dim=-1)
        return 0.5 * self.epsilon_0 * E_squared


# ============================================================================
# Unit Tests / Sanity Checks
# ============================================================================

def test_boundary_conditions():
    """
    Test boundary condition handling for differential operators.

    Tests:
        1. Periodic BC: sin wave should have correct gradient
        2. Dirichlet BC: zero at boundaries should propagate correctly
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dtype = torch.float64

    # Create test grid
    nx, ny, nz = 32, 32, 32
    Lx = 1.0  # Domain size
    h = Lx / nx

    # Create coordinate arrays
    x = torch.linspace(0, Lx - h, nx, device=device, dtype=dtype)
    y = torch.linspace(0, Lx - h, ny, device=device, dtype=dtype)
    z = torch.linspace(0, Lx - h, nz, device=device, dtype=dtype)

    X, Y, Z = torch.meshgrid(x, y, z, indexing='ij')

    print("=" * 70)
    print("Testing Boundary Conditions")
    print("=" * 70)

    # Test 1: Periodic BC with sine wave
    print("\nTest 1: Periodic BC with sin(2πx/L)")
    X4_periodic = torch.sin(2.0 * torch.pi * X / Lx)

    # Use fully periodic BC for this test (since field is truly periodic in x)
    bc_periodic = ("periodic", "periodic", "periodic")
    E_periodic = -grad_mixed(X4_periodic, h, bc_periodic)

    # Analytical gradient: -∂/∂x[sin(2πx/L)] = -(2π/L)cos(2πx/L)
    Ex_analytical = -(2.0 * torch.pi / Lx) * torch.cos(2.0 * torch.pi * X / Lx)

    error_periodic = torch.abs(E_periodic[..., 0] - Ex_analytical).max()
    print(f"  Max error in Ex: {error_periodic:.6e}")
    print(f"  Ey near zero (interior): {torch.abs(E_periodic[1:-1, 1:-1, 1:-1, 1]).max():.6e}")
    print(f"  Ez near zero (interior): {torch.abs(E_periodic[1:-1, 1:-1, 1:-1, 2]).max():.6e}")

    # Test 2: Dirichlet BC with Gaussian
    print("\nTest 2: Dirichlet BC with Gaussian bump")
    center = Lx / 2
    sigma = Lx / 8
    X4_dirichlet = torch.exp(-((X - center)**2 + (Y - center)**2 + (Z - center)**2) / (2 * sigma**2))

    bc_dirichlet = ("dirichlet0", "dirichlet0", "dirichlet0")
    E_dirichlet = -grad_mixed(X4_dirichlet, h, bc_dirichlet)

    # Divergence of E should be close to Laplacian of X4
    div_E = div_mixed(E_dirichlet, h, bc_dirichlet)
    lap_X4 = laplacian_mixed(X4_dirichlet, h, bc_dirichlet)

    error_div = torch.abs(div_E - lap_X4).max()
    print(f"  Max error in ∇·(-∇X4) = ∇²X4: {error_div:.6e}")

    # Test 3: Laplacian identity ∇²φ = ∇·(∇φ)
    print("\nTest 3: Laplacian identity ∇²φ = ∇·(∇φ)")
    lap_direct = laplacian_mixed(X4_dirichlet, h, bc_dirichlet)
    grad_X4 = grad_mixed(X4_dirichlet, h, bc_dirichlet)
    div_grad = div_mixed(grad_X4, h, bc_dirichlet)

    error_laplacian = torch.abs(lap_direct - div_grad).max()
    print(f"  Max error: {error_laplacian:.6e}")

    print("\n" + "=" * 70)
    print("Boundary condition tests complete!")
    print("=" * 70)


if __name__ == "__main__":
    # Run tests if executed as script
    test_boundary_conditions()