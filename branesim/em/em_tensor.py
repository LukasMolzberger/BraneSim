"""
Electromagnetic field tensor computation.

This module computes physical E and B fields from the four-potential A^μ,
and constructs the electromagnetic tensor F_μν = ∂_μ A_ν - ∂_ν A_μ.
"""

import torch
from branesim.em.derivatives import gradient_scalar, curl, central_diff


def potentials_to_EB(A: torch.Tensor, A_dot: torch.Tensor, grid, c: float, bc: str = "periodic"):
    """
    Compute electric and magnetic fields from four-potential.

    Uses the relations:
        E = -c ∇A^0 - ∂_t A (where A^0 = Φ/c)
        B = ∇ × A

    Args:
        A: Four-potential [N, 4] with components (A^0, Ax, Ay, Az)
        A_dot: Time derivative ∂_t A [N, 4]
        grid: BraneGrid instance with spacing and dimension info
        c: Speed of light in m/s
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        E: Electric field [N, 3] in V/m
        B: Magnetic field [N, 3] in Tesla
    """
    N = A.shape[0]
    A0 = A[:, 0]
    Avec = A[:, 1:4]
    Adot_vec = A_dot[:, 1:4]

    # Reshape to grid for spatial derivatives
    A0g = A0.view(*grid.grid_shape)
    Avecg = Avec.view(*grid.grid_shape, 3)
    Adotg = Adot_vec.view(*grid.grid_shape, 3)

    # Compute E = -c ∇A^0 - ∂_t A
    gradA0 = gradient_scalar(A0g, grid, bc=bc)  # [..., 3]
    E = -c * gradA0 - Adotg  # [..., 3]

    # Compute B = ∇ × A
    B = curl(Avecg, grid, bc=bc)  # [..., 3]

    return E.reshape(N, 3), B.reshape(N, 3)


def potentials_to_Fmunu(A: torch.Tensor, A_dot: torch.Tensor, grid, c: float, bc: str = "periodic"):
    """
    Build electromagnetic field tensor F_μν from four-potential.

    Uses covariant potential A_μ = (A^0, -A^1, -A^2, -A^3) and
    computes F_μν = ∂_μ A_ν - ∂_ν A_μ with ∂_0 = (1/c) ∂_t.

    The tensor components relate to fields as:
        E_i = c F_{0i}  (i = 1, 2, 3)
        B_x = -F_{23}, B_y = -F_{31}, B_z = -F_{12}

    Args:
        A: Four-potential [N, 4] with components (A^0, Ax, Ay, Az)
        A_dot: Time derivative ∂_t A [N, 4]
        grid: BraneGrid instance with spacing and dimension info
        c: Speed of light in m/s
        bc: Boundary condition: "periodic" or "dirichlet"

    Returns:
        F: Electromagnetic tensor [N, 4, 4] (antisymmetric)
    """
    N = A.shape[0]
    ndim = grid.dimension.value

    # Reshape to grid
    Ag = A.view(*grid.grid_shape, 4)
    Vg = A_dot.view(*grid.grid_shape, 4)

    # Convert to covariant components: A_μ = (A^0, -A^1, -A^2, -A^3)
    A_cov = Ag.clone()
    A_cov[..., 1:4] = -A_cov[..., 1:4]

    V_cov = Vg.clone()
    V_cov[..., 1:4] = -V_cov[..., 1:4]

    # Build derivative tensor dA[μ, ν] = ∂_μ A_ν
    # Shape: [4, 4, *grid_shape]
    dA = torch.zeros((4, 4, *grid.grid_shape), device=A.device, dtype=A.dtype)

    # Time derivatives (μ=0): ∂_0 = (1/c) ∂_t
    for nu in range(4):
        dA[0, nu] = (1.0 / c) * V_cov[..., nu]

    # Spatial derivatives (μ=1,2,3): ∂_i A_ν
    h = grid.spacing
    for i in range(min(3, ndim)):
        for nu in range(4):
            dA[i + 1, nu] = central_diff(A_cov[..., nu], axis=i, h=h, bc=bc)

    # Build F_μν = ∂_μ A_ν - ∂_ν A_μ
    F = torch.zeros((4, 4, *grid.grid_shape), device=A.device, dtype=A.dtype)
    for mu in range(4):
        for nu in range(4):
            F[mu, nu] = dA[mu, nu] - dA[nu, mu]

    # Reshape to [N, 4, 4]
    # Move grid dimensions to front, then reshape
    F_flat = F.reshape(4, 4, -1).permute(2, 0, 1)
    return F_flat