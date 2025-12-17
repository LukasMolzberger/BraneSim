"""
Four-potential electromagnetic field solver.

This package implements a four-potential (A^μ) based EM solver using the
d'Alembertian wave equation in Lorenz gauge. It exposes the electromagnetic
tensor F_μν and derived E, B fields.

Key modules:
    - em_state: EMState container for A^μ field
    - potential_solver: Velocity Verlet solver for wave equation
    - derivatives: Finite difference operators (gradient, curl, laplacian)
    - em_tensor: Field tensor computation (E, B, F_μν)
    - initial_conditions_em: Initial condition generators
"""

from branesim.em.em_state import EMState
from branesim.em.potential_solver import FourPotentialVerletSolver
from branesim.em.em_tensor import potentials_to_EB, potentials_to_Fmunu

__all__ = [
    'EMState',
    'FourPotentialVerletSolver',
    'potentials_to_EB',
    'potentials_to_Fmunu',
]