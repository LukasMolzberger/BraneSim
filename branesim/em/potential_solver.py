"""
Four-potential electromagnetic field solver using Velocity Verlet integration.

This module implements a wave equation solver for A^μ using the d'Alembertian:
    ∂²A^μ/∂t² = c² (∇²A^μ + μ₀ J^μ)

in Lorenz gauge with optional gauge constraint enforcement.
"""

import torch
from typing import Optional, Callable
from branesim.em.derivatives import laplacian
from branesim.em.em_tensor import potentials_to_EB


class FourPotentialVerletSolver:
    """
    Evolves A^μ under the wave equation using Velocity Verlet integration.

    The solver advances the four-potential according to:
        ∂²A^μ/∂t² = c² (∇²A^μ + μ₀ J^μ)

    where:
        - A^μ = (Φ/c, Ax, Ay, Az) is the four-potential
        - c is the speed of light
        - μ₀ is the vacuum permeability
        - J^μ = (cρ, Jx, Jy, Jz) is the four-current (optional source)

    The Lorenz gauge condition ∇·A + (1/c) ∂A^0/∂t = 0 can be optionally
    enforced at each time step.

    Attributes:
        dt: Time step in seconds
        c: Speed of light in m/s
        mu0: Vacuum permeability in H/m
        grid: BraneGrid instance
        bc: Boundary condition type ("periodic" or "dirichlet")
        enforce_lorenz: Whether to enforce Lorenz gauge at each step
        current_source: Optional callable returning J^μ [N, 4] at time t
        time: Current simulation time
        step_count: Number of steps taken
    """

    def __init__(
        self,
        dt: float,
        c: float,
        mu0: float,
        grid,
        bc: str = "periodic",
        enforce_lorenz: bool = False,
        current_source: Optional[Callable[[float], torch.Tensor]] = None,
    ):
        """
        Initialize the four-potential solver.

        Args:
            dt: Time step in seconds
            c: Speed of light in m/s
            mu0: Vacuum permeability in H/m
            grid: BraneGrid instance with spacing and dimension info
            bc: Boundary condition: "periodic" or "dirichlet"
            enforce_lorenz: If True, enforce Lorenz gauge at each step
            current_source: Optional function(t) -> J^μ [N, 4] in A/m²
        """
        self.dt = dt
        self.c = c
        self.mu0 = mu0
        self.grid = grid
        self.bc = bc
        self.enforce_lorenz = enforce_lorenz
        self.current_source = current_source

        self.time = 0.0
        self.step_count = 0

    def _compute_accel(self, state) -> torch.Tensor:
        """
        Compute acceleration ∂²A^μ/∂t² = c² (∇²A^μ + μ₀ J^μ).

        Args:
            state: EMState instance

        Returns:
            Acceleration tensor [N, 4]
        """
        # Reshape to grid for spatial derivatives
        Ag = state.view_grid(state.potential)  # [..., 4]
        lapA = laplacian(Ag, self.grid, bc=self.bc)  # [..., 4]

        # Add source term if present
        if self.current_source is None:
            J = 0.0
        else:
            J = self.current_source(self.time).to(device=state.device, dtype=state.dtype)  # [N, 4]
            J = state.view_grid(J)  # [..., 4]

        # Wave equation: A_tt = c² (∇²A + μ₀ J)
        A_ddot = (self.c * self.c) * (lapA + self.mu0 * J)
        return state.flatten_grid(A_ddot)

    def initialize_accelerations(self, state) -> None:
        """
        Initialize acceleration fields before starting time evolution.

        Args:
            state: EMState instance
        """
        state.accel = self._compute_accel(state)
        state.new_accel = state.accel.clone()
        state.apply_fixed_boundaries()

    def step(self, state):
        """
        Advance the system by one time step using Velocity Verlet.

        The Velocity Verlet algorithm:
            1. A(t+dt) = A(t) + v(t)*dt + 0.5*a(t)*dt²
            2. Compute a(t+dt) from new positions
            3. v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt

        Args:
            state: EMState instance to update in-place

        Returns:
            state: Updated EMState (same object)
        """
        dt = self.dt

        # Phase 1: Update positions (potentials)
        state.potential += state.velocity * dt + 0.5 * state.accel * dt * dt
        state.apply_fixed_boundaries()

        # Phase 2: Compute new acceleration at updated positions
        state.new_accel = self._compute_accel(state)
        state.apply_fixed_boundaries()

        # Phase 3: Update velocities using average acceleration
        state.velocity += 0.5 * (state.accel + state.new_accel) * dt
        state.apply_fixed_boundaries()

        # Swap accelerations
        state.accel = state.new_accel.clone()

        # Optional Lorenz gauge enforcement
        if self.enforce_lorenz:
            state.enforce_lorenz_gauge(self.grid, c=self.c, bc=self.bc)

        self.time += dt
        self.step_count += 1
        return state

    def compute_fields(self, state):
        """
        Compute E and B fields from current four-potential.

        Args:
            state: EMState instance

        Returns:
            E: Electric field [N, 3] in V/m
            B: Magnetic field [N, 3] in Tesla
        """
        return potentials_to_EB(state.potential, state.velocity, self.grid, c=self.c, bc=self.bc)

    def compute_dalembert_residual(self, state) -> torch.Tensor:
        """
        Compute d'Alembertian residual for diagnostics.

        The residual is: R^μ = (1/c²) ∂²A^μ/∂t² - ∇²A^μ - μ₀ J^μ
        This should be near zero if the wave equation is satisfied.

        Args:
            state: EMState instance

        Returns:
            Residual tensor [N, 4]
        """
        # Get acceleration
        A_tt = state.accel  # [N, 4]

        # Compute Laplacian
        Ag = state.view_grid(state.potential)
        lapA = laplacian(Ag, self.grid, bc=self.bc)
        lapA_flat = state.flatten_grid(lapA)

        # Get source
        if self.current_source is None:
            J = torch.zeros_like(state.potential)
        else:
            J = self.current_source(self.time).to(device=state.device, dtype=state.dtype)

        # Residual
        residual = (1.0 / (self.c * self.c)) * A_tt - lapA_flat - self.mu0 * J
        return residual

    def cfl_ok(self) -> bool:
        """
        Check if CFL condition is satisfied.

        Conservative CFL condition: dt ≤ h / (c * √D)
        where D is the spatial dimensionality.

        Returns:
            True if CFL condition is satisfied
        """
        D = self.grid.dimension.value
        cfl_limit = self.grid.spacing / (self.c * (D ** 0.5))
        return self.dt <= cfl_limit

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"FourPotentialVerletSolver(dt={self.dt:.2e}s, "
            f"c={self.c:.2e}m/s, t={self.time:.2e}s, "
            f"steps={self.step_count}, CFL_ok={self.cfl_ok()})"
        )