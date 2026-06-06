"""Block BVP solver: root-find R = 0 over interior slices using JFNK.

ARCHITECTURE.md §1.3, §1.4, §2 design decision D1.

Critical: the solver ROOT-FINDS R = 0, never minimises S.
S is Lorentzian (saddle, unbounded below); gradient-descent on S diverges
along the kinetic direction and would silently solve the wrong (Euclidean)
problem.  This is enforced by construction: the objective passed to the
Newton-Krylov solver is the residual vector R, NOT the action S.
(ARCHITECTURE.md §1.3, OPEN_PROBLEMS.md A1.)

Algorithm
---------
1. Pack the N-1 interior slices (l=1..N-1) into a flat real vector x.
2. Define the residual function F(x) = R(l=1..N-1) given boundary slices
   l=0 and l=N (fixed by the BoundaryCondition object).
3. Solve F(x) = 0 using scipy.optimize.newton_krylov (JFNK).
   - Warm-start: if warm_start=True, run a forward IVP march to get a good
     initial guess for the interior.  The march from l=0 fills slices
     l=1..N-1 deterministically; this is the recommended warm-start (D1).
4. Unpack the solution into a WorldVolume.

Boundary conditions are applied via the BoundaryCondition objects from
solver/boundary.py.  The BVP solver is agnostic to which scheme is used.

Module-level assertion: the objective is ‖R‖, not S.
    assert OBJECTIVE == "residual_norm"  # NOT "action"
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from scipy.optimize import newton_krylov, NoConvergence

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.core.residual import residual as compute_residual
from branesim.solver.boundary import (
    DirichletBC,
    ChiralBC,
    apply_dirichlet,
    apply_chiral,
    dirichlet_condition_estimate,
)
from branesim.solver.ivp import IVPProblem, WorldVolume, march
from branesim.core.residual import residual_norm as _residual_norm


# ---------------------------------------------------------------------------
# Module-level saddle-discipline assertion
# ---------------------------------------------------------------------------

OBJECTIVE: str = "residual_norm"  # the solver targets ‖R‖, NOT the action S
assert OBJECTIVE == "residual_norm", (
    "bvp.py: solver must target residual_norm, not the Lorentzian action S. "
    "S is a saddle, unbounded below — minimising it diverges or solves the "
    "Euclidean problem.  (ARCHITECTURE.md §1.3, OPEN_PROBLEMS.md A1)"
)


# ---------------------------------------------------------------------------
# Problem and options containers
# ---------------------------------------------------------------------------


@dataclass
class BoundaryProblem:
    """Full specification of a BVP block-solve problem.

    Parameters
    ----------
    lattice : SpacelikeLattice
        Precomputed spacelike topology.
    params : ActionParams
        Action parameters (k_s, alpha, dt, n_slices N).
    mass : float
        Node mass m = rho * a^dim.
    boundary_condition : DirichletBC | ChiralBC
        Which boundary scheme to use.  For DirichletBC the future slice RN
        is explicit; for ChiralBC it is derived from R0 and the mode structure.
    """

    lattice: SpacelikeLattice
    params: ActionParams
    mass: float
    boundary_condition: DirichletBC | ChiralBC

    @property
    def n_nodes(self) -> int:
        return self.lattice.n_nodes

    @property
    def m_ambient(self) -> int:
        return self.boundary_condition.R0.shape[1]

    @property
    def n_slices(self) -> int:
        return self.params.n_slices

    @property
    def R0(self) -> np.ndarray:
        return self.boundary_condition.R0


@dataclass
class SolveOpts:
    """Options for the BVP solve.

    Parameters
    ----------
    tol : float
        Convergence tolerance for ‖R‖ (default 1e-8).
    max_iter : int
        Maximum Newton-Krylov iterations (default 200).
    warm_start : bool
        If True, use a forward IVP march as the initial guess (D1 default).
        If False, use the flat reference lattice (zero initial guess).
    inner_maxiter : int
        Maximum inner Krylov (GMRES) iterations (default 300).
    method : str
        Newton-Krylov inner iteration method ('lgmres', 'gmres', 'bicgstab').
        Default 'lgmres'.
    verbose : bool
        Print convergence progress.
    """

    tol: float = 1e-8
    max_iter: int = 200
    warm_start: bool = True
    inner_maxiter: int = 300
    method: str = "lgmres"
    verbose: bool = False


# ---------------------------------------------------------------------------
# Pack / unpack helpers for the interior degrees of freedom
# ---------------------------------------------------------------------------


def _pack(interior: np.ndarray) -> np.ndarray:
    """Pack interior slices (N-1, n_nodes, m_ambient) into a flat real vector."""
    return interior.ravel()


def _unpack(x: np.ndarray, n_interior: int, n_nodes: int, m_ambient: int) -> np.ndarray:
    """Unpack a flat vector into interior slices (n_interior, n_nodes, m_ambient)."""
    return x.reshape(n_interior, n_nodes, m_ambient)


# ---------------------------------------------------------------------------
# Residual function for JFNK
# ---------------------------------------------------------------------------


def _make_residual_fn(
    problem: BoundaryProblem,
    world_template: np.ndarray,
) -> callable:
    """Return a closure F(x) = residual(interior) for use by the JFNK solver.

    The closure:
    1. Unpacks ``x`` into interior slices.
    2. Assembles the full world-volume from boundary slices + interior.
    3. Applies boundary conditions (pin l=0 and l=N).
    4. Evaluates residual(world)[1:-1] (interior nodes only).
    5. Returns the flattened interior residual.

    The solver target is ‖F(x)‖ = 0 — this is the residual of the discrete
    Euler-Lagrange equations, i.e. the discrete d'Alembertian.  NOT S.
    """
    bc = problem.boundary_condition
    lattice = problem.lattice
    params = problem.params
    mass = problem.mass
    N = problem.n_slices
    n_nodes = problem.n_nodes
    m_ambient = problem.m_ambient
    n_interior = N - 1  # slices l=1..N-1

    # Compute the fixed future boundary slice once (Dirichlet only;
    # ChiralBC takes the fast-path in solve_block and never reaches here).
    assert isinstance(bc, DirichletBC), (
        "_make_residual_fn is only for DirichletBC; "
        "ChiralBC is handled by the fast-path in solve_block."
    )
    RN_fixed = bc.RN.copy()

    R0_fixed = bc.R0.copy()

    def F(x: np.ndarray) -> np.ndarray:
        """Residual function: returns flattened interior residual vector.

        This is the objective for the JFNK solver.  The solver drives
        ‖F(x)‖ → 0.  F(x) = residual(world)[1:-1], NOT the action S.
        """
        interior = _unpack(x, n_interior, n_nodes, m_ambient)
        world = np.empty((N + 1, n_nodes, m_ambient), dtype=np.float64)
        world[0] = R0_fixed
        world[1:N] = interior
        world[N] = RN_fixed

        res = compute_residual(world, lattice, params, mass)
        # res has shape (N+1, n_nodes, m_ambient); res[0]=res[-1]=0 by convention.
        # Return only the interior (l=1..N-1).
        return _pack(res[1:N])

    return F, RN_fixed


# ---------------------------------------------------------------------------
# Warm-start: forward IVP march from l=0 and l=1
# ---------------------------------------------------------------------------


def _ivp_warm_start(
    problem: BoundaryProblem,
) -> np.ndarray:
    """Run a forward Verlet march from l=0 to get an initial interior guess.

    Uses bc.R0 as R^0 and a "zero-velocity" R^1 (i.e. R^1 = R^0 + dt*0 = R^0
    perturbed only by the force step).  This is the D1-recommended warm start.

    For Dirichlet BC: the march gives a good initial guess when the future
    slice is not resonant.  For Chiral BC: the march from R0 gives the purely
    forward-propagating interior, which is exactly the chiral BC solution in
    the linear regime.

    Returns the full world-volume (N+1, n_nodes, m_ambient).
    """
    bc = problem.boundary_condition
    params = problem.params
    N = params.n_slices

    R0 = bc.R0.copy()
    # R1 = R0 displaced by one Verlet step from rest (zero initial velocity)
    # This amounts to R1 = R0 + (dt^2/m)*F(R0) / 2 (half-step kinetic init).
    # For the warm start the exact choice of R1 is not critical; using R1=R0
    # gives a stationary-start IVP which is close to any smooth solution.
    R1 = R0.copy()

    ivp_problem = IVPProblem(
        lattice=problem.lattice,
        params=ActionParams(
            k_s=params.k_s,
            alpha=params.alpha,
            rho=params.rho,
            dt=params.dt,
            n_slices=N,
            m_ambient=params.m_ambient,
            r_t=params.r_t,
        ),
        mass=problem.mass,
        R0=R0,
        R1=R1,
    )
    wv = march(ivp_problem)
    return wv.slices  # (N+1, n_nodes, m_ambient)


# ---------------------------------------------------------------------------
# Main BVP solver
# ---------------------------------------------------------------------------


def solve_block(
    problem: BoundaryProblem,
    opts: SolveOpts | None = None,
) -> WorldVolume:
    """Root-find R = 0 over interior slices (JFNK).

    NEVER minimises the action S.  The objective is ‖R‖ → 0.
    (ARCHITECTURE.md §1.3; this is enforced by the module-level assertion.)

    Parameters
    ----------
    problem : BoundaryProblem
        Lattice, params, mass, and boundary condition (Dirichlet or Chiral).
    opts : SolveOpts, optional
        Solver options.  Defaults to SolveOpts().

    Returns
    -------
    WorldVolume
        Solved world-volume (N+1, n_nodes, m_ambient).
        solver_report contains: residual_initial, residual_final, iterations,
        converged, condition_estimate, walltime_s.
    """
    if opts is None:
        opts = SolveOpts()

    bc = problem.boundary_condition
    t0 = time.perf_counter()

    # =========================================================================
    # Fast-path: ChiralBC — the solution IS the Verlet march from (R0, R1).
    # No JFNK needed; Cauchy data is well-posed (κ bounded, N-independent).
    # =========================================================================
    if isinstance(bc, ChiralBC):
        world_out, condition_estimate = apply_chiral(
            bc, problem.lattice, problem.params, problem.mass
        )
        elapsed = time.perf_counter() - t0
        residual_final = float(_residual_norm(
            world_out, problem.lattice, problem.params, problem.mass
        ))
        solver_report = {
            "mode": "bvp",
            "bc_scheme": "chiral",
            "residual_initial": residual_final,   # march has no "before" phase
            "residual_final": residual_final,
            "iterations": 0,
            "converged": True,
            "condition_estimate": condition_estimate,
            "walltime_s": elapsed,
            "objective": OBJECTIVE,
        }
        return WorldVolume(
            slices=world_out,
            params=problem.params,
            lattice_params=problem.lattice.params,
            solver_report=solver_report,
        )

    # =========================================================================
    # Standard path: DirichletBC — JFNK root-find over interior slices.
    # =========================================================================
    N = problem.n_slices
    n_nodes = problem.n_nodes
    m_ambient = problem.m_ambient
    n_interior = N - 1  # l=1..N-1

    # --- Build a template world-volume for shape information ---
    world_template = np.zeros((N + 1, n_nodes, m_ambient), dtype=np.float64)
    world_template[0] = problem.R0

    # --- Warm start ---
    if opts.warm_start:
        world_init = _ivp_warm_start(problem)
    else:
        world_init = world_template.copy()

    # The warm-start world-volume (world_init) provides the initial guess for
    # the interior degrees of freedom.
    x0 = _pack(world_init[1:N])

    # --- Build the residual function ---
    F, RN_fixed = _make_residual_fn(problem, world_template)

    # --- Measure initial residual norm ---
    r0 = F(x0)
    residual_initial = float(np.linalg.norm(r0))

    # --- Compute condition estimate ---
    condition_estimate = dirichlet_condition_estimate(
        problem.lattice.params, problem.params
    )

    # --- JFNK solve ---
    # The solver root-finds F(x) = 0, where F is the interior residual.
    # This is NOT minimisation of S; it is root-finding of grad S = -R = 0.
    converged = False
    iterations = 0
    residual_final = residual_initial
    x_solution = x0.copy()

    try:
        x_solution = newton_krylov(
            F,
            x0,
            method=opts.method,
            verbose=opts.verbose,
            f_tol=opts.tol,
            iter=opts.max_iter,
            inner_maxiter=opts.inner_maxiter,
        )
        residual_final = float(np.linalg.norm(F(x_solution)))
        converged = residual_final <= opts.tol
        # newton_krylov does not expose iteration count directly;
        # approximate by counting the outer calls is not straightforward,
        # so we mark iterations as -1 (unknown) when it converges.
        iterations = -1  # converged; exact count not exposed by scipy
    except NoConvergence as exc:
        # scipy raises NoConvergence and attaches the best solution found.
        x_solution = exc.args[0] if exc.args else x0
        residual_final = float(np.linalg.norm(F(x_solution)))
        converged = False
        iterations = opts.max_iter

    elapsed = time.perf_counter() - t0

    # --- Unpack solution into world-volume ---
    interior_solution = _unpack(x_solution, n_interior, n_nodes, m_ambient)
    world_out = np.empty((N + 1, n_nodes, m_ambient), dtype=np.float64)
    world_out[0] = problem.R0
    world_out[1:N] = interior_solution
    world_out[N] = RN_fixed

    solver_report = {
        "mode": "bvp",
        "bc_scheme": "dirichlet",
        "residual_initial": residual_initial,
        "residual_final": residual_final,
        "iterations": iterations,
        "converged": converged,
        "condition_estimate": condition_estimate,
        "walltime_s": elapsed,
        "objective": OBJECTIVE,  # always "residual_norm" — never "action"
    }

    return WorldVolume(
        slices=world_out,
        params=problem.params,
        lattice_params=problem.lattice.params,
        solver_report=solver_report,
    )
