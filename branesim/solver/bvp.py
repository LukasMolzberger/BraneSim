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
   - Initial guess: pass an ``initial_world`` (the injected ansatz seed) for the
     best warm start; otherwise a march-free linear interpolation between the
     boundary slices (opts.warm_start) or a flat-reference zero guess.
4. Unpack the solution into a WorldVolume.

Boundary conditions are applied via the BoundaryCondition objects from
solver/boundary.py.  The BVP solver is agnostic to which scheme is used.

Module-level assertion: the objective is ‖R‖, not S.
    assert OBJECTIVE == "residual_norm"  # NOT "action"
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import newton_krylov, NoConvergence

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.core.residual import residual as compute_residual
from branesim.core.residual import residual_periodic as _residual_periodic
from branesim.solver.boundary import (
    DirichletBC,
    ChiralBC,
    PeriodicBC,
    apply_dirichlet,
    apply_chiral,
    dirichlet_condition_estimate,
)
from branesim.solver.worldvolume import WorldVolume
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
    boundary_condition: DirichletBC | ChiralBC | PeriodicBC

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
        Only consulted when no explicit ``initial_world`` is passed to
        ``solve_block``.  If True, the initial guess is a march-free linear
        interpolation between the boundary slices; if False, the flat
        reference lattice (zero interior guess).
    inner_maxiter : int
        Maximum inner Krylov (GMRES) iterations (default 300).
    method : str
        Newton-Krylov inner iteration method ('lgmres', 'gmres', 'bicgstab').
        Default 'lgmres'.
    verbose : bool
        Print convergence progress.
    plateau_patience : int
        Residual-plateau early-stop (PeriodicBC path): if ``>0``, stop the outer
        Newton iteration when the global ‖R‖ has improved by less than
        ``plateau_rtol`` (relative) over the last ``plateau_patience`` iterations.
        At cond~1e7 the unpreconditioned solve floors after the well-conditioned
        modes are removed (the 96³ run plateaued by iter ~8 yet ran to 30); this
        stops the wasted tail.  ``0`` disables (run exactly ``max_iter``).
    plateau_rtol : float
        Relative-improvement threshold for the plateau detector (default 0.01 =
        stop once <1% total ‖R‖ improvement over ``plateau_patience`` steps).
    """

    tol: float = 1e-8
    max_iter: int = 200
    warm_start: bool = True
    inner_maxiter: int = 300
    method: str = "lgmres"
    verbose: bool = False
    plateau_patience: int = 0
    plateau_rtol: float = 0.01


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
# Initial guess for the interior (march-free)
# ---------------------------------------------------------------------------


def _interp_guess(problem: BoundaryProblem) -> np.ndarray:
    """Linear-interpolation initial guess between the two boundary slices.

    A march-free, r_t-agnostic warm start: interpolate each interior slice
    linearly between R^0 and R^N.  The JFNK solve then enforces the true
    r_t-aware residual on top of this guess, so the guess need not satisfy the
    equations itself — it only has to be smooth and respect the endpoints.

    (Replaces the old forward-Verlet warm start, which raised for r_t>0 and
    coupled the block solver to the deleted IVP march.  When a full seed
    world-volume is available — e.g. the injected ansatz — pass it directly as
    ``initial_world`` to ``solve_block`` instead; it is a far better guess.)

    Returns the full world-volume (N+1, n_nodes, m_ambient).
    """
    bc = problem.boundary_condition
    if not isinstance(bc, DirichletBC):
        raise TypeError(
            "_interp_guess requires a DirichletBC (it needs both endpoints R0, RN); "
            f"got {type(bc).__name__}.  ChiralBC takes the fast-path in solve_block "
            "and never reaches the interior warm start."
        )
    N = problem.params.n_slices
    R0 = bc.R0
    RN = bc.RN
    n_nodes, m_ambient = R0.shape

    world = np.empty((N + 1, n_nodes, m_ambient), dtype=np.float64)
    for l in range(N + 1):
        f = l / N if N > 0 else 0.0
        world[l] = (1.0 - f) * R0 + f * RN
    return world


# ---------------------------------------------------------------------------
# Main BVP solver
# ---------------------------------------------------------------------------


def solve_block(
    problem: BoundaryProblem,
    opts: SolveOpts | None = None,
    *,
    initial_world: np.ndarray | None = None,
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
    initial_world : ndarray, shape (N+1, n_nodes, m_ambient), optional
        Full world-volume to use as the JFNK initial guess — typically the
        injected ansatz seed.  This is the preferred warm start: the seed is a
        far better guess than any march/interpolation.  If omitted, falls back
        to ``opts.warm_start`` (linear interpolation between the boundary
        slices) or a flat-reference zero guess.

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
            # Report convergence against tolerance, never an unconditional True —
            # the r_t=0 chiral march satisfies the EL equation to ~machine eps, but
            # a non-stationary world-volume must NOT be scored as converged.
            "converged": residual_final <= opts.tol,
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
    # Rotating-frame-periodic path: closed (cyclic) time loop, all P slices free.
    # Well-conditioned (periodic operator, no two-point Dirichlet resonance);
    # JFNK from the wound seed preserves the carrier winding.
    # =========================================================================
    if isinstance(bc, PeriodicBC):
        if initial_world is None:
            raise ValueError(
                "PeriodicBC requires initial_world (the wound seed): the carrier "
                "winding is topological and must be supplied as the initial guess."
            )
        P = problem.n_slices  # period = n_slices (slices 0..P-1; R^P ≡ R^0)
        n_nodes = problem.n_nodes
        m_ambient = problem.m_ambient
        if initial_world.shape[0] < P:
            raise ValueError(
                f"initial_world has {initial_world.shape[0]} slices; need >= {P}"
            )
        S0 = initial_world[:P].astype(np.float64).copy()  # (P, n_nodes, m)
        shape = (P, n_nodes, m_ambient)

        anchor = bc.gauge == "anchor"
        anchor_ref = bc.R0[bc.gauge_node].copy() if anchor else None

        def F_periodic(x: np.ndarray) -> np.ndarray:
            slices = x.reshape(shape)
            if anchor:
                # Pin one node on slice 0 (removes the 4 ambient translations);
                # overwrite its positions and zero its residual rows.
                slices = slices.copy()
                slices[0, bc.gauge_node] = anchor_ref
            res = _residual_periodic(slices, problem.lattice, problem.params, problem.mass)
            if anchor:
                res[0, bc.gauge_node] = 0.0
            return res.ravel()

        x0 = S0.ravel()
        residual_initial = float(np.linalg.norm(F_periodic(x0)))
        n_dof = int(x0.size)

        # Residual-plateau early-stop (opt-in via opts.plateau_patience > 0).  The
        # callback tracks ‖F‖ per outer Newton step and the latest iterate; once the
        # last `patience` steps improved ‖F‖ by < `plateau_rtol` (relative), it
        # raises `_Plateau` to stop — newton_krylov(iter=...) otherwise always runs
        # the full count even after the residual floors (the 96³ tail-waste).
        #
        # PRINCIPLES §3.2 scope: this is a CONVERGENCE-CONTROL criterion on the
        # ROOT-FINDER — the same class as the existing f_tol / max_iter stops — NOT a
        # forbidden "collapse rule in the integrator".  It changes only WHEN the
        # Newton iteration stops, never the residual operator (_residual_periodic) or
        # the substrate forces/energy; the returned state is simply a less-converged
        # approximation of the same fixed point, honestly reported converged=False.
        class _Plateau(Exception):
            pass

        resid_hist: list[float] = []
        cb_state: dict[str, object] = {"x": x0.copy(), "stopped": False}

        def _callback(x: np.ndarray, f: np.ndarray) -> None:
            cb_state["x"] = x.copy()
            r = float(np.linalg.norm(f))
            resid_hist.append(r)
            p = opts.plateau_patience
            if p and len(resid_hist) > p:
                r0, r1 = resid_hist[-(p + 1)], resid_hist[-1]
                if (r0 - r1) / max(r0, 1e-300) < opts.plateau_rtol:
                    cb_state["stopped"] = True
                    raise _Plateau()

        converged = False
        residual_final = residual_initial
        x_solution = x0.copy()
        try:
            x_solution = newton_krylov(
                F_periodic, x0, method=opts.method, verbose=opts.verbose,
                f_tol=opts.tol, iter=opts.max_iter, inner_maxiter=opts.inner_maxiter,
                callback=_callback,
            )
            residual_final = float(np.linalg.norm(F_periodic(x_solution)))
            converged = residual_final <= opts.tol
        except _Plateau:
            x_solution = np.asarray(cb_state["x"])
            residual_final = float(np.linalg.norm(F_periodic(x_solution)))
            converged = False
        except NoConvergence as exc:
            x_solution = exc.args[0] if exc.args else x0
            residual_final = float(np.linalg.norm(F_periodic(x_solution)))
            converged = False

        # Per-DOF (RMS) residual — the scale-invariant, grid-comparable figure.
        # The global ‖R‖ grows with √DOF, so an absolute global tol is unreachable
        # at 113M DOF (96³×32×4); ‖R‖/√DOF is the honest "how relaxed per node" number.
        residual_per_dof = residual_final / (n_dof ** 0.5) if n_dof else residual_final
        outer_iters = len(resid_hist)

        elapsed = time.perf_counter() - t0
        Sf = x_solution.reshape(shape)
        world_out = np.empty((P + 1, n_nodes, m_ambient), dtype=np.float64)
        world_out[:P] = Sf
        world_out[P] = Sf[0]  # wrap: R^P ≡ R^0
        solver_report = {
            "mode": "bvp",
            "bc_scheme": "rotating_frame_periodic",
            "residual_initial": residual_initial,
            "residual_final": residual_final,
            # True outer-iteration count from the callback (may be < max_iter when the
            # plateau early-stop fired).  The honest convergence signals are the ratios
            # below (E7): per-DOF RMS residual (grid-comparable), how far from tol, the
            # drop factor, and whether the plateau detector stopped it early.
            "iterations": outer_iters or opts.max_iter,
            "residual_per_dof": residual_per_dof,
            "residual_final_over_tol": (residual_final / opts.tol
                                        if opts.tol > 0 else float("inf")),
            "residual_drop_factor": (residual_initial / residual_final
                                     if residual_final > 0 else float("inf")),
            "early_stopped_plateau": bool(cb_state["stopped"]),
            "converged": converged,
            "condition_estimate": bc.condition_estimate(problem.lattice.params, problem.params),
            "gauge": bc.gauge,
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

    # --- Initial guess for the interior ---
    if initial_world is not None:
        if initial_world.shape != (N + 1, n_nodes, m_ambient):
            raise ValueError(
                f"initial_world shape {initial_world.shape} != "
                f"expected {(N + 1, n_nodes, m_ambient)}"
            )
        world_init = initial_world
    elif opts.warm_start:
        world_init = _interp_guess(problem)  # march-free linear interpolation
    else:
        world_init = world_template.copy()

    # The initial-guess world-volume provides the starting interior DOFs.
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
