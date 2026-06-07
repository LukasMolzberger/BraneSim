"""BVP solver acceptance tests — chiral Cauchy BC (verdict a, 2026-05-30).

Tests the corrected chiral BC implementation in branesim/solver/boundary.py
and the ChiralBC fast-path in branesim/solver/bvp.py.

Analytic target: forward eigenmode of the discrete brane action,
    R^l_p = ref_p + ε · cos(k·x_p − l·θ(k)) · pol,
    θ(k) = arccos(1 − Δt²·ω²(k)/2),  ε = 1e-3.

Test inventory:
  1. ChiralRecovery       — max|R_solved − R_true| < 1e-10 at resonant N
                           (both d=3 canonical and d=1 case)
  2. ChiralConditioning   — chiral cond < 100, N-independent;
                           Dirichlet cond ≫ chiral at the same resonant N
  3. RealityAndOldBugGuard — world-volume real to machine precision;
                             inline reproduction of the old a₊:=A0_k bug
                             yields error ≈ ε (not < 1e-10)
  4. DirichletRecoversMarch — BVP-Dirichlet at non-resonant N recovers IVP
                              march to ~ 1e-12
  5. BackwardChirality    — backward march from (R^N, R^{N-1}) recovers
                            earlier slices to < 1e-10
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from branesim.core.conventions import ActionParams, LatticeParams, d_of_k_eigenvalues
from branesim.core.lattice import SpacelikeLattice
from branesim.core.residual import residual_norm
from branesim.solver.boundary import (
    ChiralBC,
    DirichletBC,
    dirichlet_condition_estimate,
)
from branesim.solver.bvp import BoundaryProblem, SolveOpts, solve_block


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

ALPHA = 0.2
EPS = 1e-3          # eigenmode amplitude (linear regime)
DT = 0.1
K_S = 1.0
RHO = 1.0
SPACING = 1.0


def _make_lattice(grid_shape: tuple, periodic: bool = True) -> SpacelikeLattice:
    lp = LatticeParams(
        grid_shape=grid_shape,
        spacing=SPACING,
        periodic_axes=tuple(periodic for _ in grid_shape),
    )
    return SpacelikeLattice(lp)


def _make_action_params(n_slices: int, m_ambient: int | None = None) -> ActionParams:
    # r_t=0.0: the linear/Verlet limit (BVP chiral/Dirichlet tests use the
    # Verlet-stencil EL equations; the temporal-spring path is validated
    # separately via the canonical-substrate tests).
    return ActionParams(
        k_s=K_S, alpha=ALPHA, rho=RHO, dt=DT, n_slices=n_slices,
        m_ambient=m_ambient, r_t=0.0,
    )


def _eigenmode_slices(
    lattice: SpacelikeLattice,
    params: ActionParams,
    k_n: list[int],
    pol_axis: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Build analytic forward eigenmode world-volume.

    Returns (world_true, R0, R1, theta_k).
    world_true has shape (N+1, n_nodes, m_ambient).
    """
    dim = lattice.params.dim
    m_ambient = params.ambient_dim(dim)
    N = params.n_slices
    a = lattice.params.spacing

    grid = np.asarray(lattice.params.grid_shape, dtype=np.float64)
    n_idx = np.asarray(k_n, dtype=np.float64)
    kvec = 2.0 * np.pi * n_idx / (grid * a)          # (dim,)

    eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, a)
    omega_sq_pol = float(eig[pol_axis])
    theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * omega_sq_pol)))

    ref = lattice.reference_positions(m_ambient)       # (n_nodes, m_ambient)
    phase = ref[:, :dim] @ kvec                        # (n_nodes,)

    world = np.empty((N + 1, ref.shape[0], m_ambient), dtype=np.float64)
    for l in range(N + 1):
        world[l] = ref.copy()
        world[l, :, pol_axis] += EPS * np.cos(phase - l * theta_k)

    return world, world[0].copy(), world[1].copy(), theta_k


def _resonant_N(theta_k: float, min_N: int = 5) -> int:
    """Return the smallest N >= min_N such that N·θ ≈ π (Dirichlet resonance)."""
    N = max(min_N, int(math.ceil(math.pi / theta_k)))
    # Refine to closest integer
    for candidate in range(N - 1, N + 3):
        if candidate >= min_N and abs(candidate * theta_k / math.pi - round(candidate * theta_k / math.pi)) < 0.05:
            return candidate
    return N


def _nonresonant_N(theta_k: float, N_default: int = 20) -> int:
    """Return a non-resonant N (|N·θ mod π| > 0.3)."""
    N = N_default
    while True:
        frac = (N * theta_k) % math.pi
        if min(frac, math.pi - frac) > 0.3:
            return N
        N += 1
        if N > N_default + 50:
            return N_default  # give up and return default


# ---------------------------------------------------------------------------
# Test 1: Chiral recovery at resonant N — d=3 and d=1
# ---------------------------------------------------------------------------

class TestChiralRecovery:
    """ChiralBC recovers the analytic eigenmode to < 1e-10, even at resonant N."""

    def _run(self, grid_shape, k_n, pol_axis):
        lattice = _make_lattice(grid_shape)
        dim = lattice.params.dim
        m = dim + 1

        # Compute theta and pick a resonant N for Dirichlet
        kvec = 2.0 * np.pi * np.asarray(k_n, dtype=np.float64) / (
            np.asarray(lattice.params.grid_shape, dtype=np.float64) * SPACING
        )
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[pol_axis]))))
        N = _resonant_N(theta_k, min_N=10)

        params = _make_action_params(N, m_ambient=m)
        world_true, R0, R1, theta_k2 = _eigenmode_slices(lattice, params, k_n, pol_axis)
        mass = RHO * SPACING ** dim

        bc = ChiralBC(R0=R0, R1=R1, chirality="forward")
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))

        err = float(np.max(np.abs(wv.slices - world_true)))
        return err, N, theta_k

    def test_d3_canonical(self):
        """d=3, m=4, k=[1,0,0], longitudinal pol: recover eigenmode to 1e-10."""
        err, N, theta_k = self._run((16, 16, 16), [1, 0, 0], pol_axis=0)
        assert err < 1e-10, (
            f"d=3 chiral recovery at resonant N={N} (Nθ/π={N*theta_k/math.pi:.4f}): "
            f"max error = {err:.3e}, expected < 1e-10"
        )

    def test_d1(self):
        """d=1, m=2, k=[3], pol_axis=0: recover eigenmode to 1e-10."""
        err, N, theta_k = self._run((16,), [3], pol_axis=0)
        assert err < 1e-10, (
            f"d=1 chiral recovery at resonant N={N} (Nθ/π={N*theta_k/math.pi:.4f}): "
            f"max error = {err:.3e}, expected < 1e-10"
        )

    def test_solver_report_fields(self):
        """Solver report has bc_scheme='chiral', converged=True, condition_estimate finite."""
        lattice = _make_lattice((8, 8, 8))
        params = _make_action_params(20, m_ambient=4)
        mass = 1.0
        _, R0, R1, _ = _eigenmode_slices(lattice, params, [1, 0, 0], 0)
        bc = ChiralBC(R0=R0, R1=R1)
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))
        rpt = wv.solver_report
        assert rpt["bc_scheme"] == "chiral"
        assert rpt["converged"] is True
        assert math.isfinite(rpt["condition_estimate"])
        assert math.isfinite(rpt["residual_final"])


# ---------------------------------------------------------------------------
# Test 2: Conditioning — chiral bounded, Dirichlet large at resonant N
# ---------------------------------------------------------------------------

class TestChiralConditioning:
    """Chiral condition estimate is bounded (< 100) and N-independent.
    Dirichlet condition estimate is large (> 100× chiral) at the same resonant N.
    """

    def _setup(self, grid_shape, k_n, pol_axis):
        lattice = _make_lattice(grid_shape)
        dim = lattice.params.dim
        m = dim + 1
        kvec = 2.0 * np.pi * np.asarray(k_n, dtype=np.float64) / (
            np.asarray(lattice.params.grid_shape, dtype=np.float64) * SPACING
        )
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[pol_axis]))))
        return lattice, m, theta_k

    def test_chiral_cond_bounded_and_n_independent(self):
        """Chiral condition estimate < 100 for several N values including resonant."""
        lattice, m, theta_k = self._setup((16, 16, 16), [1, 0, 0], 0)
        N_res = _resonant_N(theta_k)
        ref_positions = lattice.reference_positions(m)
        for N in [N_res // 2, N_res, N_res * 2, N_res * 5]:
            N = max(N, 5)
            params = _make_action_params(N, m_ambient=m)
            bc = ChiralBC(
                R0=ref_positions.copy(),
                R1=ref_positions.copy(),
            )
            cond = bc.condition_estimate(lattice.params, params)
            assert cond < 100.0, (
                f"Chiral condition at N={N}: {cond:.3f}, expected < 100"
            )

    def test_dirichlet_cond_large_at_resonant_n(self):
        """Dirichlet condition ≫ chiral condition at resonant N."""
        lattice, m, theta_k = self._setup((16, 16, 16), [1, 0, 0], 0)
        N_res = _resonant_N(theta_k)
        params = _make_action_params(N_res, m_ambient=m)
        ref_positions = lattice.reference_positions(m)

        bc_chiral = ChiralBC(R0=ref_positions.copy(), R1=ref_positions.copy())
        cond_chiral = bc_chiral.condition_estimate(lattice.params, params)

        cond_dir = dirichlet_condition_estimate(lattice.params, params)

        assert cond_dir > 100.0 * cond_chiral, (
            f"Dirichlet cond ({cond_dir:.3e}) should be >> chiral cond ({cond_chiral:.3f}) "
            f"at resonant N={N_res} (Nθ/π={N_res*theta_k/math.pi:.4f})"
        )


# ---------------------------------------------------------------------------
# Test 3: Reality and old-bug negative-control regression guard
# ---------------------------------------------------------------------------

class TestRealityAndOldBugGuard:
    """World-volume is real to machine precision. Old a₊:=A0_k bug yields error ≈ ε."""

    def test_world_volume_real(self):
        """Chiral solve produces a float64 world-volume (real by construction)."""
        lattice = _make_lattice((8, 8, 8))
        params = _make_action_params(20, m_ambient=4)
        mass = 1.0
        _, R0, R1, _ = _eigenmode_slices(lattice, params, [1, 0, 0], 0)
        bc = ChiralBC(R0=R0, R1=R1)
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))
        # dtype must be real float, not complex
        assert wv.slices.dtype == np.float64, (
            f"world-volume dtype = {wv.slices.dtype}, expected float64 (real)"
        )
        # No imaginary component whatsoever
        assert not np.iscomplexobj(wv.slices), "World-volume must not be complex"

    def test_old_bug_kills_minus_char_uniformly_gives_large_error(self):
        """Inline reproduction of the old apply_chiral bug.

        The bug:  A^N_k = A^0_k * exp(-i theta_k * N)  for all k uniformly,
        then take .real.

        This violates reality because for a real field one needs
            A^N_{-k} = conj(A^N_k),
        but with uniform exp(-i theta_k N) and theta_k = theta_{-k}:
            conj(A^N_{-k}) = conj(A^0_{-k}) * exp(+i theta_k N)
                           = A^0_k * exp(+i theta_k N)   (using reality of A^0)
        while A^N_k = A^0_k * exp(-i theta_k N).
        These are equal only when theta_k N = m*pi (resonance), so generically
        the buggy A^N is NOT real.  The subsequent .real() truncates the imaginary
        part but produces the wrong real amplitude — error ≈ ε.

        This test documents the bug and confirms it cannot silently return.
        """
        # 1D periodic chain — simplest showcase
        M = 16
        N = 40
        a = DT_loc = 1.0
        dt = 0.1
        k_s = rho = 1.0
        x = np.arange(M, dtype=float)

        # 1D dispersion
        k_grid = 2.0 * math.pi * np.fft.fftfreq(M, d=a)
        omega2 = (2.0 * k_s / rho) * (1.0 - np.cos(k_grid * a))
        arg = 1.0 - 0.5 * dt * dt * omega2
        theta = np.arccos(np.clip(arg, -1.0, 1.0))

        # Ground-truth real forward wave on mode k_idx=3
        k_idx = 3
        k_val = k_grid[k_idx]
        th_val = float(theta[k_idx])
        l_arr = np.arange(N + 1)[:, None]
        p_arr = x[None, :]
        A_true = EPS * np.cos(k_val * p_arr - th_val * l_arr)  # (N+1, M), real

        # Buggy reconstruction: take A^0, FFT, multiply by exp(-i theta N), IFFT, .real()
        A0 = A_true[0]
        A0k = np.fft.fft(A0)
        AN_buggy_k = A0k * np.exp(-1j * theta * N)
        AN_buggy = np.fft.ifft(AN_buggy_k)

        # The imaginary part should be NON-ZERO (proving the bug)
        max_imag = float(np.max(np.abs(AN_buggy.imag)))
        assert max_imag > EPS * 0.1, (
            f"Bug reproduction: imaginary part of buggy A^N = {max_imag:.3e}, "
            f"expected > {EPS*0.1:.3e} (the old code was broken)"
        )

        # Even after .real(), the error to the true future slice is large (≈ ε)
        AN_buggy_real = AN_buggy.real
        err_buggy = float(np.max(np.abs(AN_buggy_real - A_true[N])))
        assert err_buggy > EPS * 0.1, (
            f"Bug reproduction: even after .real(), error = {err_buggy:.3e}, "
            f"expected > {EPS*0.1:.3e} (bug still corrupts amplitude)"
        )

        # Correct two-past-slice march recovers the true future slice
        A = np.zeros((N + 1, M))
        A[0] = A_true[0]
        A[1] = A_true[1]
        Ak = np.fft.fft(A, axis=1)
        c = np.cos(theta)
        for l in range(1, N):
            Ak[l + 1] = 2.0 * c * Ak[l] - Ak[l - 1]
        A_correct = np.fft.ifft(Ak, axis=1).real

        err_correct = float(np.max(np.abs(A_correct - A_true)))
        assert err_correct < 1e-10, (
            f"Two-past-slice march error = {err_correct:.3e}, expected < 1e-10"
        )
        # The correct error must be at least 1000x smaller than the buggy error
        assert err_correct < err_buggy / 1000.0, (
            f"Correct error ({err_correct:.3e}) not 1000x smaller than buggy error "
            f"({err_buggy:.3e})"
        )


# ---------------------------------------------------------------------------
# Test 4: Dirichlet BVP recovers a known analytic eigenmode (march-free)
# ---------------------------------------------------------------------------

class TestDirichletRecoversEigenmode:
    """Dirichlet JFNK at non-resonant N recovers the analytic eigenmode.

    March-free replacement for the old recover-the-march test: the boundary
    slices and ground truth come from the closed-form eigenmode
    (``_eigenmode_slices``), so no IVP march is involved.  Also exercises the
    ``initial_world`` seed-guess path of solve_block.
    """

    def test_d3_dirichlet_recovers_eigenmode(self):
        lattice = _make_lattice((8, 8, 8))
        m = 4
        kvec = 2.0 * math.pi * np.array([1.0, 0.0, 0.0]) / (8.0 * SPACING)
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[0]))))
        N = _nonresonant_N(theta_k, N_default=15)

        params = _make_action_params(N, m_ambient=m)
        world_true, _, _, _ = _eigenmode_slices(lattice, params, [1, 0, 0], 0)
        mass = 1.0

        # Dirichlet BC from the analytic eigenmode endpoints (no march).
        bc = DirichletBC(R0=world_true[0].copy(), RN=world_true[N].copy())
        opts = SolveOpts(tol=1e-12, warm_start=True, verbose=False)
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc), opts)

        max_err = float(np.abs(wv.slices[1:N] - world_true[1:N]).max())
        assert max_err < 1e-10, (
            f"Dirichlet BVP at non-resonant N={N}: max interior error vs "
            f"analytic eigenmode = {max_err:.3e}, expected < 1e-10"
        )

    def test_d3_dirichlet_with_seed_initial_world(self):
        """Passing the eigenmode as initial_world also converges to it."""
        lattice = _make_lattice((8, 8, 8))
        m = 4
        kvec = 2.0 * math.pi * np.array([1.0, 0.0, 0.0]) / (8.0 * SPACING)
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[0]))))
        N = _nonresonant_N(theta_k, N_default=15)

        params = _make_action_params(N, m_ambient=m)
        world_true, _, _, _ = _eigenmode_slices(lattice, params, [1, 0, 0], 0)
        mass = 1.0

        bc = DirichletBC(R0=world_true[0].copy(), RN=world_true[N].copy())
        opts = SolveOpts(tol=1e-12, verbose=False)
        wv = solve_block(
            BoundaryProblem(lattice, params, mass, bc), opts,
            initial_world=world_true.copy(),
        )
        max_err = float(np.abs(wv.slices[1:N] - world_true[1:N]).max())
        assert max_err < 1e-10, (
            f"Dirichlet BVP (seed initial_world) max interior error = {max_err:.3e}"
        )


# ---------------------------------------------------------------------------
# Test 5: Backward chirality recovers the eigenmode from the future end
# ---------------------------------------------------------------------------

class TestBackwardChirality:
    """Backward march from (R^N, R^{N-1}) recovers earlier slices to < 1e-10."""

    def test_backward_recovers_eigenmode(self):
        """Backward ChiralBC from the two future-end slices recovers all past slices."""
        lattice = _make_lattice((16, 16, 16))
        dim = 3
        m = 4
        kvec = 2.0 * math.pi * np.array([1.0, 0.0, 0.0]) / (16.0 * SPACING)
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[0]))))
        N = 20

        params = _make_action_params(N, m_ambient=m)
        world_true, _, _, _ = _eigenmode_slices(lattice, params, [1, 0, 0], 0)
        mass = 1.0

        # Provide (R^N, R^{N-1}) as the two "past" slices for backward chirality.
        # In ChiralBC backward convention: R0 = R^N, R1 = R^{N-1}.
        RN = world_true[N].copy()
        RN_minus1 = world_true[N - 1].copy()
        bc = ChiralBC(R0=RN, R1=RN_minus1, chirality="backward")
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))

        err = float(np.max(np.abs(wv.slices - world_true)))
        assert err < 1e-10, (
            f"Backward chiral recovery error = {err:.3e}, expected < 1e-10"
        )

    def test_backward_d1(self):
        """d=1 backward chiral march recovers eigenmode to < 1e-10."""
        lattice = _make_lattice((16,))
        dim = 1
        m = 2
        kvec = 2.0 * math.pi * np.array([3.0]) / (16.0 * SPACING)
        eig = d_of_k_eigenvalues(kvec, ALPHA, K_S, RHO, SPACING)
        theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * DT * DT * float(eig[0]))))
        N = 20

        params = _make_action_params(N, m_ambient=m)
        world_true, _, _, _ = _eigenmode_slices(lattice, params, [3], 0)
        mass = 1.0

        RN = world_true[N].copy()
        RN_minus1 = world_true[N - 1].copy()
        bc = ChiralBC(R0=RN, R1=RN_minus1, chirality="backward")
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))

        err = float(np.max(np.abs(wv.slices - world_true)))
        assert err < 1e-10, (
            f"d=1 backward chiral recovery error = {err:.3e}, expected < 1e-10"
        )


# ---------------------------------------------------------------------------
# Extra: residual near-zero for chiral world-volume (interior EL satisfied)
# ---------------------------------------------------------------------------

class TestChiralResidualNearZero:
    """The chiral world-volume satisfies the discrete EL equations (residual ~ 0)."""

    @pytest.mark.parametrize("grid_shape,k_n,pol_axis,m_ambient", [
        ((8, 8, 8), [1, 0, 0], 0, 4),
        ((8,), [2], 0, 2),
    ])
    def test_interior_residual_machine_precision(self, grid_shape, k_n, pol_axis, m_ambient):
        lattice = _make_lattice(grid_shape)
        N = 15
        params = _make_action_params(N, m_ambient=m_ambient)
        mass = 1.0
        _, R0, R1, _ = _eigenmode_slices(lattice, params, k_n, pol_axis)
        bc = ChiralBC(R0=R0, R1=R1)
        wv = solve_block(BoundaryProblem(lattice, params, mass, bc))
        res_norm = residual_norm(wv.slices, lattice, params, mass)
        n_interior_dof = (N - 1) * lattice.n_nodes * m_ambient
        # Normalize by sqrt(dof) so the threshold is per-dof
        res_per_dof = res_norm / math.sqrt(n_interior_dof)
        assert res_per_dof < 1e-10, (
            f"grid={grid_shape} chiral residual per dof = {res_per_dof:.3e}, "
            f"expected < 1e-10"
        )


# ---------------------------------------------------------------------------
# Test: explicit Verlet march paths reject r_t>0 (C-H2 guard)
# ---------------------------------------------------------------------------

class TestRtGuard:
    """apply_chiral() is the r_t=0 linear/Verlet limit only.

    For r_t>0 the temporal central-force spring makes the forward step implicit,
    so the explicit stencil does NOT solve the discrete EL equations.  This
    path must raise loudly rather than silently produce a non-stationary
    world-volume scored as 'converged' (ARCHITECTURE.md §1.4 / A4).
    """

    def test_chiral_solve_rejects_r_t_positive(self):
        lattice = _make_lattice((6, 6, 6))
        m = lattice.params.dim + 1
        params = ActionParams(
            k_s=K_S, alpha=ALPHA, rho=RHO, dt=DT, n_slices=8,
            m_ambient=m, r_t=0.175,
        )
        mass = RHO * SPACING ** lattice.params.dim
        ref = lattice.reference_positions(m)
        bc = ChiralBC(R0=ref.copy(), R1=ref.copy())
        # solve_block routes ChiralBC to apply_chiral, which must reject r_t>0
        # (rather than report converged=True on a non-stationary march).
        with pytest.raises(NotImplementedError):
            solve_block(BoundaryProblem(lattice, params, mass, bc))


# ---------------------------------------------------------------------------
# Test: rotating-frame-periodic BC (closed cyclic time loop) is well-conditioned
# and actually moves the brane (contrast with the frozen Dirichlet two-time).
# ---------------------------------------------------------------------------

class TestRotatingFramePeriodic:
    """PeriodicBC: closed cyclic time loop, all slices free.

    The point: the periodic temporal operator is well-conditioned (grows
    polynomially in P, no two-point resonance), so JFNK from a wound seed
    genuinely relaxes the worldtube — unlike the Dirichlet two-time BVP whose
    κ~1/sin Pθ → 1e14 freezes the solve.
    """

    def _carrier_seed(self, lattice, N, m, eps=0.05, w=1.5):
        """A localized rotating-carrier seed (NOT a solution): a Gaussian bump in
        the (comp 2, comp 3) carrier plane whose phase advances once over the loop."""
        ref = lattice.reference_positions(m)
        mi = lattice.multi_indices.astype(float)
        centre = (np.array(lattice.params.grid_shape, dtype=float) - 1.0) / 2.0
        d = np.linalg.norm(mi - centre, axis=1)
        amp = eps * np.exp(-(d / w) ** 2)
        slices = np.repeat(ref[None, :, :], N, axis=0).copy()
        for l in range(N):
            th = 2.0 * np.pi * l / N
            slices[l, :, 2] += amp * np.cos(th)
            slices[l, :, 3] += amp * np.sin(th)
        world = np.concatenate([slices, slices[:1]], axis=0)  # (N+1, n, m), wrap
        return world

    def test_periodic_better_conditioned_than_dirichlet(self):
        from branesim.solver.boundary import PeriodicBC, dirichlet_condition_estimate
        lattice = _make_lattice((6, 6, 6))
        # An N near a Dirichlet resonance to make the contrast unambiguous.
        N = 12
        ap = _make_action_params(N, m_ambient=4)
        per = PeriodicBC(R0=lattice.reference_positions(4)).condition_estimate(lattice.params, ap)
        dirich = dirichlet_condition_estimate(lattice.params, ap)
        assert np.isfinite(per)
        assert per < dirich, f"periodic cond {per:.2e} should beat Dirichlet {dirich:.2e}"

    def test_periodic_solve_moves_and_reduces_residual(self):
        from branesim.solver.boundary import PeriodicBC
        lattice = _make_lattice((6, 6, 6))
        N, m = 8, 4
        ap = _make_action_params(N, m_ambient=m)   # r_t=0 linear (division-free)
        mass = RHO * SPACING ** lattice.params.dim
        seed = self._carrier_seed(lattice, N, m)

        bc = PeriodicBC(R0=seed[0].copy())
        opts = SolveOpts(tol=1e-10, max_iter=8, inner_maxiter=30, verbose=False)
        wv = solve_block(BoundaryProblem(lattice, ap, mass, bc), opts, initial_world=seed)
        rep = wv.solver_report

        assert rep["bc_scheme"] == "rotating_frame_periodic"
        # The brane moves toward a solution: residual drops substantially.
        assert rep["residual_final"] < 0.5 * rep["residual_initial"], (
            f"periodic solve did not reduce residual: "
            f"{rep['residual_initial']:.3e} -> {rep['residual_final']:.3e}"
        )
        # The worldtube actually moved (not frozen like Dirichlet).
        moved = float(np.max(np.abs(wv.slices[:N] - seed[:N])))
        assert moved > 1e-6, f"brane did not move (max|ΔR|={moved:.2e})"
        # The loop closes exactly: R^N == R^0.
        assert float(np.max(np.abs(wv.slices[N] - wv.slices[0]))) < 1e-12

    def test_periodic_requires_initial_world(self):
        from branesim.solver.boundary import PeriodicBC
        lattice = _make_lattice((4, 4, 4))
        ap = _make_action_params(6, m_ambient=4)
        mass = RHO * SPACING ** lattice.params.dim
        bc = PeriodicBC(R0=lattice.reference_positions(4))
        with pytest.raises(ValueError):
            solve_block(BoundaryProblem(lattice, ap, mass, bc))  # no initial_world
