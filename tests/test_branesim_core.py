"""Acceptance tests for branesim core physics.

Covers the analytic / solver-independent acceptance criteria:

  2. Dispersion regression: c_L = 1 +/- 5%, c_T/c_L = sqrt(1-alpha) +/- 5%
     (closed-form D(k); the march-propagation measurements were removed with
     the IVP solver, 2026-06-06 — the analytic conventions remain the check).
  3. D(k) diagonal / closed form: eigenvalues match analytic formula.
  4. Flat prestressed lattice: zero net force on periodic flat lattice.
  6. Saddle guard: S is a saddle along the kinetic direction.

Physics targets (alpha=0.2, dimensionless k_s=a=rho=1):
  c_L = 1.0,  c_T = sqrt(0.8) = 0.8944
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from branesim.core.conventions import (
    ActionParams,
    LatticeParams,
    c_longitudinal,
    c_transverse,
    d_of_k_eigenvalues,
    speed_ratio,
)
from branesim.core.lattice import SpacelikeLattice
from branesim.core.action import spacelike_force, spacelike_potential, action


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

ALPHA = 0.2
ATOL_FORCE = 1e-12          # flat-lattice zero-force tolerance
ATOL_RESIDUAL = 1e-10       # residual ~ 0 tolerance (machine precision target)
DISPERSION_TOL = 0.05       # 5% tolerance for c_L and c_T/c_L


def make_lattice(grid_shape, periodic=True, spacing=1.0):
    periodic_axes = tuple(periodic for _ in grid_shape)
    lp = LatticeParams(
        grid_shape=grid_shape,
        spacing=spacing,
        periodic_axes=periodic_axes,
    )
    return SpacelikeLattice(lp)


def flat_positions(lattice: SpacelikeLattice, m_ambient: int) -> np.ndarray:
    """Return the flat (unstressed) reference positions."""
    return lattice.reference_positions(m_ambient)


def make_params(n_slices=20, dt=0.1, alpha=ALPHA):
    # r_t=0.0: the linear/Verlet limit (dispersion and residual validation tests).
    return ActionParams(
        k_s=1.0, alpha=alpha, rho=1.0, dt=dt, n_slices=n_slices, r_t=0.0,
    )


# ---------------------------------------------------------------------------
# Test 4: Flat prestressed lattice -- zero net force
# ---------------------------------------------------------------------------

class TestFlatLatticeZeroForce:
    """Acceptance criterion 4: no net force on a fully periodic flat lattice.

    On the flat lattice all springs are at the same length (= a = spacing).
    By the inversion symmetry of the axial stencil (every +delta offset is
    paired with a -delta offset), the forces from opposite neighbors cancel
    exactly, giving zero net force at every node regardless of alpha.
    """

    @pytest.mark.parametrize("grid_shape", [(4, 4, 4), (6, 6, 6)])
    def test_d3_periodic(self, grid_shape):
        lattice = make_lattice(grid_shape, periodic=True)
        params = make_params()
        pos = flat_positions(lattice, m_ambient=4)

        forces = spacelike_force(pos, lattice, params)

        assert np.max(np.abs(forces)) < ATOL_FORCE, (
            f"Flat lattice d=3 periodic: max |F| = {np.max(np.abs(forces)):.3e}, "
            f"expected < {ATOL_FORCE}"
        )

    @pytest.mark.parametrize("grid_shape", [(4, 4), (8, 8)])
    def test_d2_periodic(self, grid_shape):
        lattice = make_lattice(grid_shape, periodic=True)
        params = make_params()
        pos = flat_positions(lattice, m_ambient=3)

        forces = spacelike_force(pos, lattice, params)

        assert np.max(np.abs(forces)) < ATOL_FORCE, (
            f"Flat lattice d=2 periodic: max |F| = {np.max(np.abs(forces)):.3e}"
        )

    def test_d1_periodic(self):
        lattice = make_lattice((8,), periodic=True)
        params = make_params()
        pos = flat_positions(lattice, m_ambient=2)

        forces = spacelike_force(pos, lattice, params)

        assert np.max(np.abs(forces)) < ATOL_FORCE, (
            f"Flat lattice d=1 periodic: max |F| = {np.max(np.abs(forces)):.3e}"
        )


# ---------------------------------------------------------------------------
# Test 2: Dispersion regression -- analytic c_L and c_T/c_L from closed-form D(k)
# ---------------------------------------------------------------------------

class TestDispersionRegression:
    """Acceptance criterion 2: c_L = 1 +/- 5%, c_T/c_L = sqrt(1-alpha) +/- 5%.

    These are analytic (closed-form D(k)) regressions of the continuum wave
    speeds — Level-1 of the emergence hierarchy.  (The march-propagation
    measurements that previously fed these checks were removed with the IVP
    solver, 2026-06-06; the closed-form conventions remain the regression.)
    """

    def test_c_longitudinal_continuum_limit(self):
        """At small ka=0.1, closed-form D(k) gives c_L within 5% of 1.0 (ka->0)."""
        k_phys = 0.1
        eig = d_of_k_eigenvalues(np.array([k_phys, 0.0, 0.0]), ALPHA)
        c_L_discrete = math.sqrt(eig[0]) / k_phys
        c_L_target = c_longitudinal()  # = 1.0
        assert abs(c_L_discrete - c_L_target) / c_L_target < DISPERSION_TOL, (
            f"D(k) c_L at ka=0.1: {c_L_discrete:.4f} vs target {c_L_target:.4f}"
        )

    def test_c_transverse_ratio_continuum_limit(self):
        """At small ka, closed-form D(k) gives c_T/c_L = sqrt(1-alpha) within 5%.

        Longitudinal e_x-polarized vs transverse e_y-polarized eigenvalues of
        D(k e_x) give the wave speeds; their ratio -> sqrt(1-alpha) as ka->0.
        """
        k_phys = 0.1
        eig = d_of_k_eigenvalues(np.array([k_phys, 0.0, 0.0]), ALPHA)
        c_L = math.sqrt(eig[0]) / k_phys   # longitudinal (along k)
        c_T = math.sqrt(eig[1]) / k_phys   # transverse (perp to k)
        ratio_meas = c_T / c_L
        ratio_target = speed_ratio(ALPHA)  # sqrt(0.8) = 0.8944
        assert abs(ratio_meas - ratio_target) / ratio_target < DISPERSION_TOL, (
            f"c_T/c_L = {ratio_meas:.4f}, target = {ratio_target:.4f}, "
            f"error = {abs(ratio_meas - ratio_target)/ratio_target*100:.1f}%"
        )


# ---------------------------------------------------------------------------
# Test 3: D(k) diagonal / closed form
# ---------------------------------------------------------------------------

class TestDynamicalMatrix:
    """Acceptance criterion 3: D(k) eigenvalues match the analytic formula."""

    def test_eigenvalues_match_analytic_formula(self):
        """Closed form agrees with the manual formula at several k."""
        k_vectors = [
            np.array([0.1, 0.0, 0.0]),
            np.array([0.2, 0.1, 0.0]),
            np.array([0.1, 0.1, 0.1]),
        ]
        k_s, rho, a, alpha = 1.0, 1.0, 1.0, ALPHA

        for k in k_vectors:
            eig = d_of_k_eigenvalues(k, alpha, k_s, rho, a)
            h = 1.0 - np.cos(k * a)
            h_sum = h.sum()
            expected = (2.0 * k_s / rho) * (alpha * h + (1.0 - alpha) * h_sum)
            np.testing.assert_allclose(
                eig, expected, rtol=1e-12, err_msg=f"D(k) mismatch at k={k}"
            )

    def test_longitudinal_speed_from_dofk(self):
        """D(k) longitudinal eigenvalue / k^2 matches continuum c_L^2."""
        for ka in [0.05, 0.1, 0.2]:
            k = np.array([ka, 0.0, 0.0])
            eig = d_of_k_eigenvalues(k, ALPHA)
            # omega^2 = (2 k_s/rho)*[alpha*h_x + (1-alpha)*h_x] = 2*h_x
            # c^2 = omega^2/k^2 = 2*(1-cos(ka))/ka^2 -> 1 as ka->0
            c2_L_exact = 2.0 * (1.0 - math.cos(ka)) / ka ** 2
            assert abs(eig[0] - c2_L_exact * ka ** 2) < 1e-12, (
                f"D_L(k) at ka={ka}: eig={eig[0]:.10f}, expected={c2_L_exact*ka**2:.10f}"
            )


# ---------------------------------------------------------------------------
# Test 6: Saddle guard -- gradient descent on S diverges
# ---------------------------------------------------------------------------

class TestSaddleGuard:
    """Acceptance criterion 6: S is Lorentzian -- a saddle.

    A positive shift to an interior slice increases the kinetic term T
    (which enters S with + sign) and hence S itself.  Successive shifts
    in the same direction keep increasing S -- showing the saddle nature
    of S along the kinetic (temporal) direction.

    The corollary is that gradient descent (which moves in the -gradS
    direction) would diverge along this direction, while root-finding
    (targetting gradS = R = 0) converges.
    """

    def test_lorentzian_unboundedness(self):
        """S is Lorentzian (saddle): unbounded above along the kinetic direction
        and unbounded below along the potential direction.

        This test would FAIL for a Euclidean (all-positive) action because a
        Euclidean action is bounded below and therefore cannot decrease without
        bound when potential-direction displacements grow.

        Kinetic direction (T term, enters S with + sign):
            Shifting an interior slice away from its neighbours grows the
            velocity squared, so T and hence S grow.  We verify S(A) - S(0)
            is strictly increasing and reaches >> 1 as A spans many decades.
            S - S(0) scales as n_nodes * A^2 / dt^2 (purely kinetic), so
            it is analytically unbounded above.

        Potential direction (V term, enters S with - sign):
            Radially stretching the open-boundary lattice grows all spring
            strains monotonically, so V grows and S = T - V decreases without
            bound.  T = 0 (all slices identical) throughout this sweep.
        """
        # --- Kinetic direction: S - S(0) must grow without bound (periodic OK) ---
        lattice_per = make_lattice((4, 4, 4), periodic=True)
        params = make_params(n_slices=4, dt=0.1)
        mass = 1.0
        R0_per = flat_positions(lattice_per, m_ambient=4)
        world_ref = np.stack([R0_per] * (params.n_slices + 1))
        S0 = action(world_ref, lattice_per, params, mass)

        # Displace interior slice l=2 uniformly; this increases T^{3/2} and T^{5/2}.
        amplitudes_kinetic = [1e-2, 1e-1, 1.0, 10.0, 100.0]
        delta_S_prev = None
        for A in amplitudes_kinetic:
            world = world_ref.copy()
            world[2, :, 0] += A
            delta_S = action(world, lattice_per, params, mass) - S0
            assert delta_S > 0.0, (
                f"Kinetic direction: S - S(0) not positive at A={A:.3e}: "
                f"delta_S={delta_S:.6g}"
            )
            if delta_S_prev is not None:
                assert delta_S > delta_S_prev, (
                    f"Kinetic direction: S - S(0) did not grow at A={A:.3e}: "
                    f"prev={delta_S_prev:.6g}, curr={delta_S:.6g}"
                )
            delta_S_prev = delta_S

        # Must grow without bound: S(A=100) - S(0) >> 1
        assert delta_S > 1e5, (
            f"Kinetic direction: S(A=100) - S(0) = {delta_S:.3e}, "
            f"expected > 1e5 (Lorentzian unbounded above; scales as n_nodes*A^2/dt^2)"
        )

        # --- Potential direction: S must decrease without bound (open BC) ---
        # Open-boundary lattice: radial stretch grows all spring strains monotonically.
        # T = 0 throughout (all slices identical), so S = -dt * n_slices * V strictly.
        lattice_open = make_lattice((4, 4, 4), periodic=False)
        R0_open = flat_positions(lattice_open, m_ambient=4)
        center = R0_open[:, :3].mean(axis=0)
        S_prev = None
        amplitudes_potential = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
        for A in amplitudes_potential:
            pos_stretched = R0_open.copy()
            pos_stretched[:, :3] = center + (1.0 + A) * (R0_open[:, :3] - center)
            world = np.stack([pos_stretched] * (params.n_slices + 1))
            S = action(world, lattice_open, params, mass)
            if S_prev is not None:
                assert S < S_prev, (
                    f"Potential direction: S did not decrease at A={A:.2f}: "
                    f"S_prev={S_prev:.6g}, S={S:.6g}"
                )
            S_prev = S

        # Must decrease without bound: S(A=10) << S(A=0)
        S_at_zero = action(
            np.stack([R0_open] * (params.n_slices + 1)), lattice_open, params, mass
        )
        assert S_prev < S_at_zero - 1e3, (
            f"Potential direction: S(A=10)={S_prev:.6g} should be << "
            f"S(A=0)={S_at_zero:.6g} by > 1e3 (Lorentzian unbounded below)"
        )

    def test_flat_world_S_is_finite_and_negative(self):
        """S = N*dt*(-V) < 0 for a stationary flat world (T=0, V>0 due to prestress)."""
        lattice = make_lattice((4, 4, 4), periodic=True)
        params = make_params(n_slices=4, dt=0.1)
        mass = 1.0
        R0 = flat_positions(lattice, m_ambient=4)

        # Stationary: all slices identical; T=0 everywhere
        world_static = np.stack([R0] * (params.n_slices + 1))
        S_static = action(world_static, lattice, params, mass)

        # V > 0 on each slice (prestressed lattice: springs not at rest length)
        V = spacelike_potential(R0, lattice, params)
        assert V > 0.0, f"Prestressed flat lattice should have V > 0; got {V:.6f}"

        # S = sum_l dt*(0 - V) = -N*dt*V < 0 (Lorentzian sign: T enters +, V enters -)
        S_expected = -params.n_slices * params.dt * V
        assert abs(S_static - S_expected) < 1e-10, (
            f"S = {S_static:.6f}, expected {S_expected:.6f}"
        )
        assert S_static < 0.0, f"Stationary flat world S should be < 0; got {S_static}"


# ---------------------------------------------------------------------------
# Conventions: unit tests for closed-form speed helpers
# ---------------------------------------------------------------------------

class TestConventions:
    def test_c_longitudinal_dimensionless(self):
        assert abs(c_longitudinal() - 1.0) < 1e-15

    def test_c_transverse_dimensionless(self):
        assert abs(c_transverse(alpha=ALPHA) - math.sqrt(1.0 - ALPHA)) < 1e-15

    def test_speed_ratio(self):
        assert abs(speed_ratio(ALPHA) - math.sqrt(1.0 - ALPHA)) < 1e-15

    def test_speed_ratio_alpha0(self):
        assert abs(speed_ratio(0.0) - 1.0) < 1e-15

    def test_speed_ratio_alpha1(self):
        assert abs(speed_ratio(1.0) - 0.0) < 1e-15

    def test_d_of_k_zero_wavevector(self):
        """At k=0 all eigenvalues are 0 (acoustic modes)."""
        k = np.zeros(3)
        eig = d_of_k_eigenvalues(k, ALPHA)
        np.testing.assert_allclose(eig, 0.0, atol=1e-15)

    def test_d_of_k_symmetry(self):
        """Eigenvalues are symmetric: omega_a^2(k) = omega_a^2(-k)."""
        k = np.array([0.3, 0.1, 0.2])
        eig_pos = d_of_k_eigenvalues(k, ALPHA)
        eig_neg = d_of_k_eigenvalues(-k, ALPHA)
        np.testing.assert_allclose(eig_pos, eig_neg, rtol=1e-14)

    @pytest.mark.parametrize("dim", [1, 2, 3])
    def test_d_of_k_dimension_agnostic(self, dim):
        """d_of_k_eigenvalues works for any dimension."""
        k = 0.1 * np.ones(dim)
        eig = d_of_k_eigenvalues(k, ALPHA)
        assert eig.shape == (dim,)
        assert np.all(eig >= 0.0)


# ---------------------------------------------------------------------------
# Lattice topology: unit tests
# ---------------------------------------------------------------------------

class TestLatticTopology:
    def test_n_neighbors_d1(self):
        lattice = make_lattice((8,), periodic=True)
        assert lattice.n_neighbors == 2

    def test_n_neighbors_d2(self):
        lattice = make_lattice((4, 4), periodic=True)
        assert lattice.n_neighbors == 4

    def test_n_neighbors_d3(self):
        lattice = make_lattice((4, 4, 4), periodic=True)
        assert lattice.n_neighbors == 6

    def test_periodic_all_valid(self):
        """Every node has all neighbors valid in a fully periodic lattice."""
        lattice = make_lattice((4, 4, 4), periodic=True)
        assert np.all(lattice.neighbors >= 0)

    def test_open_boundary_missing_neighbors(self):
        """Corner nodes have fewer than 2*dim valid neighbors in open BC."""
        lattice = make_lattice((4, 4, 4), periodic=False)
        has_invalid = np.any(lattice.neighbors == -1)
        assert has_invalid, "Open BC lattice should have some -1 entries"

    def test_reference_positions_shape(self):
        lattice = make_lattice((4, 4, 4), periodic=True)
        ref = lattice.reference_positions(m_ambient=4)
        assert ref.shape == (64, 4)

    def test_reference_positions_spacing(self):
        """Adjacent nodes (multi-index distance 1) are spacing apart."""
        lattice = make_lattice((4, 4, 4), periodic=True)
        ref = lattice.reference_positions(m_ambient=4)
        mi = lattice.multi_indices
        for node_a in range(8):
            for node_b in range(8):
                if np.sum(np.abs(mi[node_a] - mi[node_b])) == 1:
                    dist = np.linalg.norm(ref[node_a] - ref[node_b])
                    assert abs(dist - 1.0) < 1e-14, f"Spacing error: {dist}"


# ---------------------------------------------------------------------------
# I/O contracts: round-trip test
# ---------------------------------------------------------------------------

class TestIOContracts:
    def test_boundary_problem_roundtrip(self, tmp_path):
        from branesim.io.contracts import save_boundary_problem, load_boundary_problem

        lattice = make_lattice((4, 4, 4), periodic=True)
        params = make_params(n_slices=10, dt=0.1)

        R0 = flat_positions(lattice, m_ambient=4)
        R1 = R0.copy()
        R1[:, 0] += 1e-3

        path = tmp_path / "bp.npz"
        lattice_dict = {
            "grid_shape": list(lattice.params.grid_shape),
            "spacing": lattice.params.spacing,
            "periodic_axes": list(lattice.params.periodic_axes),
            "axial_weight": lattice.params.axial_weight,
            "dim": lattice.dim,
        }
        action_dict = {
            "k_s": params.k_s,
            "alpha": params.alpha,
            "rho": params.rho,
            "dt": params.dt,
            "n_slices": params.n_slices,
            "r_t": params.r_t,
        }

        save_boundary_problem(
            path,
            ref_positions=R0,
            boundary_slices=np.stack([R0, R1]),
            boundary_indices=np.array([0, 1]),
            lattice=lattice_dict,
            action=action_dict,
        )

        loaded = load_boundary_problem(path)
        np.testing.assert_array_equal(loaded["boundary_indices"], [0, 1])
        np.testing.assert_allclose(loaded["boundary_slices"][0], R0)
        np.testing.assert_allclose(loaded["boundary_slices"][1], R1)
        assert loaded["lattice"]["dim"] == 3

    def test_worldvolume_roundtrip(self, tmp_path):
        from branesim.io.contracts import WorldVolumeWriter, load_manifest, iter_slices

        lattice = make_lattice((4, 4, 4), periodic=True)
        params = make_params(n_slices=3, dt=0.1)
        R0 = flat_positions(lattice, m_ambient=4)

        # Synthetic world-volume (the IO roundtrip needs an array, not dynamics):
        # stack R0 with a small per-slice perturbation.
        slices = np.stack([R0] * (params.n_slices + 1))
        for l in range(1, params.n_slices + 1):
            slices[l, :, 0] += 1e-3 * l

        path = tmp_path / "wv.zip"
        with WorldVolumeWriter(path, {"mode": "synthetic"}) as writer:
            for l in range(params.n_slices + 1):
                writer.write_slice(l, l * params.dt, slices[l])
            writer.write_npy("aux/ref_positions.npy", R0)

        manifest = load_manifest(path)
        assert manifest["format_version"] == "branesim-block-v1"
        assert len(manifest["slices"]) == params.n_slices + 1

        slices_loaded = list(iter_slices(path))
        assert len(slices_loaded) == params.n_slices + 1
        idx, t, pos = slices_loaded[0]
        np.testing.assert_allclose(pos, R0)


# ---------------------------------------------------------------------------
# Test W3: Force–energy consistency (finite-difference gradient check)
# ---------------------------------------------------------------------------

class TestForceEnergyConsistency:
    """F_i = -dV/dR_i verified by central finite differences on random configurations.

    For a random, non-flat, non-symmetric configuration we check:
        F_i(R) ≈ -(V(R + eps*e_i) - V(R - eps*e_i)) / (2*eps)

    over a sample of random (node, component) pairs.  Agreement should be
    ~1e-6 or better at eps=1e-5 (central-difference truncation O(eps^2)).

    This is the primary physics-core consistency check: if force is not the
    exact gradient of the potential, the Verlet integrator violates
    Substrate-Only Evolution (principles §1.2, non-negotiable #1).
    """

    RNG_SEED = 42
    EPS = 1e-5
    ATOL_FD = 1e-6        # central-difference agreement tolerance

    def _random_displaced_positions(
        self,
        lattice: SpacelikeLattice,
        m_ambient: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Flat reference + random node-wise displacements of O(0.1)."""
        pos = flat_positions(lattice, m_ambient).astype(np.float64)
        # Random perturbation: amplitude ~0.1 (non-linear but not collapsing)
        pos += 0.1 * rng.standard_normal(pos.shape)
        return pos

    def _fd_gradient(
        self,
        pos: np.ndarray,
        node: int,
        comp: int,
        lattice: SpacelikeLattice,
        params: ActionParams,
    ) -> float:
        """Central-difference estimate of -dV/dR[node, comp]."""
        pos_p = pos.copy()
        pos_p[node, comp] += self.EPS
        pos_m = pos.copy()
        pos_m[node, comp] -= self.EPS
        V_p = spacelike_potential(pos_p, lattice, params)
        V_m = spacelike_potential(pos_m, lattice, params)
        return -(V_p - V_m) / (2.0 * self.EPS)

    @pytest.mark.parametrize("grid_shape,m_ambient", [
        ((6, 6, 6), 4),
        ((4, 4), 3),
        ((8,), 2),
    ])
    def test_force_equals_negative_grad_potential(self, grid_shape, m_ambient):
        """F_i = -dV/dR_i for 20 random (node, component) pairs."""
        rng = np.random.default_rng(self.RNG_SEED)
        lattice = make_lattice(grid_shape, periodic=True)
        params = make_params()

        pos = self._random_displaced_positions(lattice, m_ambient, rng)

        # Analytic forces
        forces = spacelike_force(pos, lattice, params)

        # Sample 20 random (node, component) pairs
        n_samples = 20
        nodes = rng.integers(0, lattice.n_nodes, size=n_samples)
        comps = rng.integers(0, m_ambient, size=n_samples)

        max_abs_err = 0.0
        max_rel_err = 0.0
        worst = None

        for node, comp in zip(nodes, comps):
            F_analytic = forces[node, comp]
            F_fd = self._fd_gradient(pos, int(node), int(comp), lattice, params)
            abs_err = abs(F_analytic - F_fd)
            ref = max(abs(F_fd), 1e-12)
            rel_err = abs_err / ref
            if abs_err > max_abs_err:
                max_abs_err = abs_err
                max_rel_err = rel_err
                worst = (node, comp, F_analytic, F_fd)

        assert max_abs_err < self.ATOL_FD, (
            f"Force–gradient mismatch on {grid_shape}: "
            f"max |F_analytic - F_fd| = {max_abs_err:.3e} "
            f"(node={worst[0]}, comp={worst[1]}, "
            f"F_analytic={worst[2]:.6g}, F_fd={worst[3]:.6g}), "
            f"rel_err={max_rel_err:.3e}, eps={self.EPS:.0e}"
        )
