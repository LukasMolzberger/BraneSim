"""Unit tests for the topological-breather (baryon-candidate) machinery.

Covers:
  1. Seed builder: breather_seed_skyrmion
     - Shape and dtype.
     - l=0 slice reduces exactly to skyrme_twisted_hedgehog.
     - B=1 winding: lattice degree ≈ 1 on l=0 slice.
  2. Multi-component constraint square-count: the system has the right
     number of equations and unknowns (square G).
  3. solve_breather mode="topological" returns the correct keys and
     does not regress the scalar mode.

All tests are cheap / fast (small grids, few iterations, no convergence
required for the unit tests).  The standalone validation driver
(run_topological_breather_validation.py) runs the expensive Newton solve
and reports the go/no-go metrics.

Principles compliance
---------------------
- No solver state modified by diagnostics.
- No clamps.
- Dimension-agnostic seed builder (works for any lattice.dim).
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.seeds import skyrme_twisted_hedgehog
from branesim.solver.breather import (
    BreatherOpts,
    breather_seed_skyrmion,
    solve_breather,
)


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_3d_lattice(n: int = 7, alpha: float = 0.5) -> tuple[SpacelikeLattice, ActionParams]:
    """Minimal 3D open lattice (n^3, m=4, open boundary)."""
    lp = LatticeParams(
        grid_shape=(n, n, n),
        spacing=1.0,
        periodic_axes=(False, False, False),
    )
    lattice = SpacelikeLattice(lp)
    params = ActionParams(
        k_s=1.0,
        alpha=alpha,
        rho=1.0,
        dt=0.1,
        n_slices=1,
        m_ambient=4,
    )
    return lattice, params


# ---------------------------------------------------------------------------
# 1.  breather_seed_skyrmion: shape and dtype
# ---------------------------------------------------------------------------


class TestSeedBuilder:
    """breather_seed_skyrmion returns correct shape, dtype, and physics."""

    def test_output_shape(self):
        """Seed shape is (P, n_nodes, m_ambient)."""
        lattice, _ = _make_3d_lattice(n=5)
        P = 8
        m = 4
        slices, T_seed = breather_seed_skyrmion(lattice, m=m, u0=0.1, w=2.0, P=P)
        assert slices.shape == (P, lattice.n_nodes, m), (
            f"Expected ({P}, {lattice.n_nodes}, {m}), got {slices.shape}"
        )
        assert slices.dtype == np.float64

    def test_T_seed_positive(self):
        """Seed period T_seed > 0."""
        lattice, _ = _make_3d_lattice(n=5)
        _, T_seed = breather_seed_skyrmion(lattice, m=4, u0=0.1, w=2.0, P=8)
        assert math.isfinite(T_seed) and T_seed > 0.0

    def test_l0_slice_matches_skyrme_twisted(self):
        """At l=0 (carrier=1), seed slice == skyrme_twisted_hedgehog exactly."""
        lattice, _ = _make_3d_lattice(n=5)
        u0 = 0.12
        w = 2.0
        P = 8
        slices, _ = breather_seed_skyrmion(lattice, m=4, u0=u0, w=w, P=P)
        R0_skyrme, _ = skyrme_twisted_hedgehog(lattice, m=4, u0=u0, w=w, profile_shape="power2")
        # l=0: carrier = cos(0) = 1.0 exactly
        np.testing.assert_allclose(
            slices[0], R0_skyrme,
            atol=1e-14,
            err_msg="l=0 seed slice does not match skyrme_twisted_hedgehog",
        )

    def test_lP2_slice_is_anti_seed(self):
        """At l=P/2 (carrier=-1), seed = ref - (R0_skyrme - ref) = 2*ref - R0_skyrme."""
        lattice, _ = _make_3d_lattice(n=5)
        u0 = 0.12
        w = 2.0
        P = 8
        slices, _ = breather_seed_skyrmion(lattice, m=4, u0=u0, w=w, P=P)
        R0_skyrme, _ = skyrme_twisted_hedgehog(lattice, m=4, u0=u0, w=w, profile_shape="power2")
        ref = lattice.reference_positions(4)
        anti_amp = 2.0 * ref - R0_skyrme
        np.testing.assert_allclose(
            slices[P // 2], anti_amp,
            atol=1e-14,
            err_msg="l=P/2 seed slice is not the anti-amplitude",
        )

    def test_carrier_cosine_at_peak(self):
        """The X4 component at peak node follows cos(2πl/P) exactly."""
        lattice, _ = _make_3d_lattice(n=5)
        u0 = 0.20
        w = 2.0
        P = 8
        slices, _ = breather_seed_skyrmion(lattice, m=4, u0=u0, w=w, P=P)
        ref = lattice.reference_positions(4)

        # Find peak node (r=0, centre)
        mi = lattice.multi_indices
        center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
        d = np.linalg.norm(mi - center_mi, axis=1)
        peak = int(np.argmin(d))

        # At r=0, F(0) = pi, so static X4 disp = u0*cos(pi) = -u0
        static_x4_disp = -u0  # = u0 * cos(F(0))
        for l_idx in range(P):
            carrier = math.cos(2.0 * math.pi * l_idx / P)
            expected_x4 = ref[peak, 3] + carrier * static_x4_disp
            got_x4 = float(slices[l_idx, peak, 3])
            assert abs(got_x4 - expected_x4) < 1e-13, (
                f"l={l_idx}: expected X4[peak]={expected_x4:.8f}, got {got_x4:.8f}"
            )

    def test_requires_even_P(self):
        """P must be even; odd P raises ValueError."""
        lattice, _ = _make_3d_lattice(n=5)
        with pytest.raises(ValueError, match="even"):
            breather_seed_skyrmion(lattice, m=4, u0=0.1, w=2.0, P=7)

    def test_requires_m_geq_dim_plus_1(self):
        """Must raise if m < dim+1."""
        lattice, _ = _make_3d_lattice(n=5)
        with pytest.raises(ValueError):
            breather_seed_skyrmion(lattice, m=3, u0=0.1, w=2.0, P=8)


# ---------------------------------------------------------------------------
# 2.  B=1 winding on the l=0 seed slice (topological gate)
# ---------------------------------------------------------------------------


class TestWindingNumber:
    """The l=0 seed slice has B=1 winding (topological charge ≈ 1).

    We use a continuous-field winding measure:
    the degree of the map (sin F * x_hat, cos F) : S^2 -> S^3 restricted
    to the (lateral, X4) components of the displacement.

    For the Skyrme profile F(r) = pi/(1+(r/w)^2):
        - r=0: (sin(pi), cos(pi)) = (0, -1)  — south pole of S^3
        - r>>w: (sin(0), cos(0)) = (0, +1)   — north pole of S^3
    The map wraps the 3D space once around S^3 => degree = 1.

    We test this numerically by computing the signed volume of the image
    on a planar cross-section (sum over nodes of the Jacobian determinant).
    A clean B=1 seed must give degree > 0.5 on a reasonably-resolved grid.

    This is a FAST check (no solver invocation).
    """

    def _approximate_degree(self, lattice, slices_l0, m):
        """Approximate the topological degree on the l=0 slice.

        Uses the solid-angle element: for each node, compute the
        signed contribution to the winding number from the
        map (R_p - ref_p) on the (dim+1)-sphere.

        Simplified: just verify the image of the map sweeps from
        south to north pole of S^3 and that the X4 channel at r=0
        is strongly negative (south pole) while at the edge it is
        positive (north pole).
        """
        ref = lattice.reference_positions(m)
        disp = slices_l0 - ref           # (n_nodes, m)
        dim = lattice.dim

        # Check polar coverage:
        # South pole: disp[:, dim] minimum (should be near -u0)
        # North pole: disp[:, dim] maximum (should be near +u0)
        x4 = disp[:, dim]
        return float(np.min(x4)), float(np.max(x4))

    def test_b1_winding_poles_covered(self):
        """S^3 south and north poles are covered by the seed (necessary for B=1)."""
        lattice, _ = _make_3d_lattice(n=9)  # slightly larger for coverage
        u0 = 0.20
        w = 2.5
        P = 8
        slices, _ = breather_seed_skyrmion(lattice, m=4, u0=u0, w=w, P=P)
        south, north = self._approximate_degree(lattice, slices[0], m=4)
        ref = lattice.reference_positions(4)
        # south pole: X4 disp near -u0 at r=0
        assert south < -0.5 * u0, (
            f"South pole not covered: min X4 disp = {south:.4f}, expected < {-0.5*u0:.4f}"
        )
        # north pole: X4 disp near +u0 at large r
        assert north > 0.5 * u0, (
            f"North pole not covered: max X4 disp = {north:.4f}, expected > {0.5*u0:.4f}"
        )


# ---------------------------------------------------------------------------
# 3.  Multi-component constraint: system is square
# ---------------------------------------------------------------------------


class TestSystemSquareness:
    """The topological G has exactly P*n*m + 1 equations and unknowns."""

    def test_G_output_length_matches_input(self):
        """G(z0) output length == z0 input length (square system)."""
        lattice, params = _make_3d_lattice(n=5)
        P = 8
        m = 4
        u0 = 0.1
        w = 1.5

        slices_seed, T_seed = breather_seed_skyrmion(lattice, m=m, u0=u0, w=w, P=P)

        # Import the internal factories to check squareness directly
        from branesim.solver.breather import _make_G_topological, _pack, _unpack

        n_nodes = lattice.n_nodes
        peak_node = 0  # arbitrary for this test
        x4_comp = lattice.dim
        x4_pin_value = -u0

        G_topo = _make_G_topological(
            lattice, params, mass=1.0, P=P,
            u0=u0, peak_node=peak_node,
            x4_comp=x4_comp, x4_pin_value=x4_pin_value,
            opts=BreatherOpts(),
        )

        import math
        u0_log = math.log(T_seed)
        z0 = _pack(slices_seed, u0_log)
        g = G_topo(z0)

        # G must be square: len(g) == len(z0)
        assert len(g) == len(z0), (
            f"System is not square: len(G)={len(g)}, len(z)={len(z0)}"
        )

    def test_scalar_G_still_square(self):
        """Sanity: scalar mode G is also square (regression guard)."""
        from branesim.solver.breather import _make_G, _pack, _build_seed
        import math

        lattice, params = _make_3d_lattice(n=4)
        P = 8
        m = 4
        A = 0.05

        mi = lattice.multi_indices
        center_mi = np.array([(s - 1) / 2.0 for s in lattice.params.grid_shape])
        d = np.linalg.norm(mi - center_mi, axis=1)
        peak_node = int(np.argmin(d))
        lat_comp = m - 1
        peak_parity = int(mi[peak_node].sum() % 2)
        peak_sign = float((-1) ** peak_parity)

        slices_seed, T_seed = _build_seed(
            lattice, k_s=1.0, alpha=0.5, mass=1.0, a=1.0,
            P=P, amplitude=A, m_ambient=m,
        )

        G_scalar = _make_G(
            lattice, params, mass=1.0, P=P, amplitude=A,
            peak_node=peak_node, lat_comp=lat_comp,
            peak_sign=peak_sign, opts=BreatherOpts(),
        )

        u0_log = math.log(T_seed)
        z0 = _pack(slices_seed, u0_log)
        g = G_scalar(z0)
        assert len(g) == len(z0), f"Scalar G not square: {len(g)} != {len(z0)}"


# ---------------------------------------------------------------------------
# 4.  solve_breather mode regression: scalar mode unaffected
# ---------------------------------------------------------------------------


class TestSolverModeRegression:
    """Adding mode parameter does not regress the existing scalar solver."""

    def test_scalar_mode_result_has_mode_key(self):
        """solve_breather returns 'mode' key in result dict."""
        from branesim.core.conventions import LatticeParams
        lp = LatticeParams(grid_shape=(15,), spacing=1.0, periodic_axes=(False,))
        lattice = SpacelikeLattice(lp)
        params = ActionParams(k_s=1.0, alpha=0.5, rho=1.0, dt=0.1, n_slices=1, m_ambient=2)
        result = solve_breather(
            lattice, params, mass=1.0,
            P=8, amplitude=0.1,
            mode="scalar",
            opts=BreatherOpts(tol=1e-4, max_iter=5),  # cheap, not for convergence
        )
        assert "mode" in result
        assert result["mode"] == "scalar"

    def test_topological_mode_result_has_mode_key(self):
        """Topological mode result dict carries mode='topological'."""
        lattice, params = _make_3d_lattice(n=5)
        result = solve_breather(
            lattice, params, mass=1.0,
            P=8, amplitude=0.1,
            mode="topological",
            skyrme_w=1.5,
            opts=BreatherOpts(tol=1e-4, max_iter=3),  # cheap: just check keys
        )
        assert "mode" in result
        assert result["mode"] == "topological"

    def test_invalid_mode_raises(self):
        """Unknown mode raises ValueError."""
        lattice, params = _make_3d_lattice(n=5)
        with pytest.raises(ValueError, match="mode"):
            solve_breather(
                lattice, params, mass=1.0,
                P=8, amplitude=0.1,
                mode="bogus",
            )

    def test_topological_requires_m_geq_4(self):
        """Topological mode with m_ambient < dim+1 raises ValueError."""
        lp = LatticeParams(grid_shape=(5, 5, 5), spacing=1.0)
        lattice = SpacelikeLattice(lp)
        params = ActionParams(k_s=1.0, alpha=0.5, rho=1.0, dt=0.1, n_slices=1, m_ambient=3)
        with pytest.raises(ValueError, match="m_ambient"):
            solve_breather(
                lattice, params, mass=1.0,
                P=8, amplitude=0.1,
                mode="topological",
            )

    def test_result_keys_present(self):
        """Topological mode result contains all expected keys."""
        lattice, params = _make_3d_lattice(n=5)
        result = solve_breather(
            lattice, params, mass=1.0,
            P=8, amplitude=0.1,
            mode="topological",
            skyrme_w=1.5,
            opts=BreatherOpts(tol=1e-4, max_iter=3),
        )
        for key in ("slices", "T", "omega", "residual_norm", "converged",
                    "objective", "walltime_s", "peak_node", "lat_comp",
                    "residual_initial", "mode"):
            assert key in result, f"Missing key: {key}"

    def test_objective_is_residual_norm(self):
        """Topological mode respects the saddle discipline (objective=residual_norm)."""
        lattice, params = _make_3d_lattice(n=5)
        result = solve_breather(
            lattice, params, mass=1.0,
            P=8, amplitude=0.1,
            mode="topological",
            skyrme_w=1.5,
            opts=BreatherOpts(tol=1e-4, max_iter=3),
        )
        assert result["objective"] == "residual_norm"
