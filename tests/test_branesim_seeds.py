"""Tests for branesim.initialization.seeds — VSH soliton seed menu.

Canonical setup: d=3, m=4 (three spacelike + one amplitude channel).
All seeds are pure initial-condition generators; no solver is invoked.

Falsification discriminator
---------------------------
The radial-lock metric distinguishes locked (VSH) seeds from the negative
control:

    radial_lock = sum_active (xi . x_hat)^2 / sum_active |xi|^2

For perfectly colour-locked seeds (hedgehog, Skyrme-twisted lateral part)
radial_lock = 1.0.  For the axis_triplet (independent scalars, isotropic
weights) the expected value is 1/dim ≈ 0.333 in d=3.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from branesim.core.conventions import LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.seeds import (
    hedgehog,
    skyrme_twisted_hedgehog,
    axis_triplet,
)


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def make_lattice(grid_shape, periodic=True, spacing=1.0) -> SpacelikeLattice:
    lp = LatticeParams(
        grid_shape=grid_shape,
        spacing=spacing,
        periodic_axes=tuple(periodic for _ in grid_shape),
    )
    return SpacelikeLattice(lp)


def radial_geometry(lattice: SpacelikeLattice, m: int):
    """Return dx, r, x_hat for a lattice."""
    dim = lattice.params.dim
    ref = lattice.reference_positions(m)
    coords = ref[:, :dim]
    centre = coords.mean(axis=0)
    dx = coords - centre
    r = np.linalg.norm(dx, axis=1)
    r_safe = np.where(r > 1e-30, r, 1.0)
    x_hat = dx / r_safe[:, None]
    x_hat[r <= 1e-30] = 0.0
    return dx, r, x_hat


def radial_lock(disp_lateral: np.ndarray, x_hat: np.ndarray, r: np.ndarray, eps_r: float = 1e-6):
    """Compute sum (xi . x_hat)^2 / sum |xi|^2 over nodes with r > eps_r."""
    active = r > eps_r
    xi = disp_lateral[active]                   # (n_active, dim)
    xh = x_hat[active]                          # (n_active, dim)
    proj = np.einsum("ni,ni->n", xi, xh)       # (n_active,)
    num = float(np.sum(proj ** 2))
    den = float(np.sum(xi ** 2))
    if den < 1e-30:
        return float("nan")
    return num / den


# ---------------------------------------------------------------------------
# d=3, m=4 canonical fixtures
# ---------------------------------------------------------------------------

GRID_D3 = (16, 16, 16)
M_D3 = 4
U0 = 0.05
W = 4.0    # half-width: a few lattice spacings inside a 16^3 box


# ---------------------------------------------------------------------------
# 1.  Hedgehog: colour-lock and amplitude channel
# ---------------------------------------------------------------------------

class TestHedgehog:
    """Colour-lock ≈ 1.0; amplitude channel exactly 0."""

    @pytest.mark.parametrize("profile_shape", ["gaussian", "sech", "power2"])
    def test_colour_lock_d3(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, meta = hedgehog(lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape)
        assert meta["ansatz"] == "hedgehog"
        ref = lattice.reference_positions(M_D3)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, M_D3)
        lock = radial_lock(disp[:, :3], x_hat, r)
        assert abs(lock - 1.0) < 1e-10, (
            f"hedgehog ({profile_shape}) radial_lock = {lock:.8f}, expected 1.0"
        )

    @pytest.mark.parametrize("profile_shape", ["gaussian", "sech", "power2"])
    def test_amplitude_channel_zero(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, _ = hedgehog(lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape)
        ref = lattice.reference_positions(M_D3)
        amp_disp = R0[:, 3] - ref[:, 3]
        assert float(np.max(np.abs(amp_disp))) < 1e-15, (
            f"hedgehog amplitude channel not zero: max = {np.max(np.abs(amp_disp)):.3e}"
        )

    def test_output_shape(self):
        lattice = make_lattice(GRID_D3)
        R0, _ = hedgehog(lattice, M_D3, u0=U0, w=W)
        assert R0.shape == (lattice.n_nodes, M_D3)
        assert R0.dtype == np.float64

    def test_metadata(self):
        lattice = make_lattice(GRID_D3)
        _, meta = hedgehog(lattice, M_D3, u0=U0, w=W)
        assert meta["J"] == 0
        assert meta["L"] == 1
        assert meta["B_winding"] == 0
        assert meta["colour_structure"] == "locked_to_x_hat"


# ---------------------------------------------------------------------------
# 2.  Skyrme-twisted hedgehog: winding, lateral lock, amplitude nonzero
# ---------------------------------------------------------------------------

class TestSkyrmeTwistedHedgehog:
    """Lateral colour-lock ≈ 1.0; amplitude channel encodes cos F; S^3 wrap."""

    @pytest.mark.parametrize("profile_shape", ["power2", "tanh"])
    def test_lateral_colour_lock(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, meta = skyrme_twisted_hedgehog(
            lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape
        )
        assert meta["B_winding"] == 1
        ref = lattice.reference_positions(M_D3)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, M_D3)
        lock = radial_lock(disp[:, :3], x_hat, r)
        assert abs(lock - 1.0) < 1e-10, (
            f"skyrme_twisted ({profile_shape}) lateral radial_lock = {lock:.8f}"
        )

    @pytest.mark.parametrize("profile_shape", ["power2", "tanh"])
    def test_amplitude_channel_nonzero(self, profile_shape):
        """Amplitude channel must be non-trivially nonzero (S^3 winding)."""
        lattice = make_lattice(GRID_D3)
        R0, _ = skyrme_twisted_hedgehog(
            lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape
        )
        ref = lattice.reference_positions(M_D3)
        amp_disp = R0[:, 3] - ref[:, 3]
        # At r -> 0: cos(F) = cos(pi) = -1 -> amp_disp ~ -u0
        # At large r: cos(F) = cos(0) = +1  -> amp_disp ~ +u0
        # So the range must span close to 2*u0.
        assert float(np.max(np.abs(amp_disp))) > 0.01 * U0, (
            f"Amplitude channel too small: max |amp| = {np.max(np.abs(amp_disp)):.3e}"
        )

    def test_s3_poles_power2(self):
        """cos(F) sweeps from -1 (r=0) to +1 (r=inf); verify at extremes."""
        lattice = make_lattice(GRID_D3)
        R0, _ = skyrme_twisted_hedgehog(
            lattice, M_D3, u0=U0, w=W, profile_shape="power2"
        )
        ref = lattice.reference_positions(M_D3)
        dim = 3
        amp_disp = R0[:, dim] - ref[:, dim]          # = u0 * cos(F)
        # Minimum ~ -u0 (south pole, r near 0) — should exist somewhere
        assert float(np.min(amp_disp)) < -0.5 * U0, (
            f"S^3 south pole not reached: min amp = {np.min(amp_disp):.4f}"
        )
        # Maximum ~ +u0 (north pole, r large) — should exist somewhere
        assert float(np.max(amp_disp)) > 0.5 * U0, (
            f"S^3 north pole not reached: max amp = {np.max(amp_disp):.4f}"
        )

    def test_requires_m_geq_dim_plus_1(self):
        """Must raise if m < dim+1."""
        lattice = make_lattice(GRID_D3)
        with pytest.raises(ValueError, match="dim\\+1"):
            skyrme_twisted_hedgehog(lattice, m=3, u0=U0, w=W)

    def test_metadata(self):
        lattice = make_lattice(GRID_D3)
        _, meta = skyrme_twisted_hedgehog(lattice, M_D3, u0=U0, w=W)
        assert meta["J"] == 0
        assert meta["L"] == 1
        assert meta["B_winding"] == 1
        assert meta["colour_structure"] == "locked_to_x_hat"


# ---------------------------------------------------------------------------
# 3.  Axis-triplet negative control: radial_lock ≈ 1/dim
# ---------------------------------------------------------------------------

class TestAxisTriplet:
    """Radial lock ≈ 0.333 (≠ 1.0) — falsification discriminator."""

    @pytest.mark.parametrize("profile_shape", ["gaussian", "sech", "power2"])
    def test_radial_lock_approx_one_over_dim(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, meta = axis_triplet(
            lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape
        )
        assert meta["role"] == "negative_control"
        ref = lattice.reference_positions(M_D3)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, M_D3)
        lock = radial_lock(disp[:, :3], x_hat, r)
        dim = 3
        expected = 1.0 / dim
        # Tolerance 5%: the lock should be close to 1/3, not close to 1.
        assert abs(lock - expected) < 0.05, (
            f"axis_triplet ({profile_shape}) radial_lock = {lock:.4f}, "
            f"expected ~ 1/3 = {expected:.4f}"
        )

    def test_amplitude_channel_zero(self):
        lattice = make_lattice(GRID_D3)
        R0, _ = axis_triplet(lattice, M_D3, u0=U0, w=W)
        ref = lattice.reference_positions(M_D3)
        amp_disp = R0[:, 3] - ref[:, 3]
        assert float(np.max(np.abs(amp_disp))) < 1e-15

    def test_metadata(self):
        lattice = make_lattice(GRID_D3)
        _, meta = axis_triplet(lattice, M_D3, u0=U0, w=W)
        assert meta["J"] == 1
        assert meta["L"] == 0
        assert meta["B_winding"] == 0
        assert meta["colour_structure"] == "unlocked_axis_scalars"
        assert meta["role"] == "negative_control"

    def test_distinguishes_from_hedgehog(self):
        """Hedgehog lock 1.0; axis_triplet lock 1/3: gap > 0.5."""
        lattice = make_lattice(GRID_D3)
        _, r, x_hat = radial_geometry(lattice, M_D3)
        ref = lattice.reference_positions(M_D3)

        R_hh, _ = hedgehog(lattice, M_D3, u0=U0, w=W)
        lock_hh = radial_lock((R_hh - ref)[:, :3], x_hat, r)

        R_at, _ = axis_triplet(lattice, M_D3, u0=U0, w=W)
        lock_at = radial_lock((R_at - ref)[:, :3], x_hat, r)

        assert lock_hh - lock_at > 0.5, (
            f"lock gap too small: hedgehog={lock_hh:.4f}, axis_triplet={lock_at:.4f}"
        )


# ---------------------------------------------------------------------------
# 4.  Dimension-agnostic: hedgehog at d=2
# ---------------------------------------------------------------------------

class TestDimensionAgnostic:
    """Seeds work for d=2 (and d=1 for seeds that don't require amplitude channel)."""

    def test_hedgehog_d2(self):
        lattice = make_lattice((16, 16))
        m = 3   # dim+1
        R0, meta = hedgehog(lattice, m, u0=U0, w=W)
        assert R0.shape == (lattice.n_nodes, m)
        ref = lattice.reference_positions(m)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, m)
        lock = radial_lock(disp[:, :2], x_hat, r)
        assert abs(lock - 1.0) < 1e-10, f"d=2 hedgehog lock = {lock:.8f}"

    def test_hedgehog_d1(self):
        lattice = make_lattice((16,))
        m = 2
        R0, meta = hedgehog(lattice, m, u0=U0, w=W)
        assert R0.shape == (lattice.n_nodes, m)
        # d=1: x_hat is ±1; xi = f(r)*x_hat; projection onto x_hat is f(r) -> lock=1
        ref = lattice.reference_positions(m)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, m)
        lock = radial_lock(disp[:, :1], x_hat, r)
        assert abs(lock - 1.0) < 1e-10, f"d=1 hedgehog lock = {lock:.8f}"

    def test_axis_triplet_d2_lock(self):
        """d=2 axis_triplet: lock ≈ 1/2."""
        lattice = make_lattice((16, 16))
        m = 3
        R0, _ = axis_triplet(lattice, m, u0=U0, w=W)
        ref = lattice.reference_positions(m)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, m)
        lock = radial_lock(disp[:, :2], x_hat, r)
        expected = 0.5   # 1/dim = 1/2 in d=2
        assert abs(lock - expected) < 0.05, (
            f"d=2 axis_triplet lock = {lock:.4f}, expected ≈ 0.5"
        )

    def test_skyrme_twisted_d2(self):
        lattice = make_lattice((16, 16))
        m = 3   # dim+1 = 3
        R0, meta = skyrme_twisted_hedgehog(lattice, m, u0=U0, w=W)
        assert R0.shape == (lattice.n_nodes, m)
        ref = lattice.reference_positions(m)
        disp = R0 - ref
        _, r, x_hat = radial_geometry(lattice, m)
        lock = radial_lock(disp[:, :2], x_hat, r)
        assert abs(lock - 1.0) < 1e-10, f"d=2 skyrme_twisted lock = {lock:.8f}"


# ---------------------------------------------------------------------------
# 5.  Profile shape properties: f(0)=1, monotone decay, finite everywhere
# ---------------------------------------------------------------------------

class TestProfileShapes:
    """Each profile gives f(0)=1, decays monotonically, is finite everywhere."""

    @pytest.mark.parametrize("profile_shape", ["gaussian", "sech", "power2"])
    def test_hedgehog_profile_decay(self, profile_shape):
        """Radial displacement amplitude decays from centre outward."""
        lattice = make_lattice(GRID_D3)
        R0, _ = hedgehog(lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape)
        ref = lattice.reference_positions(M_D3)
        disp = R0 - ref
        _, r, _ = radial_geometry(lattice, M_D3)
        # |xi| = u0 * f(r); should decrease as r increases
        xi_norm = np.linalg.norm(disp[:, :3], axis=1)   # (n_nodes,)

        # Bin by radius and check mean is decreasing
        r_min, r_max = r.min(), r.max()
        n_bins = 6
        bin_edges = np.linspace(r_min, r_max, n_bins + 1)
        bin_means = []
        for b in range(n_bins):
            mask = (r >= bin_edges[b]) & (r < bin_edges[b + 1])
            if mask.sum() > 0:
                bin_means.append(float(np.mean(xi_norm[mask])))

        # Must be monotone decreasing (allow small floating-point tolerance)
        for i in range(len(bin_means) - 1):
            assert bin_means[i] >= bin_means[i + 1] - 1e-10, (
                f"{profile_shape}: bin {i}={bin_means[i]:.6f} < bin {i+1}={bin_means[i+1]:.6f}"
            )

    @pytest.mark.parametrize("profile_shape", ["gaussian", "sech", "power2"])
    def test_hedgehog_profile_finite(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, _ = hedgehog(lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape)
        assert np.all(np.isfinite(R0)), f"{profile_shape}: non-finite values in R0"

    @pytest.mark.parametrize("profile_shape", ["power2", "tanh"])
    def test_skyrme_profile_finite(self, profile_shape):
        lattice = make_lattice(GRID_D3)
        R0, _ = skyrme_twisted_hedgehog(
            lattice, M_D3, u0=U0, w=W, profile_shape=profile_shape
        )
        assert np.all(np.isfinite(R0)), f"skyrme {profile_shape}: non-finite values"


# ---------------------------------------------------------------------------
# 6.  Amplitude scaling: max|displacement| proportional to u0
# ---------------------------------------------------------------------------

class TestAmplitudeScaling:
    """max|displacement| scales linearly with u0 (seeds are linear in u0)."""

    def _max_disp(self, R0, lattice, m):
        ref = lattice.reference_positions(m)
        return float(np.max(np.abs(R0 - ref)))

    @pytest.mark.parametrize("u0_a,u0_b", [(0.01, 0.05), (0.001, 0.02)])
    def test_hedgehog_scales_with_u0(self, u0_a, u0_b):
        lattice = make_lattice(GRID_D3)
        R0_a, _ = hedgehog(lattice, M_D3, u0=u0_a, w=W)
        R0_b, _ = hedgehog(lattice, M_D3, u0=u0_b, w=W)
        d_a = self._max_disp(R0_a, lattice, M_D3)
        d_b = self._max_disp(R0_b, lattice, M_D3)
        ratio = d_b / d_a
        expected = u0_b / u0_a
        assert abs(ratio - expected) / expected < 1e-10, (
            f"hedgehog: ratio={ratio:.6f}, expected {expected:.6f}"
        )

    @pytest.mark.parametrize("u0_a,u0_b", [(0.01, 0.05), (0.001, 0.02)])
    def test_skyrme_scales_with_u0(self, u0_a, u0_b):
        lattice = make_lattice(GRID_D3)
        R0_a, _ = skyrme_twisted_hedgehog(lattice, M_D3, u0=u0_a, w=W)
        R0_b, _ = skyrme_twisted_hedgehog(lattice, M_D3, u0=u0_b, w=W)
        d_a = self._max_disp(R0_a, lattice, M_D3)
        d_b = self._max_disp(R0_b, lattice, M_D3)
        ratio = d_b / d_a
        expected = u0_b / u0_a
        assert abs(ratio - expected) / expected < 1e-10, (
            f"skyrme: ratio={ratio:.6f}, expected {expected:.6f}"
        )

    @pytest.mark.parametrize("u0_a,u0_b", [(0.01, 0.05), (0.001, 0.02)])
    def test_axis_triplet_scales_with_u0(self, u0_a, u0_b):
        lattice = make_lattice(GRID_D3)
        R0_a, _ = axis_triplet(lattice, M_D3, u0=u0_a, w=W)
        R0_b, _ = axis_triplet(lattice, M_D3, u0=u0_b, w=W)
        d_a = self._max_disp(R0_a, lattice, M_D3)
        d_b = self._max_disp(R0_b, lattice, M_D3)
        ratio = d_b / d_a
        expected = u0_b / u0_a
        assert abs(ratio - expected) / expected < 1e-10, (
            f"axis_triplet: ratio={ratio:.6f}, expected {expected:.6f}"
        )


# ---------------------------------------------------------------------------
# 7.  run_experiment._build_seed integration
# ---------------------------------------------------------------------------

class TestBuildSeedIntegration:
    """Seeds are reachable from the experiment entrypoint _build_seed."""

    def test_hedgehog_via_build_seed(self):
        from branesim.run_experiment import _build_seed
        lattice = make_lattice(GRID_D3)
        cfg = {"kind": "hedgehog", "u0": 0.01, "w": 3.0, "profile_shape": "gaussian"}
        R0 = _build_seed(cfg, lattice, m=M_D3)
        assert R0.shape == (lattice.n_nodes, M_D3)

    def test_skyrme_twisted_via_build_seed(self):
        from branesim.run_experiment import _build_seed
        lattice = make_lattice(GRID_D3)
        # Skyrme-twisted profiles: "power2" or "tanh" (not "gaussian")
        cfg = {"kind": "skyrme_twisted", "u0": 0.01, "w": 3.0, "profile_shape": "power2"}
        R0 = _build_seed(cfg, lattice, m=M_D3)
        assert R0.shape == (lattice.n_nodes, M_D3)

    def test_axis_triplet_via_build_seed(self):
        from branesim.run_experiment import _build_seed
        lattice = make_lattice(GRID_D3)
        cfg = {"kind": "axis_triplet", "u0": 0.01, "w": 3.0, "profile_shape": "gaussian"}
        R0 = _build_seed(cfg, lattice, m=M_D3)
        assert R0.shape == (lattice.n_nodes, M_D3)
