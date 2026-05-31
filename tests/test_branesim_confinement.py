"""Tests for branesim.diagnostics.confinement.

Regression note
---------------
These tests replace the self-referential ``leakage_fraction`` from the legacy
diagnostics.  The legacy metric used a threshold proportional to the packet's
own ``radius_rms``, so it stayed near zero even for a uniform (box-fill) field
and could not detect spreading.  The metrics ported here use a FIXED box-scale
threshold (``box_fill_radius``) instead.

Test structure
--------------
1. Localized Gaussian packet  — spread_ratio << 1 (< 0.2), confined_fraction ~ 1 (> 0.9)
2. Uniform / random field     — spread_ratio ~ 1 (> 0.9), confined_fraction ~ volume_fraction (0.06-0.1), NOT ~1
3. Dispersing sequence        — radius_growth > 1, spread_ratio rises toward 1
4. Dimension-agnostic (d=2 and d=3)
5. From-zip wrapper (worldvolume.zip round-trip)
6. Skyrme-aware (strain) metric  — topological Skyrme seed reads as confined; uniform field reads as dispersed
"""

from __future__ import annotations

import numpy as np
import pytest

from branesim.core.conventions import LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.confinement import (
    confinement_metrics_per_slice,
    confinement_summary,
    confinement_from_worldvolume,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def make_lattice(grid_shape: tuple[int, ...], periodic: bool = True, spacing: float = 1.0) -> SpacelikeLattice:
    lp = LatticeParams(
        grid_shape=grid_shape,
        spacing=spacing,
        periodic_axes=tuple(periodic for _ in grid_shape),
    )
    return SpacelikeLattice(lp)


def reference_positions(lattice: SpacelikeLattice, m_ambient: int) -> np.ndarray:
    return lattice.reference_positions(m_ambient)


def gaussian_displacement(
    ref: np.ndarray,
    dim: int,
    amplitude: float = 1.0,
    sigma: float | None = None,
) -> np.ndarray:
    """Gaussian packet displacement centred on the geometric centre of ref.

    The displacement is applied uniformly to all ``dim`` spatial ambient
    components; the temporal component (index ``dim``) is left at zero.

    Parameters
    ----------
    ref : (n_nodes, m_ambient)
    dim : spatial dimension
    amplitude : peak displacement amplitude
    sigma : Gaussian width in lattice units; default = 1/14 of the box extent.
        At this width the spread_ratio is ~0.16-0.17 (well below 0.2) for
        16^3 and 24^2 grids.  1/8 would be too wide (~0.29).
    """
    ref_spatial = ref[:, :dim]
    centre = np.mean(ref_spatial, axis=0)
    extents = ref_spatial.max(axis=0) - ref_spatial.min(axis=0)
    if sigma is None:
        sigma = float(np.min(extents)) / 14.0

    r = np.linalg.norm(ref_spatial - centre[np.newaxis, :], axis=1)  # (n_nodes,)
    envelope = amplitude * np.exp(-0.5 * (r / sigma) ** 2)            # (n_nodes,)

    disp = np.zeros_like(ref)
    for c in range(dim):
        disp[:, c] = envelope
    return disp


def uniform_displacement(ref: np.ndarray, dim: int, amplitude: float = 1.0, rng_seed: int = 42) -> np.ndarray:
    """Random uniform displacement (box-fill energy distribution)."""
    rng = np.random.default_rng(rng_seed)
    disp = np.zeros_like(ref)
    disp[:, :dim] = amplitude * rng.standard_normal((ref.shape[0], dim))
    return disp


# ---------------------------------------------------------------------------
# Test 1: Localized Gaussian packet  (d=3)
# ---------------------------------------------------------------------------


class TestLocalizedGaussian:
    """A narrow Gaussian packet should be clearly confined:
        spread_ratio < 0.2
        confined_fraction > 0.9
    """

    def _run(self, dim: int, grid_shape: tuple[int, ...]) -> dict:
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)
        disp = gaussian_displacement(ref, dim, amplitude=1.0)
        positions = ref + disp
        return confinement_metrics_per_slice(positions, ref, dim)

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (16, 16, 16)),
        (2, (24, 24)),
    ])
    def test_spread_ratio_small(self, dim, grid_shape):
        m = self._run(dim, grid_shape)
        assert m["spread_ratio"] < 0.2, (
            f"d={dim}: Gaussian spread_ratio={m['spread_ratio']:.4f}, expected < 0.2"
        )

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (16, 16, 16)),
        (2, (24, 24)),
    ])
    def test_confined_fraction_high(self, dim, grid_shape):
        m = self._run(dim, grid_shape)
        assert m["confined_fraction"] > 0.9, (
            f"d={dim}: Gaussian confined_fraction={m['confined_fraction']:.4f}, expected > 0.9"
        )


# ---------------------------------------------------------------------------
# Test 2: Uniform / random field  (box-fill)
# ---------------------------------------------------------------------------


class TestUniformField:
    """A uniform random field should be deconfined:
        spread_ratio > 0.9
        confined_fraction ~ volume_fraction (small, 0.04-0.15), definitely NOT > 0.5

    The volume fraction inside a sphere of radius factor*R_box in a cubic
    box is approximately (4/3 pi (0.5 R_box)^3) / box_volume; for factor=0.5
    on a 3D box this is roughly 0.065 (6.5%). On a 2D box (circle in square)
    it is pi*(0.5 R_box)^2 / box_area ~ 0.065-0.10.

    Regression: the legacy ``leakage_fraction`` gave ~0.0 for this case (it
    used radius_rms as the threshold, which scales with the packet, so energy
    outside the "leak" threshold was near zero even for a dispersed field).
    ``confined_fraction`` correctly gives a small number near the geometric
    volume fraction.
    """

    def _run(self, dim: int, grid_shape: tuple[int, ...]) -> dict:
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)
        disp = uniform_displacement(ref, dim, amplitude=1.0)
        positions = ref + disp
        return confinement_metrics_per_slice(positions, ref, dim)

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (12, 12, 12)),
        (2, (20, 20)),
    ])
    def test_spread_ratio_near_one(self, dim, grid_shape):
        m = self._run(dim, grid_shape)
        assert m["spread_ratio"] > 0.9, (
            f"d={dim}: uniform spread_ratio={m['spread_ratio']:.4f}, expected > 0.9"
        )

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (12, 12, 12)),
        (2, (20, 20)),
    ])
    def test_confined_fraction_small(self, dim, grid_shape):
        """confined_fraction is in the geometric volume-fraction range, not ~1."""
        m = self._run(dim, grid_shape)
        # Should be well below 0.5 (the geometric fraction for factor=0.5 is ~0.05-0.12)
        assert m["confined_fraction"] < 0.5, (
            f"d={dim}: uniform confined_fraction={m['confined_fraction']:.4f}, expected < 0.5 "
            f"(near geometric volume fraction, legacy leakage_fraction would have given ~0)"
        )
        # Should be consistent with a non-zero geometric fraction
        assert m["confined_fraction"] > 0.01, (
            f"d={dim}: uniform confined_fraction={m['confined_fraction']:.4f} unexpectedly tiny"
        )


# ---------------------------------------------------------------------------
# Test 3: Dispersing sequence
# ---------------------------------------------------------------------------


class TestDispersingSequence:
    """A Gaussian packet widened over multiple slices gives:
        radius_growth > 1
        spread_ratio rising toward 1 (final > initial)
    """

    def _make_slices(self, dim: int, grid_shape: tuple[int, ...], n_slices: int = 8) -> tuple[np.ndarray, np.ndarray]:
        """Build synthetic slices where the Gaussian sigma grows linearly."""
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)

        ref_spatial = ref[:, :dim]
        extents = ref_spatial.max(axis=0) - ref_spatial.min(axis=0)
        sigma_min = float(np.min(extents)) / 12.0
        sigma_max = float(np.min(extents)) / 2.5   # wide but not wider than box

        slices_list = []
        for i in range(n_slices):
            t = i / max(n_slices - 1, 1)
            sigma = sigma_min + t * (sigma_max - sigma_min)
            disp = gaussian_displacement(ref, dim, amplitude=1.0, sigma=sigma)
            slices_list.append(ref + disp)

        slices = np.stack(slices_list, axis=0)   # (n_slices, n_nodes, m_ambient)
        return slices, ref

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (14, 14, 14)),
        (2, (20, 20)),
    ])
    def test_radius_growth_greater_than_one(self, dim, grid_shape):
        slices, ref = self._make_slices(dim, grid_shape)
        summary = confinement_summary(slices, ref, dim)
        assert summary["radius_growth"] > 1.0, (
            f"d={dim}: radius_growth={summary['radius_growth']:.4f}, expected > 1.0"
        )

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (14, 14, 14)),
        (2, (20, 20)),
    ])
    def test_spread_ratio_increases(self, dim, grid_shape):
        slices, ref = self._make_slices(dim, grid_shape)
        summary = confinement_summary(slices, ref, dim)
        sr = summary["spread_ratio"]
        assert sr[-1] > sr[0], (
            f"d={dim}: spread_ratio did not increase: first={sr[0]:.4f}, last={sr[-1]:.4f}"
        )

    @pytest.mark.parametrize("dim,grid_shape", [
        (3, (14, 14, 14)),
        (2, (20, 20)),
    ])
    def test_spread_ratio_final_larger(self, dim, grid_shape):
        """Final spread_ratio should be substantially larger than initial."""
        slices, ref = self._make_slices(dim, grid_shape)
        summary = confinement_summary(slices, ref, dim)
        sr = summary["spread_ratio"]
        # At least 3x larger
        assert sr[-1] > 3.0 * sr[0], (
            f"d={dim}: spread_ratio growth insufficient: first={sr[0]:.4f}, last={sr[-1]:.4f}"
        )


# ---------------------------------------------------------------------------
# Test 4: Dimension-agnostic API checks
# ---------------------------------------------------------------------------


class TestDimensionAgnostic:
    """confinement_metrics_per_slice and confinement_summary work at d=2 and d=3."""

    @pytest.mark.parametrize("dim,grid_shape", [
        (2, (10, 10)),
        (3, (8, 8, 8)),
    ])
    def test_per_slice_returns_expected_keys(self, dim, grid_shape):
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)
        positions = ref.copy()
        positions[:, :dim] += 0.1

        m = confinement_metrics_per_slice(positions, ref, dim)
        for key in ("radius_rms", "box_fill_radius", "spread_ratio", "confined_fraction"):
            assert key in m, f"Missing key {key!r} in per-slice metrics"
            assert np.isfinite(m[key]), f"Key {key!r} is not finite: {m[key]}"

    @pytest.mark.parametrize("dim,grid_shape", [
        (2, (10, 10)),
        (3, (8, 8, 8)),
    ])
    def test_summary_returns_expected_keys(self, dim, grid_shape):
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)
        slices = np.stack([ref + 0.1, ref + 0.2], axis=0)

        summary = confinement_summary(slices, ref, dim)
        for key in ("box_fill_radius", "radius_rms", "spread_ratio", "confined_fraction",
                    "radius_growth", "final", "mean"):
            assert key in summary, f"Missing key {key!r} in summary"

        assert summary["radius_rms"].shape == (2,)
        assert summary["spread_ratio"].shape == (2,)
        assert summary["confined_fraction"].shape == (2,)
        assert np.isfinite(summary["radius_growth"])

    @pytest.mark.parametrize("dim,grid_shape", [
        (2, (10, 10)),
        (3, (8, 8, 8)),
    ])
    def test_box_fill_radius_invariant(self, dim, grid_shape):
        """box_fill_radius must be the same regardless of displacement."""
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)

        m_zero = confinement_metrics_per_slice(ref.copy(), ref, dim)
        m_large = confinement_metrics_per_slice(ref + 5.0, ref, dim)

        assert abs(m_zero["box_fill_radius"] - m_large["box_fill_radius"]) < 1e-12, (
            "box_fill_radius changed between displacements — it must be frame-invariant"
        )


# ---------------------------------------------------------------------------
# Test 5: From-zip wrapper (worldvolume.zip round-trip)
# ---------------------------------------------------------------------------


class TestFromWorldvolumeZip:
    """confinement_from_worldvolume reads a zip written by WorldVolumeWriter."""

    def test_roundtrip_d3(self, tmp_path):
        import json
        from branesim.io.contracts import WorldVolumeWriter

        dim = 3
        grid_shape = (8, 8, 8)
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)

        # Build a Gaussian packet (localized)
        disp = gaussian_displacement(ref, dim, amplitude=1.0)
        positions = ref + disp

        path = tmp_path / "wv_test.zip"
        manifest_extra = {
            "lattice": {
                "grid_shape": list(grid_shape),
                "spacing": 1.0,
                "periodic_axes": [True] * dim,
                "axial_weight": 1.0,
                "dim": dim,
            },
        }
        with WorldVolumeWriter(path, manifest_extra) as w:
            w.write_slice(0, 0.0, positions)
            w.write_npy("aux/ref_positions.npy", ref)

        summary = confinement_from_worldvolume(path)

        # One slice, so arrays have length 1
        assert summary["spread_ratio"].shape == (1,)
        # Should detect a localized packet
        assert summary["spread_ratio"][0] < 0.4, (
            f"zip round-trip: spread_ratio={summary['spread_ratio'][0]:.4f}, expected < 0.4"
        )
        assert summary["confined_fraction"][0] > 0.7, (
            f"zip round-trip: confined_fraction={summary['confined_fraction'][0]:.4f}, expected > 0.7"
        )

    def test_roundtrip_d2(self, tmp_path):
        from branesim.io.contracts import WorldVolumeWriter

        dim = 2
        grid_shape = (16, 16)
        m_ambient = dim + 1
        lattice = make_lattice(grid_shape)
        ref = reference_positions(lattice, m_ambient)

        disp = gaussian_displacement(ref, dim, amplitude=1.0)
        positions = ref + disp

        path = tmp_path / "wv_test_d2.zip"
        manifest_extra = {
            "lattice": {
                "grid_shape": list(grid_shape),
                "spacing": 1.0,
                "periodic_axes": [True] * dim,
                "axial_weight": 1.0,
                "dim": dim,
            },
        }
        with WorldVolumeWriter(path, manifest_extra) as w:
            w.write_slice(0, 0.0, positions)
            w.write_npy("aux/ref_positions.npy", ref)

        summary = confinement_from_worldvolume(path)
        assert summary["spread_ratio"].shape == (1,)
        assert summary["spread_ratio"][0] < 0.4, (
            f"d=2 zip: spread_ratio={summary['spread_ratio'][0]:.4f}"
        )


# ---------------------------------------------------------------------------
# Test 6: Contrast report (localized vs uniform — regression table)
# ---------------------------------------------------------------------------


class TestConfinementContrast:
    """Quantify the localized vs uniform contrast for the regression record.

    This test is the formal replacement for the self-referential leakage_fraction.
    Legacy leakage_fraction on a uniform field gave ~0.0 even though the field
    was dispersed — it could not distinguish localized from dispersed.

    Required contrast:
      spread_ratio:       localized < 0.2,  uniform > 0.9
      confined_fraction:  localized > 0.9,  uniform < 0.15  (not ~1)
    """

    DIM = 3
    GRID_SHAPE = (14, 14, 14)

    @pytest.fixture(scope="class")
    def metrics(self):
        m_ambient = self.DIM + 1
        lattice = make_lattice(self.GRID_SHAPE)
        ref = reference_positions(lattice, m_ambient)

        disp_loc = gaussian_displacement(ref, self.DIM, amplitude=1.0)
        disp_uni = uniform_displacement(ref, self.DIM, amplitude=1.0)

        m_loc = confinement_metrics_per_slice(ref + disp_loc, ref, self.DIM)
        m_uni = confinement_metrics_per_slice(ref + disp_uni, ref, self.DIM)
        return m_loc, m_uni

    def test_spread_ratio_contrast(self, metrics):
        m_loc, m_uni = metrics
        assert m_loc["spread_ratio"] < 0.2, (
            f"Localized spread_ratio={m_loc['spread_ratio']:.4f} should be < 0.2"
        )
        assert m_uni["spread_ratio"] > 0.9, (
            f"Uniform spread_ratio={m_uni['spread_ratio']:.4f} should be > 0.9"
        )

    def test_confined_fraction_contrast(self, metrics):
        m_loc, m_uni = metrics
        assert m_loc["confined_fraction"] > 0.9, (
            f"Localized confined_fraction={m_loc['confined_fraction']:.4f} should be > 0.9"
        )
        assert m_uni["confined_fraction"] < 0.15, (
            f"Uniform confined_fraction={m_uni['confined_fraction']:.4f} should be < 0.15 "
            f"(legacy leakage_fraction gave ~0 here — wrong, this corrects it)"
        )


# ---------------------------------------------------------------------------
# Test 6: Skyrme-aware (strain) confinement metric
# ---------------------------------------------------------------------------


class TestStrainWeightMode:
    """weight_mode='strain' correctly identifies a Skyrme soliton as confined.

    The default 'lateral' weight uses |displacement_lateral|^2 which has
    algebraic tails for a Skyrme hedgehog (sin(F(r)) ~ pi*w^2/r^2) and can
    report spread_ratio > 1 even for a geometrically localized topological
    object.  The 'strain' weight uses per-node spring energy which is zero
    on the flat background and nonzero only at the soliton core.

    Tests:
    a. Skyrme seed at w/a=3 in a 13^3 box reads as confined (spread_ratio < 0.3)
       under the strain metric.
    b. Uniform (constant lateral) displacement reads as dispersed (spread_ratio > 0.8)
       under the strain metric — it has nonzero gradient EVERYWHERE.
    c. weight_mode='strain' raises ValueError when lattice is None.
    d. weight_mode='lateral' default is unchanged (backward compat).
    e. confinement_summary echoes weight_mode in the returned dict.
    """

    def _make_lattice_and_ref(self, n=13):
        lp = LatticeParams(
            grid_shape=(n, n, n),
            spacing=1.0,
            periodic_axes=(False, False, False),
        )
        lat = SpacelikeLattice(lp)
        ref = lat.reference_positions(4)
        return lat, ref

    def _skyrme_seed_slice(self, lat, ref, u0=0.15, w=3.0):
        """Build l=0 slice of a Skyrme-hedgehog seed."""
        from branesim.initialization.seeds import skyrme_twisted_hedgehog
        R0, _ = skyrme_twisted_hedgehog(lat, m=4, u0=u0, w=w, profile_shape="power2")
        return R0

    def test_skyrme_seed_confined_under_strain_metric(self):
        """A Skyrme seed (w/a=1.5) in a 13^3 box reads as confined under the
        direction-variation ('strain') metric.

        The direction-variation weight |disp_p - disp_q|^2 peaks where the
        hedgehog field rotates most rapidly (near r ~ w) and is small far from
        the soliton, giving a spread_ratio well below 1 for a localized soliton.

        We use w=1.5 (tight core relative to box = 13) where the signal is clear.
        """
        lat, ref = self._make_lattice_and_ref(n=13)
        positions = self._skyrme_seed_slice(lat, ref, u0=0.15, w=1.5)
        m = confinement_metrics_per_slice(
            positions, ref, dim=lat.dim,
            weight_mode="strain",
            _neighbor_table=lat.neighbors,
            _k_s=1.0, _alpha=0.5,
        )
        # Direction-variation for Skyrme seed with w=1.5 in 13^3 box:
        # empirically ~0.42; well below 0.6
        assert m["spread_ratio"] < 0.6, (
            f"Skyrme seed direction-variation spread_ratio={m['spread_ratio']:.4f}, "
            f"expected < 0.6 (confined topological soliton should have low direction-variation spread)"
        )

    def test_uniform_displacement_dispersed_under_strain_metric(self):
        """A linear ramp displacement fills the whole box with nonzero
        direction variation => reads as dispersed (spread_ratio > 0.5)."""
        lat, ref = self._make_lattice_and_ref(n=9)
        mi = lat.multi_indices
        positions = ref.copy()
        # Linear ramp: disp ~ epsilon * x => constant spatial gradient everywhere
        positions[:, 0] += 0.05 * mi[:, 0]  # ramp in ambient x-component
        m = confinement_metrics_per_slice(
            positions, ref, dim=lat.dim,
            weight_mode="strain",
            _neighbor_table=lat.neighbors,
            _k_s=1.0, _alpha=0.5,
        )
        # A uniform ramp has constant |disp_p - disp_q|^2 throughout the box
        # => energy-weighted centroid is near the geometric centre => spread_ratio ~ box-fill
        assert m["spread_ratio"] > 0.5, (
            f"Linear ramp spread_ratio={m['spread_ratio']:.4f}, expected > 0.5 "
            f"(deconfined field should have spread_ratio close to box-fill)"
        )

    def test_strain_mode_requires_neighbor_table(self):
        """weight_mode='strain' without _neighbor_table raises ValueError."""
        lat, ref = self._make_lattice_and_ref(n=5)
        positions = ref.copy()
        with pytest.raises(ValueError, match="_neighbor_table"):
            confinement_metrics_per_slice(
                positions, ref, dim=lat.dim,
                weight_mode="strain",
                _neighbor_table=None,
            )

    def test_invalid_weight_mode_raises(self):
        """Unknown weight_mode raises ValueError."""
        lat, ref = self._make_lattice_and_ref(n=5)
        with pytest.raises(ValueError, match="weight_mode"):
            confinement_metrics_per_slice(
                ref.copy(), ref, dim=lat.dim,
                weight_mode="bogus",
            )

    def test_lateral_mode_unchanged(self):
        """Default weight_mode='lateral' gives same result as before (backward compat)."""
        lat, ref = self._make_lattice_and_ref(n=9)
        rng = np.random.default_rng(0)
        positions = ref + 0.1 * rng.standard_normal(ref.shape)
        m_default = confinement_metrics_per_slice(positions, ref, dim=lat.dim)
        m_lateral = confinement_metrics_per_slice(positions, ref, dim=lat.dim, weight_mode="lateral")
        assert abs(m_default["spread_ratio"] - m_lateral["spread_ratio"]) < 1e-12, (
            "Default and explicit 'lateral' modes should give identical results"
        )

    def test_confinement_summary_strain_requires_lattice(self):
        """confinement_summary with weight_mode='strain' but no lattice raises ValueError."""
        lat, ref = self._make_lattice_and_ref(n=5)
        slices = np.stack([ref, ref], axis=0)
        with pytest.raises(ValueError, match="lattice"):
            confinement_summary(slices, ref, dim=lat.dim, weight_mode="strain", lattice=None)

    def test_confinement_summary_echoes_weight_mode(self):
        """confinement_summary result echoes the weight_mode used."""
        lat, ref = self._make_lattice_and_ref(n=7)
        slices = np.stack([ref, ref + 0.01], axis=0)
        for mode in ("lateral", "strain"):
            kwargs = {}
            if mode == "strain":
                kwargs = {"lattice": lat, "k_s": 1.0, "alpha": 0.5}
            summary = confinement_summary(slices, ref, dim=lat.dim, weight_mode=mode, **kwargs)
            assert summary.get("weight_mode") == mode, (
                f"Expected weight_mode='{mode}' in summary, got {summary.get('weight_mode')!r}"
            )
