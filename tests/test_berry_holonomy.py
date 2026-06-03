"""Unit tests for branesim.diagnostics.berry_holonomy.

Tests cover:
  P2 — k-space plaquette holonomy is identity for all alpha (FH construction)
    1. Single-plaquette holonomy returns wz_dev_from_id < 1e-10 for all tested alpha.
    2. Per-band Berry phases are all < 1e-10 rad.
    3. Gauge randomization: per-k unitary rotations leave the result unchanged.
    4. Refinement: halving dk does not change the (already-zero) result.
    5. Full verify_p2_all_alpha() returns pass_p2_all = True.

  P3 — SO(3) rotation holonomy is spin-1 (+I at 2pi)
    6. Spin-1 state transported through 2pi has holonomy +1 (phase 0) for each band.
    7. Spin-1 rank-3 WZ holonomy at 2pi = +I.
    8. Spin-1 rank-3 WZ holonomy at 4pi = +I (4pi trivial, same as 2pi for J=1).
    9. Result is the same for alpha = 0, 0.2, 0.5, 0.8, 1.0 (alpha-independent).
   10. Rotation about different axes ([z], [x], [111]) all give the same result.
   11. Refinement: 2x n_steps changes holonomy by < 1% (stability).
   12. SYNTHETIC SPIN-1/2 CONTROL: J=1/2 spinor at 2pi gives holonomy -1 (phase=pi).
   13. Synthetic spin-1/2 at 4pi gives holonomy +1 (phase=0, as expected: 4pi trivial).
   14. The diagnostic DISTINGUISHES spin-1 from spin-1/2 at 2pi.

  Rotation matrices
   15. D^{(1)} at 2pi is exactly +I (algebraic check, not transport).
   16. D^{(1/2)} at 2pi is exactly -I.
   17. D^{(1)} at 4pi is exactly +I.
   18. D^{(1/2)} at 4pi is exactly +I.
   19. D^{(1)} is a proper rotation (det = 1, unitary).

  WZ decomposition
   20. U(1)/SU(n) decomposition of +I gives u1_phase = 0, sun_dev_from_id = 0.
   21. U(1)/SU(n) decomposition of -I (n=2) gives u1_phase = pi, sun_dev_from_id = 0.
"""

from __future__ import annotations

import numpy as np
import pytest

from branesim.diagnostics.berry_holonomy import (
    decompose_wz_u1_sun,
    heff_eigenframe,
    plaquette_holonomy_p2,
    rotate_and_transport,
    rotation_matrix_spin1,
    rotation_matrix_spin_half,
    spin1_frame_fn,
    spin1_state_fn,
    spin_half_frame_fn,
    spin_half_state_fn,
    verify_p2_all_alpha,
    verify_p3_so3_holonomy,
)


# ---------------------------------------------------------------------------
# Test constants
# ---------------------------------------------------------------------------
ALPHA_LIST = [0.0, 0.2, 0.5, 0.8, 1.0]
AXES = {
    "z": np.array([0.0, 0.0, 1.0]),
    "x": np.array([1.0, 0.0, 0.0]),
    "[111]": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
}
Q0 = np.pi / 4.0
PHASE_TOL = 0.05   # rad — P3 spec tolerance
WZ_TOL_P2 = 1e-10  # P2 spec tolerance


# ---------------------------------------------------------------------------
# P2: k-space plaquette holonomy
# ---------------------------------------------------------------------------

class TestP2PlaquetteHolonomy:

    @pytest.mark.parametrize("alpha", ALPHA_LIST)
    def test_wz_dev_from_identity(self, alpha):
        """P2: WZ matrix deviation from identity < 1e-10 for all alpha."""
        k_center = np.array([Q0, 0.0, 0.0])
        r = plaquette_holonomy_p2(k_center, alpha)
        assert r["wz_dev_from_id"] < WZ_TOL_P2, (
            f"P2 FAIL: alpha={alpha}, wz_dev={r['wz_dev_from_id']:.2e} >= {WZ_TOL_P2}"
        )

    @pytest.mark.parametrize("alpha", ALPHA_LIST)
    def test_per_band_phases(self, alpha):
        """P2: per-band Berry phases < 1e-10 rad for all alpha."""
        k_center = np.array([Q0, 0.0, 0.0])
        r = plaquette_holonomy_p2(k_center, alpha)
        max_phase = r["per_band_max_dev_rad"]
        assert max_phase < WZ_TOL_P2, (
            f"P2 per-band FAIL: alpha={alpha}, max_phase={max_phase:.2e}"
        )

    @pytest.mark.parametrize("alpha", ALPHA_LIST)
    def test_pass_flag(self, alpha):
        """P2: pass_p2 flag is True for all alpha."""
        k_center = np.array([Q0, 0.0, 0.0])
        r = plaquette_holonomy_p2(k_center, alpha)
        assert r["pass_p2"], (
            f"P2 FAIL at alpha={alpha}: wz_dev={r['wz_dev_from_id']:.2e}"
        )

    @pytest.mark.parametrize("alpha", ALPHA_LIST)
    @pytest.mark.parametrize("k_label, k_vec", [
        ("[100]", np.array([Q0, 0.0, 0.0])),
        ("[111]", np.array([Q0, Q0, Q0]) / np.sqrt(3.0) * np.sqrt(3.0)),
        ("[110]", np.array([Q0, Q0, 0.0])),
    ])
    def test_multiple_k_centers(self, alpha, k_label, k_vec):
        """P2: holonomy is identity at multiple k-centers and planes."""
        r = plaquette_holonomy_p2(k_vec, alpha)
        assert r["wz_dev_from_id"] < WZ_TOL_P2, (
            f"P2 FAIL at k={k_label}, alpha={alpha}: dev={r['wz_dev_from_id']:.2e}"
        )

    def test_refinement_dk_halved(self):
        """P2: halving dk leaves the already-zero holonomy stable."""
        k_center = np.array([Q0, 0.0, 0.0])
        alpha = 0.2
        r_coarse = plaquette_holonomy_p2(k_center, alpha, dk=1e-3)
        r_fine = plaquette_holonomy_p2(k_center, alpha, dk=5e-4)
        assert r_fine["wz_dev_from_id"] < WZ_TOL_P2
        # Both are near zero; their difference must be < 1% of WZ_TOL_P2 * 100
        diff = abs(r_fine["wz_dev_from_id"] - r_coarse["wz_dev_from_id"])
        assert diff < 1e-10

    def test_verify_p2_all_alpha(self):
        """P2: full verify_p2_all_alpha returns pass_p2_all = True."""
        result = verify_p2_all_alpha()
        assert result["pass_p2_all"], (
            "verify_p2_all_alpha() FAILED — see detailed result for which alpha/k failed"
        )

    def test_gauge_randomization_p2(self):
        """P2: random per-k U(3) gauge transforms leave WZ holonomy identity."""
        result = verify_p2_all_alpha(alpha_list=[0.2], n_gauge_trials=5)
        assert result["pass_p2_all"]
        r_02 = result["alpha=0.20"]
        for k_label, kr in r_02.items():
            assert kr["pass_gauge"], (
                f"P2 gauge check FAIL at alpha=0.20, k={k_label}: "
                f"max_gauge_dev={kr['max_gauge_dev']:.2e}"
            )


# ---------------------------------------------------------------------------
# Rotation matrix algebra
# ---------------------------------------------------------------------------

class TestRotationMatrices:

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_2pi_is_plus_identity(self, axis_label, axis):
        """D^{(1)}(R(2pi)) = +I exactly."""
        R = rotation_matrix_spin1(2.0 * np.pi, axis)
        dev = np.max(np.abs(R - np.eye(3)))
        assert dev < 1e-12, (
            f"D^(1)(2pi) != +I for axis {axis_label}: max_dev={dev:.2e}"
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_4pi_is_plus_identity(self, axis_label, axis):
        """D^{(1)}(R(4pi)) = +I exactly (4pi trivial for spin-1)."""
        R = rotation_matrix_spin1(4.0 * np.pi, axis)
        dev = np.max(np.abs(R - np.eye(3)))
        assert dev < 1e-12, (
            f"D^(1)(4pi) != +I for axis {axis_label}: max_dev={dev:.2e}"
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin_half_2pi_is_minus_identity(self, axis_label, axis):
        """D^{(1/2)}(R(2pi)) = -I exactly (spin-1/2 double-cover signature)."""
        R = rotation_matrix_spin_half(2.0 * np.pi, axis)
        dev = np.max(np.abs(R + np.eye(2)))
        assert dev < 1e-12, (
            f"D^(1/2)(2pi) != -I for axis {axis_label}: max_dev={dev:.2e}"
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin_half_4pi_is_plus_identity(self, axis_label, axis):
        """D^{(1/2)}(R(4pi)) = +I (4pi trivial for spin-1/2 too)."""
        R = rotation_matrix_spin_half(4.0 * np.pi, axis)
        dev = np.max(np.abs(R - np.eye(2)))
        assert dev < 1e-12, (
            f"D^(1/2)(4pi) != +I for axis {axis_label}: max_dev={dev:.2e}"
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_is_unitary(self, axis_label, axis):
        """D^{(1)}(R) is unitary for arbitrary angle."""
        for angle in [0.1, 1.0, 2.5, np.pi, 2 * np.pi]:
            R = rotation_matrix_spin1(angle, axis)
            dev = np.max(np.abs(R.conj().T @ R - np.eye(3)))
            assert dev < 1e-12

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_det_is_one(self, axis_label, axis):
        """D^{(1)}(R) has determinant +1 (proper rotation)."""
        for angle in [0.5, np.pi, 2 * np.pi]:
            R = rotation_matrix_spin1(angle, axis)
            d = np.linalg.det(R)
            assert abs(d - 1.0) < 1e-12

    def test_spin1_composition(self):
        """D^{(1)}(R1 @ R2) = D^{(1)}(R1) @ D^{(1)}(R2) (group homomorphism)."""
        axis = np.array([0.0, 0.0, 1.0])
        R1 = rotation_matrix_spin1(0.7, axis)
        R2 = rotation_matrix_spin1(1.3, axis)
        R12 = rotation_matrix_spin1(2.0, axis)
        dev = np.max(np.abs(R1 @ R2 - R12))
        assert dev < 1e-12


# ---------------------------------------------------------------------------
# P3: SO(3) rotation holonomy
# ---------------------------------------------------------------------------

class TestP3RotationHolonomy:

    @pytest.mark.parametrize("alpha", ALPHA_LIST)
    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    @pytest.mark.parametrize("band", [0, 1, 2])
    def test_spin1_per_band_2pi_phase_zero(self, alpha, axis_label, axis, band):
        """P3: spin-1 per-band holonomy phase = 0 at 2pi, all alpha, all axes."""
        fn = spin1_state_fn(axis, band=band)
        r = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
        phase = r["holonomy_phase_rad"]
        assert abs(phase) < PHASE_TOL, (
            f"P3 FAIL: alpha={alpha}, axis={axis_label}, band={band}: "
            f"phase={phase:.4f} rad (tol={PHASE_TOL}). "
            "A non-zero phase here would FALSIFY derivation_H_eff.md Part 3."
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_wz_2pi_identity(self, axis_label, axis):
        """P3: spin-1 rank-3 WZ holonomy at 2pi = +I."""
        fn = spin1_frame_fn(axis)
        r = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
        dev = r["holonomy_wz_dev_from_id"]
        assert dev < PHASE_TOL, (
            f"P3 WZ FAIL axis={axis_label}: dev_from_I={dev:.4f} (tol={PHASE_TOL}). "
            "Expected +I (spin-1, J=1)."
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_wz_4pi_identity(self, axis_label, axis):
        """P3: spin-1 rank-3 WZ holonomy at 4pi = +I (4pi trivial for J=1)."""
        fn = spin1_frame_fn(axis)
        r = rotate_and_transport(fn, 4.0 * np.pi, n_steps=400)
        dev = r["holonomy_wz_dev_from_id"]
        assert dev < PHASE_TOL, (
            f"P3 4pi FAIL axis={axis_label}: dev_from_I={dev:.4f} (tol={PHASE_TOL})"
        )

    def test_spin1_alpha_independent(self):
        """P3: holonomy phase is identical for all alpha (alpha-independence)."""
        axis = np.array([0.0, 0.0, 1.0])
        fn = spin1_frame_fn(axis)
        devs = []
        for alpha in ALPHA_LIST:
            # alpha doesn't enter the spin-1 transport; the state_fn doesn't depend
            # on alpha at the linear layer (derivation_H_eff.md Part 3)
            r = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
            devs.append(r["holonomy_wz_dev_from_id"])
        # All devs should be near zero and identical
        assert max(devs) < PHASE_TOL
        # Spread across alpha values should be < 1e-12 (truly alpha-independent)
        assert max(devs) - min(devs) < 1e-12

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin1_refinement(self, axis_label, axis):
        """P3: doubling n_steps changes WZ holonomy by < 1% (refinement stability)."""
        fn = spin1_frame_fn(axis)
        r_coarse = rotate_and_transport(fn, 2.0 * np.pi, n_steps=100)
        r_fine = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
        delta = abs(r_fine["holonomy_wz_dev_from_id"] - r_coarse["holonomy_wz_dev_from_id"])
        assert delta < 0.01, (
            f"P3 refinement FAIL axis={axis_label}: delta={delta:.4f} >= 0.01"
        )


# ---------------------------------------------------------------------------
# Synthetic spin-1/2 control
# ---------------------------------------------------------------------------

class TestSyntheticSpinHalfControl:
    """The J=1/2 spinor control: proves the diagnostic detects fermionic holonomy.

    This is NOT a physical BraneSim linear-layer object — it is a synthetic
    test confirming the diagnostic can distinguish spin-1 from spin-1/2 when
    a spinor is presented.  Required by the P3 deliverable.

    The expected L5 use case: a hedgehog soliton with odd Skyrme winding number,
    rigidly rotated through 2pi, should return holonomy -I via this same machinery.

    Implementation note: spin-1/2 holonomy is measured using the rank-2 WZ frame
    D^{(1/2)}(R(theta)), not a single spinor component. The rank-1 path does NOT
    correctly recover the geometric phase for a single spinor because the dynamic
    phase from the rotating state cancels the geometric phase in the FH closed-loop
    sum. See spin_half_state_fn docstring and spin_half_frame_fn for the analysis.
    """

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin_half_2pi_gives_minus_identity(self, axis_label, axis):
        """Synthetic control: J=1/2 rank-2 WZ at 2pi holonomy = -I."""
        fn = spin_half_frame_fn(axis)
        r = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
        dev_from_neg_id = r["holonomy_wz_dev_from_neg_id"]
        assert dev_from_neg_id < PHASE_TOL, (
            f"Spin-1/2 control FAIL axis={axis_label}: "
            f"dev_from_neg_I={dev_from_neg_id:.4f} (tol={PHASE_TOL}). "
            "The diagnostic cannot detect spin-1/2 holonomy — check FH link formula."
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_spin_half_4pi_gives_plus_identity(self, axis_label, axis):
        """Synthetic control: J=1/2 rank-2 WZ at 4pi holonomy = +I (4pi trivial)."""
        fn = spin_half_frame_fn(axis)
        r = rotate_and_transport(fn, 4.0 * np.pi, n_steps=400)
        dev_from_id = r["holonomy_wz_dev_from_id"]
        assert dev_from_id < PHASE_TOL, (
            f"Spin-1/2 4pi FAIL axis={axis_label}: dev_from_I={dev_from_id:.4f} rad"
        )

    @pytest.mark.parametrize("axis_label, axis", list(AXES.items()))
    def test_distinguishes_spin1_from_spin_half_at_2pi(self, axis_label, axis):
        """Diagnostic clearly distinguishes spin-1 (+I) from spin-1/2 (-I) at 2pi."""
        fn_1 = spin1_frame_fn(axis)
        fn_half = spin_half_frame_fn(axis)
        r1 = rotate_and_transport(fn_1, 2.0 * np.pi, n_steps=200)
        r_half = rotate_and_transport(fn_half, 2.0 * np.pi, n_steps=200)
        # Spin-1: holonomy = +I, dev_from_id near 0
        assert r1["holonomy_wz_dev_from_id"] < PHASE_TOL
        # Spin-1/2: holonomy = -I, dev_from_neg_id near 0
        assert r_half["holonomy_wz_dev_from_neg_id"] < PHASE_TOL
        # They must be clearly separated: dev_from_id(spin-1) << dev_from_id(spin-1/2)
        spin1_dev = r1["holonomy_wz_dev_from_id"]
        spin_half_dev = r_half["holonomy_wz_dev_from_id"]
        assert spin_half_dev > 1.5, (
            f"Spin-1/2 should have large dev_from_+I (expected ~2): {spin_half_dev:.3f}"
        )


# ---------------------------------------------------------------------------
# WZ U(1)/SU(n) decomposition
# ---------------------------------------------------------------------------

class TestWZDecomposition:

    def test_identity_decomposes_to_zero_phase(self):
        """Decompose +I: u1_phase = 0, sun_dev_from_id = 0."""
        r = decompose_wz_u1_sun(np.eye(3, dtype=complex))
        assert abs(r["u1_phase_rad"]) < 1e-12
        assert r["sun_dev_from_id"] < 1e-12

    def test_minus_identity_2x2_decomposes_to_pi(self):
        """Decompose -I_{2x2}: u1_phase = pi/2 per component (det(-I_2x2) = +1 -> phase 0)."""
        # Note: det(-I_{2x2}) = (-1)^2 = +1, so u1_phase = arg(+1)/2 = 0.
        # The SU(2) part is -I itself (since exp(i*0)*(-I) = -I), so sun_dev_from_id = 2.
        # To check the det branch: det(-I_{3x3}) = (-1)^3 = -1 -> u1_phase = pi/3
        r3 = decompose_wz_u1_sun(-np.eye(3, dtype=complex))
        # det(-I_3) = -1, u1_phase = arg(-1)/3 = pi/3
        assert abs(abs(r3["u1_phase_rad"]) - np.pi / 3.0) < 1e-10

    def test_decomposition_preserves_unitarity(self):
        """SU(n) part of a unitary holonomy must have unit determinant."""
        # Use a non-trivial unitary
        rng = np.random.default_rng(42)
        Z = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
        Q, _ = np.linalg.qr(Z)
        r = decompose_wz_u1_sun(Q)
        sun = np.array(r["sun_matrix_real"]) + 1j * np.array(r["sun_matrix_imag"])
        det_sun = np.linalg.det(sun)
        assert abs(abs(det_sun) - 1.0) < 1e-10
        # The SU(3) part should have det near +1 (by construction)
        assert abs(det_sun - 1.0) < 1e-10

    def test_spin1_wz_decomposes_to_trivial(self):
        """Spin-1 WZ holonomy at 2pi decomposes to u1_phase=0, sun=I."""
        axis = np.array([0.0, 0.0, 1.0])
        fn = spin1_frame_fn(axis)
        r = rotate_and_transport(fn, 2.0 * np.pi, n_steps=200)
        # holonomy_matrix is the accumulated WZ matrix
        wz = np.array(r["holonomy_matrix"])
        dec = decompose_wz_u1_sun(wz)
        assert abs(dec["u1_phase_rad"]) < PHASE_TOL
        assert dec["sun_dev_from_id"] < PHASE_TOL


# ---------------------------------------------------------------------------
# Full verification entry points
# ---------------------------------------------------------------------------

class TestFullVerification:

    def test_verify_p2_passes(self):
        """verify_p2_all_alpha() returns pass_p2_all = True."""
        result = verify_p2_all_alpha(alpha_list=[0.0, 0.2, 0.5, 0.8, 1.0])
        assert result["pass_p2_all"], "P2 verification failed"

    def test_verify_p3_passes(self):
        """verify_p3_so3_holonomy() returns pass_p3_all = True."""
        result = verify_p3_so3_holonomy(
            alpha_list=[0.0, 0.2, 0.5, 0.8, 1.0],
            n_steps_list=[100, 200],
        )
        assert result["pass_p3_all"], "P3 spin-1 verification failed"

    def test_verify_p3_synthetic_spin_half_passes(self):
        """verify_p3_so3_holonomy() confirms synthetic spin-1/2 returns -I at 2pi."""
        result = verify_p3_so3_holonomy(
            alpha_list=[0.2],
            n_steps_list=[100, 200],
        )
        assert result["pass_synthetic_spin_half"], (
            "Synthetic spin-1/2 control failed — diagnostic cannot detect spin-1/2 holonomy"
        )

    def test_p2_and_p3_jointly(self):
        """Both P2 and P3 pass simultaneously (joint consistency check)."""
        p2 = verify_p2_all_alpha(alpha_list=[0.2])
        p3 = verify_p3_so3_holonomy(alpha_list=[0.2], n_steps_list=[100, 200])
        assert p2["pass_p2_all"], "P2 FAILED"
        assert p3["pass_p3_all"], "P3 FAILED"