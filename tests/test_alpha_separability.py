"""Unit tests for branesim.diagnostics.alpha_separability.

Tests cover:
  1. Projection operators: partition of unity, traces, idempotency.
  2. g-factor: g([111]) = 0 exactly; g([100]) is maximal (within cubic symmetry).
  3. Closed-form identities:
     - traceless content is exactly zero at alpha=0
     - traceless / alpha is k-independent (linear in alpha)
     - trace sum: lambda_a sum = prefactor * H * 3 * (1 - 2*alpha/3) [= lambda_bar * 3? No.]
       Actually: sum(lambda_a) = prefactor * [alpha*H + (1-alpha)*3H] = prefactor*H*(3-2*alpha)
       And lambda_bar = prefactor * H * (1 - 2*alpha/3) = sum/3. Check this.
     - rho_SU3 = g * sqrt(3) * alpha / (3 - 2*alpha)
  4. Numerical vs closed-form: tolerance ~1e-14.
  5. P1 prediction: v_T/v_L = sqrt(1-alpha) from H_eff formula.
  6. Full verify_track_a() returns all_pass=True.
"""

from __future__ import annotations

import numpy as np
import pytest

from branesim.diagnostics.alpha_separability import (
    closed_form_observables,
    g_factor,
    group_velocity_ratio_p1,
    numerical_trace_traceless,
    projection_operators,
    verify_track_a,
)


# ---------------------------------------------------------------------------
# 1. Projection operators
# ---------------------------------------------------------------------------

class TestProjectionOperators:
    def test_partition_of_unity(self):
        P_U1, P_SU3 = projection_operators()
        assert np.allclose(P_U1 + P_SU3, np.eye(3), atol=1e-15)

    def test_traces(self):
        P_U1, P_SU3 = projection_operators()
        assert abs(np.trace(P_U1) - 1.0) < 1e-15
        assert abs(np.trace(P_SU3) - 2.0) < 1e-15

    def test_idempotency_P_U1(self):
        P_U1, _ = projection_operators()
        assert np.allclose(P_U1 @ P_U1, P_U1, atol=1e-15)

    def test_idempotency_P_SU3(self):
        _, P_SU3 = projection_operators()
        assert np.allclose(P_SU3 @ P_SU3, P_SU3, atol=1e-15)

    def test_orthogonality(self):
        P_U1, P_SU3 = projection_operators()
        assert np.allclose(P_U1 @ P_SU3, np.zeros((3, 3)), atol=1e-15)

    def test_P_U1_form(self):
        P_U1, _ = projection_operators()
        expected = np.ones((3, 3)) / 3.0
        assert np.allclose(P_U1, expected, atol=1e-15)


# ---------------------------------------------------------------------------
# 2. g-factor
# ---------------------------------------------------------------------------

class TestGFactor:
    def test_g_111_zero_at_all_magnitudes(self):
        """g([111]) = 0 exactly for any |k| (h_a all equal along [111])."""
        for mag in [0.1, 0.3, 0.5, np.pi / 4, 1.0, 1.5]:
            k = mag / np.sqrt(3.0) * np.ones(3)
            assert abs(g_factor(k)) < 1e-14, f"g([111]) != 0 at mag={mag}"

    def test_g_100_nonzero(self):
        """g([100]) is nonzero (maximal cubic anisotropy)."""
        k = np.array([np.pi / 4, 0.0, 0.0])
        g = g_factor(k)
        assert g > 0.5, f"g([100]) unexpectedly small: {g}"

    def test_g_symmetry_cubic(self):
        """g is invariant under cubic permutation symmetry of k."""
        mag = np.pi / 4
        k_100 = np.array([mag, 0.0, 0.0])
        k_010 = np.array([0.0, mag, 0.0])
        k_001 = np.array([0.0, 0.0, mag])
        g_100 = g_factor(k_100)
        g_010 = g_factor(k_010)
        g_001 = g_factor(k_001)
        assert abs(g_100 - g_010) < 1e-14
        assert abs(g_100 - g_001) < 1e-14

    def test_g_zero_k(self):
        """g at k=0 returns 0 (H=0, convention)."""
        k = np.zeros(3)
        g = g_factor(k)
        assert g == 0.0

    def test_g_nonnegative(self):
        """g is non-negative by definition (it's a ratio of norms)."""
        rng = np.random.default_rng(42)
        ks = rng.uniform(0.1, 1.5, size=(20, 3))
        for k in ks:
            assert g_factor(k) >= 0.0


# ---------------------------------------------------------------------------
# 3. Closed-form identities
# ---------------------------------------------------------------------------

class TestClosedFormIdentities:
    def test_traceless_zero_at_alpha0(self):
        """At alpha=0, the dynamical matrix is proportional to identity => traceless=0."""
        for k in [np.array([0.5, 0.0, 0.0]),
                  np.array([0.3, 0.4, 0.0]),
                  np.array([0.2, 0.3, 0.5])]:
            obs = closed_form_observables(k, alpha=0.0)
            assert np.allclose(obs["traceless"], 0.0, atol=1e-14), \
                f"traceless != 0 at alpha=0, k={k}: {obs['traceless']}"

    def test_traceless_sum_is_zero(self):
        """sum(lambda_a - lambda_bar) = 0 always (by definition of mean)."""
        rng = np.random.default_rng(7)
        ks = rng.uniform(0.1, 1.5, size=(15, 3))
        alphas = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0]
        for k in ks:
            for alpha in alphas:
                obs = closed_form_observables(k, alpha)
                assert abs(obs["traceless"].sum()) < 1e-13, \
                    f"traceless sum != 0: {obs['traceless'].sum()}"

    def test_lambda_bar_is_mean(self):
        """lambda_bar = mean(lambda_a)."""
        k = np.array([np.pi / 4, np.pi / 3, np.pi / 6])
        for alpha in [0.0, 0.2, 0.5, 1.0]:
            obs = closed_form_observables(k, alpha)
            assert abs(obs["lambda_bar"] - float(obs["lambda_a"].mean())) < 1e-13

    def test_traceless_linear_in_alpha(self):
        """traceless(alpha) / alpha = const (k-dependent, alpha-independent) for alpha > 0."""
        k = np.array([np.pi / 4, 0.0, 0.0])
        alphas = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
        ratios = []
        for alpha in alphas:
            obs = closed_form_observables(k, alpha)
            ratios.append(obs["traceless"] / alpha)
        for i in range(1, len(ratios)):
            assert np.allclose(ratios[i], ratios[0], atol=1e-13), \
                f"traceless/alpha not constant: {ratios[i]} vs {ratios[0]}"

    def test_traceless_linear_in_alpha_111(self):
        """At k along [111], traceless = 0 for all alpha (g=0 case)."""
        k = np.array([0.3, 0.3, 0.3])
        for alpha in [0.1, 0.2, 0.5, 1.0]:
            obs = closed_form_observables(k, alpha)
            assert np.allclose(obs["traceless"], 0.0, atol=1e-13), \
                f"traceless != 0 at [111], alpha={alpha}: {obs['traceless']}"

    def test_rho_SU3_formula(self):
        """rho_SU3 = g * sqrt(3) * alpha / (3 - 2*alpha)."""
        k = np.array([np.pi / 4, 0.0, 0.0])
        g = g_factor(k)
        for alpha in [0.1, 0.2, 0.5, 0.8]:
            obs = closed_form_observables(k, alpha)
            expected = g * np.sqrt(3.0) * alpha / (3.0 - 2.0 * alpha)
            assert abs(obs["rho_SU3"] - expected) < 1e-14, \
                f"rho_SU3 formula mismatch at alpha={alpha}"

    def test_rho_SU3_zero_at_alpha0(self):
        """rho_SU3 = 0 at alpha=0 (full U(3) degeneracy)."""
        k = np.array([0.5, 0.2, 0.3])
        obs = closed_form_observables(k, alpha=0.0)
        assert abs(obs["rho_SU3"]) < 1e-14

    def test_eigenvalue_sum(self):
        """sum(lambda_a) = prefactor * H * (3 - 2*alpha)."""
        k = np.array([0.4, 0.3, 0.6])
        h = 1.0 - np.cos(k)
        H = float(h.sum())
        for alpha in [0.0, 0.2, 0.5, 1.0]:
            obs = closed_form_observables(k, alpha)
            expected_sum = obs["prefactor"] * H * (3.0 - 2.0 * alpha)
            assert abs(obs["lambda_a"].sum() - expected_sum) < 1e-13, \
                f"eigenvalue sum mismatch at alpha={alpha}"


# ---------------------------------------------------------------------------
# 4. Numerical vs closed-form consistency
# ---------------------------------------------------------------------------

class TestNumericalClosedFormConsistency:
    ALPHAS = [0.0, 0.1, 0.2, 0.5, 0.8, 1.0]
    K_SAMPLES = [
        np.array([np.pi / 4, 0.0, 0.0]),
        np.array([np.pi / 4, np.pi / 4, 0.0]),
        np.array([0.3, 0.3, 0.3]),
        np.array([0.7, 0.4, 0.2]),
        np.array([1.2, 0.0, 0.0]),
    ]

    def test_lambda_bar_matches(self):
        """Closed-form lambda_bar matches numerical trace / 3."""
        tol = 1e-13
        for k in self.K_SAMPLES:
            for alpha in self.ALPHAS:
                cf = closed_form_observables(k, alpha)
                num = numerical_trace_traceless(k, alpha)
                err = abs(cf["lambda_bar"] - num["lambda_bar_num"])
                assert err < tol, \
                    f"lambda_bar mismatch at k={k}, alpha={alpha}: err={err}"

    def test_traceless_matches(self):
        """Closed-form traceless matches numerical per-axis traceless."""
        tol = 1e-13
        for k in self.K_SAMPLES:
            for alpha in self.ALPHAS:
                cf = closed_form_observables(k, alpha)
                num = numerical_trace_traceless(k, alpha)
                err = float(np.max(np.abs(cf["traceless"] - num["traceless_num"])))
                assert err < tol, \
                    f"traceless mismatch at k={k}, alpha={alpha}: err={err}"

    def test_D_is_diagonal(self):
        """Numerical D(k) has near-zero off-diagonal entries."""
        for k in self.K_SAMPLES:
            num = numerical_trace_traceless(k, alpha=0.2)
            assert num["offdiag_max"] < 1e-14, \
                f"D has off-diagonal entries: {num['offdiag_max']}"


# ---------------------------------------------------------------------------
# 5. P1 group-velocity ratio
# ---------------------------------------------------------------------------

class TestGroupVelocityRatioP1:
    def test_ratio_equals_sqrt_1_minus_alpha(self):
        """v_T/v_L = sqrt(1-alpha) from H_eff formula, to machine precision."""
        k0 = np.array([np.pi / 4, 0.0, 0.0])
        for alpha in [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]:
            result = group_velocity_ratio_p1(k0, alpha)
            predicted = np.sqrt(1.0 - alpha)
            assert abs(result["ratio"] - predicted) < 1e-12, \
                f"ratio mismatch at alpha={alpha}: {result['ratio']} vs {predicted}"

    def test_ratio_independent_of_k0_magnitude(self):
        """The ratio v_T/v_L = sqrt(1-alpha) is k0-magnitude-independent."""
        alpha = 0.2
        predicted = np.sqrt(1.0 - alpha)
        for k0a in [0.1, 0.3, np.pi / 4, 0.8, 1.0]:
            k0 = np.array([k0a, 0.0, 0.0])
            result = group_velocity_ratio_p1(k0, alpha)
            assert abs(result["ratio"] - predicted) < 1e-12, \
                f"ratio not k-independent at k0a={k0a}: {result['ratio']}"

    def test_v_L_positive(self):
        """Longitudinal group velocity is positive for non-zero k0."""
        k0 = np.array([np.pi / 4, 0.0, 0.0])
        for alpha in [0.0, 0.2, 0.5]:
            result = group_velocity_ratio_p1(k0, alpha)
            assert result["v_L"] > 0, f"v_L not positive at alpha={alpha}"

    def test_v_T_leq_v_L(self):
        """v_T <= v_L for all alpha in [0, 1) (transverse is slower)."""
        k0 = np.array([np.pi / 4, 0.0, 0.0])
        for alpha in [0.0, 0.1, 0.2, 0.5, 0.8, 0.99]:
            result = group_velocity_ratio_p1(k0, alpha)
            assert result["v_T"] <= result["v_L"] + 1e-12, \
                f"v_T > v_L at alpha={alpha}: {result['v_T']} > {result['v_L']}"

    def test_ratio_at_canonical_alpha_020(self):
        """At alpha=0.2, ratio = sqrt(0.8) = 0.8944...."""
        k0 = np.array([np.pi / 4, 0.0, 0.0])
        result = group_velocity_ratio_p1(k0, alpha=0.2)
        expected = np.sqrt(0.8)
        assert abs(result["ratio"] - expected) < 1e-12

    def test_ratio_at_alpha_050(self):
        """At alpha=0.5, ratio = sqrt(0.5) = 0.7071...."""
        k0 = np.array([np.pi / 4, 0.0, 0.0])
        result = group_velocity_ratio_p1(k0, alpha=0.5)
        expected = np.sqrt(0.5)
        assert abs(result["ratio"] - expected) < 1e-12


# ---------------------------------------------------------------------------
# 6. Full verify_track_a integration test
# ---------------------------------------------------------------------------

class TestVerifyTrackA:
    def test_all_checks_pass(self):
        """verify_track_a() must return all_pass=True."""
        results = verify_track_a()
        assert results["all_pass"], \
            f"verify_track_a failed: {results}"

    def test_check1_closed_form_vs_numerical(self):
        results = verify_track_a()
        assert results["check1_pass"]
        assert results["max_closed_form_vs_numerical_err"] < 1e-10

    def test_check2_linearity_in_alpha(self):
        results = verify_track_a()
        assert results["check2_pass"]
        assert results["max_linearity_rel_err"] < 1e-10

    def test_check3_g_111_zero(self):
        results = verify_track_a()
        assert results["check3_strict_pass"]
        assert results["max_g_111_over_magnitudes"] < 1e-14

    def test_check4_rho_SU3_coefficient(self):
        results = verify_track_a()
        assert results["check4_pass"]
        # rho_SU3 = g * coeff where coeff = sqrt(3)*0.2/(3-0.4)
        # Numerically: coeff ≈ 0.13323
        assert abs(results["expected_rho_SU3_coeff"] - 0.13323) < 1e-4