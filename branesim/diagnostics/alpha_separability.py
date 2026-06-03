"""Track-A alpha-separability diagnostic: U(1) x SU(3) split over an alpha-grid.

Read-only diagnostic: no back-reaction, no state mutation. Dimension-agnostic
where reasonable; the trace/traceless decomposition is specific to the lateral
3x3 block (dim=3) because SU(3) lives there, but the underlying lattice
operator construction is kept dim-agnostic via d_of_k_eigenvalues.

Physics background (alpha_separability/SPEC.md, derivation_H_eff.md)
----------------------------------------------------------------------
The lateral dynamical-matrix eigenvalues on the 6-neighbor axial stencil are

    omega_a^2(k) = (2 k_s / rho) * [alpha * h_a + (1-alpha) * H]
    h_a = 1 - cos(k_a * a),   H = sum_b h_b

Projection onto the 3x3 lateral block:

    P_U1  = (1/3) * 1 * 1^T   (dilatational, U(1))
    P_SU3 = I - P_U1          (shear, SU(3))

gives

    trace eigenvalue:      lambda_bar = (2 k_s / rho) * H * (1 - 2*alpha/3)
    traceless (per axis):  lambda_a - lambda_bar = (2 k_s / rho) * alpha * (h_a - H/3)

SU(3)-to-U(1) ratio:

    rho_SU3(k_hat, alpha) = g(k_hat) * sqrt(3) * alpha / (3 - 2*alpha)
    g(k_hat) = sqrt(sum_a (h_a - H/3)^2) / H

Key identities (all exact in closed form):
    * traceless content is linear in alpha (zero at alpha=0)
    * g([111]) = 0 (triplet degeneracy)
    * at alpha=0.2 along [100], rho_SU3 = g * sqrt(3)*0.2/(3-0.4) = g * 0.1323

Prediction P1 (group-velocity anisotropy)
-----------------------------------------
For a [100] carrier at k0a = pi/4:
    v_transverse / v_longitudinal = sqrt(1 - alpha)
The group velocities follow from dω/dk, and at k=(k0,0,0):
    v_L = dω_x/dk_x  ~ c_L * sin(k0*a)/(k0*a) * normalization
    v_T = dω_y/dk_x  = 0  (D is diagonal)
BUT the transverse envelope drifts because the transverse channel has its own
carrier-group velocity dω_y/dk_y evaluated at k=(k0,0,0) where k_y=0.
For a Gaussian wavepacket with carrier k0 along x:
    The y-channel is on its own branch with omega_y^2(k) at k=(k0,0,0).
    A transverse perturbation in y with carrier wavevector ALONG x does not
    propagate; instead we need to track the group velocity ratio for a
    wavepacket with carrier k0 along x, where:
        - longitudinal channel (x-pol): v_g = dω_x/dk_x|_{k=(k0,0,0)}
        - transverse channel (y-pol): v_g = dω_y/dk_x|_{k=(k0,0,0)} = 0 (y-channel insensitive to k_x shift)

Per derivation_H_eff.md the group velocities are:
    v^{(a)}_j = k_s * a * [alpha * delta_{aj} + (1-alpha)] * sin(k_{0,j} * a) / (rho * omega_a(k_0))

For carrier along [100] (j=x only), k_{0,y}=k_{0,z}=0:
    v^{(x)}_x = k_s * a * [alpha + (1-alpha)] * sin(k_{0,x}*a) / (rho * omega_x(k_0))
              = k_s * a * sin(k_{0,x}*a) / (rho * omega_x(k_0))   [longitudinal]
    v^{(y)}_x = k_s * a * (1-alpha) * sin(k_{0,x}*a) / (rho * omega_y(k_0))  [transverse]

The RATIO:
    v^{(y)}_x / v^{(x)}_x = (1-alpha) * omega_x(k_0) / omega_y(k_0)

At k=(k0,0,0):
    omega_x^2 = (2 k_s/rho) * [alpha * h_x + (1-alpha) * h_x] = (2 k_s/rho) * h_x
    omega_y^2 = (2 k_s/rho) * [alpha * 0 + (1-alpha) * h_x]  = (2 k_s/rho) * (1-alpha) * h_x

So omega_y / omega_x = sqrt(1-alpha), and the ratio becomes:
    v_T / v_L = (1-alpha) * (1/sqrt(1-alpha)) = sqrt(1-alpha)

This confirms the SPEC prediction. The physical picture: launching a y-polarized
wavepacket with x-carrier k0, the packet drifts in x at speed v^{(y)}_x, which
is sqrt(1-alpha) times the x-polarized drift speed v^{(x)}_x.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Allow standalone execution from any working directory
_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from branesim.core.conventions import d_of_k_eigenvalues


# ---------------------------------------------------------------------------
# Projection operators (lateral 3x3 block)
# ---------------------------------------------------------------------------

def projection_operators() -> tuple[np.ndarray, np.ndarray]:
    """Return the U(1) and SU(3) projection operators on the lateral 3x3 block.

    P_U1  = (1/3) * ones(3,3)   dilatational / trace channel
    P_SU3 = I_3 - P_U1          shear / traceless channel

    Returns
    -------
    P_U1, P_SU3 : ndarray, shape (3, 3) each
    """
    P_U1 = np.ones((3, 3), dtype=float) / 3.0
    P_SU3 = np.eye(3, dtype=float) - P_U1
    return P_U1, P_SU3


# ---------------------------------------------------------------------------
# Directional anisotropy factor g(k_hat)
# ---------------------------------------------------------------------------

def g_factor(k_hat: np.ndarray, a: float = 1.0) -> float:
    """Directional anisotropy factor g(k_hat) = sqrt(sum_a (h_a - H/3)^2) / H.

    Parameters
    ----------
    k_hat : array_like, shape (3,)
        Unit wavevector direction. Magnitude is used as-is; factor is
        dimensionless so a scale factor on k_hat cancels.
    a : float
        Lattice spacing (default 1.0). Only |k|*a matters.

    Returns
    -------
    float
        g(k_hat). Returns 0.0 for k_hat = [0,0,0] (H=0 limit).
    """
    k = np.asarray(k_hat, dtype=float).ravel()
    if len(k) != 3:
        raise ValueError(f"k_hat must have 3 components; got shape {k.shape}")
    h = 1.0 - np.cos(k * a)      # (3,)
    H = float(h.sum())
    if H < 1e-30:
        return 0.0
    dev = h - H / 3.0             # (3,) traceless deviation
    return float(np.sqrt(np.sum(dev**2)) / H)


# ---------------------------------------------------------------------------
# Closed-form alpha-curve observables
# ---------------------------------------------------------------------------

def closed_form_observables(
    k: np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Closed-form U(1)/SU(3) decomposition at a single (k, alpha).

    Parameters
    ----------
    k : array_like, shape (3,)
        Wavevector (must be 3D; function is specific to the lateral triplet).
    alpha : float
        Prestress ratio in [0, 1].
    k_s, rho, a : float
        Lattice parameters (dimensionless defaults: all 1.0).

    Returns
    -------
    dict with keys:
        lambda_a         : ndarray (3,)  eigenvalues = omega_a^2
        lambda_bar       : float         trace (U(1)) eigenvalue
        traceless        : ndarray (3,)  lambda_a - lambda_bar
        H                : float         sum of h_b
        g                : float         g(k_hat) directional factor
        rho_SU3          : float         SU(3)/U(1) ratio
        prefactor        : float         2*k_s/rho
    """
    k_arr = np.asarray(k, dtype=float).ravel()
    if len(k_arr) != 3:
        raise ValueError(f"k must be 3-dimensional; got {len(k_arr)}")
    h = 1.0 - np.cos(k_arr * a)
    H = float(h.sum())
    prefactor = 2.0 * k_s / rho

    lambda_a = prefactor * (alpha * h + (1.0 - alpha) * H)
    lambda_bar = prefactor * H * (1.0 - 2.0 * alpha / 3.0)
    traceless = lambda_a - lambda_bar

    g = g_factor(k_arr, a=a)
    denom = 3.0 - 2.0 * alpha
    rho_SU3 = g * np.sqrt(3.0) * alpha / denom if abs(denom) > 1e-15 else np.inf

    return {
        "lambda_a": lambda_a,
        "lambda_bar": lambda_bar,
        "traceless": traceless,
        "H": H,
        "g": g,
        "rho_SU3": rho_SU3,
        "prefactor": prefactor,
    }


def alpha_curve(
    k: np.ndarray,
    alpha_grid: np.ndarray | None = None,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Evaluate the alpha-curve (all Track-A observables vs alpha) at a fixed k.

    Parameters
    ----------
    k : array_like, shape (3,)
        Wavevector.
    alpha_grid : array_like, optional
        Alpha values; defaults to linspace(0, 1, 51).
    k_s, rho, a : float
        Lattice parameters.

    Returns
    -------
    dict with keys:
        alpha        : (N,) alpha values
        lambda_bar   : (N,) trace eigenvalue
        traceless    : (N, 3) per-axis traceless content
        rho_SU3      : (N,) SU(3)/U(1) ratio
        g            : float (direction-dependent, constant in alpha)
        H            : float
    """
    if alpha_grid is None:
        alpha_grid = np.linspace(0.0, 1.0, 51)
    alpha_grid = np.asarray(alpha_grid, dtype=float)

    lambda_bar_arr = np.empty(len(alpha_grid), dtype=float)
    traceless_arr = np.empty((len(alpha_grid), 3), dtype=float)
    rho_SU3_arr = np.empty(len(alpha_grid), dtype=float)

    k_arr = np.asarray(k, dtype=float).ravel()
    h = 1.0 - np.cos(k_arr * a)
    H = float(h.sum())
    g = g_factor(k_arr, a=a)
    prefactor = 2.0 * k_s / rho

    for i, alpha in enumerate(alpha_grid):
        obs = closed_form_observables(k_arr, alpha, k_s=k_s, rho=rho, a=a)
        lambda_bar_arr[i] = obs["lambda_bar"]
        traceless_arr[i] = obs["traceless"]
        rho_SU3_arr[i] = obs["rho_SU3"]

    return {
        "alpha": alpha_grid,
        "lambda_bar": lambda_bar_arr,
        "traceless": traceless_arr,
        "rho_SU3": rho_SU3_arr,
        "g": g,
        "H": H,
        "prefactor": prefactor,
    }


# ---------------------------------------------------------------------------
# Lattice dynamical block (numerical, from explicit bond sum)
# ---------------------------------------------------------------------------

def build_dynamical_block_3d(
    k: np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> np.ndarray:
    """Build the 3x3 dynamical matrix by explicit 6-bond summation.

    Mirrors the sprint2_subtask9 reference construction; does NOT use the
    closed form (which is what we want to verify against).

    D_ij(k) = (k_s/rho) * sum_{delta in axial_bonds} (1-cos(k.delta*a))
              * [alpha * delta_hat_i * delta_hat_j + (1-alpha) * Id_ij]
    """
    k = np.asarray(k, dtype=float).ravel()
    if len(k) != 3:
        raise ValueError(f"k must be 3-dimensional; got {len(k)}")
    D = np.zeros((3, 3), dtype=float)
    for axis in range(3):
        for sign in (+1, -1):
            delta = np.zeros(3, dtype=float)
            delta[axis] = float(sign)
            structure = 1.0 - np.cos(np.dot(k, a * delta))
            bond_tensor = alpha * np.outer(delta, delta) + (1.0 - alpha) * np.eye(3)
            D += structure * bond_tensor
    return (k_s / rho) * D


def numerical_trace_traceless(
    k: np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Apply P_U1/P_SU3 to the numerical D(k) and extract trace/traceless split.

    Builds D(k) from bond sums, projects with P_U1 and P_SU3, and returns
    the projected eigenvalues (which equal the diagonal entries since D is
    diagonal in the Cartesian basis).

    Returns
    -------
    dict with keys:
        D              : (3,3) full dynamical block
        D_U1           : (3,3) P_U1 @ D @ P_U1 (trace part, rank-1)
        D_SU3          : (3,3) P_SU3 @ D @ P_SU3 (traceless part, rank-2)
        lambda_diag    : (3,) diagonal of D (eigenvalues, Cartesian basis)
        lambda_bar_num : float  trace of D / 3 (= the U(1) eigenvalue)
        traceless_num  : (3,) diagonal of D minus lambda_bar_num
        offdiag_max    : float  max |off-diagonal| of D (should be ~1e-15)
    """
    D = build_dynamical_block_3d(k, alpha, k_s=k_s, rho=rho, a=a)
    P_U1, P_SU3 = projection_operators()

    D_U1 = P_U1 @ D @ P_U1
    D_SU3 = P_SU3 @ D @ P_SU3

    diag = np.diag(D)
    lambda_bar_num = float(np.trace(D) / 3.0)
    traceless_num = diag - lambda_bar_num

    off = D - np.diag(diag)
    offdiag_max = float(np.max(np.abs(off)))

    return {
        "D": D,
        "D_U1": D_U1,
        "D_SU3": D_SU3,
        "lambda_diag": diag,
        "lambda_bar_num": lambda_bar_num,
        "traceless_num": traceless_num,
        "offdiag_max": offdiag_max,
    }


# ---------------------------------------------------------------------------
# Verification suite (Track A numerical confirmation)
# ---------------------------------------------------------------------------

def verify_track_a(
    k_samples: dict[str, np.ndarray] | None = None,
    alpha_samples: list[float] | None = None,
    tol: float = 1e-10,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Run all Track-A numerical confirmation checks.

    Checks:
      1. Closed-form vs numerical (bond-sum) dynamical block: should match to
         near machine precision (tol default 1e-10, expected ~1e-15).
      2. Traceless content is linear in alpha: at each k, the traceless vector
         must be exactly proportional to alpha.  Test: traceless(alpha) / alpha
         is constant to tol (for alpha > 0).
      3. g([111]) = 0 exactly.
      4. At alpha=0.2 along [100], rho_SU3 = 0.1323... * g (sqrt(3)*0.2/2.6).

    Returns
    -------
    dict with pass/fail flags and measured deviations.
    """
    if k_samples is None:
        k0 = np.pi / 4.0
        k_samples = {
            "[100]": np.array([k0, 0.0, 0.0]),
            "[110]": np.array([k0, k0, 0.0]),
            "[111]": np.array([k0, k0, k0]),
            "[210]": np.array([2 * k0, k0, 0.0]),
            "[321]": np.array([3 * k0 / 4, 2 * k0 / 4, k0 / 4]),
        }
    if alpha_samples is None:
        alpha_samples = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0]

    results = {}

    # --- Check 1: closed-form vs numerical ---
    max_cf_err = 0.0
    for k_label, k_vec in k_samples.items():
        for alpha in alpha_samples:
            cf = closed_form_observables(k_vec, alpha, k_s=k_s, rho=rho, a=a)
            num = numerical_trace_traceless(k_vec, alpha, k_s=k_s, rho=rho, a=a)

            # Compare lambda_bar
            err_bar = abs(cf["lambda_bar"] - num["lambda_bar_num"])
            # Compare traceless per axis
            err_tl = float(np.max(np.abs(cf["traceless"] - num["traceless_num"])))
            max_cf_err = max(max_cf_err, err_bar, err_tl)

    results["max_closed_form_vs_numerical_err"] = max_cf_err
    results["check1_pass"] = bool(max_cf_err < tol)

    # --- Check 2: traceless content linear in alpha ---
    # For alpha > 0, traceless(alpha) / alpha must be k-dependent but
    # alpha-independent. Verify by computing the ratio at several alpha
    # values and checking consistency.
    #
    # Special case: at k along [111], g([111])=0 so the traceless content is
    # identically zero at all alpha. The "ratio" traceless/alpha = 0 for all
    # alpha > 0, which IS linear (= 0*alpha). We verify absolute smallness.
    max_linearity_err = 0.0
    for k_label, k_vec in k_samples.items():
        ref_ratio = None
        for alpha in alpha_samples:
            if alpha == 0.0:
                # At alpha=0, traceless must be exactly zero
                cf = closed_form_observables(k_vec, 0.0, k_s=k_s, rho=rho, a=a)
                err_zero = float(np.max(np.abs(cf["traceless"])))
                max_linearity_err = max(max_linearity_err, err_zero)
                continue
            cf = closed_form_observables(k_vec, alpha, k_s=k_s, rho=rho, a=a)
            ratio = cf["traceless"] / alpha
            if ref_ratio is None:
                ref_ratio = ratio.copy()
            else:
                # All ratios should be equal (same k-dependent prefactor).
                err = float(np.max(np.abs(ratio - ref_ratio)))
                # Normalize by the prefactor scale: (2*k_s/rho) * max(h_a - H/3)
                # Use the closed-form scale, not ref_ratio, to avoid 0/0 at [111].
                h = 1.0 - np.cos(k_vec * a)
                H = float(h.sum())
                scale = float(2.0 * k_s / rho * max(np.max(np.abs(h - H / 3.0)), 1e-30))
                if scale < 1e-20:
                    # [111]-like: traceless is identically 0; check absolute
                    max_linearity_err = max(max_linearity_err, err)
                else:
                    max_linearity_err = max(max_linearity_err, err / scale)

    results["max_linearity_rel_err"] = max_linearity_err
    results["check2_pass"] = bool(max_linearity_err < tol)

    # --- Check 3: g([111]) = 0 ---
    k_111 = np.array([np.pi / 4, np.pi / 4, np.pi / 4])
    g_111 = g_factor(k_111, a=a)
    results["g_111"] = float(g_111)
    results["check3_pass"] = bool(abs(g_111) < 1e-14)

    # Verify for multiple k magnitudes along [111]
    max_g_111 = abs(g_111)
    for mag in [0.1, 0.3, 0.5, 0.8, 1.2]:
        k_mag_111 = mag / np.sqrt(3.0) * np.ones(3)
        max_g_111 = max(max_g_111, abs(g_factor(k_mag_111, a=a)))
    results["max_g_111_over_magnitudes"] = float(max_g_111)
    results["check3_strict_pass"] = bool(max_g_111 < 1e-14)

    # --- Check 4: rho_SU3 at alpha=0.2 along [100] ---
    k_100 = np.array([np.pi / 4, 0.0, 0.0])
    cf_020 = closed_form_observables(k_100, 0.2, k_s=k_s, rho=rho, a=a)
    rho_SU3_020 = cf_020["rho_SU3"]
    g_020 = cf_020["g"]
    # Expected: sqrt(3)*0.2 / (3-0.4) = sqrt(3)*0.2/2.6
    expected_coeff = np.sqrt(3.0) * 0.2 / (3.0 - 0.4)
    expected_rho_SU3 = expected_coeff * g_020
    err_rho = abs(rho_SU3_020 - expected_rho_SU3)
    results["rho_SU3_alpha020_100"] = float(rho_SU3_020)
    results["g_100_alpha020"] = float(g_020)
    results["expected_rho_SU3_coeff"] = float(expected_coeff)
    results["check4_pass"] = bool(err_rho < tol)

    # Summary
    results["all_pass"] = (
        results["check1_pass"]
        and results["check2_pass"]
        and results["check3_strict_pass"]
        and results["check4_pass"]
    )
    results["tolerance_used"] = tol

    return results


# ---------------------------------------------------------------------------
# Group-velocity ratio (P1 closed form)
# ---------------------------------------------------------------------------

def group_velocity_ratio_p1(
    k0: np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Closed-form group-velocity ratio for a [100] carrier (Prediction P1).

    For carrier k0 along [100] = (k0, 0, 0), the per-channel group velocity
    in the x-direction (the drift direction for both channels) is:

        v^{(a)}_x = k_s * a * [alpha * delta_{ax} + (1-alpha)] * sin(k0*a) / (rho * omega_a(k0))

    Ratio: v_T / v_L = (1-alpha) * omega_L / omega_T
                     = (1-alpha) / sqrt(1-alpha) = sqrt(1-alpha)

    This is valid for ANY k0a along [100] since the sin(k0*a) factor cancels.

    Parameters
    ----------
    k0 : array_like, shape (3,)
        Carrier wavevector. Must be along [100] (k_y = k_z = 0) for the
        prediction to apply cleanly.
    alpha : float
    k_s, rho, a : float

    Returns
    -------
    dict with:
        v_L      : float  longitudinal group velocity (x-channel, x-direction)
        v_T      : float  transverse group velocity (y-channel, x-direction)
        ratio    : float  v_T / v_L (predicted = sqrt(1-alpha))
        predicted: float  sqrt(1-alpha)
        rel_err  : float  (ratio - predicted) / predicted
        omega_L  : float  omega of x-channel at k0
        omega_T  : float  omega of y-channel at k0
    """
    k0_arr = np.asarray(k0, dtype=float).ravel()
    omega_sq = d_of_k_eigenvalues(k0_arr, alpha, k_s=k_s, rho=rho, a=a)
    omega = np.sqrt(np.maximum(omega_sq, 0.0))

    omega_L = float(omega[0])   # x-channel at k=(k0,0,0)
    omega_T = float(omega[1])   # y-channel at k=(k0,0,0): same h_x but (1-alpha) weight

    # sin(k0_x * a) prefactor is the same for both channels: cancels in ratio
    sin_k0a = np.sin(float(k0_arr[0]) * a)

    # Group velocities in x-direction
    if abs(omega_L) < 1e-15 or abs(omega_T) < 1e-15:
        v_L = v_T = 0.0
        ratio = np.nan
    else:
        alpha_fac_L = alpha + (1.0 - alpha)      # = 1.0 for a=x, j=x
        alpha_fac_T = 0.0 + (1.0 - alpha)        # = 1-alpha for a=y, j=x
        v_L = float(k_s * a * alpha_fac_L * sin_k0a / (rho * omega_L))
        v_T = float(k_s * a * alpha_fac_T * sin_k0a / (rho * omega_T))
        ratio = v_T / v_L if abs(v_L) > 1e-15 else np.nan

    predicted = float(np.sqrt(max(1.0 - alpha, 0.0)))
    rel_err = float((ratio - predicted) / predicted) if not np.isnan(ratio) else np.nan

    return {
        "v_L": v_L,
        "v_T": v_T,
        "ratio": ratio,
        "predicted": predicted,
        "rel_err": rel_err,
        "omega_L": omega_L,
        "omega_T": omega_T,
        "alpha": alpha,
        "k0": k0_arr.tolist(),
    }


# ---------------------------------------------------------------------------
# Alpha-grid survey across multiple k-directions
# ---------------------------------------------------------------------------

def survey_alpha_grid(
    alpha_grid: np.ndarray | None = None,
    k_magnitude: float = np.pi / 4.0,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Survey alpha-curves across canonical k-directions.

    Evaluates closed-form and g-factor for [100], [110], [111], and a few
    off-axis directions over the supplied alpha_grid.

    Returns a dict keyed by direction label, each containing the output of
    alpha_curve() plus the g-factor.
    """
    if alpha_grid is None:
        alpha_grid = np.linspace(0.0, 1.0, 51)

    s2 = 1.0 / np.sqrt(2.0)
    s3 = 1.0 / np.sqrt(3.0)
    directions = {
        "[100]": np.array([1.0, 0.0, 0.0]) * k_magnitude,
        "[010]": np.array([0.0, 1.0, 0.0]) * k_magnitude,
        "[110]": np.array([s2, s2, 0.0]) * k_magnitude,
        "[111]": np.array([s3, s3, s3]) * k_magnitude,
        "[210]": np.array([2.0 / np.sqrt(5.0), 1.0 / np.sqrt(5.0), 0.0]) * k_magnitude,
        "[321]": np.array([3.0, 2.0, 1.0]) / np.sqrt(14.0) * k_magnitude,
    }

    output = {}
    for label, k_vec in directions.items():
        curve = alpha_curve(k_vec, alpha_grid=alpha_grid, k_s=k_s, rho=rho, a=a)
        output[label] = curve

    return output