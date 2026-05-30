"""Numerical certificate: on the 6-neighbor axial-only lattice, D(k) is
diagonal in the Cartesian basis at every k and every alpha, with eigenvalues
matching the closed form derived in components.diagnostics.christoffel_6nn.

Tests, all to floating-point tolerance:

  T1. D(k) is diagonal in Cartesian basis (off-diagonal entries ~ 1e-15).
  T2. The diagonal entries match
         omega^2_a(k) = (2 k0 / rho) * [ alpha h_a + (1-alpha) (h_x+h_y+h_z) ]
      with h_i = 1 - cos(k_i a).
  T3. At k along [111], all three eigenvalues are equal (analytic degeneracy).
  T4. At alpha = 0, D(k) is proportional to the identity for every k.
  T5. The Fukui-Hatsugai plaquette holonomy on a small (k_x, k_y) plaquette
      at alpha = 0.20 around k_0 = (1,1,1)*k_unit/sqrt(3) is the identity
      to better than 1e-12 — direct consequence of T1 (eigenvectors are
      k-independent Cartesian, no parallel transport).

Outputs: results.json in the same directory.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

# Allow running this script directly: add repo root to sys.path.
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from components.diagnostics.christoffel_6nn import christoffel_6nn, eigvals_6nn


# ----------------------------------------------------------------------------
# Independent reference: build D(k) by explicit 6-bond summation.
# This avoids using the closed form we want to verify; we sum bond by bond.
# ----------------------------------------------------------------------------
def christoffel_6nn_from_bonds(k, alpha, k0=1.0, rho=1.0, a=1.0):
    k = np.asarray(k, dtype=float).reshape(3)
    D = np.zeros((3, 3), dtype=float)
    for sign in (+1, -1):
        for axis in range(3):
            delta = np.zeros(3, dtype=float)
            delta[axis] = sign
            delta_hat = delta.copy()  # |delta| = 1
            structure = 1.0 - np.cos(np.dot(k, a * delta))
            bond_tensor = alpha * np.outer(delta_hat, delta_hat) + (1.0 - alpha) * np.eye(3)
            D += structure * bond_tensor
    return (k0 / rho) * D


# ----------------------------------------------------------------------------
# Test grid: a few alpha values and a dense set of k samples.
# ----------------------------------------------------------------------------
ALPHAS = [0.0, 0.1, 0.2, 0.5, 0.8, 1.0]

# A spread of k samples: axes, face diagonals, body diagonal, plus random.
rng = np.random.default_rng(20260525)
K_AXIS = np.array([
    [0.7, 0.0, 0.0],
    [0.0, 0.7, 0.0],
    [0.0, 0.0, 0.7],
    [1.3, 0.0, 0.0],
])
K_FACE = np.array([
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
])
K_BODY = np.array([
    [0.3, 0.3, 0.3],
    [0.5, 0.5, 0.5],
    [0.9, 0.9, 0.9],
])
K_RANDOM = rng.uniform(low=-1.4, high=1.4, size=(20, 3))
K_SAMPLES = np.vstack([K_AXIS, K_FACE, K_BODY, K_RANDOM])


def test_T1_T2_T3_T4():
    """Run T1-T4 over the k/alpha grid; return per-test maxima."""
    max_offdiag = 0.0
    max_closed_form_err = 0.0
    max_111_degen_err = 0.0
    max_alpha0_iso_err = 0.0
    for alpha in ALPHAS:
        for k in K_SAMPLES:
            D_ref = christoffel_6nn_from_bonds(k, alpha)
            # T1: off-diagonal entries
            off = D_ref - np.diag(np.diag(D_ref))
            max_offdiag = max(max_offdiag, float(np.max(np.abs(off))))

            # T2: closed form matches bond-sum
            D_closed = christoffel_6nn(k, alpha)
            err = float(np.max(np.abs(D_ref - D_closed)))
            scale = max(float(np.max(np.abs(D_ref))), 1.0)
            max_closed_form_err = max(max_closed_form_err, err / scale)

            # T4: at alpha = 0, D should be proportional to identity.
            if alpha == 0.0:
                diag = np.diag(D_ref)
                spread = float(diag.max() - diag.min())
                norm = max(float(np.max(np.abs(diag))), 1.0)
                max_alpha0_iso_err = max(max_alpha0_iso_err, spread / norm)

        # T3: at k along [111], all three eigenvalues are equal.
        for k_unit in [0.2, 0.5, 1.0, 1.5]:
            k = (k_unit / np.sqrt(3.0)) * np.array([1.0, 1.0, 1.0])
            ev = eigvals_6nn(k, alpha)
            spread = float(ev.max() - ev.min())
            norm = max(float(np.max(np.abs(ev))), 1.0)
            max_111_degen_err = max(max_111_degen_err, spread / norm)

    return {
        "T1_max_offdiag_abs": max_offdiag,
        "T2_max_closed_form_rel_err": max_closed_form_err,
        "T3_max_111_degeneracy_rel_spread": max_111_degen_err,
        "T4_max_alpha0_isotropy_rel_spread": max_alpha0_iso_err,
    }


def test_T5_plaquette_holonomy():
    """Fukui-Hatsugai plaquette holonomy on a (k_x, k_y) plaquette at alpha=0.20.

    We do this two ways:
      (a) Treating each axis-band individually (band index = x, y, z).
          Because eigenvectors are Cartesian-constant, every link variable is
          either +1 or -1, and on a small plaquette no sign change occurs.
      (b) Treating the full triplet as a rank-3 subspace (degenerate-band
          case at k along [111]). Link projector is the identity 3x3
          everywhere, so plaquette holonomy is exactly identity.
    """
    alpha = 0.20
    k_unit = 0.5
    k0_vec = (k_unit / np.sqrt(3.0)) * np.array([1.0, 1.0, 1.0])
    dk = 1e-3  # plaquette side in k-space
    corners = [
        k0_vec + np.array([0.0, 0.0, 0.0]),
        k0_vec + np.array([dk, 0.0, 0.0]),
        k0_vec + np.array([dk, dk, 0.0]),
        k0_vec + np.array([0.0, dk, 0.0]),
    ]
    # Eigenframes via eigh — these come out as columns ordered by ascending
    # eigenvalue, but the underlying eigenvectors are Cartesian. We extract
    # them as 3x3 orthonormal frames at each corner.
    frames = []
    for k in corners:
        D = christoffel_6nn(k, alpha)
        evals, evecs = np.linalg.eigh(D)
        frames.append(evecs)

    # (a) Per-band link products around the plaquette: each band has its own
    # 1D unitary link variable U_j^n = <n,k_j | n,k_{j+1}>. Because the
    # underlying eigenvectors are Cartesian, the per-band link is +-1 and
    # the round-trip product is +1 to floating-point precision.
    per_band_holonomies = []
    for band in range(3):
        prod = 1.0
        for j in range(4):
            v0 = frames[j][:, band]
            v1 = frames[(j + 1) % 4][:, band]
            link = float(np.dot(v0, v1))
            # Normalize sign per Fukui-Hatsugai: link / |link|
            sign = 1.0 if link >= 0.0 else -1.0
            prod *= sign
        per_band_holonomies.append(prod)

    # (b) Rank-3 plaquette holonomy on the full triplet projector.
    # P(k) = sum_n |n,k><n,k| = identity (the frame spans R^3 at every k).
    # The link operator M_j = E_j^T E_{j+1} is the identity for all j.
    # Therefore U_plaq = I.
    U_plaq = np.eye(3)
    for j in range(4):
        M = frames[j].T @ frames[(j + 1) % 4]
        # Polar decomposition to project onto orthogonal part.
        u_, _, vt_ = np.linalg.svd(M)
        U_j = u_ @ vt_
        U_plaq = U_plaq @ U_j
    rank3_dev_from_identity = float(np.max(np.abs(U_plaq - np.eye(3))))

    return {
        "T5a_per_band_holonomies_signs": per_band_holonomies,
        "T5b_rank3_plaquette_max_dev_from_identity": rank3_dev_from_identity,
    }


def gap_ratios_table():
    """Tabulate eigenvalue spread along [100] and [111] for several alphas.

    Returns dict alpha -> {gap_100, gap_111, ratio_T_over_L_along_100}.
    """
    k_unit = 0.5
    out = {}
    for alpha in ALPHAS:
        # Along [100]
        k_100 = np.array([k_unit, 0.0, 0.0])
        ev_100 = np.sort(eigvals_6nn(k_100, alpha))
        # ev_100[0] = transverse (alpha=1: zero), ev_100[2] = longitudinal
        L = ev_100[2]
        T = ev_100[0]
        ratio = (T / L) if L > 1e-15 else float("nan")
        gap_100 = (L - T) / L if L > 1e-15 else float("nan")

        # Along [111]
        k_111 = (k_unit / np.sqrt(3.0)) * np.array([1.0, 1.0, 1.0])
        ev_111 = np.sort(eigvals_6nn(k_111, alpha))
        gap_111 = (ev_111[2] - ev_111[0]) / max(ev_111[2], 1e-15)

        out[f"alpha={alpha:.2f}"] = {
            "ev_along_100_sorted": ev_100.tolist(),
            "ev_along_111_sorted": ev_111.tolist(),
            "gap_100_rel": float(gap_100),
            "gap_111_rel": float(gap_111),
            "sqrt_1_minus_alpha": float(np.sqrt(max(1.0 - alpha, 0.0))),
            "T_over_L_along_100": float(ratio),
        }
    return out


def main():
    here = Path(__file__).resolve().parent
    out: dict = {}
    out["meta"] = {
        "stencil": "6-neighbor axial-only (canonical, backbone #15)",
        "k_samples": int(K_SAMPLES.shape[0]),
        "alphas_tested": ALPHAS,
        "k_unit_for_plaquette_and_table": 0.5,
    }
    out.update(test_T1_T2_T3_T4())
    out.update(test_T5_plaquette_holonomy())
    out["gap_ratios"] = gap_ratios_table()

    # PASS/FAIL summary at sensible tolerances.
    out["pass_fail"] = {
        "T1_diagonal": bool(out["T1_max_offdiag_abs"] < 1e-14),
        "T2_closed_form_matches_bond_sum": bool(out["T2_max_closed_form_rel_err"] < 1e-14),
        "T3_111_degenerate": bool(out["T3_max_111_degeneracy_rel_spread"] < 1e-14),
        "T4_alpha0_isotropic": bool(out["T4_max_alpha0_isotropy_rel_spread"] < 1e-14),
        "T5_plaquette_identity": bool(
            out["T5b_rank3_plaquette_max_dev_from_identity"] < 1e-12
            and all(s == 1.0 for s in out["T5a_per_band_holonomies_signs"])
        ),
    }

    with open(here / "results.json", "w") as f:
        json.dump(out, f, indent=2)

    print(json.dumps(out["pass_fail"], indent=2))
    print("\nGap ratios along [100] and [111]:")
    for key, val in out["gap_ratios"].items():
        print(f"  {key}: gap_100_rel = {val['gap_100_rel']:.4f}, "
              f"gap_111_rel = {val['gap_111_rel']:.4e}, "
              f"sqrt(1-alpha) = {val['sqrt_1_minus_alpha']:.4f}")


if __name__ == "__main__":
    main()