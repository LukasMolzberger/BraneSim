#!/usr/bin/env python3
"""
Option-1 (Cauchy–Born) parameter estimation for the brane-electron test mode.

Model:
  X(q,t) = (q, w(q,t)) in R^4, q in R^3
  W(grad w) = (1/Vcell) sum_a 0.5*k_a * ( sqrt(|a|^2 + (grad w·a)^2 ) - alpha|a| )^2

Electron test ansatz:
  w(t,r) = A * exp(-r^2/(2 a^2)) * cos(omega t)

Pipeline:
  - Choose h, alpha, bond weights
  - Calibrate k0 from wave speed c (linearization about grad w = 0)
  - Either:
      Mode A (recommended): choose eps = A/a, solve rho from E = m_e c^2
      Mode B: fix rho, solve A from E = m_e c^2
  - Compute a self-consistency residual norm for the ansatz at t=0
"""

from __future__ import annotations
import math
import numpy as np

# ----------------------------
# Physical constants (SI)
# ----------------------------
c = 299_792_458.0
hbar = 1.054_571_817e-34
m_e = 9.109_383_7015e-31

E0 = m_e * c**2
omega_C = E0 / hbar
lambda_C = hbar / (m_e * c)

# ----------------------------
# Bond set (unique, no double-count)
# ----------------------------
def bond_set(h: float):
    """
    Return list of bonds a (3-vector) and category labels: 'axis','face','body'
    Unique half-set sufficient for energy density (no double counting).
    """
    bonds = []

    # 3 axis bonds
    bonds += [(np.array([h, 0, 0], dtype=float), "axis")]
    bonds += [(np.array([0, h, 0], dtype=float), "axis")]
    bonds += [(np.array([0, 0, h], dtype=float), "axis")]

    # 6 face diagonals: (h, ±h, 0) etc
    bonds += [(np.array([h,  h, 0], dtype=float), "face")]
    bonds += [(np.array([h, -h, 0], dtype=float), "face")]
    bonds += [(np.array([h, 0,  h], dtype=float), "face")]
    bonds += [(np.array([h, 0, -h], dtype=float), "face")]
    bonds += [(np.array([0,  h,  h], dtype=float), "face")]
    bonds += [(np.array([0,  h, -h], dtype=float), "face")]

    # 4 body diagonals with x positive: (h, ±h, ±h)
    for sy in (+1, -1):
        for sz in (+1, -1):
            bonds += [(np.array([h, sy*h, sz*h], dtype=float), "body")]

    return bonds

# ----------------------------
# Linear calibration: k0 from c
# ----------------------------
def k0_from_c(rho: float, h: float, alpha: float, w_axis: float, w_face: float, w_body: float) -> float:
    """
    Compute base spring constant k0 such that small-slope wave speed is c.
    Uses linearization:
      c^2 = K/rho,  K = (1-alpha)/Vcell * (1/3) sum_a k_a |a|^2
      k_a = k0 * w_cat
    """
    if not (alpha < 1.0):
        raise ValueError("Need alpha < 1 for linear (pre-tension) stiffness around grad w = 0.")

    Vcell = h**3
    bonds = bond_set(h)

    # S = sum_a w_cat |a|^2
    S = 0.0
    for a, cat in bonds:
        L2 = float(a @ a)
        w = w_axis if cat == "axis" else (w_face if cat == "face" else w_body)
        S += w * L2

    # k0
    k0 = (3.0 * rho * c**2 * Vcell) / ((1.0 - alpha) * S)
    return k0

# ----------------------------
# Energy density W for given scalar radial gradient wr = dw/dr and time phase
# ----------------------------
def W_density_from_wr(
    wr: float,
    h: float,
    alpha: float,
    k0: float,
    w_axis: float,
    w_face: float,
    w_body: float,
    n_mu: int = 48,
) -> float:
    """
    Compute W for spherical symmetry where grad w = wr * r_hat.
    We must average over the orientation between r_hat and each bond direction.
    For a given bond vector a with length L, r_hat·a = L * mu, mu in [-1,1] uniform.
    Then:
      ell = sqrt(L^2 + (wr * L * mu)^2) = L*sqrt(1 + (wr*mu)^2)
      delta = ell - alpha*L
      U = 0.5*k_a*delta^2
    W = (1/Vcell) * sum_a <U>_mu
    """
    Vcell = h**3
    bonds = bond_set(h)

    # Gauss-Legendre quadrature for mu in [-1,1]
    mu, w_mu = np.polynomial.legendre.leggauss(n_mu)

    W_sum = 0.0
    for a, cat in bonds:
        L = float(np.linalg.norm(a))
        wcat = w_axis if cat == "axis" else (w_face if cat == "face" else w_body)
        k_a = k0 * wcat

        # ell(mu) = L*sqrt(1 + (wr*mu)^2)
        s = wr * mu
        ell = L * np.sqrt(1.0 + s*s)
        delta = ell - alpha * L
        U = 0.5 * k_a * (delta*delta)

        # average over mu (uniform): (1/2)∫_{-1}^1 U dmu
        U_avg = 0.5 * float(np.sum(w_mu * U))
        W_sum += U_avg

    return W_sum / Vcell

# ----------------------------
# Electron ansatz and energy integration
# ----------------------------
def gaussian_profile(r: np.ndarray, a: float) -> np.ndarray:
    return np.exp(-0.5 * (r/a)**2)

def gaussian_profile_dr(r: np.ndarray, a: float) -> np.ndarray:
    # d/dr exp(-r^2/(2a^2)) = -(r/a^2)*exp(...)
    return -(r/(a*a)) * np.exp(-0.5 * (r/a)**2)

def total_energy_time_averaged(
    rho: float,
    h: float,
    alpha: float,
    k0: float,
    w_axis: float,
    w_face: float,
    w_body: float,
    A: float,
    a: float,
    omega: float,
    Rmax_factor: float = 10.0,
    nr: int = 4000,
    nt: int = 48,
    n_mu: int = 48,
) -> float:
    """
    Compute time-averaged total energy for w(t,r) = A*f(r)*cos(omega t), f Gaussian.

    E = ∫ 4π r^2 [ <0.5 rho (wt)^2> + <W(grad w)> ] dr
    Kinetic average done analytically. Potential average sampled in time.
    """
    Rmax = Rmax_factor * a
    r = np.linspace(0.0, Rmax, nr)
    dr = r[1] - r[0]

    f = gaussian_profile(r, a)
    fp = gaussian_profile_dr(r, a)

    # Kinetic average: <0.5 rho (wt)^2> = 0.25 rho omega^2 A^2 f^2
    kin = 0.25 * rho * (omega**2) * (A**2) * (f*f)

    # Potential average: sample time phases
    phases = np.linspace(0.0, 2.0*math.pi, nt, endpoint=False)
    W_accum = np.zeros_like(r)
    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph  # radial derivative
        # compute W at each r point
        # (vectorize over r by looping; keep it simple/robust)
        for i in range(nr):
            W_accum[i] += W_density_from_wr(
                wr=float(wr[i]),
                h=h,
                alpha=alpha,
                k0=k0,
                w_axis=w_axis,
                w_face=w_face,
                w_body=w_body,
                n_mu=n_mu,
            )

    W_avg = W_accum / nt

    integrand = 4.0*math.pi * r*r * (kin + W_avg)
    # trapezoidal integration
    E = float(np.trapz(integrand, r))
    return E

# ----------------------------
# Residual norm for "is this an approximate solution?"
# ----------------------------
def p_from_wr_numeric(
    wr: float,
    h: float,
    alpha: float,
    k0: float,
    w_axis: float,
    w_face: float,
    w_body: float,
    n_mu: int = 48,
    dwr: float = 1e-6,
) -> float:
    """
    p(wr) = dW/dwr (instantaneous, no time average) for spherical symmetry.
    Finite difference derivative of W_density_from_wr.
    """
    Wp = W_density_from_wr(wr + dwr, h, alpha, k0, w_axis, w_face, w_body, n_mu=n_mu)
    Wm = W_density_from_wr(wr - dwr, h, alpha, k0, w_axis, w_face, w_body, n_mu=n_mu)
    return (Wp - Wm) / (2.0*dwr)

def residual_norm_at_t0(
    rho: float,
    h: float,
    alpha: float,
    k0: float,
    w_axis: float,
    w_face: float,
    w_body: float,
    A: float,
    a: float,
    omega: float,
    Rmax_factor: float = 10.0,
    nr: int = 2000,
    n_mu: int = 48,
) -> float:
    """
    Compute an L2-like residual norm for the ansatz at t=0 (cos=1):
      res(r) = rho*(-omega^2 w) - (1/r^2) d/dr (r^2 p(wr))
    where p(wr) = dW/dwr (instantaneous).

    Returns sqrt( ∫ res(r)^2 4π r^2 dr / ∫ 4π r^2 dr ) for normalization.
    """
    Rmax = Rmax_factor * a
    r = np.linspace(1e-12, Rmax, nr)  # avoid r=0 singularity
    f = gaussian_profile(r, a)
    fp = gaussian_profile_dr(r, a)

    w = A * f
    wr = A * fp

    # compute p(wr) pointwise
    p = np.zeros_like(r)
    for i in range(nr):
        p[i] = p_from_wr_numeric(
            wr=float(wr[i]),
            h=h,
            alpha=alpha,
            k0=k0,
            w_axis=w_axis,
            w_face=w_face,
            w_body=w_body,
            n_mu=n_mu,
        )

    # divergence in spherical: (1/r^2) d/dr (r^2 p)
    rp2p = (r*r) * p
    dr = r[1] - r[0]
    d_rp2p = np.gradient(rp2p, dr)
    div = d_rp2p / (r*r)

    res = rho * (-omega**2 * w) - div

    wgt = 4.0*math.pi * r*r
    num = float(np.trapz((res*res)*wgt, r))
    den = float(np.trapz(wgt, r))
    return math.sqrt(num/den)

# ----------------------------
# Main estimation routines
# ----------------------------
def estimate_mode_A(
    alpha: float,
    eps: float,
    N: int,
    w_axis: float = 1.0,
    w_face: float = 1.0,
    w_body: float = 1.0,
):
    """
    Mode A: Choose eps=A/a and solve rho from E=E0, with a=lambda_C and omega=omega_C.
    """
    a = lambda_C
    omega = omega_C
    h = a / N
    A = eps * a

    # Compute k0 for rho=1
    rho1 = 1.0
    k0_1 = k0_from_c(rho=rho1, h=h, alpha=alpha, w_axis=w_axis, w_face=w_face, w_body=w_body)

    # Energy for rho=1
    E1 = total_energy_time_averaged(
        rho=rho1, h=h, alpha=alpha, k0=k0_1,
        w_axis=w_axis, w_face=w_face, w_body=w_body,
        A=A, a=a, omega=omega,
        nr=1500, nt=32, n_mu=32, Rmax_factor=10.0
    )

    rho = E0 / E1
    k0 = rho * k0_1

    res = residual_norm_at_t0(
        rho=rho, h=h, alpha=alpha, k0=k0,
        w_axis=w_axis, w_face=w_face, w_body=w_body,
        A=A, a=a, omega=omega,
        nr=1200, n_mu=24, Rmax_factor=10.0
    )

    return {
        "mode": "A",
        "alpha": alpha,
        "eps": eps,
        "N": N,
        "h": h,
        "a": a,
        "A": A,
        "omega": omega,
        "rho": rho,
        "k0": k0,
        "E_check": E0,
        "residual_norm": res,
    }

def estimate_sweep_alpha(
    alphas: np.ndarray,
    eps: float,
    N: int,
    w_axis: float = 1.0,
    w_face: float = 1.0,
    w_body: float = 1.0,
):
    out = []
    for alpha in alphas:
        out.append(estimate_mode_A(alpha=float(alpha), eps=eps, N=N, w_axis=w_axis, w_face=w_face, w_body=w_body))
        print(f"alpha={alpha:.4f} -> rho={out[-1]['rho']:.3e}, k0={out[-1]['k0']:.3e}, residual={out[-1]['residual_norm']:.3e}")
    return out

if __name__ == "__main__":
    # Example run
    # You should sweep alpha and eps. alpha controls pre-tension/coupling; eps controls amplitude scale.
    alphas = np.linspace(0.90, 0.995, 6)   # keep < 1
    eps = 0.25                             # A/a
    N = 30                                 # h = lambda_C/N

    results = estimate_sweep_alpha(alphas, eps=eps, N=N, w_axis=1.0, w_face=1.0, w_body=1.0)

    best = min(results, key=lambda d: d["residual_norm"])
    print("\nBEST (min residual):")
    for k, v in best.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6e}")
        else:
            print(f"  {k}: {v}")
