#!/usr/bin/env python3
"""
Option 2 (St. Venant–Kirchhoff) parameter estimator for a brane-electron test ansatz.

No microstructure: no h, no bond sets.

Geometry:
  X(q,t) = (q, w(q,t)) in R^4, q in R^3
  g = I + grad w (grad w)^T
Reference metric:
  g0 = alpha^2 I,   0 < alpha < 1   (pre-stress / pre-tension analogue)
Strain:
  C = g0^{-1} g = alpha^{-2} g
  E = 1/2 (C - I)
StVK energy:
  W = mu tr(E^2) + (lambda/2)(tr E)^2
with lambda = 2 mu nu / (1 - 2 nu) for isotropic 3D material.

Ansatz:
  w(t,r) = A exp(-r^2/(2a^2)) cos(omega t)
with a=lambda_C, omega=omega_C.

Calibration:
  small-slope transverse wave speed:  c^2 = K/rho  with K = mu * F(alpha,nu)
  => for rho=1: mu1 = c^2 / F(alpha,nu)

Energy normalization (Mode A):
  choose eps = A/a
  compute E1 with rho=1 (and mu1, lambda1)
  set rho = E0 / E1, and then mu = rho*mu1, lambda = rho*lambda1

Reports:
  - rho, mu, lambda
  - dimensionless relative residual at t=0
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
# StVK helpers
# ----------------------------
def lame_from_mu_nu(mu: float, nu: float) -> float:
    if not (0.0 < nu < 0.5):
        raise ValueError("Poisson ratio nu must be in (0, 0.5).")
    return 2.0 * mu * nu / (1.0 - 2.0 * nu)

def F_alpha_nu(alpha: float, nu: float) -> float:
    """
    Factor such that K = mu * F(alpha,nu) is the quadratic coefficient in W ~ W0 + 0.5 K |grad w|^2
    for the Monge-gauge transverse field with reference metric g0=alpha^2 I.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0,1) to produce pre-stress and a finite linear wave speed.")
    ar = alpha**(-2)  # alpha^{-2}
    return (ar - 1.0) * ar * ((1.0 + nu) / (1.0 - 2.0 * nu))

def W_stvk_from_wr(wr: float, alpha: float, mu: float, lam: float) -> float:
    """
    For spherical symmetry grad w = wr r_hat, the induced metric has eigenvalues:
      g: (1+wr^2, 1, 1)
    With g0=alpha^2 I => C = alpha^{-2} g has eigenvalues:
      C: (ar(1+wr^2), ar, ar), ar=alpha^{-2}
    Strain eigenvalues:
      e_r = 0.5*(ar(1+wr^2) - 1)
      e_t = 0.5*(ar - 1)   (twice)
    Then
      trE = e_r + 2 e_t
      trE2 = e_r^2 + 2 e_t^2
      W = mu trE2 + 0.5 lam (trE)^2
    """
    ar = alpha**(-2)
    et = 0.5 * (ar - 1.0)
    er = et + 0.5 * ar * (wr*wr)

    trE = er + 2.0 * et
    trE2 = er*er + 2.0 * et*et
    return mu * trE2 + 0.5 * lam * (trE*trE)

def p_stvk_from_wr(wr: float, alpha: float, mu: float, lam: float) -> float:
    """
    p(wr) = dW/dwr for spherical symmetry.
    Uses analytic derivative.

    Using x = wr^2:
      dW/dwr = 2 wr * dW/dx
    """
    ar = alpha**(-2)
    et = 0.5 * (ar - 1.0)
    x = wr*wr

    # trE = 3 et + 0.5 ar x
    trE = 3.0*et + 0.5*ar*x

    # d(trE)/dx = 0.5 ar
    # d(trE^2)/dx = 2 trE * dtrE/dx = trE * ar
    # trE2 = 3 et^2 + et ar x + 0.25 ar^2 x^2
    # d(trE2)/dx = et ar + 0.5 ar^2 x
    dWdx = mu * (et*ar + 0.5*(ar*ar)*x) + 0.5*lam * (trE * ar)

    return 2.0 * wr * dWdx

# ----------------------------
# Ansatz profile
# ----------------------------
def f_gauss(r: np.ndarray, a: float) -> np.ndarray:
    return np.exp(-0.5*(r/a)**2)

def f_gauss_dr(r: np.ndarray, a: float) -> np.ndarray:
    return -(r/(a*a)) * np.exp(-0.5*(r/a)**2)

# ----------------------------
# Energy integral (time-averaged), using EXCESS W - W(0)
# ----------------------------
def total_energy_time_averaged(
    rho: float,
    alpha: float,
    mu: float,
    lam: float,
    A: float,
    a: float,
    omega: float,
    Rmax_factor: float = 10.0,
    nr: int = 4000,
    nt: int = 64,
) -> float:
    Rmax = Rmax_factor * a
    r = np.linspace(0.0, Rmax, nr)

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    # kinetic average: <0.5 rho wt^2> = 0.25 rho omega^2 A^2 f^2
    kin = 0.25 * rho * (omega**2) * (A**2) * (f*f)

    # baseline W0 at wr=0
    W0 = W_stvk_from_wr(0.0, alpha, mu, lam)

    # time-average potential by sampling phase
    phases = np.linspace(0.0, 2.0*math.pi, nt, endpoint=False)
    Wacc = np.zeros_like(r)
    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph
        # vectorize W(wr)
        Wvals = np.array([W_stvk_from_wr(float(wri), alpha, mu, lam) for wri in wr], dtype=float)
        Wacc += Wvals

    Wavg = (Wacc / nt) - W0
    Wavg = np.maximum(Wavg, 0.0)

    integrand = 4.0*math.pi * r*r * (kin + Wavg)
    return float(np.trapezoid(integrand, r))

# ----------------------------
# Dimensionless relative residual at t=0
# ----------------------------
def relative_residual_rms_at_t0(
    rho: float,
    alpha: float,
    mu: float,
    lam: float,
    A: float,
    a: float,
    omega: float,
    Rmax_factor: float = 10.0,
    nr: int = 2500,
) -> float:
    Rmax = Rmax_factor * a
    r = np.linspace(1e-12, Rmax, nr)

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    w = A * f
    wr = A * fp

    # p(wr) pointwise
    p = np.array([p_stvk_from_wr(float(wri), alpha, mu, lam) for wri in wr], dtype=float)

    # divergence in spherical: (1/r^2) d/dr (r^2 p)
    rp2p = (r*r) * p
    dr = r[1] - r[0]
    d_rp2p = np.gradient(rp2p, dr)
    div = d_rp2p / (r*r)

    res = rho * (-omega**2 * w) - div

    scale = rho * (omega**2) * np.maximum(np.abs(w), 1e-30)
    rel = res / scale

    wgt = 4.0*math.pi * r*r
    num = float(np.trapezoid((rel*rel)*wgt, r))
    den = float(np.trapezoid(wgt, r))
    return math.sqrt(num/den)

# ----------------------------
# Mode A estimator
# ----------------------------
def estimate_mode_A(alpha: float, nu: float, eps: float, Rmax_factor: float = 10.0):
    a = lambda_C
    omega = omega_C
    A = eps * a

    # Calibrate mu for rho=1 so that small-slope wave speed is c:
    # K = mu * F(alpha,nu), c^2 = K/rho
    F = F_alpha_nu(alpha, nu)
    rho1 = 1.0
    mu1 = rho1 * c**2 / F
    lam1 = lame_from_mu_nu(mu1, nu)

    # Energy for rho=1
    E1 = total_energy_time_averaged(
        rho=rho1, alpha=alpha, mu=mu1, lam=lam1,
        A=A, a=a, omega=omega,
        Rmax_factor=Rmax_factor, nr=3500, nt=64
    )
    if E1 <= 0.0:
        raise RuntimeError(f"E1 <= 0 (E1={E1}). Try larger eps or different alpha/nu.")

    # Scale rho to hit E0, and scale mu,lam linearly with rho (keeps c fixed because mu/rho stays constant)
    rho = E0 / E1
    mu = rho * mu1
    lam = rho * lam1

    relres = relative_residual_rms_at_t0(
        rho=rho, alpha=alpha, mu=mu, lam=lam,
        A=A, a=a, omega=omega,
        Rmax_factor=Rmax_factor, nr=2000
    )

    return {
        "alpha": alpha,
        "nu": nu,
        "eps": eps,
        "a": a,
        "A": A,
        "omega": omega,
        "rho": rho,
        "mu": mu,
        "lambda": lam,
        "E_target": E0,
        "E1_rho_eq_1": E1,
        "relative_residual_rms": relres,
        "F(alpha,nu)": F,
    }

if __name__ == "__main__":
    # Example sweep
    alphas = np.linspace(0.90, 0.99, 6)  # must be < 1
    nu = 0.25
    eps = 0.25

    results = []
    for alpha in alphas:
        d = estimate_mode_A(float(alpha), nu=nu, eps=eps, Rmax_factor=10.0)
        results.append(d)
        print(
            f"alpha={d['alpha']:.6f} nu={d['nu']:.3f} eps={d['eps']:.3f} -> "
            f"rho={d['rho']:.3e}, mu={d['mu']:.3e} Pa, lambda={d['lambda']:.3e} Pa, "
            f"E1(rho=1)={d['E1_rho_eq_1']:.3e}, rel_res={d['relative_residual_rms']:.3e}"
        )

    best = min(results, key=lambda x: x["relative_residual_rms"])
    print("\nBEST (min relative residual):")
    for k, v in best.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.6e}")
        else:
            print(f"  {k}: {v}")
