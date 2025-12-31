#!/usr/bin/env python3
"""
Option 2 (St. Venant–Kirchhoff) estimator WITHOUT scalar Monge reduction.

We keep a full spherical embedding in R^4:
  X(r,theta,phi,t) = ( R(r,t) * rhat(theta,phi),  w(r,t) )

This allows lateral contraction/expansion (via R) to couple to X^4 waves (via w)
through the induced metric, while remaining rotation invariant.

Reference metric:
  g0 = alpha^2 * diag(1, r^2, r^2 sin^2 theta),   0 < alpha < 1

Constitutive law:
  E = 0.5*(C - I) with C = g0^{-1} g
  W = mu tr(E^2) + 0.5*lambda (tr E)^2

Ansatz family (search parameters):
  w(t,r) = A * f(r;a) * cos(omega t)
  R(r)   = r * (1 - beta * f(r;a))     (static radial contraction profile)
where f is Gaussian f=exp(-r^2/(2a^2)), so R(0)=0.

We enforce omega = omega_C. We do NOT force a=lambda_C during optimization; we let it vary.

Mode A scaling:
  - choose (alpha, nu, eps=A/a, beta, a)
  - calibrate mu/rho so that small-slope transverse w-waves have speed c in vacuum
  - compute E1 for rho=1 and scale rho = E0/E1, and mu,lambda scale with rho
  - report residuals for w equation at t=0 and for R equilibrium (time-averaged)

Dependencies: numpy, scipy
"""

from __future__ import annotations
import math
import numpy as np
from scipy.optimize import minimize

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
# Elastic helpers
# ----------------------------
def lame_from_mu_nu(mu: float, nu: float) -> float:
    if not (0.0 < nu < 0.5):
        raise ValueError("nu must be in (0,0.5).")
    return 2.0 * mu * nu / (1.0 - 2.0 * nu)

def F_alpha_nu_for_w_vacuum(alpha: float, nu: float) -> float:
    """
    Same factor as in the w-only vacuum linearization:
      K_w = mu * F(alpha,nu),  c^2 = K_w / rho
    This is a calibration choice: we set the small-slope w-wave speed to c in the far field.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0,1).")
    ar = alpha**(-2)
    return (ar - 1.0) * ar * ((1.0 + nu) / (1.0 - 2.0 * nu))

# ----------------------------
# Ansatz shapes
# ----------------------------
def f_gauss(r: np.ndarray, a: float) -> np.ndarray:
    return np.exp(-0.5 * (r/a)**2)

def f_gauss_dr(r: np.ndarray, a: float) -> np.ndarray:
    return -(r/(a*a)) * np.exp(-0.5 * (r/a)**2)

def R_profile(r: np.ndarray, a: float, beta: float) -> np.ndarray:
    # R = r (1 - beta f)
    f = f_gauss(r, a)
    return r * (1.0 - beta * f)

def R_profile_dr(r: np.ndarray, a: float, beta: float) -> np.ndarray:
    # R' = (1 - beta f) + r * (-beta f')
    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)
    return (1.0 - beta * f) + r * (-beta * fp)

# ----------------------------
# Strain invariants in spherical symmetry
# ----------------------------
def strain_components(
    r: float, R: float, Rp: float, wr: float, alpha: float
):
    """
    Returns (E_rr, E_t, trE) where E_t is the tangential strain eigenvalue (double).
    """
    ar = alpha**(-2)
    # avoid division by zero at r=0 by using the limit R/r -> Rp at r->0 for smooth R~Rp*r
    if r < 1e-20:
        Rt_over_r = Rp
    else:
        Rt_over_r = R / r

    Crr = ar * (Rp*Rp + wr*wr)
    Ctt = ar * (Rt_over_r * Rt_over_r)

    Err = 0.5 * (Crr - 1.0)
    Et  = 0.5 * (Ctt - 1.0)
    trE = Err + 2.0 * Et
    return Err, Et, trE

def W_stvk(Err: float, Et: float, trE: float, mu: float, lam: float) -> float:
    trE2 = Err*Err + 2.0*Et*Et
    return mu * trE2 + 0.5 * lam * (trE*trE)

def dW_dwr(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂wr = ar * wr * (2 mu Err + lam trE)
    ar = alpha**(-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    S = (2.0 * mu * Err + lam * trE)
    return ar * wr * S

def dW_dRp(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂Rp = ar * Rp * (2 mu Err + lam trE)
    ar = alpha**(-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    S = (2.0 * mu * Err + lam * trE)
    return ar * Rp * S

def dW_dR(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂R = (2 ar R / r^2) * (2 mu Et + lam trE)
    ar = alpha**(-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    T = (2.0 * mu * Et + lam * trE)
    if r < 1e-20:
        # use limit: R/r -> Rp and (R/r^2) ~ Rp/r diverges, but equilibrium uses combined terms.
        # numerically, skip r=0 by starting at small r in residual computations.
        return 0.0
    return (2.0 * ar * R / (r*r)) * T

# ----------------------------
# Energy (time averaged), subtract baseline W0 (vacuum)
# ----------------------------
def total_energy_time_averaged(
    rho: float, alpha: float, mu: float, lam: float,
    eps: float, beta: float, a: float,
    Rmax_factor: float = 10.0, nr: int = 4000, nt: int = 64
) -> float:
    omega = omega_C
    A = eps * a

    Rmax = Rmax_factor * a
    r = np.linspace(0.0, Rmax, nr)

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    R = R_profile(r, a, beta)
    Rp = R_profile_dr(r, a, beta)

    # kinetic average for w: <0.5 rho wt^2> = 0.25 rho omega^2 A^2 f^2
    kin = 0.25 * rho * (omega**2) * (A**2) * (f*f)

    # baseline (vacuum): R=r, Rp=1, wr=0
    # In the far field u->0, this is the correct baseline.
    W0 = W_stvk(*strain_components(1.0, 1.0, 1.0, 0.0, alpha), mu, lam)  # r cancels in vacuum

    phases = np.linspace(0.0, 2.0*math.pi, nt, endpoint=False)
    Wacc = np.zeros_like(r)

    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph
        for i in range(nr):
            Err, Et, trE = strain_components(float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha)
            Wacc[i] += W_stvk(Err, Et, trE, mu, lam)

    Wavg = (Wacc / nt) - W0
    Wavg = np.maximum(Wavg, 0.0)

    integrand = 4.0*math.pi * r*r * (kin + Wavg)
    return float(np.trapezoid(integrand, r))

# ----------------------------
# Residuals:
#   w equation at t=0:  rho (-omega^2 w) - (1/r^2) d/dr (r^2 dW/dwr) = 0
#   R equilibrium (time-averaged): (1/r^2) d/dr (r^2 dW/dRp) - dW/dR = 0
# ----------------------------
def residuals(
    rho: float, alpha: float, mu: float, lam: float,
    eps: float, beta: float, a: float,
    Rmax_factor: float = 10.0, nr: int = 2500, nt_avg: int = 32
):
    omega = omega_C
    A = eps * a

    Rmax = Rmax_factor * a
    r = np.linspace(1e-12, Rmax, nr)  # skip 0 for the R-derivative term

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    R = R_profile(r, a, beta)
    Rp = R_profile_dr(r, a, beta)

    # --- w residual at t=0 (cos=1) ---
    w = A * f
    wr0 = A * fp

    pw = np.array([dW_dwr(float(r[i]), float(R[i]), float(Rp[i]), float(wr0[i]), alpha, mu, lam) for i in range(nr)])
    rp2pw = (r*r) * pw
    dr = r[1] - r[0]
    d_rp2pw = np.gradient(rp2pw, dr)
    div_pw = d_rp2pw / (r*r)

    res_w = rho * (-omega**2 * w) - div_pw
    scale_w = rho * (omega**2) * np.maximum(np.abs(w), 1e-30)
    rel_w = res_w / scale_w
    rel_w_rms = math.sqrt(float(np.trapezoid((rel_w*rel_w)*4*math.pi*r*r, r)) /
                          float(np.trapezoid(4*math.pi*r*r, r)))

    # --- R equilibrium residual, time-averaged over phase ---
    phases = np.linspace(0.0, 2.0*math.pi, nt_avg, endpoint=False)
    pRp_acc = np.zeros_like(r)
    dR_acc = np.zeros_like(r)
    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph
        for i in range(nr):
            pRp_acc[i] += dW_dRp(float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha, mu, lam)
            dR_acc[i]  += dW_dR (float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha, mu, lam)

    pRp = pRp_acc / nt_avg
    dR  = dR_acc / nt_avg

    rp2pRp = (r*r) * pRp
    d_rp2pRp = np.gradient(rp2pRp, dr)
    div_pRp = d_rp2pRp / (r*r)

    res_R = div_pRp - dR
    scale_R = np.maximum(np.abs(div_pRp) + np.abs(dR), 1e-30)
    rel_R = res_R / scale_R
    rel_R_rms = math.sqrt(float(np.trapezoid((rel_R*rel_R)*4*math.pi*r*r, r)) /
                          float(np.trapezoid(4*math.pi*r*r, r)))

    return rel_w_rms, rel_R_rms

# ----------------------------
# Mode A scaling: enforce c and E0 exactly
# ----------------------------
def mode_A_scale(alpha: float, nu: float, eps: float, beta: float, a: float):
    # calibrate mu/rho so that small-slope w-wave speed in vacuum is c
    F = F_alpha_nu_for_w_vacuum(alpha, nu)
    rho1 = 1.0
    mu1 = rho1 * c**2 / F
    lam1 = lame_from_mu_nu(mu1, nu)

    # energy with rho=1
    E1 = total_energy_time_averaged(
        rho=rho1, alpha=alpha, mu=mu1, lam=lam1,
        eps=eps, beta=beta, a=a
    )
    if E1 <= 0.0:
        return None

    rho = E0 / E1
    mu = rho * mu1
    lam = rho * lam1
    return rho, mu, lam, E1

# ----------------------------
# Optimization wrapper
# ----------------------------
def pack(alpha, nu, eps, beta, a):
    return np.array([
        math.log(alpha/(1.0-alpha)),                 # alpha in (0,1)
        math.log((nu-0.01)/(0.49-nu)),               # nu in (0.01,0.49)
        math.log(eps),                               # eps > 0
        math.asinh(beta),                            # beta unbounded
        math.log(a / lambda_C),                      # a > 0 relative
    ], dtype=float)

def unpack(x):
    # alpha in (1e-4, 0.999) (allow near 0)
    amin, amax = 1e-4, 0.999
    sig_a = 1.0/(1.0+math.exp(-x[0]))
    alpha = amin + (amax-amin)*sig_a

    # nu in (0.01, 0.49)
    nmin, nmax = 0.01, 0.49
    sig_n = 1.0/(1.0+math.exp(-x[1]))
    nu = nmin + (nmax-nmin)*sig_n

    eps = math.exp(x[2])
    beta = math.sinh(x[3])
    a = lambda_C * math.exp(x[4])
    return alpha, nu, eps, beta, a

def objective(x):
    alpha, nu, eps, beta, a = unpack(x)

    scaled = mode_A_scale(alpha, nu, eps, beta, a)
    if scaled is None:
        return 1e9
    rho, mu, lam, E1 = scaled

    rel_w, rel_R = residuals(rho, alpha, mu, lam, eps, beta, a)

    # penalties to keep things sane
    pen = 0.0
    # keep a not too extreme (soft)
    pen += 0.05 * (math.log(a/lambda_C)**2)
    # discourage absurd eps (soft)
    pen += 0.05 * (max(0.0, eps-1.0)**2)
    # discourage huge beta (soft)
    pen += 0.02 * (beta**2)

    # combine residuals
    return (rel_w + rel_R) + pen

if __name__ == "__main__":
    # Start closer to your intuition: alpha small-ish, not near 1
    x0 = pack(alpha=0.2, nu=0.25, eps=0.25, beta=0.1, a=lambda_C)

    res = minimize(objective, x0, method="Powell", options={"maxiter": 80, "disp": True})

    alpha, nu, eps, beta, a = unpack(res.x)
    scaled = mode_A_scale(alpha, nu, eps, beta, a)
    rho, mu, lam, E1 = scaled
    rel_w, rel_R = residuals(rho, alpha, mu, lam, eps, beta, a)

    print("\nOPTIMUM:")
    print(f"  alpha={alpha:.6e}   (reference metric scale)")
    print(f"  nu={nu:.6e}        (Poisson ratio)")
    print(f"  eps={eps:.6e}       (A/a)")
    print(f"  beta={beta:.6e}     (radial contraction strength)")
    print(f"  a={a:.6e}           (width, a/lambda_C={a/lambda_C:.6e})")
    print(f"  rho={rho:.6e} kg/m^3")
    print(f"  mu={mu:.6e} Pa")
    print(f"  lambda={lam:.6e} Pa")
    print(f"  E1(rho=1)={E1:.6e} J")
    print(f"  rel_w_rms={rel_w:.6e}")
    print(f"  rel_R_rms={rel_R:.6e}")
    print(f"  objective={objective(res.x):.6e}")
