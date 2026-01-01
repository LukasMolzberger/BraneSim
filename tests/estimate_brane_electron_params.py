#!/usr/bin/env python3
"""
Brane + electron parameter estimator (v2)

Key changes vs. v1:
  - Fix brane density rho (or optionally optimize it); do NOT rescale rho to hit E0.
  - Enforce the electron rest-energy by solving for eps (=A/a) inside the objective.
  - Add a selectable wave-speed calibration policy:
        S: shear speed c_T = c  (rigid substrate)
        W: w-branch speed c_w = c (emergent U(1) phase sector)
        L: longitudinal speed c_L = c
  - Enforce that the geometric (Pythagorean) nonlinearity is actually in play via a
    slope constraint |w_r|_max >= s_min (otherwise confinement mechanisms cannot act).
  - Use bounded optimization so alpha, nu, a, beta cannot run into nonsensical corners.

Model (as in v1):
  X(r,theta,phi,t) = ( R(r) * rhat(theta,phi),  w(r,t) )
  w(r,t) = A f(r;a) cos(omega_C t),  A = eps*a
  R(r)   = r * exp( - beta * r^2 / (2 a^2) )  (simple contraction family)

Reference metric:
  g0 = alpha^2 * diag(1, r^2, r^2 sin^2 theta)

Constitutive law (St. Venant–Kirchhoff):
  E = 0.5 (C - I), C = g0^{-1} g
  W = mu tr(E^2) + 0.5 lambda (tr E)^2

Residuals (as in v1):
  - w equation at phase=0 (cos=1):  rho (-omega^2 w) - (1/r^2) d/dr (r^2 dW/dwr) = 0
  - R equilibrium (time/angle averaged): (1/r^2) d/dr (r^2 dW/dRp) - dW/dR = 0

Energy:
  time-averaged total energy integral with vacuum subtraction (as in v1)

Usage example:
  python estimate_brane_electron_params_v2.py --policy S --alpha-min 0.9 --alpha-max 0.999 \\
      --nu-min 0.05 --nu-max 0.35 --s-min 0.3 --maxiter 60

The script prints the optimum and key diagnostics (including c_T, c_L, c_w).
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Dict, Tuple, Optional

import numpy as np
from scipy.optimize import minimize


# numpy integration helper: supports both numpy<2.0 (trapz) and numpy>=2.0 (trapezoid)
def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(_trapz(y, x))
    return float(np.trapz(y, x))
# ----------------------------
# Physical constants (SI)
# ----------------------------
c = 299_792_458.0  # m/s
hbar = 1.054_571_817e-34  # J s
m_e = 9.109_383_7015e-31  # kg
E0 = m_e * c * c  # J

lambda_C = hbar / (m_e * c)  # reduced Compton wavelength
omega_C = m_e * c * c / hbar  # Compton angular frequency


# ----------------------------
# Helpers: isotropic elastic parameters and speed policies
# ----------------------------
def lame_from_mu_nu(mu: float, nu: float) -> float:
    """lambda in terms of mu, nu for isotropic elasticity."""
    denom = (1.0 - 2.0 * nu)
    if denom <= 0.0:
        return float("nan")
    return 2.0 * mu * nu / denom

def speeds(mu: float, lam: float, rho: float, alpha: float, nu: float) -> Tuple[float, float, float]:
    """(cT, cL, cw) where cw is the vacuum small-slope w-branch speed from the pre-tension factor."""
    cT = math.sqrt(max(mu / rho, 0.0))
    cL = math.sqrt(max((lam + 2.0 * mu) / rho, 0.0))
    F = F_alpha_nu_for_w_vacuum(alpha, nu)
    cw = math.sqrt(max(mu * F / rho, 0.0))
    return cT, cL, cw

def mu_lam_from_policy(policy: str, rho: float, alpha: float, nu: float) -> Tuple[float, float]:
    """
    Return (mu, lambda) from a chosen calibration policy, with rho fixed.

      S: c_T = c  -> mu = rho c^2
      W: c_w = c  -> mu = rho c^2 / F(alpha,nu)
      L: c_L = c  -> lambda+2mu = rho c^2 and nu fixes split
    """
    policy = policy.upper()
    if policy == "S":
        mu = rho * c * c
        lam = lame_from_mu_nu(mu, nu)
        return mu, lam

    if policy == "W":
        F = F_alpha_nu_for_w_vacuum(alpha, nu)
        if not np.isfinite(F) or F <= 0.0:
            return float("nan"), float("nan")
        mu = rho * c * c / F
        lam = lame_from_mu_nu(mu, nu)
        return mu, lam

    if policy == "L":
        # nu => lam = 2 mu nu/(1-2nu), and lam+2mu = 2mu(1-nu)/(1-2nu)
        denom = (1.0 - nu)
        if denom <= 0.0:
            return float("nan"), float("nan")
        mu = (rho * c * c) * 0.5 * (1.0 - 2.0 * nu) / denom
        lam = lame_from_mu_nu(mu, nu)
        return mu, lam

    raise ValueError(f"Unknown policy '{policy}'. Use S, W, or L.")


# ----------------------------
# Brane ansatz functions
# ----------------------------
def f_gauss(r: np.ndarray, a: float) -> np.ndarray:
    return np.exp(-0.5 * (r / a) ** 2)

def f_gauss_dr(r: np.ndarray, a: float) -> np.ndarray:
    # d/dr exp(-r^2/(2a^2)) = -(r/a^2) exp(...)
    return -(r / (a * a)) * np.exp(-0.5 * (r / a) ** 2)

def R_profile(r: np.ndarray, a: float, beta: float) -> np.ndarray:
    # R = r exp(-beta r^2/(2 a^2))
    return r * np.exp(-0.5 * beta * (r / a) ** 2)

def R_profile_dr(r: np.ndarray, a: float, beta: float) -> np.ndarray:
    # derivative: exp(s)*(1 + r*s') with s = -beta r^2/(2a^2), s' = -beta r/a^2
    s = -0.5 * beta * (r / a) ** 2
    return np.exp(s) * (1.0 - beta * (r / a) ** 2)

# ----------------------------
# StVK strain components in the spherical embedding
# ----------------------------
def strain_components(r: float, R: float, Rp: float, wr: float, alpha: float) -> Tuple[float, float, float]:
    """
    C = g0^{-1} g in spherical symmetry:
      C_rr = alpha^{-2}(Rp^2 + wr^2)
      C_tt = alpha^{-2}(R/r)^2  (and same for phi)
    E = 0.5(C - I)
    Return (Err, Et, trE) with Et = E_theta = E_phi.
    """
    ar = alpha ** (-2)
    Crr = ar * (Rp * Rp + wr * wr)
    if r < 1e-30:
        # near origin use limit R/r -> Rp
        Ctt = ar * (Rp * Rp)
    else:
        Ctt = ar * (R / r) ** 2
    Err = 0.5 * (Crr - 1.0)
    Et = 0.5 * (Ctt - 1.0)
    trE = Err + 2.0 * Et
    return Err, Et, trE

def W_stvk(Err: float, Et: float, trE: float, mu: float, lam: float) -> float:
    trE2 = Err * Err + 2.0 * Et * Et
    return mu * trE2 + 0.5 * lam * (trE * trE)

def dW_dwr(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂wr = ar * wr * (2 mu Err + lam trE)
    ar = alpha ** (-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    S = (2.0 * mu * Err + lam * trE)
    return ar * wr * S

def dW_dRp(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂Rp = ar * Rp * (2 mu Err + lam trE)
    ar = alpha ** (-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    S = (2.0 * mu * Err + lam * trE)
    return ar * Rp * S

def dW_dR(r: float, R: float, Rp: float, wr: float, alpha: float, mu: float, lam: float) -> float:
    # ∂W/∂R = (2 ar R / r^2) * (2 mu Et + lam trE)
    ar = alpha ** (-2)
    Err, Et, trE = strain_components(r, R, Rp, wr, alpha)
    T = (2.0 * mu * Et + lam * trE)
    if r < 1e-20:
        return 0.0
    return (2.0 * ar * R / (r * r)) * T


# ----------------------------
# Pre-tension factor for the w-branch speed (vacuum, small slope)
# ----------------------------
def F_alpha_nu_for_w_vacuum(alpha: float, nu: float) -> float:
    """
    Factor F so that small-slope w satisfies rho w_tt = (mu F) Δ w in the far field.

    From v1: F = (alpha^{-2}-1) alpha^{-2} * (1+nu)/(1-2nu)

    NOTE: this diverges as alpha -> 0 or nu -> 0.5.
    """
    if alpha <= 0.0:
        return float("nan")
    if (1.0 - 2.0 * nu) <= 0.0:
        return float("nan")
    ar2 = alpha ** (-2)
    return (ar2 - 1.0) * ar2 * (1.0 + nu) / (1.0 - 2.0 * nu)


# ----------------------------
# Energy and residuals (copied from v1 with only light cleanup)
# ----------------------------
def total_energy_time_averaged(
    rho: float, alpha: float, mu: float, lam: float,
    eps: float, beta: float, a: float,
    Rmax_factor: float = 10.0, nr: int = 3000, nt: int = 48
) -> float:
    omega = omega_C
    A = eps * a

    Rmax = Rmax_factor * a
    r = np.linspace(0.0, Rmax, nr)

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    R = R_profile(r, a, beta)
    Rp = R_profile_dr(r, a, beta)

    # kinetic avg for w: <0.5 rho wt^2> = 0.25 rho omega^2 A^2 f^2
    kin = 0.25 * rho * (omega**2) * (A**2) * (f * f)

    # baseline (vacuum): R=r, Rp=1, wr=0  (r cancels for vacuum)
    W0 = W_stvk(*strain_components(1.0, 1.0, 1.0, 0.0, alpha), mu, lam)

    phases = np.linspace(0.0, 2.0 * math.pi, nt, endpoint=False)
    Wacc = np.zeros_like(r)

    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph
        # per-point W
        for i in range(nr):
            Err, Et, trE = strain_components(float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha)
            Wacc[i] += W_stvk(Err, Et, trE, mu, lam)

    # vacuum-subtracted time average (clipped to avoid tiny negative numerical noise)
    Wavg = (Wacc / nt) - W0
    Wavg = np.maximum(Wavg, 0.0)

    integrand = 4.0 * math.pi * r * r * (kin + Wavg)
    return float(_trapz(integrand, r))


def residuals(
    rho: float, alpha: float, mu: float, lam: float,
    eps: float, beta: float, a: float,
    Rmax_factor: float = 10.0, nr: int = 2200, nt_avg: int = 28
) -> Tuple[float, float]:
    omega = omega_C
    A = eps * a

    Rmax = Rmax_factor * a
    r = np.linspace(1e-12, Rmax, nr)  # skip 0 for R-derivative term
    dr = float(r[1] - r[0])

    f = f_gauss(r, a)
    fp = f_gauss_dr(r, a)

    R = R_profile(r, a, beta)
    Rp = R_profile_dr(r, a, beta)

    # --- w residual at phase=0 (cos=1) ---
    w = A * f
    wr0 = A * fp

    pw = np.array([dW_dwr(float(r[i]), float(R[i]), float(Rp[i]), float(wr0[i]), alpha, mu, lam)
                   for i in range(nr)], dtype=float)
    rp2pw = (r * r) * pw
    d_rp2pw = np.gradient(rp2pw, dr)
    div_pw = d_rp2pw / (r * r)

    res_w = rho * (-omega * omega) * w - div_pw
    scale_w = np.maximum(np.abs(rho * omega * omega * w) + np.abs(div_pw), 1e-30)
    rel_w = res_w / scale_w
    rel_w_rms = math.sqrt(float(_trapz((rel_w * rel_w) * 4 * math.pi * r * r, r)) /
                          float(_trapz(4 * math.pi * r * r, r)))

    # --- R residual: time/angle avg approximated by phase average only (as in v1) ---
    # The contraction ansatz is purely radial; the wr term carries the phase.
    phases = np.linspace(0.0, 2.0 * math.pi, nt_avg, endpoint=False)
    pRp_acc = np.zeros_like(r)
    dR_acc = np.zeros_like(r)

    for ph in phases:
        cosph = math.cos(ph)
        wr = (A * fp) * cosph
        for i in range(nr):
            pRp_acc[i] += dW_dRp(float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha, mu, lam)
            dR_acc[i] += dW_dR(float(r[i]), float(R[i]), float(Rp[i]), float(wr[i]), alpha, mu, lam)

    pRp = pRp_acc / nt_avg
    dR = dR_acc / nt_avg

    rp2pRp = (r * r) * pRp
    d_rp2pRp = np.gradient(rp2pRp, dr)
    div_pRp = d_rp2pRp / (r * r)

    res_R = div_pRp - dR
    scale_R = np.maximum(np.abs(div_pRp) + np.abs(dR), 1e-30)
    rel_R = res_R / scale_R
    rel_R_rms = math.sqrt(float(_trapz((rel_R * rel_R) * 4 * math.pi * r * r, r)) /
                          float(_trapz(4 * math.pi * r * r, r)))

    return rel_w_rms, rel_R_rms


# ----------------------------
# Solve eps from the energy constraint (fast quasi-quadratic solver)
# ----------------------------
@dataclass(frozen=True)
class EpsSolveResult:
    eps: float
    E_beta: float
    E_coeff: float
    E_final: float

def solve_eps_for_energy(
    rho: float, alpha: float, mu: float, lam: float,
    beta: float, a: float,
    target_E: float = E0,
    eps_max: float = 8.0,
    iters: int = 2,
    energy_kwargs: Optional[dict] = None,
) -> Optional[EpsSolveResult]:
    """
    Treat E(eps) approximately as:
        E(eps) ≈ E_beta + eps^2 * E_coeff
    where E_beta = E(0) includes the cost of the radial contraction ansatz (beta).
    Then refine with a couple of multiplicative updates.

    Returns None if no positive eps can satisfy the target.
    """
    if energy_kwargs is None:
        energy_kwargs = {}

    E_beta = total_energy_time_averaged(rho, alpha, mu, lam, eps=0.0, beta=beta, a=a, **energy_kwargs)
    if not np.isfinite(E_beta):
        return None
    if target_E <= E_beta:
        return None

    E1 = total_energy_time_averaged(rho, alpha, mu, lam, eps=1.0, beta=beta, a=a, **energy_kwargs)
    if not np.isfinite(E1):
        return None
    E_coeff = E1 - E_beta
    if E_coeff <= 0.0:
        return None

    eps = math.sqrt((target_E - E_beta) / E_coeff)
    if not np.isfinite(eps) or eps <= 0.0:
        return None
    eps = min(eps, eps_max)

    # refine a bit for non-quadratic regimes
    for _ in range(iters):
        E_eps = total_energy_time_averaged(rho, alpha, mu, lam, eps=eps, beta=beta, a=a, **energy_kwargs)
        if not np.isfinite(E_eps) or E_eps <= E_beta:
            return None
        eps = eps * math.sqrt((target_E - E_beta) / (E_eps - E_beta))
        if eps <= 0.0 or not np.isfinite(eps):
            return None
        eps = min(eps, eps_max)

    E_final = total_energy_time_averaged(rho, alpha, mu, lam, eps=eps, beta=beta, a=a, **energy_kwargs)
    if not np.isfinite(E_final):
        return None

    return EpsSolveResult(eps=eps, E_beta=E_beta, E_coeff=E_coeff, E_final=E_final)


# ----------------------------
# Objective wrapper
# ----------------------------
@dataclass
class EstimatorConfig:
    policy: str
    rho: float
    s_min: float
    slope_penalty: float
    slope_hard: bool
    min_inplane_speed: float
    max_internal_speed: float
    speed_penalty: float
    beta_penalty: float
    alpha_prior: float
    alpha_prior_w: float
    energy_Rmax_factor: float
    energy_nr: int
    energy_nt: int
    res_Rmax_factor: float
    res_nr: int
    res_nt: int
    eps_max: float

def objective_from_params(alpha: float, nu: float, beta: float, a: float, cfg: EstimatorConfig) -> Tuple[float, Dict[str, float]]:
    """
    Returns (objective_value, diagnostics_dict).
    """
    diag: Dict[str, float] = {
        "alpha": alpha, "nu": nu, "beta": beta, "a": a,
    }

    mu, lam = mu_lam_from_policy(cfg.policy, cfg.rho, alpha, nu)
    if not np.isfinite(mu) or not np.isfinite(lam) or mu <= 0.0 or lam <= 0.0:
        return 1e9, {**diag, "fail": 1.0}

    cT, cL, cw = speeds(mu, lam, cfg.rho, alpha, nu)
    diag.update({"rho": cfg.rho, "mu": mu, "lambda": lam, "cT": cT, "cL": cL, "cw": cw})

    # speed sanity (hard-ish penalties)
    pen = 0.0
    if min(cT, cL) < cfg.min_inplane_speed * c:
        pen += 50.0 * (cfg.min_inplane_speed * c / max(min(cT, cL), 1e-30) - 1.0) ** 2
    if max(cT, cL, cw) > cfg.max_internal_speed * c:
        pen += 10.0 * (max(cT, cL, cw) / (cfg.max_internal_speed * c) - 1.0) ** 2

    # optional soft single-cone penalty
    if cfg.speed_penalty > 0.0:
        pen += cfg.speed_penalty * (
            math.log(max(cT / c, 1e-30)) ** 2 +
            math.log(max(cL / c, 1e-30)) ** 2 +
            math.log(max(cw / c, 1e-30)) ** 2
        )

    # solve eps from energy constraint
    energy_kwargs = dict(Rmax_factor=cfg.energy_Rmax_factor, nr=cfg.energy_nr, nt=cfg.energy_nt)
    eps_sol = solve_eps_for_energy(
        rho=cfg.rho, alpha=alpha, mu=mu, lam=lam,
        beta=beta, a=a, target_E=E0,
        eps_max=cfg.eps_max,
        iters=2,
        energy_kwargs=energy_kwargs,
    )
    if eps_sol is None:
        # Return a graded penalty (not a flat 1e9) so optimizers can move.
        E_beta = total_energy_time_averaged(cfg.rho, alpha, mu, lam, eps=0.0, beta=beta, a=a, **energy_kwargs)
        diag.update({"E_beta": E_beta})
        pen_local = 0.0
        if np.isfinite(E_beta) and E_beta >= E0:
            pen_local += 200.0 * ((E_beta / E0) - 1.0) ** 2
            obj = 1e3 + pen + pen_local
            diag.update({"penalty": pen + pen_local, "objective": obj, "fail": 2.1})
            return obj, diag
        obj = 1e6 + pen
        diag.update({"penalty": pen, "objective": obj, "fail": 2.0})
        return obj, diag

    eps = eps_sol.eps
    diag.update({"eps": eps, "E": eps_sol.E_final, "E_beta": eps_sol.E_beta, "E_coeff": eps_sol.E_coeff})

    # enforce that geometric nonlinearity is "on"
    # For Gaussian: max |w_r| = eps * e^{-1/2}
    slope_max = eps * math.exp(-0.5)
    diag["slope_max"] = slope_max
    if slope_max < cfg.s_min:
        if cfg.slope_hard:
            obj = 1e9
            diag.update({"penalty": pen, "objective": obj, "fail": 3.0})
            return obj, diag
        if cfg.slope_penalty > 0.0 and cfg.s_min > 0.0:
            pen += cfg.slope_penalty * ((cfg.s_min - slope_max) / cfg.s_min) ** 2

    # residuals
    res_kwargs = dict(Rmax_factor=cfg.res_Rmax_factor, nr=cfg.res_nr, nt_avg=cfg.res_nt)
    rel_w, rel_R = residuals(cfg.rho, alpha, mu, lam, eps, beta, a, **res_kwargs)
    diag.update({"rel_w_rms": rel_w, "rel_R_rms": rel_R})

    # mild priors/regularizers
    if cfg.beta_penalty > 0.0:
        pen += cfg.beta_penalty * (beta ** 2)

    if cfg.alpha_prior_w > 0.0:
        # keep alpha near alpha_prior (typically close to 1.0), but do not force it.
        pen += cfg.alpha_prior_w * (math.log(alpha / cfg.alpha_prior) ** 2)

    obj = (rel_w + rel_R) + pen
    diag["objective"] = obj
    diag["penalty"] = pen
    return obj, diag


# ----------------------------
# Optimization driver
# ----------------------------
def run(cfg: EstimatorConfig, bounds: Dict[str, Tuple[float, float]], x0: Dict[str, float], maxiter: int, method: str) -> Dict[str, float]:
    """
    Optimize over (alpha, nu, beta, log_a) with bounds.
    """
    b_alpha = bounds["alpha"]
    b_nu = bounds["nu"]
    b_beta = bounds["beta"]
    b_loga = bounds["log_a"]

    x0_vec = np.array([x0["alpha"], x0["nu"], x0["beta"], math.log(x0["a"])], dtype=float)

    opt_bounds = [b_alpha, b_nu, b_beta, b_loga]

    def fun(xvec: np.ndarray) -> float:
        alpha, nu, beta, loga = float(xvec[0]), float(xvec[1]), float(xvec[2]), float(xvec[3])
        a = math.exp(loga)
        obj, _ = objective_from_params(alpha, nu, beta, a, cfg)
        return obj

    res = minimize(
        fun,
        x0_vec,
        method=method,
        bounds=opt_bounds,
        options={"maxiter": int(maxiter), "disp": True},
    )

    alpha, nu, beta, loga = float(res.x[0]), float(res.x[1]), float(res.x[2]), float(res.x[3])
    a = math.exp(loga)
    obj, diag = objective_from_params(alpha, nu, beta, a, cfg)
    diag["scipy_success"] = float(bool(res.success))
    diag["scipy_status"] = float(res.status)
    diag["scipy_nfev"] = float(getattr(res, "nfev", -1))
    diag["scipy_nit"] = float(getattr(res, "nit", -1))
    return diag


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--policy", choices=["S", "W", "L"], default="S", help="speed calibration policy")
    ap.add_argument("--rho", type=float, default=float(m_e / (lambda_C ** 3)),
                    help="fixed brane density rho in kg/m^3 (default: m_e/lambda_C^3)")
    ap.add_argument("--s-min", type=float, default=0.0,
                    help="target minimum slope_max (~eps*e^-1/2). Used with --slope-penalty or --slope-hard.")
    ap.add_argument("--slope-penalty", type=float, default=50.0,
                    help="soft penalty weight for slope_max falling below s-min (0 disables)")
    ap.add_argument("--slope-hard", action="store_true",
                    help="treat slope_max < s-min as hard infeasibility")
    ap.add_argument("--min-inplane-speed", type=float, default=0.1, help="minimum allowed min(cT,cL) as fraction of c")
    ap.add_argument("--max-internal-speed", type=float, default=30.0, help="maximum allowed max(cT,cL,cw) as multiple of c")
    ap.add_argument("--speed-penalty", type=float, default=0.0, help="soft penalty weight to keep all speeds near c (0 disables)")
    ap.add_argument("--beta-penalty", type=float, default=0.02, help="soft penalty weight on beta^2")
    ap.add_argument("--alpha-prior", type=float, default=0.98, help="alpha prior (for soft regularization)")
    ap.add_argument("--alpha-prior-w", type=float, default=0.0, help="alpha prior weight (0 disables)")
    ap.add_argument("--eps-max", type=float, default=8.0, help="cap for eps during energy solve")

    # numerical settings
    ap.add_argument("--energy-Rmax-factor", type=float, default=10.0)
    ap.add_argument("--energy-nr", type=int, default=2200)
    ap.add_argument("--energy-nt", type=int, default=40)
    ap.add_argument("--res-Rmax-factor", type=float, default=10.0)
    ap.add_argument("--res-nr", type=int, default=1800)
    ap.add_argument("--res-nt", type=int, default=24)

    # bounds
    ap.add_argument("--alpha-min", type=float, default=0.90)
    ap.add_argument("--alpha-max", type=float, default=0.999)
    ap.add_argument("--nu-min", type=float, default=0.05)
    ap.add_argument("--nu-max", type=float, default=0.35)
    ap.add_argument("--beta-min", type=float, default=0.0)
    ap.add_argument("--beta-max", type=float, default=1e-2)
    ap.add_argument("--a-min", type=float, default=0.5, help="min a in units of lambda_C")
    ap.add_argument("--a-max", type=float, default=20.0, help="max a in units of lambda_C")

    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--method", type=str, default="Powell", help="scipy minimize method supporting bounds (Powell or L-BFGS-B)")
    ap.add_argument("--seed", type=int, default=0)

    args = ap.parse_args()

    cfg = EstimatorConfig(
        policy=args.policy,
        rho=float(args.rho),
        s_min=float(args.s_min),
        slope_penalty=float(args.slope_penalty),
        slope_hard=bool(args.slope_hard),
        min_inplane_speed=float(args.min_inplane_speed),
        max_internal_speed=float(args.max_internal_speed),
        speed_penalty=float(args.speed_penalty),
        beta_penalty=float(args.beta_penalty),
        alpha_prior=float(args.alpha_prior),
        alpha_prior_w=float(args.alpha_prior_w),
        energy_Rmax_factor=float(args.energy_Rmax_factor),
        energy_nr=int(args.energy_nr),
        energy_nt=int(args.energy_nt),
        res_Rmax_factor=float(args.res_Rmax_factor),
        res_nr=int(args.res_nr),
        res_nt=int(args.res_nt),
        eps_max=float(args.eps_max),
    )

    bounds = {
        "alpha": (float(args.alpha_min), float(args.alpha_max)),
        "nu": (float(args.nu_min), float(args.nu_max)),
        "beta": (float(args.beta_min), float(args.beta_max)),
        "log_a": (math.log(float(args.a_min) * lambda_C), math.log(float(args.a_max) * lambda_C)),
    }

    # starting point (middle of bounds)
    rng = np.random.default_rng(args.seed)
    x0 = {
        "alpha": 0.5 * (bounds["alpha"][0] + bounds["alpha"][1]),
        "nu": 0.5 * (bounds["nu"][0] + bounds["nu"][1]),
        "beta": max(bounds["beta"][0], 1e-10),
        "a": math.exp(0.5 * (bounds["log_a"][0] + bounds["log_a"][1])),
    }

    # small random jitter to avoid symmetry traps (optional)
    x0["alpha"] = float(np.clip(x0["alpha"] * (1.0 + 0.02 * rng.standard_normal()), *bounds["alpha"]))
    x0["nu"] = float(np.clip(x0["nu"] * (1.0 + 0.05 * rng.standard_normal()), *bounds["nu"]))
    x0["beta"] = float(np.clip(x0["beta"], *bounds["beta"]))
    x0["a"] = float(np.clip(x0["a"] * (1.0 + 0.10 * rng.standard_normal()),
                            math.exp(bounds["log_a"][0]), math.exp(bounds["log_a"][1])))

    diag = run(cfg, bounds=bounds, x0=x0, maxiter=args.maxiter, method=args.method)

    # pretty print
    print("\nOPTIMUM (v2):")
    print(f"  policy = {cfg.policy}")
    print(f"  alpha  = {diag.get('alpha'):.6e}")
    print(f"  nu     = {diag.get('nu'):.6e}")
    print(f"  beta   = {diag.get('beta'):.6e}")
    print(f"  a      = {diag.get('a'):.6e}   (a/lambda_C = {diag.get('a')/lambda_C:.6e})")
    print(f"  eps    = {diag.get('eps', float('nan')):.6e}   (A/a)")
    if 'fail' in diag:
        print(f"  fail   = {diag.get('fail')}")
    print(f"  slope_max = {diag.get('slope_max', float('nan')):.6e}")
    print(f"  rho    = {diag.get('rho'):.6e} kg/m^3")
    print(f"  mu     = {diag.get('mu'):.6e} Pa")
    print(f"  lambda = {diag.get('lambda'):.6e} Pa")
    print(f"  speeds: cT={diag.get('cT'):.6e}, cL={diag.get('cL'):.6e}, cw={diag.get('cw'):.6e} (m/s)")
    print(f"  E      = {diag.get('E', float('nan')):.6e} J   (target E0={E0:.6e} J)")
    print(f"  rel_w_rms = {diag.get('rel_w_rms', float('nan')):.6e}")
    print(f"  rel_R_rms = {diag.get('rel_R_rms', float('nan')):.6e}")
    print(f"  penalty   = {diag.get('penalty', float('nan')):.6e}")
    print(f"  objective = {diag.get('objective', float('nan')):.6e}")
    print(f"  scipy: success={bool(diag.get('scipy_success'))}, nit={int(diag.get('scipy_nit'))}, nfev={int(diag.get('scipy_nfev'))}")

if __name__ == "__main__":
    main()
