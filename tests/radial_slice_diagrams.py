#!/usr/bin/env python3
"""
radial_slice_diagrams.py

Standalone “radial slice” diagram generator for the electron confinement
discussion (static, 1D cross-section through the tube).

This script does NOT read a simulation snapshot. Instead, it generates
*design-target* profiles that satisfy the two key confinement constraints we
derived:

(1) Bending / guiding requirement (order-of-magnitude):
        Δn / n  ~  a / R

(2) Inner/outer rail coherence (path-length compensation at fixed ω):
        k_in (R - w0) ≈ k_out (R + w0)
    ⇒  k_in / k_out ≈ (R + w0) / (R - w0)
    and since k ∝ n_eff at fixed ω, we enforce:
        n_eff(-w0) / n_eff(+w0) = (R + w0) / (R - w0)

We then derive:
- effective speed: c_eff(r) = c_ref / n_eff(r)
- axial “tension/stress proxy”: σ_s(r) ∝ ρ c_eff(r)^2  (tension-dominated intuition)
- a simple radial “tension” from radial stretch in 4D: T_r(r) ∝ (|∂_r X| - 1)

Plots are styled similarly to your 1D brane experiment plots: clean line plots,
grid, bold titles, annotated parameter boxes, and marked reference positions.

Outputs:
- radial_slice_geometry.png
- radial_slice_propagation.png
- radial_slice_stresses.png

Dependencies:
- numpy
- matplotlib
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Utilities
# -----------------------------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def safe_div(a: float, b: float, eps: float = 1e-12) -> float:
    return a / (b if abs(b) > eps else (eps if b >= 0 else -eps))

def solve_delta_n_for_gaussian(target_ratio: float) -> float:
    """
    Choose Δn_core for n_sym(r)=1+Δn_core*exp(-(r/a)^2) such that
    (n_sym(0)-n_sym(a))/n_sym(0) ≈ target_ratio = a/R.

    n_sym(0)=1+Δ
    n_sym(a)=1+Δ*exp(-1)

    Δn/n ≈ (Δ*(1-exp(-1))) / (1+Δ)
    => Δ = t / ((1-exp(-1)) - t)  , requires t < (1-exp(-1)) ~ 0.632.
    """
    one_minus_e = 1.0 - math.exp(-1.0)  # ~0.6321
    t = target_ratio
    if t <= 0:
        return 0.0
    if t >= 0.95 * one_minus_e:
        # clamp to keep finite
        t = 0.95 * one_minus_e
    return t / (one_minus_e - t)

def format_units_scale(units: str) -> tuple[str, float]:
    """
    Return (axis_label, scale) mapping meters -> chosen unit.
    """
    u = units.lower().strip()
    if u in ("m", "meter", "meters"):
        return "m", 1.0
    if u in ("nm", "nanometer", "nanometers"):
        return "nm", 1e9
    if u in ("pm", "picometer", "picometers"):
        return "pm", 1e12
    if u in ("au", "arb", "arbitrary", "arb."):
        return "arb. units", 1.0
    raise ValueError(f"Unknown units '{units}'. Use m, nm, pm, or au.")


# -----------------------------
# Model (design-target profiles)
# -----------------------------

@dataclass
class SliceParams:
    # Geometry
    R: float            # torus major radius
    a: float            # tube (mode) radius / width scale in the slice
    w0: float           # rail offset (inner/outer sampling position)

    # Wave / medium references (dimensionless-friendly)
    c_ref: float = 1.0
    rho: float = 1.0
    omega: float = 1.0

    # Displacement amplitudes (purely for plotting geometry/proxies)
    x4_amp: float = 1.0e-12    # meters (or au if you use au)
    ur_amp: float = 2.0e-13    # meters (or au)

    # Radial “constitutive” proxy strength
    Kr: float = 1.0            # radial stiffness proxy for T_r


def build_profiles(p: SliceParams, r: np.ndarray) -> dict[str, np.ndarray | float]:
    """
    Build target profiles that satisfy the confinement constraints.

    Returns a dict with arrays for:
      - n_eff, c_eff, k_eff
      - x4 (X^4 displacement), u_r (radial), eps_s (tangential strain proxy)
      - sigma_s (axial stress proxy)
      - T_r (radial tension proxy), and component projections (sigma_r_comp, sigma_4_comp)
    plus some scalars for achieved constraint checks.
    """
    R, a, w0 = p.R, p.a, p.w0

    # --- Constraint (1): Δn/n ~ a/R via symmetric Gaussian n_sym
    target_dn_over_n = safe_div(a, R)
    delta_n_core = solve_delta_n_for_gaussian(target_dn_over_n)
    n_sym = 1.0 + delta_n_core * np.exp(-(r / a) ** 2)

    # --- Constraint (2): enforce inner/outer ratio at ±w0 via linear asymmetry term
    # target ratio for k (and hence n_eff at fixed ω)
    target_ratio = (R + w0) / (R - w0)  # > 1
    # Use n_eff(r) = n_sym(r) - gamma * r  (so n_eff(-w0) > n_eff(+w0))
    n_at_w0 = float(1.0 + delta_n_core * math.exp(-(w0 / a) ** 2))
    # Solve exactly:
    # (n_at_w0 + gamma*w0) / (n_at_w0 - gamma*w0) = target_ratio
    gamma = (n_at_w0 * (target_ratio - 1.0)) / (w0 * (target_ratio + 1.0))

    n_eff = n_sym - gamma * r

    # Keep n_eff positive and not crazy small
    n_eff = np.clip(n_eff, 1e-6, None)

    # Propagation proxies
    c_eff = p.c_ref / n_eff
    k_eff = (n_eff * p.omega) / p.c_ref

    # --- Geometry proxies:
    # X^4 bulge: symmetric (even), tied to n_sym so it is “even around the ring”
    # (This is a *proxy*: in the real model n_eff comes from linearization of W(g).)
    x4 = p.x4_amp * np.sqrt(np.maximum(n_sym - 1.0, 0.0) / max(delta_n_core, 1e-12))

    # radial displacement (even, small contraction in the tube core)
    u_r = -p.ur_amp * np.exp(-(r / a) ** 2)

    # tangential strain proxy from torus curvature (geometry):
    # For a straight reference fiber, bending into a circle gives ε_s ≈ -r/R.
    eps_s = -r / R
    stretch_s = 1.0 + eps_s  # = 1 - r/R

    # --- Stresses / tensions:
    # Axial stress proxy from tension-dominated wave speed (σ_s ∝ ρ c_eff^2)
    sigma_s = p.rho * (c_eff ** 2)

    # Radial “stretch” in 4D from ∂_r X components:
    # ∂_r X ≈ (1 + u_r') e_r + (x4') e_4   (keep simple)
    du_r = np.gradient(u_r, r)
    dx4 = np.gradient(x4, r)
    drX_norm = np.sqrt((1.0 + du_r) ** 2 + (dx4) ** 2)
    # Radial tension proxy from a simple quadratic stretch energy
    T_r = p.Kr * (drX_norm - 1.0)

    # Component “stress contributions” along e_r and e_4
    # (Interpretation: projection of radial tension along those derivative components)
    sigma_r_comp = T_r * (1.0 + du_r) / np.maximum(drX_norm, 1e-12)
    sigma_4_comp = T_r * (dx4) / np.maximum(drX_norm, 1e-12)

    # --- Achieved constraint checks
    # (1) Δn/n across [0,a] using center and r=a
    n0 = float(1.0 + delta_n_core)
    na = float(1.0 + delta_n_core * math.exp(-1.0))
    achieved_dn_over_n = (n0 - na) / n0

    # (2) inner/outer ratio at ±w0
    # Evaluate by interpolation
    n_in = float(np.interp(-w0, r, n_eff))
    n_out = float(np.interp(+w0, r, n_eff))
    achieved_ratio = n_in / n_out

    return {
        "n_eff": n_eff,
        "n_sym": n_sym,
        "c_eff": c_eff,
        "k_eff": k_eff,
        "x4": x4,
        "u_r": u_r,
        "eps_s": eps_s,
        "stretch_s": stretch_s,
        "sigma_s": sigma_s,
        "T_r": T_r,
        "sigma_r_comp": sigma_r_comp,
        "sigma_4_comp": sigma_4_comp,
        "gamma": gamma,
        "delta_n_core": delta_n_core,
        "target_dn_over_n": target_dn_over_n,
        "achieved_dn_over_n": achieved_dn_over_n,
        "target_ratio": target_ratio,
        "achieved_ratio": achieved_ratio,
    }


# -----------------------------
# Plotting (1D brane-style)
# -----------------------------

def add_param_box(ax, text: str) -> None:
    ax.text(
        0.02, 0.95, text,
        transform=ax.transAxes,
        fontsize=11, va="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.85),
    )

def add_markers(ax, xvals, labels=None):
    for i, xv in enumerate(xvals):
        ax.axvline(x=xv, linestyle="--", linewidth=1.0, alpha=0.5)
        if labels and i < len(labels):
            ax.text(
                xv, 0.98, labels[i],
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.7),
            )

def save_fig(fig, outpath: Path) -> None:
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)

def plot_geometry(r_plot, unit_label, outdir: Path, profiles: dict, p: SliceParams) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Geometry Proxies", fontsize=16, fontweight="bold")

    # X^4 bulge
    axes[0].plot(r_plot, profiles["x4"], linewidth=2)
    axes[0].set_ylabel(f"X⁴ displacement ξ₄ [{unit_label}]" if unit_label != "arb. units" else "X⁴ displacement ξ₄ [au]")
    axes[0].grid(True, alpha=0.3)

    # Radial displacement
    axes[1].plot(r_plot, profiles["u_r"], linewidth=2)
    axes[1].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[1].set_ylabel(f"Radial displacement u_r [{unit_label}]" if unit_label != "arb. units" else "Radial displacement u_r [au]")
    axes[1].grid(True, alpha=0.3)

    # Tangential strain / stretch factor
    axes[2].plot(r_plot, profiles["eps_s"], linewidth=2, label="ε_s ≈ -r/R")
    axes[2].plot(r_plot, profiles["stretch_s"], linewidth=2, alpha=0.85, label="stretch_s = 1 + ε_s")
    axes[2].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[2].set_ylabel("Tangential strain / stretch [–]")
    axes[2].set_xlabel(f"Signed radial coordinate r [{unit_label}] (inner < 0 < outer)")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right", fontsize=10)

    # markers at ±w0
    w0_plot = (p.w0 * (r_plot[-1] / (r_plot[-1] if r_plot[-1] != 0 else 1)))  # no-op, r_plot already scaled
    add_markers(axes[0], [-w0_plot, +w0_plot], labels=["inner (-w0)", "outer (+w0)"])
    add_markers(axes[1], [-w0_plot, +w0_plot])
    add_markers(axes[2], [-w0_plot, +w0_plot])

    param_text = (
        f"R={p.R:.3e}  a={p.a:.3e}  w0={p.w0:.3e}\n"
        f"Target Δn/n≈a/R={profiles['target_dn_over_n']:.3f}, achieved={profiles['achieved_dn_over_n']:.3f}\n"
        f"Target n(-w0)/n(+w0)={(profiles['target_ratio']):.3f}, achieved={profiles['achieved_ratio']:.3f}"
    )
    add_param_box(axes[0], param_text)

    save_fig(fig, outdir / "radial_slice_geometry.png")


def plot_propagation(r_plot, unit_label, outdir: Path, profiles: dict, p: SliceParams) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Effective Propagation Proxies", fontsize=16, fontweight="bold")

    axes[0].plot(r_plot, profiles["n_eff"], linewidth=2, label="n_eff(r)")
    axes[0].plot(r_plot, profiles["n_sym"], linewidth=2, alpha=0.6, label="n_sym(r) (even part)")
    axes[0].set_ylabel("Effective index n_eff [–]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=10)

    axes[1].plot(r_plot, profiles["c_eff"], linewidth=2)
    axes[1].set_ylabel("Effective speed c_eff [–]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(r_plot, profiles["k_eff"], linewidth=2)
    axes[2].set_ylabel("Effective wavenumber k_eff [–]")
    axes[2].set_xlabel(f"Signed radial coordinate r [{unit_label}] (inner < 0 < outer)")
    axes[2].grid(True, alpha=0.3)

    w0_plot = (p.w0 * (r_plot[-1] / (r_plot[-1] if r_plot[-1] != 0 else 1)))
    add_markers(axes[0], [-w0_plot, +w0_plot], labels=["inner (-w0)", "outer (+w0)"])
    add_markers(axes[1], [-w0_plot, +w0_plot])
    add_markers(axes[2], [-w0_plot, +w0_plot])

    param_text = (
        f"c_ref={p.c_ref:.3g}  ρ={p.rho:.3g}  ω={p.omega:.3g}\n"
        f"Δn_core={profiles['delta_n_core']:.3f}  asymmetry slope γ={profiles['gamma']:.3e}\n"
        f"Target Δn/n≈a/R={profiles['target_dn_over_n']:.3f}, achieved={profiles['achieved_dn_over_n']:.3f}\n"
        f"Target ratio={(profiles['target_ratio']):.3f}, achieved={profiles['achieved_ratio']:.3f}"
    )
    add_param_box(axes[0], param_text)

    save_fig(fig, outdir / "radial_slice_propagation.png")


def plot_stresses(r_plot, unit_label, outdir: Path, profiles: dict, p: SliceParams) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Stress/Tension Proxies", fontsize=16, fontweight="bold")

    axes[0].plot(r_plot, profiles["sigma_s"], linewidth=2)
    axes[0].set_ylabel("Axial stress proxy σ_s ∝ ρ c_eff² [–]")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(r_plot, profiles["T_r"], linewidth=2)
    axes[1].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[1].set_ylabel("Radial tension proxy T_r [–]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(r_plot, profiles["sigma_r_comp"], linewidth=2, label="radial component")
    axes[2].plot(r_plot, profiles["sigma_4_comp"], linewidth=2, label="X⁴ component")
    axes[2].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[2].set_ylabel("Projected radial-tension components [–]")
    axes[2].set_xlabel(f"Signed radial coordinate r [{unit_label}] (inner < 0 < outer)")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right", fontsize=10)

    w0_plot = (p.w0 * (r_plot[-1] / (r_plot[-1] if r_plot[-1] != 0 else 1)))
    add_markers(axes[0], [-w0_plot, +w0_plot], labels=["inner (-w0)", "outer (+w0)"])
    add_markers(axes[1], [-w0_plot, +w0_plot])
    add_markers(axes[2], [-w0_plot, +w0_plot])

    param_text = (
        f"Kr={p.Kr:.3g}\n"
        f"σ_s uses tension-dominated proxy: σ_s=ρ(c_ref/n_eff)²\n"
        f"Radial proxy uses |∂_r X| with X⁴ + u_r only (static slice)."
    )
    add_param_box(axes[0], param_text)

    save_fig(fig, outdir / "radial_slice_stresses.png")


# -----------------------------
# Main
# -----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Generate radial-slice confinement diagrams (design targets).")
    ap.add_argument("--outdir", type=str, default="radial_slice_plots", help="Output directory for PNGs.")
    ap.add_argument("--units", type=str, default="pm", help="Axis units for r, x4, u_r: pm|nm|m|au")

    # If you want a Compton-ish default, leave these as None and use --use-compton
    ap.add_argument("--R", type=float, default=None, help="Torus major radius R (meters if units != au).")
    ap.add_argument("--a", type=float, default=None, help="Tube width scale a (meters if units != au).")
    ap.add_argument("--w0", type=float, default=None, help="Rail offset w0 (meters if units != au).")

    ap.add_argument("--use-compton", action="store_true",
                    help="Use R = λC/(4π), a = 0.3 R, w0 = 0.3 a with λC = 3.86159e-13 m.")
    ap.add_argument("--r-mult", type=float, default=3.0, help="Plot range is r ∈ [-r_mult*a, +r_mult*a].")

    ap.add_argument("--c-ref", type=float, default=1.0, help="Reference speed (dimensionless ok).")
    ap.add_argument("--rho", type=float, default=1.0, help="Density proxy (dimensionless ok).")
    ap.add_argument("--omega", type=float, default=1.0, help="Angular frequency (dimensionless ok).")

    ap.add_argument("--x4-amp", type=float, default=None, help="X^4 bulge amplitude (same length units as R).")
    ap.add_argument("--ur-amp", type=float, default=None, help="Radial displacement amplitude (same length units as R).")
    ap.add_argument("--Kr", type=float, default=1.0, help="Radial stiffness proxy for radial tension.")

    args = ap.parse_args()

    unit_label, scale = format_units_scale(args.units)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    if args.use_compton:
        lambda_C = 3.8615926796e-13  # m
        R = lambda_C / (4.0 * math.pi)
        a = 0.30 * R
        w0 = 0.30 * a
    else:
        if args.R is None or args.a is None or args.w0 is None:
            ap.error("Provide --R --a --w0 (or use --use-compton).")
        R, a, w0 = args.R, args.a, args.w0

    # Default amplitudes if not provided
    x4_amp = args.x4_amp if args.x4_amp is not None else 0.8 * a
    ur_amp = args.ur_amp if args.ur_amp is not None else 0.2 * a

    p = SliceParams(
        R=R, a=a, w0=w0,
        c_ref=args.c_ref, rho=args.rho, omega=args.omega,
        x4_amp=x4_amp, ur_amp=ur_amp, Kr=args.Kr
    )

    # Signed radial coordinate array (meters unless units=au)
    r = np.linspace(-args.r_mult * a, +args.r_mult * a, 2000)

    profiles = build_profiles(p, r)

    # r_plot is in chosen units
    r_plot = r * scale
    x4_plot = np.asarray(profiles["x4"]) * scale if unit_label != "arb. units" else np.asarray(profiles["x4"])
    ur_plot = np.asarray(profiles["u_r"]) * scale if unit_label != "arb. units" else np.asarray(profiles["u_r"])

    # Replace arrays in dict for plotting (so y-axes show chosen units)
    profiles_plot = dict(profiles)
    profiles_plot["x4"] = x4_plot
    profiles_plot["u_r"] = ur_plot

    # Print a quick constraint report
    print("\nRadial slice constraint report")
    print("------------------------------")
    print(f"R   = {R:.6e} m")
    print(f"a   = {a:.6e} m")
    print(f"w0  = {w0:.6e} m")
    print()
    print(f"Target Δn/n ≈ a/R = {profiles['target_dn_over_n']:.6f}")
    print(f"Achieved Δn/n      = {profiles['achieved_dn_over_n']:.6f}")
    print()
    print(f"Target ratio n(-w0)/n(+w0) = (R+w0)/(R-w0) = {profiles['target_ratio']:.6f}")
    print(f"Achieved ratio            = {profiles['achieved_ratio']:.6f}")
    print(f"Δn_core used              = {profiles['delta_n_core']:.6f}")
    print(f"Asymmetry slope gamma      = {profiles['gamma']:.6e}")
    print()

    # Plots
    plot_geometry(r_plot, unit_label, outdir, profiles_plot, p)
    plot_propagation(r_plot, unit_label, outdir, profiles, p)  # propagation is dimensionless
    plot_stresses(r_plot, unit_label, outdir, profiles, p)      # stresses are dimensionless

    print(f"Saved PNGs to: {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
