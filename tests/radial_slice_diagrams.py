#!/usr/bin/env python3
"""
radial_slice_diagrams_v2.py

Standalone “radial slice” diagram generator for the electron confinement discussion
(static 1D cross-section through the tube).

Fixes vs v1:
- Coordinate convention: default is an outward slice starting at the centerline (r>=0).
  Optional signed slice (inner<0<outer) via --domain signed.
- Default x-axis is global radius ρ = R + r (so plots are not centered at 0).
- Correct scaling of marker positions (±w0) to plotted units.
- Tangential strain sign: outer side has longer circumference => stretch_s = 1 + r/R.
- Remove the unphysical linear asymmetry term in n_eff (which made n_eff cross 0 and
  blew up c_eff and σ_s). Replace by a localized odd function
    f_odd(r) = (r/a) * exp(-(r/a)^2).
- Improved label placement to reduce overlap.

Important: these are *design-target proxies*. In the full theory, n_eff(r) comes from
linearizing the dynamics around a deformed background defined by X(r) and W(g).

Outputs (PNG):
- radial_slice_geometry.png
- radial_slice_propagation.png
- radial_slice_stresses.png
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def safe_div(a: float, b: float, eps: float = 1e-12) -> float:
    return a / (b if abs(b) > eps else (eps if b >= 0 else -eps))


def solve_delta_n_for_gaussian(target_ratio: float) -> float:
    """
    Choose Δn_core for n_sym(r)=1+Δn_core*exp(-(r/a)^2) such that
    (n_sym(0)-n_sym(a))/n_sym(0) ≈ target_ratio = a/R.
    """
    one_minus_e = 1.0 - math.exp(-1.0)  # ~0.6321
    t = float(target_ratio)
    if t <= 0:
        return 0.0
    if t >= 0.95 * one_minus_e:
        t = 0.95 * one_minus_e
    return t / (one_minus_e - t)


def format_units_scale(units: str) -> tuple[str, float]:
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


@dataclass
class SliceParams:
    R: float          # torus major radius (m or au)
    a: float          # tube (mode) radius scale (m or au)
    w0: float         # rail offset (m or au)

    # reference constants (dimensionless-friendly)
    c_ref: float = 1.0
    rho: float = 1.0
    omega: float = 1.0

    # geometry proxy amplitudes (m or au)
    x4_amp: float = 1.0
    ur_amp: float = 0.2

    # radial tension proxy stiffness
    Kr: float = 1.0


def f_odd_localized(r: np.ndarray, a: float) -> np.ndarray:
    """Odd, localized shape: r/a * exp(-(r/a)^2)."""
    x = r / a
    return x * np.exp(-(x ** 2))


def build_profiles(p: SliceParams, r: np.ndarray) -> dict[str, np.ndarray | float]:
    """
    Build design-target profiles satisfying:
      (1) Δn/n ~ a/R  (using n_sym Gaussian)
      (2) n(-w0)/n(+w0) = (R+w0)/(R-w0)  (using localized odd modulation)

    Note: constraint (2) is evaluated on the symmetric extension even if r>=0 is plotted.
    """
    R, a, w0 = p.R, p.a, p.w0

    # (1) Δn/n ~ a/R
    target_dn_over_n = safe_div(a, R)
    delta_n_core = solve_delta_n_for_gaussian(target_dn_over_n)
    n_sym = 1.0 + delta_n_core * np.exp(-(r / a) ** 2)

    # (2) coherence ratio across inner/outer rails at ±w0
    target_ratio = (R + w0) / (R - w0)  # k_in/k_out, hence n_in/n_out at fixed ω
    f_w = float(f_odd_localized(np.array([w0], dtype=float), a=a)[0])
    f_w = max(abs(f_w), 1e-12)  # avoid divide-by-zero

    # Choose alpha so that: (1 + alpha f_w)/(1 - alpha f_w) = target_ratio
    # using n_eff(r) = n_sym(r) * (1 - alpha f_odd(r))
    alpha = (target_ratio - 1.0) / (f_w * (target_ratio + 1.0))

    f = f_odd_localized(r, a=a)
    n_eff = n_sym * (1.0 - alpha * f)

    # Keep strictly positive
    n_eff = np.clip(n_eff, 1e-6, None)

    # propagation proxies
    c_eff = p.c_ref / n_eff
    k_eff = (n_eff * p.omega) / p.c_ref

    # geometry proxies
    if delta_n_core > 1e-12:
        x4 = p.x4_amp * np.sqrt(np.maximum(n_sym - 1.0, 0.0) / delta_n_core)
    else:
        x4 = np.zeros_like(r)

    u_r = -p.ur_amp * np.exp(-(r / a) ** 2)

    # Correct torus stretch: ds = (R+r)dθ  => stretch = (R+r)/R = 1 + r/R
    eps_s = r / R
    stretch_s = 1.0 + eps_s

    # stresses/tensions (proxies)
    sigma_s = p.rho * (c_eff ** 2)

    # radial stretch proxy in 4D from |∂_r X|
    du_r = np.gradient(u_r, r)
    dx4 = np.gradient(x4, r)
    drX_norm = np.sqrt((1.0 + du_r) ** 2 + (dx4) ** 2)
    T_r = p.Kr * (drX_norm - 1.0)

    sigma_r_comp = T_r * (1.0 + du_r) / np.maximum(drX_norm, 1e-12)
    sigma_4_comp = T_r * (dx4) / np.maximum(drX_norm, 1e-12)

    # achieved checks
    n0 = float(1.0 + delta_n_core)
    na = float(1.0 + delta_n_core * math.exp(-1.0))
    achieved_dn_over_n = (n0 - na) / n0

    # ratio at ±w0 analytically
    n_sym_w0 = float(1.0 + delta_n_core * math.exp(-(w0 / a) ** 2))
    n_in = n_sym_w0 * (1.0 + alpha * f_w)   # r=-w0 => increases
    n_out = n_sym_w0 * (1.0 - alpha * f_w)  # r=+w0 => decreases
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
        "alpha": alpha,
        "delta_n_core": delta_n_core,
        "target_dn_over_n": target_dn_over_n,
        "achieved_dn_over_n": achieved_dn_over_n,
        "target_ratio": target_ratio,
        "achieved_ratio": achieved_ratio,
    }


def add_param_box(ax, text: str, loc: str = "tl") -> None:
    xy = (0.02, 0.95) if loc == "tl" else (0.98, 0.95)
    ha = "left" if loc == "tl" else "right"
    ax.text(
        xy[0], xy[1], text,
        transform=ax.transAxes,
        fontsize=10, va="top", ha=ha,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.85),
    )


def add_vertical_markers(ax, xs, labels=None, y=0.02):
    for i, xv in enumerate(xs):
        ax.axvline(x=xv, linestyle="--", linewidth=1.0, alpha=0.5)
        if labels and i < len(labels) and labels[i]:
            ax.text(
                xv, y, labels[i],
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.7),
            )


def save_fig(fig, outpath):
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_geometry(x_plot, x_label, unit_label, outdir, prof, p, marker_xs, marker_labels):
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Geometry Proxies", fontsize=16, fontweight="bold")

    axes[0].plot(x_plot, prof["x4"], linewidth=2)
    axes[0].set_ylabel(f"X⁴ displacement ξ₄ [{unit_label}]" if unit_label != "arb. units" else "X⁴ displacement ξ₄ [au]")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(x_plot, prof["u_r"], linewidth=2)
    axes[1].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[1].set_ylabel(f"Radial displacement u_r [{unit_label}]" if unit_label != "arb. units" else "Radial displacement u_r [au]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(x_plot, prof["eps_s"], linewidth=2, label="ε_s ≈ r/R")
    axes[2].plot(x_plot, prof["stretch_s"], linewidth=2, alpha=0.85, label="stretch_s = 1 + ε_s")
    axes[2].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[2].set_ylabel("Tangential strain / stretch [–]")
    axes[2].set_xlabel(x_label)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right", fontsize=10)

    for ax in axes:
        add_vertical_markers(ax, marker_xs, marker_labels)

    param_text = (
        f"R={p.R:.3e}  a={p.a:.3e}  w0={p.w0:.3e}\n"
        f"Target Δn/n≈a/R={prof['target_dn_over_n']:.3f}, achieved={prof['achieved_dn_over_n']:.3f}\n"
        f"Target n(-w0)/n(+w0)={prof['target_ratio']:.3f}, achieved={prof['achieved_ratio']:.3f}"
    )
    add_param_box(axes[0], param_text, loc="tl")

    save_fig(fig, Path(outdir) / "radial_slice_geometry.png")


def plot_propagation(x_plot, x_label, outdir, prof, p, marker_xs, marker_labels):
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Effective Propagation Proxies", fontsize=16, fontweight="bold")

    axes[0].plot(x_plot, prof["n_eff"], linewidth=2, label="n_eff(r)")
    axes[0].plot(x_plot, prof["n_sym"], linewidth=2, alpha=0.6, label="n_sym(r) (even part)")
    axes[0].set_ylabel("Effective index n_eff [–]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=10)

    axes[1].plot(x_plot, prof["c_eff"], linewidth=2)
    axes[1].set_ylabel("Effective speed c_eff [–]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(x_plot, prof["k_eff"], linewidth=2)
    axes[2].set_ylabel("Effective wavenumber k_eff [–]")
    axes[2].set_xlabel(x_label)
    axes[2].grid(True, alpha=0.3)

    for ax in axes:
        add_vertical_markers(ax, marker_xs, marker_labels)

    param_text = (
        f"c_ref={p.c_ref:.3g}  ρ={p.rho:.3g}  ω={p.omega:.3g}\n"
        f"Δn_core={prof['delta_n_core']:.3f}  odd-mod amplitude α={prof['alpha']:.3g}\n"
        f"Target Δn/n≈a/R={prof['target_dn_over_n']:.3f}, achieved={prof['achieved_dn_over_n']:.3f}\n"
        f"Target ratio={prof['target_ratio']:.3f}, achieved={prof['achieved_ratio']:.3f}"
    )
    add_param_box(axes[0], param_text, loc="tl")

    save_fig(fig, Path(outdir) / "radial_slice_propagation.png")


def plot_stresses(x_plot, x_label, outdir, prof, p, marker_xs, marker_labels):
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    fig.suptitle("Radial Slice – Stress/Tension Proxies", fontsize=16, fontweight="bold")

    axes[0].plot(x_plot, prof["sigma_s"], linewidth=2)
    axes[0].set_ylabel("Axial stress proxy σ_s ∝ ρ c_eff² [–]")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(x_plot, prof["T_r"], linewidth=2)
    axes[1].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[1].set_ylabel("Radial tension proxy T_r [–]")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(x_plot, prof["sigma_r_comp"], linewidth=2, label="radial component")
    axes[2].plot(x_plot, prof["sigma_4_comp"], linewidth=2, label="X⁴ component")
    axes[2].axhline(0, linestyle="--", linewidth=0.8, alpha=0.4)
    axes[2].set_ylabel("Projected radial-tension components [–]")
    axes[2].set_xlabel(x_label)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right", fontsize=10)

    for ax in axes:
        add_vertical_markers(ax, marker_xs, marker_labels)

    param_text = (
        f"Kr={p.Kr:.3g}\n"
        f"σ_s uses proxy: σ_s=ρ(c_ref/n_eff)²\n"
        f"Radial proxy uses |∂_r X| with X⁴ + u_r only (static slice)."
    )
    add_param_box(axes[0], param_text, loc="tl")

    save_fig(fig, Path(outdir) / "radial_slice_stresses.png")


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate radial-slice confinement diagrams (design targets).")
    ap.add_argument("--outdir", type=str, default="radial_slice_plots_v2", help="Output directory for PNGs.")
    ap.add_argument("--units", type=str, default="pm", help="Axis units: pm|nm|m|au")
    ap.add_argument("--domain", choices=["outward", "signed"], default="outward",
                    help="outward: r∈[0,rmax]; signed: r∈[-rmax,+rmax]")
    ap.add_argument("--x-axis", choices=["local", "global"], default="global",
                    help="local: x=r; global: x=R+r (distance from torus center).")

    ap.add_argument("--R", type=float, default=None, help="Torus major radius R (meters if units != au).")
    ap.add_argument("--a", type=float, default=None, help="Tube width scale a (meters if units != au).")
    ap.add_argument("--w0", type=float, default=None, help="Rail offset w0 (meters if units != au).")

    ap.add_argument("--use-compton", action="store_true",
                    help="Use R = λC/(4π), a = 0.30 R, w0 = 0.30 a with λC = 3.86159e-13 m.")
    ap.add_argument("--r-mult", type=float, default=3.0, help="rmax = r_mult*a.")

    ap.add_argument("--c-ref", type=float, default=1.0)
    ap.add_argument("--rho", type=float, default=1.0)
    ap.add_argument("--omega", type=float, default=1.0)

    ap.add_argument("--x4-amp", type=float, default=None, help="X^4 bulge amplitude (same length unit as R).")
    ap.add_argument("--ur-amp", type=float, default=None, help="Radial displacement amplitude (same length unit as R).")
    ap.add_argument("--Kr", type=float, default=1.0)

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

    x4_amp = args.x4_amp if args.x4_amp is not None else 0.8 * a
    ur_amp = args.ur_amp if args.ur_amp is not None else 0.2 * a

    p = SliceParams(R=R, a=a, w0=w0, c_ref=args.c_ref, rho=args.rho, omega=args.omega,
                    x4_amp=x4_amp, ur_amp=ur_amp, Kr=args.Kr)

    rmax = args.r_mult * a
    r = np.linspace(-rmax, +rmax, 2000) if args.domain == "signed" else np.linspace(0.0, rmax, 2000)

    prof = build_profiles(p, r)

    # x coordinate for plotting
    x = (R + r) if args.x_axis == "global" else r
    x_plot = x * scale
    x_label = (f"Global radial coordinate ρ = R + r [{unit_label}]"
               if args.x_axis == "global" else
               f"Local radial coordinate r [{unit_label}] (from centerline outward)")

    # scale geometry arrays to chosen units
    prof_plot = dict(prof)
    if unit_label != "arb. units":
        prof_plot["x4"] = np.asarray(prof["x4"]) * scale
        prof_plot["u_r"] = np.asarray(prof["u_r"]) * scale

    # marker positions
    if args.domain == "signed":
        marker_rs = [-w0, +w0]
        marker_labels = ["inner (-w0)", "outer (+w0)"]
    else:
        marker_rs = [w0]
        marker_labels = ["outer (+w0)"]

    marker_x = [(R + rr) if args.x_axis == "global" else rr for rr in marker_rs]
    marker_x_plot = [mx * scale for mx in marker_x]

    print("\nRadial slice constraint report (v2)")
    print("----------------------------------")
    print(f"R={R:.6e}, a={a:.6e}, w0={w0:.6e}  domain={args.domain}  x-axis={args.x_axis}")
    print(f"Target Δn/n≈a/R = {prof['target_dn_over_n']:.6f}   achieved = {prof['achieved_dn_over_n']:.6f}")
    print(f"Target ratio n(-w0)/n(+w0) = {prof['target_ratio']:.6f}   achieved = {prof['achieved_ratio']:.6f}")
    print(f"Δn_core = {prof['delta_n_core']:.6f}   alpha = {prof['alpha']:.6g}\n")

    plot_geometry(x_plot, x_label, unit_label, outdir, prof_plot, p, marker_x_plot, marker_labels)
    plot_propagation(x_plot, x_label, outdir, prof, p, marker_x_plot, marker_labels)
    plot_stresses(x_plot, x_label, outdir, prof, p, marker_x_plot, marker_labels)

    print(f"Saved PNGs to: {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
