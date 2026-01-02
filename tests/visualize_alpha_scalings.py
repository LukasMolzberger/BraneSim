# tests/visualize_alpha_scalings.py
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Helpers: robust imports
# -----------------------------
def _get_attr_any(obj, names: Sequence[str], default=None):
    for n in names:
        if hasattr(obj, n):
            return getattr(obj, n)
    return default


def load_physical_constants() -> Tuple[float, float]:
    """
    Returns (lambda_C, c) from your project's physical_constants.py if possible,
    otherwise returns sensible SI defaults.

    Tries several import paths + several common constant names.
    """
    lambda_c = None
    c = None

    # Try a few common module paths.
    candidates = [
        "physical_constants",
        "branesim.physical_constants",
        "branesim.physics.physical_constants",
        "branesim.constants.physical_constants",
    ]

    for mod in candidates:
        try:
            pc = __import__(mod, fromlist=["*"])
            lambda_c = _get_attr_any(
                pc,
                ["COMPTON_WAVELENGTH", "compton_wavelength", "LAMBDA_C", "lambda_C", "LAMBDA_COMPTON"],
                None,
            )
            c = _get_attr_any(
                pc,
                ["SPEED_OF_LIGHT", "speed_of_light", "C", "c"],
                None,
            )
            if lambda_c is not None and c is not None:
                break
        except Exception:
            continue

    # Fallbacks (SI)
    if lambda_c is None:
        # electron Compton wavelength (reduced): 3.8615926796e-13 m (CODATA-ish)
        lambda_c = 3.8615926796e-13
    if c is None:
        c = 299_792_458.0

    return float(lambda_c), float(c)


# -----------------------------
# TestRunManager integration (best-effort)
# -----------------------------
class _FallbackRunManager:
    """
    Minimal replacement if your TestRunManager isn't importable.

    Creates a run folder in: tests/_runs/<script>/<timestamp>/
    Provides .run_dir and .savefig(fig, name).
    """
    def __init__(self, script_name: str, root: Optional[Path] = None):
        self.script_name = script_name
        self.root = root or Path(__file__).resolve().parent / "_runs"
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_dir = self.root / script_name / ts
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def savefig(self, fig, name: str, dpi: int = 160):
        out = self.run_dir / f"{name}.png"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"[saved] {out}")


def get_run_manager(script_name: str):
    """
    Try to import your real TestRunManager, else fallback.

    Adjust/extend the import list if your class lives elsewhere.
    """
    import_paths = [
        ("branesim.testing.test_run_manager", "TestRunManager"),
        ("branesim.tests.test_run_manager", "TestRunManager"),
        ("branesim.test_run_manager", "TestRunManager"),
        ("tests.test_run_manager", "TestRunManager"),
    ]

    for mod, cls in import_paths:
        try:
            m = __import__(mod, fromlist=[cls])
            TRM = getattr(m, cls)
            return TRM(script_name=script_name)  # type: ignore
        except Exception:
            continue

    return _FallbackRunManager(script_name=script_name)


# -----------------------------
# Model formulas used for plots
# -----------------------------
def amplitude_threshold_discrete(alpha: np.ndarray, lam: float) -> np.ndarray:
    """
    Discrete-derived 'onset' estimate (from earlier expansion):
        A_th ~ (λ/π) * sqrt((1-α)/α)

    Interpretation: amplitude at wavelength λ required such that the
    quartic geometric term is comparable to the quadratic tension term
    (within the small-slope expansion logic).
    """
    return (lam / math.pi) * np.sqrt((1.0 - alpha) / alpha)


def amplitude_threshold_continuum(alpha: np.ndarray, lam: float) -> np.ndarray:
    """
    Continuum SVK incremental estimate (pre-tension dominates):
    Using scaling W4/W2 ~ α^2 |∇w|^2, onset ~ |∇w| ~ 1/α.
    For sinusoid w=A sin(2πx/λ): |∇w| ~ 2πA/λ => A_th ~ λ/(2π α).
    """
    return lam / (2.0 * math.pi * alpha)


def wave_speeds(alpha: np.ndarray, c_ref: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Minimal scaling model (matches the 'α→0 => equal speeds' intuition):

      c_L(α) = c_ref
      c_T(α) = c_ref * sqrt(1 - α)

    Here c_ref is the 'longitudinal reference speed' (often you set that to c).

    Notes:
    - This is the clean 1D-chain / tension-shear scaling you discussed.
    - In your full 3D + diagonals implementation, numerical prefactors can differ,
      but the sqrt(1-α) trend is the main point.
    """
    c_L = np.full_like(alpha, c_ref, dtype=float)
    c_T = c_ref * np.sqrt(np.maximum(0.0, 1.0 - alpha))
    return c_L, c_T


def lattice_dispersion(q: np.ndarray, h: float, c_branch: float) -> np.ndarray:
    """
    1D nearest-neighbor lattice dispersion that matches the continuum slope at small q:

        ω(q) = (2 c / h) * sin(q h / 2)

    where c is the continuum-limit wave speed for that branch.
    """
    return (2.0 * c_branch / h) * np.sin(0.5 * q * h)


def lattice_group_velocity(q: np.ndarray, h: float, c_branch: float) -> np.ndarray:
    """
    Group velocity for ω(q) = (2 c / h) sin(qh/2):

        v_g(q) = dω/dq = c * cos(qh/2)
    """
    return c_branch * np.cos(0.5 * q * h)


# -----------------------------
# Config + plotting
# -----------------------------
@dataclass
class PlotConfig:
    alpha_min: float = 1e-4
    alpha_max: float = 0.999
    alpha_count: int = 600

    # Lattice spacing: choose relative to λ_C to check representability.
    # If h > λ_C, then λ_C is sub-grid and q_C will exceed Brillouin zone.
    h_over_lambda_c: float = 0.25

    # Dispersion sampling
    q_count: int = 800
    dispersion_alphas: Tuple[float, ...] = (1e-4, 0.2, 0.5, 0.9, 0.99)

    # If True, set c_ref = speed_of_light for the longitudinal branch
    use_si_speed: bool = True


def plot_amplitude_vs_alpha(cfg: PlotConfig, lam_c: float, run):
    alpha = np.linspace(cfg.alpha_min, cfg.alpha_max, cfg.alpha_count)

    A_disc = amplitude_threshold_discrete(alpha, lam=lam_c)
    A_cont = amplitude_threshold_continuum(alpha, lam=lam_c)

    fig = plt.figure(figsize=(8.2, 5.0))
    ax = fig.add_subplot(111)

    ax.plot(alpha, A_disc / lam_c, label="A_th / λ_C (discrete onset estimate)")
    ax.plot(alpha, A_cont / lam_c, label="A_th / λ_C (continuum SVK scaling)")

    ax.set_xlabel("alpha (pre-stretch parameter)")
    ax.set_ylabel("threshold amplitude  A_th  in units of λ_C")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()

    ax.set_title("Estimated amplitude threshold vs alpha at λ = λ_C")

    _save(run, fig, "alpha_vs_A_threshold")


def plot_speeds_vs_alpha(cfg: PlotConfig, c: float, run):
    alpha = np.linspace(cfg.alpha_min, cfg.alpha_max, cfg.alpha_count)

    c_ref = c if cfg.use_si_speed else 1.0
    c_L, c_T = wave_speeds(alpha, c_ref=c_ref)

    fig = plt.figure(figsize=(8.2, 5.0))
    ax = fig.add_subplot(111)

    ax.plot(alpha, c_L, label="c_L(α) longitudinal (model)")
    ax.plot(alpha, c_T, label="c_T(α) transverse (model)")

    ax.set_xlabel("alpha (pre-stretch parameter)")
    ax.set_ylabel("wave speed (m/s)" if cfg.use_si_speed else "wave speed (arb. units)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()
    ax.set_title("Branch wave speeds vs alpha (simple scaling model)")

    _save(run, fig, "alpha_vs_branch_speeds")


def plot_dispersion(cfg: PlotConfig, lam_c: float, c: float, run):
    h = cfg.h_over_lambda_c * lam_c

    # Brillouin zone (1D NN lattice): q in [0, π/h]
    q_max = math.pi / h
    q = np.linspace(0.0, q_max, cfg.q_count)

    q_C = 2.0 * math.pi / lam_c

    c_ref = c if cfg.use_si_speed else 1.0

    # ω(q) plots
    fig1 = plt.figure(figsize=(8.6, 5.2))
    ax1 = fig1.add_subplot(111)

    # Also plot continuum reference lines ω = c q for comparison (non-dispersive)
    ax1.plot(q, c_ref * q, label="continuum ω = c_L q", linewidth=1.2)

    for a in cfg.dispersion_alphas:
        a = float(a)
        if not (0.0 < a < 1.0):
            continue
        # Branch speeds for that alpha
        c_L = c_ref
        c_T = c_ref * math.sqrt(max(0.0, 1.0 - a))
        wL = lattice_dispersion(q, h=h, c_branch=c_L)
        wT = lattice_dispersion(q, h=h, c_branch=c_T)

        ax1.plot(q, wL, label=f"ω_L(q), α={a:g}")
        ax1.plot(q, wT, label=f"ω_T(q), α={a:g}")

    # Mark Compton wavenumber if within zone
    if q_C <= q_max:
        ax1.axvline(q_C, linestyle="--", linewidth=1.0, label="q_C = 2π/λ_C")
    else:
        ax1.axvline(q_max, linestyle="--", linewidth=1.0, label="q_max = π/h (BZ edge)")
        ax1.text(
            0.02, 0.04,
            "Note: q_C exceeds Brillouin zone for this h.\nλ_C is sub-grid here.",
            transform=ax1.transAxes
        )

    ax1.set_xlabel("wavenumber q (1/m)")
    ax1.set_ylabel("angular frequency ω (1/s)" if cfg.use_si_speed else "ω (arb. units)")
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax1.legend(fontsize=8)
    ax1.set_title("Dispersion relations (1D NN lattice proxy)")

    _save(run, fig1, "dispersion_omega_vs_q")

    # Group velocity plot
    fig2 = plt.figure(figsize=(8.6, 5.2))
    ax2 = fig2.add_subplot(111)

    ax2.plot(q, np.full_like(q, c_ref), label="continuum v_g = c_L", linewidth=1.2)

    for a in cfg.dispersion_alphas:
        a = float(a)
        if not (0.0 < a < 1.0):
            continue
        c_L = c_ref
        c_T = c_ref * math.sqrt(max(0.0, 1.0 - a))
        vgL = lattice_group_velocity(q, h=h, c_branch=c_L)
        vgT = lattice_group_velocity(q, h=h, c_branch=c_T)
        ax2.plot(q, vgL, label=f"v_g,L(q), α={a:g}")
        ax2.plot(q, vgT, label=f"v_g,T(q), α={a:g}")

    if q_C <= q_max:
        ax2.axvline(q_C, linestyle="--", linewidth=1.0, label="q_C = 2π/λ_C")
    else:
        ax2.axvline(q_max, linestyle="--", linewidth=1.0, label="q_max = π/h (BZ edge)")

    ax2.set_xlabel("wavenumber q (1/m)")
    ax2.set_ylabel("group velocity v_g (m/s)" if cfg.use_si_speed else "v_g (arb. units)")
    ax2.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax2.legend(fontsize=8)
    ax2.set_title("Group velocity vs q (lattice dispersion proxy)")

    _save(run, fig2, "dispersion_group_velocity_vs_q")


def _save(run, fig, name: str):
    """
    Save via your run manager if it supports a method, else fallback to .savefig.
    """
    # Common method names in projects
    for meth in ["save_figure", "savefig", "add_figure", "write_figure"]:
        if hasattr(run, meth):
            try:
                getattr(run, meth)(fig, name)  # type: ignore
                plt.close(fig)
                return
            except TypeError:
                # Some managers expect (name, fig) order
                try:
                    getattr(run, meth)(name, fig)  # type: ignore
                    plt.close(fig)
                    return
                except Exception:
                    pass

    # Fallback
    if hasattr(run, "savefig"):
        run.savefig(fig, name)  # type: ignore
    else:
        out = Path(getattr(run, "run_dir", Path.cwd())) / f"{name}.png"
        fig.savefig(out, dpi=160, bbox_inches="tight")
        print(f"[saved] {out}")
    plt.close(fig)


def main():
    script_name = Path(__file__).stem
    cfg = PlotConfig()

    lam_c, c = load_physical_constants()
    h = cfg.h_over_lambda_c * lam_c
    q_C = 2.0 * math.pi / lam_c
    q_max = math.pi / h

    with get_run_manager(script_name) as run:
        # Write config + derived scalars
        run_dir = Path(getattr(run, "run_dir", Path.cwd()))
        meta = {
            "config": asdict(cfg),
            "lambda_c_m": lam_c,
            "speed_of_light_m_per_s": c,
            "h_m": h,
            "q_C_1_per_m": q_C,
            "q_max_1_per_m": q_max,
            "qC_within_BZ": bool(q_C <= q_max),
        }
        (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
        print(f"[run] {run_dir}")
        print(f"[meta] λ_C={lam_c:.6e} m, h={h:.6e} m, q_C={q_C:.6e} 1/m, q_max={q_max:.6e} 1/m")

        plot_amplitude_vs_alpha(cfg, lam_c=lam_c, run=run)
        plot_speeds_vs_alpha(cfg, c=c, run=run)
        plot_dispersion(cfg, lam_c=lam_c, c=c, run=run)


if __name__ == "__main__":
    main()
