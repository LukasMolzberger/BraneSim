"""Parameter search for a *dynamical* (undamped) electron-like initial state.

Standalone test (no modifications to existing modules):
- Builds a 3D brane-in-4D state
- Initializes an electron mode using ElectronModeSpec/initialize_electron_mode_3d
- Sweeps parameters (including rest_length) and runs short conservative simulations
- Scores candidates by how well the envelope stays localized (low leakage, low radius growth)

Outputs (via TestRunManager):
- CSV of all candidates + score
- Plots: score vs rest_length, best-candidate leakage/radius curves
- JSON of best parameters

Run:
    python tests/test_electron_equilibrium_parameter_search.py

Tuning knobs are at the top of `main()`.
"""

from __future__ import annotations

import os
import sys
import json
import math
import time
from dataclasses import dataclass
from typing import Dict, Any, List, Tuple

import numpy as np

# Ensure project root is on sys.path
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)
sys.path.append(PROJECT_ROOT)

import torch
import matplotlib.pyplot as plt

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.core.dimensions import MassModel
from branesim.physics.forces import SpringForceComputer
from branesim.utils import TestRunManager

from branesim.electron.electron_initialization import ElectronModeSpec, initialize_electron_mode_3d


@dataclass
class Candidate:
    rest_length_sim: float
    amplitude: float
    radius_sim: float
    containment_depth: float
    containment_sigma: float
    polarization: str = "spatial_x4"
    l: int = 1
    m: int = 1
    n: int = 1


def _pick_device_and_dtype() -> Tuple[torch.device, torch.dtype]:
    # Mirror your experiments: prefer MPS, then CUDA, then CPU.
    if torch.backends.mps.is_available():
        return torch.device("mps"), torch.float32
    if torch.cuda.is_available():
        return torch.device("cuda"), torch.float32
    return torch.device("cpu"), torch.float64



def _energy_to_float(E: Any) -> float:
    """Convert solver.compute_energy output to a scalar float.

    The codebase sometimes returns a dict with components (e.g. kinetic/potential),
    or a torch scalar, depending on solver/energy implementation.
    This helper keeps the test robust without modifying core modules.
    """
    if E is None:
        return float("nan")

    # torch / numpy scalar
    if hasattr(E, "item") and callable(getattr(E, "item")):
        try:
            return float(E.item())
        except Exception:
            pass

    # plain numbers
    if isinstance(E, (float, int, np.floating, np.integer)):
        return float(E)

    # dict energies (common patterns)
    if isinstance(E, dict):
        for k in ("total", "total_energy", "energy", "E", "sum", "total_J"):
            if k in E:
                return _energy_to_float(E[k])

        if "kinetic" in E and "potential" in E:
            return _energy_to_float(E["kinetic"]) + _energy_to_float(E["potential"])
        if "kinetic_energy" in E and "potential_energy" in E:
            return _energy_to_float(E["kinetic_energy"]) + _energy_to_float(E["potential_energy"])

        s = 0.0
        found = False
        for v in E.values():
            try:
                s += _energy_to_float(v)
                found = True
            except Exception:
                continue
        return float(s) if found else float("nan")

    # sequences: sum what we can
    if isinstance(E, (list, tuple)):
        s = 0.0
        found = False
        for v in E:
            try:
                s += _energy_to_float(v)
                found = True
            except Exception:
                continue
        return float(s) if found else float("nan")

    try:
        return float(E)
    except Exception:
        return float("nan")

def _precompute_radial_distance(grid: BraneGrid, center_xyz: Tuple[float, float, float], device: torch.device) -> torch.Tensor:
    coords = grid.get_spatial_coordinates()  # (N,3) tensor
    cx, cy, cz = center_xyz
    dx = coords[:, 0] - cx
    dy = coords[:, 1] - cy
    dz = coords[:, 2] - cz
    r = torch.sqrt(dx * dx + dy * dy + dz * dz)
    return r.to(device)


def _project_quadratures(dX: torch.Tensor, p1: torch.Tensor, p2: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """a = dX·p1, b = dX·p2 (both are (N,) tensors)."""
    a = torch.sum(dX * p1, dim=-1)
    b = torch.sum(dX * p2, dim=-1)
    # remove global DC offsets (mean deformation) to focus on the oscillatory envelope
    a = a - torch.mean(a)
    b = b - torch.mean(b)
    return a, b


def _leakage_and_rms_radius(weight: torch.Tensor, r: torch.Tensor, radius: float) -> Tuple[float, float]:
    """Return (leakage_fraction_outside_R, rms_radius_over_R)."""
    w = torch.clamp(weight, min=0.0)
    wsum = torch.sum(w)
    if float(wsum.item()) <= 0.0:
        return 1.0, float("inf")

    outside = (r > radius)
    w_out = torch.sum(w[outside])
    leak = float((w_out / wsum).item())

    r2 = r * r
    rms = torch.sqrt(torch.sum(w * r2) / wsum)
    rms_over_R = float((rms / radius).item())
    return leak, rms_over_R


def _score_series(leak: List[float], rms: List[float]) -> float:
    """Lower is better. Emphasize boundedness and low leakage."""
    leak0, leakN = leak[0], leak[-1]
    rms0, rmsN = rms[0], rms[-1]

    leak_max = max(leak)
    rms_max = max(rms)

    leak_drift = leakN - leak0
    rms_drift = rmsN - rms0

    # Penalize: high leakage, growth, and excursions beyond R.
    score = (
        3.0 * leak_max
        + 1.0 * abs(leak_drift)
        + 1.5 * max(0.0, rms_max - 1.0)
        + 0.5 * abs(rms_drift)
    )
    return float(score)


def evaluate_candidate(
    cand: Candidate,
    *,
    grid_shape: Tuple[int, int, int],
    h_sim: float,
    dt_sim: float,
    steps: int,
    sample_stride: int,
    device: torch.device,
    dtype: torch.dtype,
) -> Dict[str, Any]:
    """Build state/solver, initialize electron, run undamped sim, compute metrics."""

    nx, ny, nz = grid_shape

    # Build state
    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
    state.initialize_flat_configuration(h_sim)
    X_ref = state.positions.clone()

    # Fixed boundary like your experiment
    state.set_fixed_boundaries()

    # Build grid
    grid = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h_sim, device)

    # Physics
    k_sim = 1.0
    m_sim = 1.0
    rest_length_sim = float(cand.rest_length_sim)
    physics = SpringForceComputer(k_sim, rest_length_sim)

    rho_sim = m_sim / (h_sim ** 3)
    mass_model = MassModel.from_density(density=rho_sim, intrinsic_dim=3, spacing=h_sim)
    solver = VelocityVerletSolver(dt_sim, mass_model, physics, grid)

    # Center in sim coordinates (domain is (n-1)*h)
    Lx = (nx - 1) * h_sim
    Ly = (ny - 1) * h_sim
    Lz = (nz - 1) * h_sim
    center = (0.5 * Lx, 0.5 * Ly, 0.5 * Lz)

    # Init electron
    spec = ElectronModeSpec(
        l=cand.l,
        m=cand.m,
        n=cand.n,
        radius=float(cand.radius_sim),
        center=center,
        amplitude=float(cand.amplitude),
        wave_speed=1.0,
        polarization=cand.polarization,
        containment_component=3,
        containment_depth=float(cand.containment_depth),
        containment_sigma=float(cand.containment_sigma),
        smooth_edge=2.0,
    )

    debug = initialize_electron_mode_3d(state, grid, spec, set_velocities=True, normalize=True, return_debug=True)

    # Prepare diagnostics
    p1 = torch.tensor(debug["p1"], device=device, dtype=state.positions.dtype)
    p2 = torch.tensor(debug["p2"], device=device, dtype=state.positions.dtype)

    r = _precompute_radial_distance(grid, center, device=device)

    times = []
    leak_env = []
    rms_env = []
    leak_x4 = []
    rms_x4 = []
    energies = []

    # Evaluate at t=0 too
    def sample(t_sim: float):
        dX = state.positions - X_ref
        a, b = _project_quadratures(dX, p1, p2)
        w_env = a * a + b * b
        leak1, rms1 = _leakage_and_rms_radius(w_env, r, float(cand.radius_sim))

        x4 = dX[:, 3] - torch.mean(dX[:, 3])
        w4 = x4 * x4
        leak2, rms2 = _leakage_and_rms_radius(w4, r, float(cand.radius_sim))

        E = _energy_to_float(solver.compute_energy(state))

        times.append(float(t_sim))
        leak_env.append(leak1)
        rms_env.append(rms1)
        leak_x4.append(leak2)
        rms_x4.append(rms2)
        energies.append(E)

    sample(0.0)

    t_sim = 0.0
    for step in range(1, steps + 1):
        solver.step(state)
        t_sim += dt_sim
        if (step % sample_stride) == 0 or step == steps:
            sample(t_sim)

    # Score
    score_env = _score_series(leak_env, rms_env)
    score_x4 = _score_series(leak_x4, rms_x4)

    # Keep X^4 in check, but prioritize envelope localization
    score = score_env + 0.5 * score_x4

    # Energy drift metric
    E0, EN = energies[0], energies[-1]
    edrift = 0.0
    if abs(E0) > 0:
        edrift = (EN - E0) / E0
    score += 0.2 * abs(edrift)

    return {
        "score": float(score),
        "score_env": float(score_env),
        "score_x4": float(score_x4),
        "energy_initial": float(E0),
        "energy_final": float(EN),
        "energy_drift_frac": float(edrift),
        "rest_length_sim": float(cand.rest_length_sim),
        "amplitude": float(cand.amplitude),
        "radius_sim": float(cand.radius_sim),
        "containment_depth": float(cand.containment_depth),
        "containment_sigma": float(cand.containment_sigma),
        "polarization": cand.polarization,
        "l": cand.l,
        "m": cand.m,
        "n": cand.n,
        "omega": float(debug.get("omega", float("nan"))),
        "k": float(debug.get("k", float("nan"))),
        "p1": debug.get("p1"),
        "p2": debug.get("p2"),
        "times": times,
        "leak_env": leak_env,
        "rms_env": rms_env,
        "leak_x4": leak_x4,
        "rms_x4": rms_x4,
        "energies": energies,
    }


def _latin_hypercube(rng: np.random.Generator, n: int, d: int) -> np.ndarray:
    """Simple LHS in [0,1]^d."""
    u = (rng.random((n, d)) + np.arange(n)[:, None]) / n
    for j in range(d):
        rng.shuffle(u[:, j])
    return u



def _to_jsonable(obj):
    # Recursively convert numpy/torch objects to JSON-serializable Python types.
    try:
        import numpy as _np
    except Exception:
        _np = None
    try:
        import torch as _torch
    except Exception:
        _torch = None

    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj

    # numpy scalars / arrays
    if _np is not None:
        if isinstance(obj, _np.generic):
            return obj.item()
        if isinstance(obj, _np.ndarray):
            return obj.tolist()

    # torch tensors
    if _torch is not None and hasattr(_torch, 'is_tensor') and _torch.is_tensor(obj):
        t = obj.detach().cpu()
        if t.numel() == 1:
            return t.item()
        return t.tolist()

    if isinstance(obj, dict):
        return {str(k): _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]

    # Fallback: try to stringify
    return str(obj)

def main() -> None:
    run_manager = TestRunManager(experiment_name="electron_equilibrium_parameter_search")
    print(run_manager.get_summary())

    # ------------------------------
    # Search configuration
    # ------------------------------
    # Keep this modest: this is meant to be a *diagnostic* search.
    GRID_N = int(os.environ.get("EQUIL_GRID_N", "56"))  # 56^3 ~ 175k points
    STEPS = int(os.environ.get("EQUIL_STEPS", "300"))
    SAMPLE_STRIDE = int(os.environ.get("EQUIL_SAMPLE_STRIDE", "10"))

    # Candidate budget
    N_CANDIDATES = int(os.environ.get("EQUIL_CANDIDATES", "36"))
    SEED = int(os.environ.get("EQUIL_SEED", "0"))

    # Sim units
    h_sim = 1.0
    dt_sim = float(os.environ.get("EQUIL_DT", "0.1"))

    # Electron shape
    radius_sim = float(os.environ.get("EQUIL_RADIUS", str(max(8.0, GRID_N / 5.0))))

    device, dtype = _pick_device_and_dtype()
    print(f"\nDevice: {device}  dtype: {dtype}")
    print(f"Grid: {GRID_N}^3  steps: {STEPS}  stride: {SAMPLE_STRIDE}  candidates: {N_CANDIDATES}")

    grid_shape = (GRID_N, GRID_N, GRID_N)

    rng = np.random.default_rng(SEED)

    # Parameter domains
    # rest_length is in sim units, relative to spacing h=1.0.
    # Use tension fraction τ = h - L0 to sample both high and moderate pretension.
    # τ in [1e-3, 1.0] => L0 in [0.0, 0.999]
    tension_min = float(os.environ.get("EQUIL_TENSION_MIN", "1e-3"))
    tension_max = float(os.environ.get("EQUIL_TENSION_MAX", "1.0"))

    # containment depth/sigma
    depth_min = float(os.environ.get("EQUIL_DEPTH_MIN", "0.0"))
    depth_max = float(os.environ.get("EQUIL_DEPTH_MAX", "0.8"))

    sigma_min = float(os.environ.get("EQUIL_SIGMA_MIN", str(0.25 * radius_sim)))
    sigma_max = float(os.environ.get("EQUIL_SIGMA_MAX", str(1.00 * radius_sim)))

    amp_min = float(os.environ.get("EQUIL_AMP_MIN", "0.15"))
    amp_max = float(os.environ.get("EQUIL_AMP_MAX", "0.8"))

    # Build LHS samples for (tension, depth, sigma, amplitude)
    U = _latin_hypercube(rng, N_CANDIDATES, 4)

    results: List[Dict[str, Any]] = []

    t_start = time.time()
    for idx in range(N_CANDIDATES):
        u_t, u_d, u_s, u_a = U[idx]

        # Sample tension logarithmically for better coverage of small values
        log_tmin = math.log10(tension_min)
        log_tmax = math.log10(tension_max)
        tension = 10 ** (log_tmin + (log_tmax - log_tmin) * float(u_t))
        rest_length = max(1e-12, min(0.999999, 1.0 - tension))

        depth = depth_min + (depth_max - depth_min) * float(u_d)
        sigma = sigma_min + (sigma_max - sigma_min) * float(u_s)
        amp = amp_min + (amp_max - amp_min) * float(u_a)

        cand = Candidate(
            rest_length_sim=float(rest_length),
            amplitude=float(amp),
            radius_sim=float(radius_sim),
            containment_depth=float(depth),
            containment_sigma=float(sigma),
            polarization="spatial_x4",
            l=1,
            m=1,
            n=1,
        )

        print(
            f"\n[{idx+1:02d}/{N_CANDIDATES}] rest_length={cand.rest_length_sim:.3e}  "
            f"depth={cand.containment_depth:.3f}  sigma={cand.containment_sigma:.2f}  amp={cand.amplitude:.3f}"
        )

        out = evaluate_candidate(
            cand,
            grid_shape=grid_shape,
            h_sim=h_sim,
            dt_sim=dt_sim,
            steps=STEPS,
            sample_stride=SAMPLE_STRIDE,
            device=device,
            dtype=dtype,
        )
        print(
            f"    score={out['score']:.4f}  leak_end={out['leak_env'][-1]:.3f}  "
            f"rms_end={out['rms_env'][-1]:.3f}  E_drift={out['energy_drift_frac']:.3e}"
        )

        # Drop heavy time-series for non-top candidates later; keep for now.
        results.append(out)

        # Help MPS/CUDA memory fragmentation in loops
        if device.type == "mps":
            try:
                torch.mps.empty_cache()
            except Exception:
                pass
        if device.type == "cuda":
            torch.cuda.empty_cache()

    elapsed = time.time() - t_start
    print(f"\nSearch finished in {elapsed:.1f}s")

    # Sort by score
    results_sorted = sorted(results, key=lambda r: r["score"])
    best = results_sorted[0]

    # Save CSV (flattened)
    csv_path = os.path.join(run_manager.data_dir, "electron_equilibrium_candidates.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        header = [
            "score",
            "score_env",
            "score_x4",
            "rest_length_sim",
            "amplitude",
            "radius_sim",
            "containment_depth",
            "containment_sigma",
            "energy_drift_frac",
            "omega",
            "k",
        ]
        f.write(",".join(header) + "\n")
        for r in results_sorted:
            f.write(",".join(
                [
                    f"{r['score']}",
                    f"{r['score_env']}",
                    f"{r['score_x4']}",
                    f"{r['rest_length_sim']}",
                    f"{r['amplitude']}",
                    f"{r['radius_sim']}",
                    f"{r['containment_depth']}",
                    f"{r['containment_sigma']}",
                    f"{r['energy_drift_frac']}",
                    f"{r['omega']}",
                    f"{r['k']}",
                ]
            ) + "\n")

    # Save best params JSON (without the heavy time-series)
    best_light = {k: v for k, v in best.items() if k not in ("times", "leak_env", "rms_env", "leak_x4", "rms_x4", "energies")}
    json_path = os.path.join(run_manager.data_dir, "electron_equilibrium_best_params.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(_to_jsonable(best_light), f, indent=2)

    # Print top 5
    print("\nTop 5 candidates:")
    for j, r in enumerate(results_sorted[:5], 1):
        print(
            f"  #{j}: score={r['score']:.4f}  L0={r['rest_length_sim']:.3e}  "
            f"depth={r['containment_depth']:.3f}  sigma={r['containment_sigma']:.2f}  amp={r['amplitude']:.3f}"
        )

    # Plot: score vs rest_length and tension
    rest = np.array([r["rest_length_sim"] for r in results_sorted], dtype=float)
    score = np.array([r["score"] for r in results_sorted], dtype=float)
    tension = 1.0 - rest

    # Tension fraction is the directly relevant pretension control parameter: tau = h - L0 (with h=1).
    plt.figure(figsize=(8, 4))
    plt.scatter(tension, score, s=25)
    plt.xscale("log")
    plt.xlabel("tension fraction tau = (h - L0)  (h=1)")
    plt.ylabel("score (lower is better)")
    plt.title("Electron equilibrium search: score vs pretension")
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path("equilibrium_score_vs_tension.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.scatter(rest, score, s=25)
    plt.xlabel("rest_length_sim (L0)")
    plt.ylabel("score (lower is better)")
    plt.title("Electron equilibrium search: score vs rest_length")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path("equilibrium_score_vs_rest_length_linear.png"), dpi=150)
    plt.close()

    # Plot best leakage/radius curves
    t = np.array(best["times"], dtype=float)
    plt.figure(figsize=(10, 5))
    plt.plot(t, best["leak_env"], label="Leakage env |ψ|²")
    plt.plot(t, best["leak_x4"], label="Leakage X4_dyn²")
    plt.xlabel("time (sim units)")
    plt.ylabel("leakage fraction outside R")
    plt.title("Best candidate leakage")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path("equilibrium_best_leakage.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(10, 5))
    plt.plot(t, best["rms_env"], label="RMS radius / R (env)")
    plt.plot(t, best["rms_x4"], label="RMS radius / R (X4_dyn)")
    plt.axhline(1.0, linestyle="--", linewidth=1)
    plt.xlabel("time (sim units)")
    plt.ylabel("RMS radius / R")
    plt.title("Best candidate radius")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(run_manager.get_plot_path("equilibrium_best_radius.png"), dpi=150)
    plt.close()

    # Save config
    run_manager.save_config(
        {
            "grid_n": GRID_N,
            "steps": STEPS,
            "sample_stride": SAMPLE_STRIDE,
            "n_candidates": N_CANDIDATES,
            "seed": SEED,
            "dt_sim": dt_sim,
            "radius_sim": radius_sim,
            "tension_min": tension_min,
            "tension_max": tension_max,
            "depth_min": depth_min,
            "depth_max": depth_max,
            "sigma_min": sigma_min,
            "sigma_max": sigma_max,
            "amp_min": amp_min,
            "amp_max": amp_max,
            "best": best_light,
        }
    )

    print(f"\n✓ Saved results to: {run_manager.run_dir}")
    print(f"  - {csv_path}")
    print(f"  - {json_path}")


if __name__ == "__main__":
    main()
