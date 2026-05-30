"""Numerical dispersion on the 6-neighbor axial-only lattice (Sprint 1 #2 + #3).

Measures `ω(k)` for plane-wave standing-wave initial conditions along the
[100] and [111] directions at the canonical operating point
``(α=0.2, k₀=ρ=a=1)`` and extrapolates to ``k → 0`` to recover the
long-wavelength branch speeds:

  - c_L([100])   ≡ longitudinal speed along [100]              (analytic: 1.0000)
  - c_T([100])   ≡ transverse speed along [100]                (analytic: √(1-α) ≈ 0.8944)
  - c_T([111])   ≡ triplet (Cartesian-pol) speed along [111]   (analytic: √((3-2α)/3) ≈ 0.9309)
  - ratio        ≡ c_L([100]) / c_T([111])                     (analytic: √(3/(3-2α)) ≈ 1.0744)

All four numbers fall out of the closed-form `D(k)` derived in
`components/diagnostics/christoffel_6nn.py`. This script does *not* re-derive
them; it verifies that the actual nonlinear Velocity-Verlet integrator
reproduces them at finite amplitude ε = 1e-3 (well inside the linear regime,
ε·|k| ≲ 6×10⁻⁴ at the largest k tested).

Method:
  1. Initialize `u(x,0) = ε p̂ cos(k·x)`, `v(x,0) = 0` on a periodic N³ lattice.
  2. Integrate with Velocity-Verlet for ~6 analytic periods.
  3. Project the trajectory onto the standing-wave mode:
         A(t) = (2 / N³) · Σᵢ (u_i · p̂) · cos(k · x_i)
  4. Fit A(t) = ε cos(ω t) → ω.
  5. Repeat over k-modes (n = 1, 2, 3 on N=32), extrapolate c(|k|·a) → c(0)
     via a quadratic in (|k|·a)².

Output:
  - dispersion_results.json    : per-run measured ω, c, residual.
  - dispersion_summary.json    : extrapolated branch speeds + anisotropy ratio.
  - dispersion_raw.npz         : A(t) traces per run for inspection.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import curve_fit

from components.shared import BraneState3D
from components.simulation import (
    BraneGrid3D,
    NodeMassModel,
    SpringForceComputer,
    VelocityVerletSolver,
)
from components.diagnostics.christoffel_6nn import eigvals_6nn


HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------- #
# Analytic predictions (closed form)
# --------------------------------------------------------------------------- #

def analytic_omega(k_vec: np.ndarray, polarization: np.ndarray,
                   alpha: float, k0: float = 1.0, rho: float = 1.0,
                   a: float = 1.0) -> float:
    """Predict ω for a Cartesian-aligned polarization at wavevector k.

    The eigenframe of D(k) is Cartesian at every k; pick the eigenvalue
    matching whichever Cartesian axis the polarization aligns with.
    """
    eigvals = eigvals_6nn(k_vec, alpha, k0=k0, rho=rho, a=a)
    p = np.abs(polarization)
    axis = int(np.argmax(p))
    if p[axis] < 0.99:
        raise ValueError("This script only handles Cartesian-aligned polarizations; "
                         "got %r" % (polarization,))
    return float(np.sqrt(eigvals[axis]))


def analytic_c_long_wavelength(direction: str, polarization: str, alpha: float) -> float:
    """Long-wavelength branch speeds, closed form."""
    if direction == "100" and polarization == "L":
        return 1.0
    if direction == "100" and polarization == "T":
        return float(np.sqrt(1.0 - alpha))
    if direction == "111":
        # All three Cartesian polarizations degenerate along [111].
        return float(np.sqrt((3.0 - 2.0 * alpha) / 3.0))
    raise ValueError(f"Unknown branch: direction={direction!r}, pol={polarization!r}")


# --------------------------------------------------------------------------- #
# Single dispersion measurement
# --------------------------------------------------------------------------- #

@dataclass
class RunSpec:
    label: str
    direction: str           # "100" or "111"
    polarization_tag: str    # "L" or "T"
    N: int
    k_index: tuple[int, int, int]
    polarization: tuple[float, float, float]
    alpha: float
    n_periods_target: float
    dt: float
    checkpoint_stride: int


@dataclass
class RunResult:
    label: str
    direction: str
    polarization_tag: str
    N: int
    k_index: list[int]
    polarization: list[float]
    alpha: float
    k_mag: float
    k_dot_a: float
    omega_analytic: float
    omega_measured: float
    omega_err: float
    c_analytic: float
    c_measured: float
    c_rel_err: float
    num_steps: int
    num_periods_run: float
    fit_residual_max: float
    wall_time_s: float


def run_dispersion(spec: RunSpec, k0: float = 1.0, rho: float = 1.0,
                   spacing: float = 1.0, epsilon: float = 1e-3) -> tuple[RunResult, np.ndarray, np.ndarray]:
    device = torch.device("cpu")
    dtype = torch.float64

    grid_shape = (spec.N, spec.N, spec.N)
    state = BraneState3D(grid_shape, device=device, dtype=dtype)
    state.initialize_flat_configuration(spacing)
    # Free / periodic everywhere — no fixed boundaries.

    grid = BraneGrid3D(
        grid_shape=grid_shape,
        spacing=spacing,
        device=device,
        periodic_axes=(True, True, True),
        axial_weight=1.0,
    )
    physics = SpringForceComputer(spring_constant=k0, rest_length=spec.alpha * spacing)
    mass_model = NodeMassModel.from_density(density=rho, spacing=spacing)
    solver = VelocityVerletSolver(dt=spec.dt, mass_model=mass_model, physics=physics, grid=grid)

    coords = state.positions[:, :3].detach().cpu().numpy()  # rest positions, shape (N³, 3)
    L_box = np.array([spec.N * spacing] * 3)
    k_vec = 2.0 * np.pi * np.asarray(spec.k_index, dtype=np.float64) / L_box
    k_mag = float(np.linalg.norm(k_vec))

    polarization = np.asarray(spec.polarization, dtype=np.float64)
    polarization = polarization / np.linalg.norm(polarization)

    phase = coords @ k_vec
    cos_phase = np.cos(phase)
    sin_phase = np.sin(phase)

    u_np = np.zeros((coords.shape[0], 4), dtype=np.float64)
    u_np[:, 0] = polarization[0] * epsilon * cos_phase
    u_np[:, 1] = polarization[1] * epsilon * cos_phase
    u_np[:, 2] = polarization[2] * epsilon * cos_phase

    state.set_kinematics(
        torch.from_numpy(u_np).to(device=device, dtype=dtype),
        torch.zeros_like(state.velocities),
    )
    solver.initialize_accelerations(state)

    omega_analytic = analytic_omega(k_vec, polarization, spec.alpha, k0=k0, rho=rho, a=spacing)
    T_pred = 2.0 * np.pi / omega_analytic
    num_steps = max(200, int(np.ceil(spec.n_periods_target * T_pred / spec.dt)))

    # Periodic projection onto the cos(k·x) standing-wave envelope, normalised so A(0) = epsilon.
    # (Σ_i cos²(k·x_i) ≈ N³/2 for non-trivial k, so factor = 2/N³.)
    npoints = coords.shape[0]
    norm_factor = 2.0 / npoints

    # Convert basis vectors to torch tensors once.
    cos_phase_t = torch.from_numpy(cos_phase).to(device=device, dtype=dtype)
    sin_phase_t = torch.from_numpy(sin_phase).to(device=device, dtype=dtype)
    pol_t = torch.from_numpy(polarization[:3].astype(np.float64)).to(device=device, dtype=dtype)
    rest_xyz = state.rest_positions[:, :3]

    times: list[float] = []
    A_cos: list[float] = []
    A_sin: list[float] = []

    def record():
        u_xyz = state.positions[:, :3] - rest_xyz
        u_pol = (u_xyz * pol_t).sum(dim=1)
        Ac = float(norm_factor * (u_pol * cos_phase_t).sum().item())
        As = float(norm_factor * (u_pol * sin_phase_t).sum().item())
        times.append(float(solver.time))
        A_cos.append(Ac)
        A_sin.append(As)

    record()
    t0 = time.time()
    for step in range(1, num_steps + 1):
        solver.step(state)
        if step % spec.checkpoint_stride == 0:
            record()
    wall = time.time() - t0

    times_arr = np.asarray(times)
    A_cos_arr = np.asarray(A_cos)
    A_sin_arr = np.asarray(A_sin)

    # The standing-wave initial condition (u·p̂ = ε cos(k·x), v = 0) populates only the
    # cos(k·x) channel; in the linearised problem this stays at A_cos(t) = ε cos(ω t).
    # Fit that one-parameter model.
    def model(t, omega):
        return epsilon * np.cos(omega * t)

    popt, pcov = curve_fit(model, times_arr, A_cos_arr, p0=[omega_analytic])
    omega_meas = float(popt[0])
    omega_err = float(np.sqrt(pcov[0, 0]))
    fit_residual = float(np.max(np.abs(A_cos_arr - model(times_arr, omega_meas))))

    c_analytic_branch = analytic_c_long_wavelength(spec.direction, spec.polarization_tag, spec.alpha)
    c_analytic_at_k = omega_analytic / k_mag
    c_meas = omega_meas / k_mag

    result = RunResult(
        label=spec.label,
        direction=spec.direction,
        polarization_tag=spec.polarization_tag,
        N=spec.N,
        k_index=list(int(v) for v in spec.k_index),
        polarization=list(float(v) for v in polarization),
        alpha=float(spec.alpha),
        k_mag=k_mag,
        k_dot_a=float(k_mag * spacing),
        omega_analytic=omega_analytic,
        omega_measured=omega_meas,
        omega_err=omega_err,
        c_analytic=float(c_analytic_at_k),
        c_measured=float(c_meas),
        c_rel_err=float((c_meas - c_analytic_at_k) / c_analytic_at_k),
        num_steps=num_steps,
        num_periods_run=float(solver.time / T_pred),
        fit_residual_max=fit_residual,
        wall_time_s=wall,
    )
    # Tag with long-wavelength branch speed too (for k=0 extrapolation later).
    return result, times_arr, A_cos_arr


# --------------------------------------------------------------------------- #
# Sweep configuration
# --------------------------------------------------------------------------- #

ALPHA = 0.2
DT = 0.01
N_PERIODS = 6.0
CHECKPOINT_STRIDE = 20
N_GRID = 32
K_INDEX_RANGE = (1, 2, 3)


def make_sweep() -> list[RunSpec]:
    specs: list[RunSpec] = []
    for n in K_INDEX_RANGE:
        # [100] longitudinal — pol along [100] is x-axis eigenmode → c_L = 1.
        specs.append(RunSpec(
            label=f"100_L_n{n}",
            direction="100", polarization_tag="L",
            N=N_GRID, k_index=(n, 0, 0), polarization=(1.0, 0.0, 0.0),
            alpha=ALPHA, n_periods_target=N_PERIODS, dt=DT,
            checkpoint_stride=CHECKPOINT_STRIDE,
        ))
        # [100] transverse — pol perp to k. Use y-axis (eigenmode) at k·a from x-axis.
        specs.append(RunSpec(
            label=f"100_T_n{n}",
            direction="100", polarization_tag="T",
            N=N_GRID, k_index=(n, 0, 0), polarization=(0.0, 1.0, 0.0),
            alpha=ALPHA, n_periods_target=N_PERIODS, dt=DT,
            checkpoint_stride=CHECKPOINT_STRIDE,
        ))
        # [111] — all three Cartesian polarizations are eigenmodes with the same ω.
        # Run all three to verify the triplet degeneracy.
        for axis_tag, axis_vec in (("x", (1.0, 0.0, 0.0)),
                                   ("y", (0.0, 1.0, 0.0)),
                                   ("z", (0.0, 0.0, 1.0))):
            specs.append(RunSpec(
                label=f"111_p{axis_tag}_n{n}",
                direction="111", polarization_tag="T",
                N=N_GRID, k_index=(n, n, n), polarization=axis_vec,
                alpha=ALPHA, n_periods_target=N_PERIODS, dt=DT,
                checkpoint_stride=CHECKPOINT_STRIDE,
            ))
    return specs


# --------------------------------------------------------------------------- #
# k → 0 extrapolation
# --------------------------------------------------------------------------- #

def extrapolate_to_zero(c_vs_ka: list[tuple[float, float]]) -> tuple[float, float, list[float]]:
    """Fit c(|k|·a) = c0 + c2 · (|k|·a)² and return (c0, c2, residuals).

    Uses ordinary least squares on at least 2 points; with 3+ points the fit
    is over-determined and the residual is informative.
    """
    arr = np.asarray(c_vs_ka, dtype=np.float64)
    ka = arr[:, 0]
    c = arr[:, 1]
    A = np.stack([np.ones_like(ka), ka ** 2], axis=1)
    sol, *_ = np.linalg.lstsq(A, c, rcond=None)
    c0, c2 = float(sol[0]), float(sol[1])
    fit = A @ sol
    residuals = (c - fit).tolist()
    return c0, c2, residuals


def main():
    specs = make_sweep()
    print(f"Running {len(specs)} dispersion sims at α={ALPHA}, N={N_GRID}, dt={DT}, ~{N_PERIODS} periods each")

    raw_traces: dict[str, dict] = {}
    results: list[RunResult] = []
    t_start = time.time()
    for i, spec in enumerate(specs, 1):
        print(f"  [{i:2d}/{len(specs)}] {spec.label}  k_idx={spec.k_index} pol={spec.polarization}  ...", flush=True)
        result, ts, A_cos = run_dispersion(spec)
        results.append(result)
        raw_traces[spec.label] = {
            "times": ts.tolist(),
            "A_cos": A_cos.tolist(),
        }
        print(f"      ω_pred={result.omega_analytic:.6f}  ω_meas={result.omega_measured:.6f}  "
              f"c_meas={result.c_measured:.6f}  rel_err={result.c_rel_err:+.3e}  "
              f"wall={result.wall_time_s:.1f}s")
    t_total = time.time() - t_start
    print(f"Total wall time: {t_total:.1f}s")

    # --- Persist per-run results --------------------------------------------
    HERE.mkdir(parents=True, exist_ok=True)
    with (HERE / "dispersion_results.json").open("w") as fh:
        json.dump(
            {
                "alpha": ALPHA,
                "N": N_GRID,
                "dt": DT,
                "amplitude": 1e-3,
                "k0": 1.0,
                "rho": 1.0,
                "spacing": 1.0,
                "runs": [asdict(r) for r in results],
                "wall_time_s_total": float(t_total),
            },
            fh, indent=2,
        )

    np.savez_compressed(HERE / "dispersion_raw.npz",
                        **{k: np.asarray(v["A_cos"]) for k, v in raw_traces.items()},
                        **{f"_t__{k}": np.asarray(v["times"]) for k, v in raw_traces.items()})

    # --- Per-branch summary + extrapolation ---------------------------------
    branches = {
        "100_L": [r for r in results if r.direction == "100" and r.polarization_tag == "L"],
        "100_T": [r for r in results if r.direction == "100" and r.polarization_tag == "T"],
        "111_T": [r for r in results if r.direction == "111"],   # all three pol axes pooled
    }

    summary = {"alpha": ALPHA, "analytic": {}, "measured": {}, "extrapolated": {}}
    for tag, runs in branches.items():
        if not runs:
            continue
        c_vs_ka = sorted([(r.k_dot_a, r.c_measured) for r in runs])
        c0, c2, residuals = extrapolate_to_zero(c_vs_ka)
        if tag == "100_L":
            c_analytic_0 = 1.0
        elif tag == "100_T":
            c_analytic_0 = float(np.sqrt(1.0 - ALPHA))
        else:
            c_analytic_0 = float(np.sqrt((3.0 - 2.0 * ALPHA) / 3.0))
        summary["analytic"][tag] = c_analytic_0
        summary["measured"][tag] = [{"k_dot_a": ka, "c": c} for ka, c in c_vs_ka]
        summary["extrapolated"][tag] = {
            "c0": c0,
            "c2": c2,
            "residuals": residuals,
            "rel_err_vs_analytic": float((c0 - c_analytic_0) / c_analytic_0),
        }

    # Anisotropy ratio (subtask 3)
    c_100_L = summary["extrapolated"]["100_L"]["c0"]
    c_111   = summary["extrapolated"]["111_T"]["c0"]
    ratio_meas = c_100_L / c_111
    ratio_analytic = float(np.sqrt(3.0 / (3.0 - 2.0 * ALPHA)))
    summary["anisotropy"] = {
        "ratio_measured": ratio_meas,
        "ratio_analytic": ratio_analytic,
        "rel_err": float((ratio_meas - ratio_analytic) / ratio_analytic),
    }

    # Triplet degeneracy along [111] — spread across the three Cartesian polarizations
    # at the smallest k-mode tested (least nonlinear / finite-ka contamination).
    triplet = {}
    for n in K_INDEX_RANGE:
        runs_n = [r for r in results if r.direction == "111" and tuple(r.k_index) == (n, n, n)]
        if len(runs_n) == 3:
            omegas = np.asarray([r.omega_measured for r in runs_n])
            cs = np.asarray([r.c_measured for r in runs_n])
            triplet[f"n{n}"] = {
                "k_dot_a": float(runs_n[0].k_dot_a),
                "omega_per_pol": {r.label[-2:]: r.omega_measured for r in runs_n},
                "omega_mean": float(omegas.mean()),
                "omega_std": float(omegas.std()),
                "relative_spread": float(omegas.std() / omegas.mean()),
                "c_mean": float(cs.mean()),
            }
    summary["triplet_degeneracy_111"] = triplet

    with (HERE / "dispersion_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    # --- Console summary ----------------------------------------------------
    print("\nExtrapolated branch speeds (k → 0):")
    for tag in ("100_L", "100_T", "111_T"):
        c0 = summary["extrapolated"][tag]["c0"]
        ca = summary["analytic"][tag]
        err = summary["extrapolated"][tag]["rel_err_vs_analytic"]
        print(f"  {tag}:  c_meas = {c0:.6f}   c_analytic = {ca:.6f}   rel_err = {err:+.3e}")

    print(f"\nAnisotropy ratio c_L([100]) / c([111]) at α={ALPHA}:")
    print(f"  measured:  {ratio_meas:.6f}")
    print(f"  analytic:  {ratio_analytic:.6f}  (≡ √(3/(3-2α)))")
    print(f"  rel_err :  {summary['anisotropy']['rel_err']:+.3e}")

    print("\n[111] triplet degeneracy:")
    for tag, info in triplet.items():
        print(f"  {tag} (k·a≈{info['k_dot_a']:.3f}):  "
              f"σ(ω)/⟨ω⟩ = {info['relative_spread']:.3e}")


if __name__ == "__main__":
    main()