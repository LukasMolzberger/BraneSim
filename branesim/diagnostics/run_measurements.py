"""Per-run diagnostics measurement suite.

Single entry point: run_measurements(run_dir) or CLI usage.

Reads a finished worldvolume (seed_world.npz + config.json) from a run folder
and writes into <run_dir>/diagnostics/:

    energy.csv / energy.png               — D1: energy & consistency
    confinement.csv / confinement.png     — D2: confinement metrics
    winding.csv / winding.png             — D3: U(1) phase winding (per slice)
    berry.csv / berry.png                 — D4: Berry/phase envelope time series
    em_fields.png                         — D5: E_i, B_i quiver on mid-plane
    color_channels.csv / color_channels.png — D6: U(1)+SU(3) per-channel split
    spectra.csv / spectra.png             — D7: spatial-FFT energy spectrum
    report.md                             — stitched verdict document

Design constraints (PRINCIPLES.md):
  - Read-only: no modification of solver state, no back-reaction.
  - Vacuum-subtracted excess energy (LESSONS_LEARNED discipline).
  - No hard-coded winding estimator — uses actual field plaquette sum.
  - Dimension-agnostic where possible; dim=3 specialisation only in plot layout.
  - Agg backend throughout; tight_layout on every figure.
  - Residual norm is noted as "N/A — seed only" when world has no interior
    slices of a solved worldvolume.

Usage::

    python -m branesim.diagnostics.run_measurements runs/vortex_seed_YYYY-MM-DD_HHMMSS/

    # or programmatically:
    from branesim.diagnostics.run_measurements import run_measurements
    paths = run_measurements("runs/vortex_seed_2026-06-06_141527")
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.core.action import spacelike_potential
from branesim.diagnostics.confinement import confinement_summary
from branesim.diagnostics.alpha_separability import projection_operators
from branesim.initialization.vortex_worldtube import (
    CARRIER_RE_COMPONENTS,
    CARRIER_RE_WEIGHTS,
    CARRIER_IM,
    measure_winding_closure,
    project_carrier_re,
    vacuum_offsets,
)
from branesim.diagnostics.binding_probe import device_binding_probe
from branesim.diagnostics.contracts import (
    CONTRACTS,
    build_vacuum_world,
    build_ledger,
    render_ledger_md,
)
from branesim.diagnostics._plot_helpers import _apply_style, _savefig, _phase_to_rgb


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------


def _load_run(run_dir: Path) -> dict[str, Any]:
    """Load world array + config from a run folder.

    Returns a dict with keys:
        world      : (n_slices+1, n_nodes, m_ambient) float64
        config     : dict
        lattice    : SpacelikeLattice
        params     : ActionParams
        ref        : (n_nodes, m_ambient) reference positions
        grid_shape : (nx, ny, nz)
        spacing    : float
        alpha      : float
        n_slices   : int
    """
    config_path = run_dir / "config.json"
    # Accept either a solved/iterated world (world.npz) or the bare seed
    # (seed_world.npz); prefer world.npz when both are present.
    world_path = run_dir / "world.npz"
    if not world_path.exists():
        world_path = run_dir / "seed_world.npz"

    if not config_path.exists():
        raise FileNotFoundError(f"config.json not found in {run_dir}")
    if not world_path.exists():
        raise FileNotFoundError(
            f"neither world.npz nor seed_world.npz found in {run_dir}"
        )

    config = json.loads(config_path.read_text())
    grid_shape = tuple(int(v) for v in config["grid_shape"])
    spacing = float(config.get("spacing", 1.0))
    alpha = float(config.get("alpha", 0.7))
    n_slices = int(config.get("n_slices", 32))
    dt = float(config.get("dt", 0.25))
    beta = float(config.get("beta", 1.0))
    r_t = float(config.get("r_t", alpha * beta * dt))

    lattice_params = LatticeParams(
        grid_shape=grid_shape,
        spacing=spacing,
        periodic_axes=(True, True, True),
    )
    action_params = ActionParams(
        k_s=1.0,
        alpha=alpha,
        rho=1.0,
        dt=dt,
        n_slices=n_slices,
        m_ambient=4,
        r_t=r_t,
        beta=beta,
    )
    lattice = SpacelikeLattice(lattice_params)
    m_ambient = action_params.ambient_dim(lattice.dim)
    ref = lattice.reference_positions(m_ambient)

    data = np.load(str(world_path))
    world = data["world"].astype(np.float64)

    return {
        "world": world,
        "config": config,
        "lattice": lattice,
        "params": action_params,
        "ref": ref,
        "grid_shape": grid_shape,
        "spacing": spacing,
        "alpha": alpha,
        "n_slices": n_slices,
    }


# ---------------------------------------------------------------------------
# D1: Energy & consistency
# ---------------------------------------------------------------------------


def _vacuum_energy(ref: np.ndarray, lattice: SpacelikeLattice, params: ActionParams) -> float:
    """Spacelike potential of the vacuum (reference) configuration."""
    return spacelike_potential(ref, lattice, params)


def device_energy(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D1: Energy & consistency.

    Computes per-slice:
      - V_total     : total spacelike potential energy
      - V_vacuum    : vacuum reference (constant)
      - V_excess    : V_total - V_vacuum (vacuum-subtracted, LESSONS_LEARNED)
      - T_kinetic   : (m/2) * |velocity|^2 per slice (central difference)
      - E_total     : V_total + T_kinetic (Hamiltonian proxy)

    Residual norm is noted as N/A for a seed (no BVP solve).
    Vacuum stability check: at ref positions V_excess = 0 by construction.
    """
    n_slices_plus1 = world.shape[0]
    n_slices = n_slices_plus1 - 1
    dt = params.dt

    V_vac = _vacuum_energy(ref, lattice, params)
    V_total = np.empty(n_slices_plus1)
    for l in range(n_slices_plus1):
        V_total[l] = spacelike_potential(world[l], lattice, params)
    V_excess = V_total - V_vac

    # Kinetic energy of the EXCITATION (carrier), not the vacuum.  The prestressed
    # periodic vacuum is an N-gon that rotates in the carrier plane every slice, so
    # its global per-slice velocity (uniform across nodes, comps 2,3) would dominate
    # T_kinetic.  Subtract it: build the carrier field w = world - vacuum-offset and
    # take its velocity.  (V_total above is unaffected — the offset is a global
    # translation that cancels in the within-slice spatial bonds.)
    off_re, off_im = vacuum_offsets(n_slices, params.r_t)
    w = world.copy()
    # Re offset lives along the trace direction (1,1,1)/sqrt(3) on comps 0,1,2.
    w[:, :, 0:3] -= off_re[:, np.newaxis, np.newaxis] * np.asarray(CARRIER_RE_WEIGHTS)
    w[:, :, CARRIER_IM] -= off_im[:, np.newaxis]

    mass = params.mass(lattice.params)
    T_kinetic = np.zeros(n_slices_plus1)
    for l in range(n_slices_plus1):
        if l == 0:
            vel = (w[1] - w[0]) / dt
        elif l == n_slices:
            vel = (w[-1] - w[-2]) / dt
        else:
            vel = (w[l + 1] - w[l - 1]) / (2.0 * dt)
        T_kinetic[l] = 0.5 * mass * float(np.sum(vel ** 2))

    E_total = V_total + T_kinetic
    slice_indices = np.arange(n_slices_plus1)
    times = slice_indices * dt

    # Conservation check: max deviation from mean E_total
    E_mean = float(np.mean(E_total))
    E_dev = float(np.max(np.abs(E_total - E_mean))) / max(abs(E_mean), 1e-40)

    # Vacuum stability: V_excess at ref should be ~0
    V_ref_check = spacelike_potential(ref, lattice, params) - V_vac
    vacuum_stable = abs(V_ref_check) < 1e-10

    # Save CSV
    csv_path = out_dir / "energy.csv"
    header = "slice,time,V_total,V_vacuum,V_excess,T_kinetic,E_total"
    rows = np.column_stack([
        slice_indices, times, V_total,
        np.full(n_slices_plus1, V_vac), V_excess, T_kinetic, E_total,
    ])
    np.savetxt(str(csv_path), rows, delimiter=",", header=header, comments="")

    # Paper-ready PNG
    _apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("D1 — Energy & Consistency", fontweight="bold")

    ax = axes[0, 0]
    ax.plot(times, V_total, label="V_total")
    ax.plot(times, np.full_like(times, V_vac), "--", color="gray", label="V_vacuum")
    ax.set_xlabel("time"); ax.set_ylabel("V"); ax.legend(); ax.set_title("Potential energy")

    ax = axes[0, 1]
    ax.plot(times, V_excess, color="tab:orange")
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("time"); ax.set_ylabel("V_excess = V - V_vac")
    ax.set_title("Vacuum-subtracted excess energy")

    ax = axes[1, 0]
    ax.plot(times, T_kinetic, color="tab:green", label="T_kinetic")
    ax.plot(times, E_total, color="tab:red", label="E_total=V+T")
    ax.set_xlabel("time"); ax.set_ylabel("Energy"); ax.legend()
    ax.set_title("Kinetic + total energy")

    ax = axes[1, 1]
    E_norm = E_total - E_mean
    ax.plot(times, E_norm, color="tab:purple")
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("time"); ax.set_ylabel("E_total - <E>")
    ax.set_title(f"Energy conservation  (max dev = {E_dev:.2e})")

    _savefig(fig, out_dir / "energy.png")

    return {
        "V_vacuum": float(V_vac),
        "V_excess_mean": float(np.mean(V_excess)),
        "V_excess_max": float(np.max(V_excess)),
        "E_conservation_max_dev": E_dev,
        "vacuum_stable": vacuum_stable,
        "residual_norm_DOF": "N/A — seed only (no BVP solve)",
        "csv": str(csv_path),
        "png": str(out_dir / "energy.png"),
    }


# ---------------------------------------------------------------------------
# D2: Confinement
# ---------------------------------------------------------------------------


def device_confinement(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D2: Confinement metrics (spread_ratio, radius_rms, confined_fraction).

    Uses ``weight_mode="strain"`` (per-node spacelike spring strain energy):
    weighting by bond differences |disp_p - disp_q|^2 makes the metric INVARIANT
    under the prestressed-vacuum N-gon offset (a global per-slice translation
    cancels in the differences).  The ``"lateral"`` |disp|^2 weight would be
    swamped by the spatially-uniform rho-sized vacuum background and report
    spurious box-fill.
    """
    dim = lattice.dim
    result = confinement_summary(
        world, ref, dim,
        confinement_radius_factor=0.5,
        weight_mode="strain",
        lattice=lattice,
        k_s=params.k_s,
        alpha=params.alpha,
    )

    n_slices_plus1 = world.shape[0]
    dt = params.dt
    times = np.arange(n_slices_plus1) * dt

    # CSV
    csv_path = out_dir / "confinement.csv"
    rows = np.column_stack([
        times,
        result["radius_rms"],
        result["spread_ratio"],
        result["confined_fraction"],
    ])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,radius_rms,spread_ratio,confined_fraction", comments="",
    )

    # PNG
    _apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("D2 — Confinement Metrics", fontweight="bold")

    ax = axes[0]
    ax.plot(times, result["radius_rms"])
    ax.axhline(result["box_fill_radius"], color="gray", linestyle="--", label="box_fill_radius")
    ax.set_xlabel("time"); ax.set_ylabel("radius_rms"); ax.legend()
    ax.set_title("Energy-weighted RMS radius")

    ax = axes[1]
    ax.plot(times, result["spread_ratio"])
    ax.axhline(1.0, color="gray", linestyle="--", label="box-fill = 1")
    ax.set_xlabel("time"); ax.set_ylabel("spread_ratio"); ax.legend()
    ax.set_title("Spread ratio  (<<1 = confined, ~1 = dispersed)")

    ax = axes[2]
    ax.plot(times, result["confined_fraction"])
    ax.set_xlabel("time"); ax.set_ylabel("confined_fraction")
    ax.set_title("Confined fraction  (within 0.5 × box_fill_radius)")

    _savefig(fig, out_dir / "confinement.png")

    return {
        "box_fill_radius": float(result["box_fill_radius"]),
        "spread_ratio_mean": float(np.mean(result["spread_ratio"])),
        "confined_fraction_mean": float(np.mean(result["confined_fraction"])),
        "radius_growth": float(result["radius_growth"]),
        "csv": str(csv_path),
        "png": str(out_dir / "confinement.png"),
    }


# ---------------------------------------------------------------------------
# D3: Winding (U(1) phase winding per slice)
# ---------------------------------------------------------------------------


def _winding_per_slice(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    r_t: float = 0.0,
) -> dict[str, np.ndarray]:
    """Compute winding number through each pair of periodic faces per time slice.

    Uses the same discrete plaquette method as measure_winding_closure().
    ``r_t`` is forwarded so the prestressed-vacuum offset is subtracted.
    Returns dict of three arrays shaped (n_slices+1,).
    """
    n_slices_plus1 = world.shape[0]
    wz = np.zeros(n_slices_plus1)
    wy = np.zeros(n_slices_plus1)
    wx = np.zeros(n_slices_plus1)

    for l in range(n_slices_plus1):
        w = measure_winding_closure(world, lattice, slice_index=l, r_t=r_t)
        wz[l] = w["winding_through_z_normal"]
        wy[l] = w["winding_through_y_normal"]
        wx[l] = w["winding_through_x_normal"]

    return {"z": wz, "y": wy, "x": wx}


def device_winding(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
    m_expected: int | None = None,
) -> dict[str, Any]:
    """D3: U(1) phase winding per time slice.

    No hard-coded estimator: uses a discrete plaquette sum through the actual
    field on a contour *enclosing* each central axis (see
    ``measure_winding_closure``).  For the canonical ``Y_l^m`` seed the contour
    about z encloses the vortex axis, so the **target** is

        winding_z ≈ m   (the U(1) charge),   winding_x ≈ winding_y ≈ 0,

    held **constant across every slice** of the loop (the topological charge is
    conserved).  ``m_expected`` (from ``config["vortex_params"]["m"]``) sets the
    z-target; when it is ``None`` the device reports the measured winding without
    a pass/fail.  (Earlier versions targeted net winding = 0 — that was the
    contractible smoke-ring seed, and is wrong for the line-vortex ``Y_l^m``.)
    """
    n_slices_plus1 = world.shape[0]
    dt = params.dt
    times = np.arange(n_slices_plus1) * dt

    windings = _winding_per_slice(world, lattice, r_t=params.r_t)
    wz, wy, wx = windings["z"], windings["y"], windings["x"]

    csv_path = out_dir / "winding.csv"
    rows = np.column_stack([times, wz, wy, wx])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,winding_z_normal,winding_y_normal,winding_x_normal", comments="",
    )

    # Target: winding_z ≈ m, off-axis ≈ 0, all constant across the loop.
    wz_mean = float(np.mean(wz))
    wz_spread = float(np.max(wz) - np.min(wz))          # constancy over the loop
    max_off_axis = max(float(np.max(np.abs(wy))), float(np.max(np.abs(wx))))
    tol = 0.1

    if m_expected is not None:
        z_dev = float(np.max(np.abs(wz - m_expected)))
        z_ok = z_dev < tol
        off_ok = max_off_axis < tol
        const_ok = wz_spread < tol
        winding_ok = bool(z_ok and off_ok and const_ok)
        if winding_ok:
            verdict = (
                f"PASS — winding_z = {wz_mean:+.3f} ≈ m = {m_expected:+d}, "
                f"off-axis max = {max_off_axis:.2e} ≈ 0, slice spread = {wz_spread:.2e} ≈ 0: "
                "a unit U(1) line vortex about z with the charge conserved across the loop."
            )
        else:
            reasons = []
            if not z_ok:
                reasons.append(f"|winding_z - m| = {z_dev:.2e} > {tol}")
            if not off_ok:
                reasons.append(f"off-axis winding {max_off_axis:.2e} > {tol}")
            if not const_ok:
                reasons.append(f"winding varies over the loop (spread {wz_spread:.2e} > {tol})")
            verdict = (
                f"FAIL — winding_z = {wz_mean:+.3f} (target m = {m_expected:+d}): "
                + "; ".join(reasons) + "."
            )
    else:
        winding_ok = None
        verdict = (
            f"winding_z = {wz_mean:+.3f}, off-axis max = {max_off_axis:.2e}, "
            f"slice spread = {wz_spread:.2e}. No m_expected in config — "
            "reporting measured winding without a pass/fail."
        )

    # PNG
    _apply_style()
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(times, wz, label="winding_z (xy-plane)")
    ax.plot(times, wy, label="winding_y (xz-plane)", linestyle="--")
    ax.plot(times, wx, label="winding_x (yz-plane)", linestyle=":")
    ax.axhline(0, color="gray", linewidth=0.8)
    if m_expected is not None:
        ax.axhline(m_expected, color="tab:blue", linewidth=0.9, linestyle="-",
                   alpha=0.35, label=f"target winding_z = m = {m_expected}")
        status = "PASS" if winding_ok else "FAIL"
        subtitle = (f"winding_z={wz_mean:+.2f} (target m={m_expected:+d}), "
                    f"off-axis={max_off_axis:.1e}, spread={wz_spread:.1e}  [{status}]")
    else:
        subtitle = (f"winding_z={wz_mean:+.2f}, off-axis={max_off_axis:.1e}, "
                    f"spread={wz_spread:.1e}  [no m_expected]")
    ax.set_xlabel("time"); ax.set_ylabel("Net winding number")
    ax.legend()
    ax.set_title(
        f"D3 — U(1) Phase Winding (discrete plaquette, actual field)\n{subtitle}"
    )
    fig.suptitle("D3 — Winding", fontweight="bold")

    _savefig(fig, out_dir / "winding.png")

    return {
        "m_expected": m_expected,
        "winding_z_mean": wz_mean,
        "winding_z_slice_spread": wz_spread,
        "max_off_axis_winding": max_off_axis,
        "winding_ok": winding_ok,
        "verdict": verdict,
        "csv": str(csv_path),
        "png": str(out_dir / "winding.png"),
    }


# ---------------------------------------------------------------------------
# D4: Berry / phase envelope
# ---------------------------------------------------------------------------


def _build_complex_envelope(
    world: np.ndarray,
    ref: np.ndarray,
    omega0: float,
    off_re: np.ndarray | None = None,
    off_im: np.ndarray | None = None,
) -> np.ndarray:
    """Build the complex carrier field Psi from the carrier 2-plane.

    Complex displacement A = (Re-off_re) + i(R[3]-ref[3]-off_im), where Re is the
    lateral displacement projected onto the unit trace direction (1,1,1)/sqrt(3)
    (``project_carrier_re``) and ``off_re/off_im`` are the per-slice prestressed-
    vacuum N-gon offsets (so the rho-sized vacuum background is removed and only
    the physical carrier remains).

    Returns
    -------
    Psi : complex array, shape (n_slices+1, n_nodes)
    """
    # Re(Psi) = lateral displacement projected onto the unit trace direction (E1).
    re_disp = project_carrier_re(world[:, :, 0:3] - ref[np.newaxis, :, 0:3])  # (T, N)
    im_disp = world[:, :, CARRIER_IM] - ref[np.newaxis, :, CARRIER_IM]
    if off_re is not None:
        re_disp = re_disp - off_re[:, np.newaxis]
        im_disp = im_disp - off_im[:, np.newaxis]

    # Complex amplitude: A = re + i*im
    A = re_disp + 1j * im_disp  # (T, N)
    return A


def device_berry(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D4: Berry/phase envelope time series.

    Builds Psi = re_disp + i*im_disp from the carrier 2-plane (trace-direction Re
    via project_carrier_re, CARRIER_IM).  Computes per-slice:
      - amplitude |Psi| (node-averaged)
      - U(1) carrier phase arg(Psi) (node-averaged over bright region)
      - Berry connection A_t = Im(<Psi|d_t Psi>) (spatial average)
      - FH discrete link: <Psi(t)|Psi(t+1)> / |...| → accumulated phase

    Mid-plane XY slice phase maps are saved as a multi-panel PNG.
    """
    n_slices_plus1 = world.shape[0]
    n_slices = n_slices_plus1 - 1
    dt = params.dt
    times = np.arange(n_slices_plus1) * dt

    # Prestressed-vacuum N-gon offsets to subtract (else the rho background swamps).
    off_re, off_im = vacuum_offsets(n_slices, params.r_t)
    Psi = _build_complex_envelope(world, ref, params.alpha, off_re, off_im)  # (T, N) complex

    amp = np.abs(Psi)  # (T, N)
    amp_mean = np.mean(amp, axis=1)  # (T,)
    amp_max = np.max(amp, axis=1)

    # U(1) carrier phase: spatial average weighted by amplitude^2
    phase_field = np.angle(Psi)  # (T, N)
    w_sq = amp ** 2 + 1e-40
    phase_weighted_mean = np.sum(phase_field * w_sq, axis=1) / np.sum(w_sq, axis=1)

    # Berry connection: A_t(t) = Im(<Psi(t) | dPsi/dt>)
    # Discretized as Im(<Psi_t | (Psi_{t+1} - Psi_{t-1}) / 2dt>) / <Psi_t|Psi_t>
    A_t = np.zeros(n_slices_plus1)
    for l in range(n_slices_plus1):
        if l == 0:
            dPsi = (Psi[1] - Psi[0]) / dt
        elif l == n_slices:
            dPsi = (Psi[-1] - Psi[-2]) / dt
        else:
            dPsi = (Psi[l + 1] - Psi[l - 1]) / (2.0 * dt)
        norm_sq = float(np.sum(amp[l] ** 2)) + 1e-40
        A_t[l] = float(np.imag(np.sum(np.conj(Psi[l]) * dPsi))) / norm_sq

    # Accumulated Berry phase (FH discrete link)
    berry_phase_accum = np.zeros(n_slices_plus1)
    for l in range(n_slices):
        inner = np.sum(np.conj(Psi[l]) * Psi[l + 1])
        amp_inner = abs(inner)
        if amp_inner > 1e-30:
            berry_phase_accum[l + 1] = berry_phase_accum[l] + np.angle(inner)
        else:
            berry_phase_accum[l + 1] = berry_phase_accum[l]

    # CSV
    csv_path = out_dir / "berry.csv"
    rows = np.column_stack([times, amp_mean, amp_max, phase_weighted_mean, A_t, berry_phase_accum])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,amp_mean,amp_max,phase_mean,berry_connection_At,berry_phase_accum",
        comments="",
    )

    # PNG: 2-row figure
    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape
    _apply_style()

    # --- Row 1: time-series panels ---
    fig = plt.figure(figsize=(16, 10), layout="constrained")
    fig.suptitle("D4 — Berry / Phase Envelope", fontweight="bold")

    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(times, amp_mean, label="|Psi| mean")
    ax.plot(times, amp_max, label="|Psi| max", linestyle="--")
    ax.set_xlabel("time"); ax.set_ylabel("|Psi|"); ax.legend()
    ax.set_title("Carrier amplitude")

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(times, phase_weighted_mean)
    ax.set_xlabel("time"); ax.set_ylabel("arg(Psi) [rad]")
    ax.set_title("Amplitude-weighted mean carrier phase")

    ax = fig.add_subplot(gs[0, 2])
    ax.plot(times, berry_phase_accum)
    ax.set_xlabel("time"); ax.set_ylabel("Accumulated Berry phase [rad]")
    ax.set_title("FH discrete Berry phase")

    # --- Row 2: mid-plane phase maps at t=0, t_mid, t_end ---
    t_indices = [0, n_slices // 2, n_slices]
    for col, ti in enumerate(t_indices):
        pos = world[ti]
        re_d = project_carrier_re(pos[:, 0:3] - ref[:, 0:3]) - off_re[ti]
        im_d = pos[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[ti]
        amp_3d = np.sqrt(re_d ** 2 + im_d ** 2).reshape(nx, ny, nz)
        ph_3d = np.arctan2(im_d, re_d + 1e-300).reshape(nx, ny, nz)

        # XY mid-plane
        sl_amp = amp_3d[:, :, nz // 2]
        sl_ph = ph_3d[:, :, nz // 2]
        amp_max_sl = float(np.max(sl_amp)) or 1.0

        rgb = _phase_to_rgb(sl_ph)
        alpha_ch = np.clip(sl_amp / amp_max_sl, 0, 1)[..., None]
        rgba = np.concatenate([rgb, alpha_ch], axis=-1)

        ax = fig.add_subplot(gs[1, col])
        ax.imshow(
            rgba.swapaxes(0, 1), origin="lower",
            extent=[0, nx, 0, ny], aspect="equal",
        )
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(f"XY mid-plane phase, t={times[ti]:.2f}")

    _savefig(fig, out_dir / "berry.png")

    return {
        "amp_mean_t0": float(amp_mean[0]),
        "amp_max_t0": float(amp_max[0]),
        "berry_phase_accum_final": float(berry_phase_accum[-1]),
        "berry_connection_At_mean": float(np.mean(A_t)),
        "csv": str(csv_path),
        "png": str(out_dir / "berry.png"),
    }


# ---------------------------------------------------------------------------
# D5: EM fields (A_mu, F_mu_nu)
# ---------------------------------------------------------------------------


def device_em_fields(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D5: EM field quiver on mid-plane slice.

    A_mu = i<u|d_mu u>  (Berry connection in Psi space, CARRIER 2-plane)
    F_munu = d_mu A_nu - d_nu A_mu  (field strength)

    For the seed (single time slice t=0):
      - Compute A_x, A_y, A_z as spatial Berry connections from neighbouring nodes
      - Compute A_t from central-time difference (as in D4)
      - Render E_i = F_0i and B_i = (1/2)eps_ijk F_jk on XY mid-plane

    Note: on a seed worldvolume with smooth slowly-varying phase, these fields
    are dominated by the vortex topology.
    """
    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape
    a = params.dt  # spatial spacing same as lattice.params.spacing

    # Prestressed-vacuum N-gon offsets (subtract to recover the physical carrier).
    off_re, off_im = vacuum_offsets(world.shape[0] - 1, params.r_t)

    # Build complex field at t=0 on the 3D grid
    pos0 = world[0]
    re_d = project_carrier_re(pos0[:, 0:3] - ref[:, 0:3]) - off_re[0]
    im_d = pos0[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[0]
    Psi = (re_d + 1j * im_d).reshape(nx, ny, nz)

    # Normalize node-wise: psi_hat = Psi / |Psi| (or 0 where |Psi| < eps)
    amp3 = np.abs(Psi)
    amp_max = float(np.max(amp3))
    eps = 1e-6 * amp_max if amp_max > 0 else 1e-30
    psi_hat = np.where(amp3 > eps, Psi / (amp3 + 1e-300), 0.0 + 0.0j)

    # Spatial Berry connection A_i(r) = Im(<psi_hat(r) | d_i psi_hat(r)>)
    # Discretized: Im(<psi_hat(r)^* | psi_hat(r + a*e_i) - psi_hat(r) >) / a
    def _A_spatial(arr: np.ndarray, axis: int) -> np.ndarray:
        shifted = np.roll(arr, -1, axis=axis)
        conn = np.imag(arr.conj() * (shifted - arr))
        return conn / lattice.params.spacing

    A_x = _A_spatial(psi_hat, 0)
    A_y = _A_spatial(psi_hat, 1)
    A_z = _A_spatial(psi_hat, 2)

    # Temporal: use slices 0 and 1
    if world.shape[0] > 1:
        pos1 = world[1]
        re_d1 = project_carrier_re(pos1[:, 0:3] - ref[:, 0:3]) - off_re[1]
        im_d1 = pos1[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[1]
        Psi1 = (re_d1 + 1j * im_d1).reshape(nx, ny, nz)
        amp3_1 = np.abs(Psi1)
        psi_hat1 = np.where(amp3_1 > eps, Psi1 / (amp3_1 + 1e-300), 0.0 + 0.0j)
        A_t = np.imag(psi_hat.conj() * (psi_hat1 - psi_hat)) / params.dt
    else:
        A_t = np.zeros_like(A_x)

    # Field strength tensor (antisymmetric) F_mu_nu = d_mu A_nu - d_nu A_mu
    # E_i = F_0i = d_t A_i - d_i A_t   (approximated from t=0 slice only)
    # B_z = F_xy = d_x A_y - d_y A_x
    # B_y = F_zx = d_z A_x - d_x A_z
    # B_x = F_yz = d_y A_z - d_z A_y

    def _d(arr, axis):
        return (np.roll(arr, -1, axis=axis) - arr) / lattice.params.spacing

    B_z = _d(A_y, 0) - _d(A_x, 1)
    B_y = _d(A_x, 2) - _d(A_z, 0)
    B_x = _d(A_z, 1) - _d(A_y, 2)

    # E fields: need d_t A_i — use finite difference from slice 0 to 1
    # and d_i A_t
    # For seed, A_t is approximately zero so E_i ~ -d_i A_t (negligible)
    # We report E_i = -d_i phi where phi = A_t as the electrostatic piece
    E_x = -_d(A_t, 0)
    E_y = -_d(A_t, 1)
    E_z = -_d(A_t, 2)

    # Mid-plane slice: z = nz//2
    z_mid = nz // 2
    _apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("D5 — EM Fields  (A_mu Berry connection, F_munu field strength)\nt=0 mid-plane z=" + str(z_mid), fontweight="bold")

    extent = [0, nx * lattice.params.spacing, 0, ny * lattice.params.spacing]

    def _quiver_slice(ax, U, V, C, title, xlabel="x", ylabel="y"):
        xs = np.arange(nx) * lattice.params.spacing
        ys = np.arange(ny) * lattice.params.spacing
        Xg, Yg = np.meshgrid(xs, ys, indexing="ij")
        stride = max(1, nx // 16)
        Xs = Xg[::stride, ::stride]
        Ys = Yg[::stride, ::stride]
        Us = U[::stride, ::stride, z_mid]
        Vs = V[::stride, ::stride, z_mid]
        Cs = C[:, :, z_mid]
        vmax = float(np.max(np.abs(Cs))) or 1.0
        im = ax.imshow(Cs.T, origin="lower", extent=extent, aspect="equal",
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        scale = float(np.max(np.sqrt(Us ** 2 + Vs ** 2))) or 1.0
        ax.quiver(Xs, Ys, Us / scale, Vs / scale, alpha=0.6, scale=25, color="k", width=0.003)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)

    _quiver_slice(axes[0, 0], A_x, A_y, A_z, "A_z (color) + (A_x,A_y) quiver")
    _quiver_slice(axes[0, 1], E_x, E_y, E_z, "E_z (color) + (E_x,E_y) quiver")
    _quiver_slice(axes[0, 2], B_x, B_y, B_z, "B_z (color) + (B_x,B_y) quiver")

    # Row 2: amplitude overlay on XZ slice (ring cross-section)
    ax = axes[1, 0]
    amp_xz = amp3[:, ny // 2, :]
    ax.imshow(amp_xz.T, origin="lower", extent=[0, nx * lattice.params.spacing, 0, nz * lattice.params.spacing],
              aspect="equal", cmap="inferno")
    ax.set_xlabel("x"); ax.set_ylabel("z"); ax.set_title("Amplitude |Psi| XZ midplane (ring cross-section)")

    ax = axes[1, 1]
    amp_xy = amp3[:, :, nz // 2]
    ax.imshow(amp_xy.T, origin="lower", extent=extent, aspect="equal", cmap="inferno")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("Amplitude |Psi| XY midplane (ring plan view)")

    ax = axes[1, 2]
    phase_xy = np.angle(Psi[:, :, nz // 2])
    rgb_xy = _phase_to_rgb(phase_xy)
    alpha_xy = np.clip(amp_xy / (float(np.max(amp_xy)) or 1.0), 0, 1)[..., None]
    rgba_xy = np.concatenate([rgb_xy, alpha_xy], axis=-1)
    ax.imshow(rgba_xy.swapaxes(0, 1), origin="lower", extent=extent, aspect="equal")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("Phase → RGB (XY midplane)")

    _savefig(fig, out_dir / "em_fields.png")

    return {
        "A_x_max": float(np.max(np.abs(A_x))),
        "A_y_max": float(np.max(np.abs(A_y))),
        "A_z_max": float(np.max(np.abs(A_z))),
        "B_z_max_midplane": float(np.max(np.abs(B_z[:, :, z_mid]))),
        "png": str(out_dir / "em_fields.png"),
    }


# ---------------------------------------------------------------------------
# D6: Per-color-channel / QCD split
# ---------------------------------------------------------------------------


def device_color_channels(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D6: U(1) trace + SU(3) traceless split of the lateral triplet.

    The lateral triplet = displacements in components 0,1,2 (the spacelike/
    colour channels).  The vortex seed writes the real carrier along the trace
    direction (1,1,1)/sqrt(3) of the triplet (E1 fix), so the bare seed is a
    PURE U(1)/EM excitation: P_SU3 annihilates it and ``u1_fraction -> 1``.
    Any SU(3)/colour content is therefore genuinely coexcited by the relaxation,
    not injected by construction.

    Projects each node's lateral displacement (d0, d1, d2) with P_U1 and
    P_SU3 from alpha_separability.projection_operators().

    Per-channel amplitude = |P @ disp| and per-channel energy = |P @ disp|^2.
    Reports:
      - U(1) channel amplitude map (XY midplane)
      - SU(3) 8-channel breakdown (summed shear energy)
      - Per-channel time series
    """
    n_slices_plus1 = world.shape[0]
    dt = params.dt
    times = np.arange(n_slices_plus1) * dt

    P_U1, P_SU3 = projection_operators()  # (3,3) each

    # The prestressed-vacuum N-gon offset lives along the trace direction
    # (1,1,1)/sqrt(3) of the lateral triplet; subtract it so the trace/traceless
    # split reads the carrier, not the spatially-uniform vacuum background.
    # (Component 3 is timelike, not in the lateral triplet, so off_im is absent.)
    off_re, _off_im = vacuum_offsets(n_slices_plus1 - 1, params.r_t)
    _trace_w = np.asarray(CARRIER_RE_WEIGHTS)

    # Per-slice: compute lateral displacement (components 0,1,2) and project
    u1_energy = np.zeros(n_slices_plus1)
    su3_energy = np.zeros(n_slices_plus1)
    u1_amp_mean = np.zeros(n_slices_plus1)
    su3_amp_mean = np.zeros(n_slices_plus1)

    for l in range(n_slices_plus1):
        disp_lat = world[l, :, :3] - ref[np.newaxis, :, :3].reshape(-1, 3)  # (N, 3)
        disp_lat -= off_re[l] * _trace_w  # remove vacuum offset (trace direction)
        d_u1 = disp_lat @ P_U1.T   # (N, 3)  — trace/dilatational component
        d_su3 = disp_lat @ P_SU3.T  # (N, 3) — shear/traceless component

        u1_energy[l] = float(np.sum(d_u1 ** 2))
        su3_energy[l] = float(np.sum(d_su3 ** 2))
        u1_amp_mean[l] = float(np.mean(np.linalg.norm(d_u1, axis=1)))
        su3_amp_mean[l] = float(np.mean(np.linalg.norm(d_su3, axis=1)))

    total_lat_energy = u1_energy + su3_energy
    u1_fraction = u1_energy / (total_lat_energy + 1e-40)
    su3_fraction = su3_energy / (total_lat_energy + 1e-40)

    # CSV
    csv_path = out_dir / "color_channels.csv"
    rows = np.column_stack([times, u1_energy, su3_energy, u1_fraction, su3_fraction,
                            u1_amp_mean, su3_amp_mean])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,u1_energy,su3_energy,u1_fraction,su3_fraction,u1_amp_mean,su3_amp_mean",
        comments="",
    )

    # PNG: time series + spatial maps at t=0
    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape

    disp_lat0 = world[0, :, :3] - ref[:, :3]  # (N, 3)
    disp_lat0 = disp_lat0 - off_re[0] * _trace_w  # remove vacuum pedestal (trace)
    d_u1_0 = (disp_lat0 @ P_U1.T).reshape(nx, ny, nz, 3)
    d_su3_0 = (disp_lat0 @ P_SU3.T).reshape(nx, ny, nz, 3)

    amp_u1_map = np.linalg.norm(d_u1_0, axis=-1)     # (nx, ny, nz)
    amp_su3_map = np.linalg.norm(d_su3_0, axis=-1)   # (nx, ny, nz)

    _apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("D6 — Per-Color-Channel Split  (U(1) trace + SU(3) traceless)", fontweight="bold")

    extent = [0, nx, 0, ny]
    z_mid = nz // 2

    ax = axes[0, 0]
    ax.plot(times, u1_energy, label="U(1) trace", color="tab:blue")
    ax.plot(times, su3_energy, label="SU(3) shear", color="tab:orange")
    ax.set_xlabel("time"); ax.set_ylabel("Energy proxy"); ax.legend()
    ax.set_title("Per-channel lateral energy")

    ax = axes[0, 1]
    ax.plot(times, u1_fraction, label="U(1) fraction", color="tab:blue")
    ax.plot(times, su3_fraction, label="SU(3) fraction", color="tab:orange")
    ax.set_xlabel("time"); ax.set_ylabel("Fraction"); ax.legend()
    ax.set_title("Channel energy fraction (U(1) + SU(3) = 1)")

    ax = axes[0, 2]
    ax.plot(times, u1_amp_mean, label="|P_U1 disp| mean", color="tab:blue")
    ax.plot(times, su3_amp_mean, label="|P_SU3 disp| mean", color="tab:orange")
    ax.set_xlabel("time"); ax.set_ylabel("Displacement amplitude"); ax.legend()
    ax.set_title("Per-channel displacement amplitude")

    ax = axes[1, 0]
    vmax = float(np.max(amp_u1_map)) or 1.0
    im = ax.imshow(amp_u1_map[:, :, z_mid].T, origin="lower", extent=extent, aspect="equal",
                   cmap="Blues", vmin=0, vmax=vmax)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("U(1) channel amplitude  XY midplane (t=0)")

    ax = axes[1, 1]
    vmax_su3 = float(np.max(amp_su3_map)) or 1.0
    im = ax.imshow(amp_su3_map[:, :, z_mid].T, origin="lower", extent=extent, aspect="equal",
                   cmap="Oranges", vmin=0, vmax=vmax_su3)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("SU(3) channel amplitude  XY midplane (t=0)")

    # Per-component breakdown (d0, d1, d2)
    ax = axes[1, 2]
    for comp_idx, color, label in zip(
        range(3),
        ["tab:blue", "tab:orange", "tab:green"],
        ["comp-0 (trace 1/3)", "comp-1 (trace 1/3)", "comp-2 (trace 1/3)"],
    ):
        amp_comp = np.zeros(n_slices_plus1)
        for l in range(n_slices_plus1):
            amp_comp[l] = float(np.mean(np.abs(world[l, :, comp_idx] - ref[:, comp_idx])))
        ax.plot(times, amp_comp, label=label, color=color)
    ax.set_xlabel("time"); ax.set_ylabel("Mean |disp_i|")
    ax.legend(fontsize=8); ax.set_title("Per-lateral-component displacement")

    _savefig(fig, out_dir / "color_channels.png")

    return {
        "u1_fraction_mean": float(np.mean(u1_fraction)),
        "su3_fraction_mean": float(np.mean(su3_fraction)),
        "u1_energy_t0": float(u1_energy[0]),
        "su3_energy_t0": float(su3_energy[0]),
        "note": (
             "The bare seed writes the carrier along the trace direction "
            "(1,1,1)/sqrt(3) of components (0,1,2) plus the timelike component 3 "
            "(Im). P_SU3 annihilates a pure-trace displacement, so the bare seed "
            "reads u1_fraction -> 1 and su3_fraction -> 0; any nonzero SU(3) "
            "content is genuinely coexcited by the relaxation (E1 fix)."
        ),
        "csv": str(csv_path),
        "png": str(out_dir / "color_channels.png"),
    }


# ---------------------------------------------------------------------------
# D7: Dispersion / spectra
# ---------------------------------------------------------------------------


def device_spectra(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D7: Spatial FFT energy spectrum at t=0.

    Computes the 3D FFT of the complex carrier amplitude Psi on the t=0 slice.
    Bins the power spectrum by |k| and reports which k-shells carry energy.
    Also marks the analytic dispersion omega(k) from d_of_k_eigenvalues as
    reference.
    """
    from branesim.core.conventions import d_of_k_eigenvalues

    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape
    a = lattice.params.spacing

    off_re, off_im = vacuum_offsets(world.shape[0] - 1, params.r_t)
    pos0 = world[0]
    re_d = project_carrier_re(pos0[:, 0:3] - ref[:, 0:3]) - off_re[0]
    im_d = pos0[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[0]
    Psi0 = (re_d + 1j * im_d).reshape(nx, ny, nz)

    # 3D FFT
    Psi_k = np.fft.fftn(Psi0)
    power = np.abs(Psi_k) ** 2

    # Wavevectors
    kx = np.fft.fftfreq(nx, d=a) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=a) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=a) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    K_mag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)

    # Radially bin
    k_max = float(np.max(K_mag))
    n_bins = 32
    k_edges = np.linspace(0, k_max, n_bins + 1)
    k_centres = 0.5 * (k_edges[:-1] + k_edges[1:])
    power_radial = np.zeros(n_bins)
    count_radial = np.zeros(n_bins, dtype=int)

    k_flat = K_mag.ravel()
    p_flat = power.ravel()
    for i in range(n_bins):
        mask = (k_flat >= k_edges[i]) & (k_flat < k_edges[i + 1])
        power_radial[i] = float(np.sum(p_flat[mask]))
        count_radial[i] = int(np.sum(mask))

    # Analytic dispersion reference (isotropic average: omega_L along [100])
    alpha = params.alpha
    omega_ref = np.array([
        float(np.sqrt(d_of_k_eigenvalues(
            np.array([kc, 0.0, 0.0]), alpha, k_s=1.0, rho=1.0, a=a
        )[0]))
        for kc in k_centres
    ])

    # CSV
    csv_path = out_dir / "spectra.csv"
    rows = np.column_stack([k_centres, power_radial, count_radial, omega_ref])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="k_mag,power_radial,count,omega_analytic_L", comments="",
    )

    # PNG
    _apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("D7 — Spatial FFT Energy Spectrum (t=0)", fontweight="bold")

    ax = axes[0]
    ax.semilogy(k_centres, power_radial + 1e-30)
    ax.set_xlabel("|k| (rad/a)"); ax.set_ylabel("Radial power spectrum")
    ax.set_title("Carrier energy vs k-shell  (vortex ring structure)")

    ax2 = ax.twinx()
    ax2.plot(k_centres, omega_ref, "r--", alpha=0.7, label="omega_L analytic [100]")
    ax2.set_ylabel("omega analytic [rad/step]", color="r")
    ax2.legend(loc="upper right")

    ax = axes[1]
    # 2D slice of the power spectrum (kz=0 plane)
    power_2d = np.abs(np.fft.fftn(Psi0)[:, :, 0]) ** 2
    kx_plot = np.fft.fftshift(kx)
    ky_plot = np.fft.fftshift(ky)
    power_2d_shift = np.fft.fftshift(power_2d)
    vmax = float(np.percentile(power_2d_shift, 99)) or 1.0
    ax.imshow(
        np.log1p(power_2d_shift).T, origin="lower",
        extent=[kx_plot[0], kx_plot[-1], ky_plot[0], ky_plot[-1]],
        aspect="equal", cmap="inferno",
    )
    ax.set_xlabel("kx (rad/a)"); ax.set_ylabel("ky (rad/a)")
    ax.set_title("log(1+power) 2D spectrum  kz=0 plane")

    _savefig(fig, out_dir / "spectra.png")

    # Peak k
    peak_bin = int(np.argmax(power_radial))
    return {
        "peak_k_mag": float(k_centres[peak_bin]),
        "peak_power": float(power_radial[peak_bin]),
        "total_power": float(np.sum(power_radial)),
        "csv": str(csv_path),
        "png": str(out_dir / "spectra.png"),
    }


# ---------------------------------------------------------------------------
# Report assembler
# ---------------------------------------------------------------------------


def _fmt(value: Any, spec: str = ".4g") -> str:
    """Format ``value`` with ``spec`` if it is a real number, else pass it through.

    Guards report rendering against absent devices: when ``run_measurements`` is
    invoked on a device subset (e.g. ``devices=["winding"]`` to re-check one
    device on an existing run folder), the other devices' metrics are missing and
    default to the string ``"N/A"``.  Feeding ``"N/A"`` to a numeric format spec
    (``:.4g`` / ``:.2e``) raised ``ValueError`` and crashed ``_write_report``
    (E13).  Booleans pass through as-is (they are formatted as text, not numbers).
    """
    if value is None:
        return "N/A"
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return str(value)
    return format(value, spec)


def _write_report(
    run_dir: Path,
    diag_dir: Path,
    config: dict,
    results: dict[str, dict],
    ledger_md: str | None = None,
) -> Path:
    """Write diagnostics/report.md with key figures and verdicts."""
    ts = config.get("timestamp", "unknown")
    alpha = config.get("alpha", "?")
    grid = config.get("grid_shape", "?")

    lines = [
        f"# BraneSim Diagnostics Report",
        f"",
        f"**Run:** `{run_dir.name}`  ",
        f"**Timestamp:** {ts}  ",
        f"**Grid:** {grid}  alpha={alpha}",
        f"",
        f"---",
        f"",
    ]

    # Trust ledger first — what is trustworthy on this run, read bottom-up.
    if ledger_md:
        lines += [ledger_md]

    # D1
    d1 = results.get("energy", {})
    lines += [
        f"## D1 — Energy & Consistency",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| V_vacuum | {_fmt(d1.get('V_vacuum'), '.4g')} |",
        f"| V_excess mean | {_fmt(d1.get('V_excess_mean'), '.4g')} |",
        f"| V_excess max | {_fmt(d1.get('V_excess_max'), '.4g')} |",
        f"| E_conservation max dev | {_fmt(d1.get('E_conservation_max_dev'), '.2e')} |",
        f"| Vacuum stable | {d1.get('vacuum_stable', 'N/A')} |",
        f"| Residual norm/DOF | {d1.get('residual_norm_DOF', 'N/A')} |",
        f"",
        f"![Energy plot](energy.png)",
        f"",
    ]

    # D2
    d2 = results.get("confinement", {})
    lines += [
        f"## D2 — Confinement",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| box_fill_radius | {_fmt(d2.get('box_fill_radius'), '.4g')} |",
        f"| spread_ratio mean | {_fmt(d2.get('spread_ratio_mean'), '.4g')} |",
        f"| confined_fraction mean | {_fmt(d2.get('confined_fraction_mean'), '.4g')} |",
        f"| radius_growth | {_fmt(d2.get('radius_growth'), '.4g')} |",
        f"",
        f"**Verdict:** spread_ratio << 1 confirms localized ring seed. "
        f"radius_growth ~ 1 expected (seed, no dynamics evolved).",
        f"",
        f"![Confinement plot](confinement.png)",
        f"",
    ]

    # D3
    d3 = results.get("winding", {})
    _m_exp = d3.get("m_expected")
    _wok = d3.get("winding_ok")
    _wok_str = {True: "PASS", False: "FAIL", None: "n/a (no m_expected)"}.get(_wok, "?")
    lines += [
        f"## D3 — Phase Winding",
        f"",
        f"| Metric | Value | Target |",
        f"|--------|-------|--------|",
        f"| winding_z mean | {d3.get('winding_z_mean', float('nan')):.4g} | m = {_m_exp if _m_exp is not None else '?'} |",
        f"| winding_z slice spread | {d3.get('winding_z_slice_spread', float('nan')):.2e} | ~0 (const over loop) |",
        f"| max off-axis winding | {d3.get('max_off_axis_winding', float('nan')):.2e} | ~0 |",
        f"| winding check | {_wok_str} | PASS |",
        f"",
        f"**Verdict:** {d3.get('verdict', '')}",
        f"",
        f"![Winding plot](winding.png)",
        f"",
    ]

    # D4
    d4 = results.get("berry", {})
    lines += [
        f"## D4 — Berry / Phase Envelope",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| amp_mean (t=0) | {_fmt(d4.get('amp_mean_t0'), '.4g')} |",
        f"| amp_max (t=0) | {_fmt(d4.get('amp_max_t0'), '.4g')} |",
        f"| Berry phase accumulated | {_fmt(d4.get('berry_phase_accum_final'), '.4g')} rad |",
        f"| Berry connection <A_t> | {_fmt(d4.get('berry_connection_At_mean'), '.4g')} |",
        f"",
        f"![Berry plot](berry.png)",
        f"",
    ]

    # D5
    d5 = results.get("em_fields", {})
    lines += [
        f"## D5 — EM Fields",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| max |A_x| | {_fmt(d5.get('A_x_max'), '.4g')} |",
        f"| max |A_y| | {_fmt(d5.get('A_y_max'), '.4g')} |",
        f"| max |A_z| | {_fmt(d5.get('A_z_max'), '.4g')} |",
        f"| max |B_z| midplane | {_fmt(d5.get('B_z_max_midplane'), '.4g')} |",
        f"",
        f"![EM fields plot](em_fields.png)",
        f"",
    ]

    # D6
    d6 = results.get("color_channels", {})
    lines += [
        f"## D6 — Color / QCD Channel Split",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| U(1) energy fraction | {_fmt(d6.get('u1_fraction_mean'), '.4g')} |",
        f"| SU(3) energy fraction | {_fmt(d6.get('su3_fraction_mean'), '.4g')} |",
        f"| U(1) energy (t=0) | {_fmt(d6.get('u1_energy_t0'), '.4g')} |",
        f"| SU(3) energy (t=0) | {_fmt(d6.get('su3_energy_t0'), '.4g')} |",
        f"",
        f"**Note:** {d6.get('note', '')}",
        f"",
        f"![Color channels plot](color_channels.png)",
        f"",
    ]

    # D7
    d7 = results.get("spectra", {})
    lines += [
        f"## D7 — Dispersion / Spectra",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| peak k | {_fmt(d7.get('peak_k_mag'), '.4g')} rad/a |",
        f"| peak power | {_fmt(d7.get('peak_power'), '.4g')} |",
        f"| total power | {_fmt(d7.get('total_power'), '.4g')} |",
        f"",
        f"![Spectra plot](spectra.png)",
        f"",
    ]

    # D9 — ring azimuthal asymmetry
    d9 = results.get("ring_azimuth", {})
    lines += [
        f"## D9 — Ring Azimuthal Asymmetry",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| asymmetry mean (Σ_m≥1 P_m / ΣP) | {_fmt(d9.get('asymmetry_mean'), '.4g')} |",
        f"| asymmetry max | {_fmt(d9.get('asymmetry_max'), '.4g')} |",
        f"| m=1 dipole fraction mean | {_fmt(d9.get('m1_dipole_frac_mean'), '.4g')} |",
        f"| core centroid offset mean | {_fmt(d9.get('centroid_offset_mean'), '.4g')} a |",
        f"| ring radius r_peak mean | {_fmt(d9.get('r_peak_mean'), '.4g')} a |",
        f"",
        f"**Verdict:** symmetric ring → asymmetry ≈ 0 (all power in m=0). Nonzero "
        f"m≥1 power / centroid offset = the donut deformed (one side brighter) — a "
        f"progress tracker of the relaxation, not a pass/fail.",
        f"",
        f"![Ring azimuth plot](ring_azimuth.png)",
        f"",
    ]

    # D8
    d8 = results.get("binding_probe", {})
    _d8_cl = d8.get("closure_locked")
    _d8_sign = d8.get("du_par_sign", 0)
    _d8_sign_str = {-1: "NEGATIVE (binding-consistent)", 0: "ZERO / undefined", 1: "POSITIVE (repulsive)"}.get(_d8_sign, "?")
    _d8_kahler_flag = "PRESENT (R > 0.1)" if (d8.get("R_kahler", 0.0) or 0.0) > 0.1 else "absent (R < 0.1)"
    lines += [
        f"## D8 — U(1)↔SU(3) Binding Probe",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| sep_d mean (P1) | {_fmt(d8.get('sep_d_mean'), '.4g')} lattice units |",
        f"| sep_d max (P1) | {_fmt(d8.get('sep_d_max'), '.4g')} |",
        f"| Δu_∥ sign (P2) | {_d8_sign_str} |",
        f"| Δu_∥ extremum (P2) | {_fmt(d8.get('du_par_extremum'), '.4g')} |",
        f"| ρ of extremum (P2) | {_fmt(d8.get('rho_extremum'), '.4g')} |",
        f"| ⟨Δu_∥⟩ scalar (P2) | {_fmt(d8.get('du_par_scalar'), '.4g')} |",
        f"| ω mean (P3) | {_fmt(d8.get('omega_mean'), '.4g')} rad/time_unit |",
        f"| ω std (P3) | {_fmt(d8.get('omega_std'), '.4g')} |",
        f"| ω_ref closure-locked (P3) | {d8.get('omega_ref', 'N/A')} |",
        f"| closure-locked verdict (P3) | {_d8_cl} |",
        f"| γ_Γ loop holonomy (P4) | {_fmt(d8.get('gamma_Gamma'), '.4g')} rad |",
        f"| R_kahler (P5) | {_fmt(d8.get('R_kahler'), '.4g')} |",
        f"| Kähler part (P5) | {_d8_kahler_flag} |",
        f"",
        f"**One-line verdict:** "
        f"sep_d={d8.get('sep_d_mean', float('nan')):.3g} (U(1)/SU(3) co-located if < 1); "
        f"Δu_∥ {_d8_sign_str}; "
        f"ω {'closure-locked' if _d8_cl else ('NOT locked' if _d8_cl is False else 'untested')}; "
        f"γ_Γ={d8.get('gamma_Gamma', float('nan')):.3g} rad "
        f"(expected ≈ 2π·n_t={round(d8.get('gamma_Gamma', 0)/(2*3.14159), 1) if d8.get('gamma_Gamma') else '?'} turns); "
        f"Kähler R={d8.get('R_kahler', float('nan')):.3g} — {_d8_kahler_flag}.",
        f"",
        f"![Binding probe plot](binding_probe.png)",
        f"",
        f"---",
        f"",
        f"## Verdict Summary",
        f"",
        f"- **D1 Energy**: V_excess localized to ring; vacuum_stable={'OK' if d1.get('vacuum_stable') else 'FAIL'}. "
        f"Residual: {d1.get('residual_norm_DOF', 'N/A')} (seed; full check requires BVP solve).",
        f"- **D2 Confinement**: spread_ratio ~ {_fmt(d2.get('spread_ratio_mean', '?'), '.3g')}; "
        f"confined_fraction ~ {_fmt(d2.get('confined_fraction_mean', '?'), '.3g')}.",
        f"- **D3 Winding**: winding_z = {d3.get('winding_z_mean', float('nan')):.3g} "
        f"(target m={_m_exp if _m_exp is not None else '?'}), "
        f"off-axis ~ {d3.get('max_off_axis_winding', float('nan')):.2g} — check {_wok_str}.",
        f"- **D4 Berry**: carrier amplitude peak {_fmt(d4.get('amp_max_t0', '?'), '.3g')} ≈ 0.607 × A0 (donut profile). "
        f"Berry phase accumulated over all slices.",
        f"- **D5 EM**: Berry connection A_mu computed from phase gradient; B field shows vortex topology.",
        f"- **D6 Color**: pure U(1) trace channel for bare seed "
        f"(U(1) frac ~ {_fmt(d6.get('u1_fraction_mean', '?'), '.3g')}); "
        f"SU(3) content {'appears' if d6.get('su3_fraction_mean', 0) > 0.01 else 'absent'} for this seed.",
        f"- **D7 Spectra**: ring geometry imprints k-ring structure on FFT.",
        f"- **D8 Binding Probe**: sep_d ~ {_fmt(d8.get('sep_d_mean', '?'), '.3g')}; "
        f"Δu_∥ {_d8_sign_str}; "
        f"γ_Γ ~ {_fmt(d8.get('gamma_Gamma', '?'), '.3g')} rad; "
        f"Kähler {_d8_kahler_flag}.",
        f"",
        f"*Generated by `branesim.diagnostics.run_measurements`.*",
    ]

    report_path = diag_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# D9: Ring azimuthal-asymmetry (Fourier decomposition of |Psi| around the ring)
# ---------------------------------------------------------------------------


def device_ring_azimuth(
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
) -> dict[str, Any]:
    """D9: azimuthal Fourier decomposition of |Psi| around the donut ring.

    For the bare ``Y_1^1`` seed, ``|Psi| ~ |sin(theta) e^{i phi}| = sin(theta)`` is
    **azimuthally uniform** in the z-midplane (all power in m=0).  When the
    relaxation deforms the ring — the "one side brighter / other side dimmer"
    effect — the broken symmetry shows up as **m>=1 azimuthal harmonics** (m=1 =
    the dipole / off-axis core drift, m=2 = ellipticity, …).

    Per time slice, in the z-midplane:
      1. find the donut radius ``r_peak`` (peak of the azimuthally-averaged radial
         profile of |Psi|),
      2. average |Psi| into azimuthal bins over a radial band around ``r_peak`` →
         ``A(phi)``,
      3. ``rfft(A(phi))`` → harmonic power ``P_m``; report fractions ``P_m/Σ P``,
         the asymmetry ``Σ_{m>=1} P_m / Σ_{m>=0} P_m``, and the energy-weighted
         core-centroid offset from the box axis.

    A symmetric ring reads ``asymmetry ≈ 0`` and ``centroid_offset ≈ 0``; growth of
    either over the loop / across relaxation iterations quantifies the deformation.
    Read-only; no pass/fail (it is a progress tracker, not a contract).
    """
    nx, ny, nz = lattice.params.grid_shape
    n_T = world.shape[0]
    times = np.arange(n_T) * params.dt
    off_re, off_im = vacuum_offsets(n_T - 1, params.r_t)

    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    zc = nz // 2
    xs = np.arange(nx)[:, None].astype(float)
    ys = np.arange(ny)[None, :].astype(float)
    rr = np.sqrt((xs - cx) ** 2 + (ys - cy) ** 2)         # (nx, ny)
    phi = np.arctan2(ys - cy, xs - cx)                    # (nx, ny)
    n_bins = 64
    bin_idx = (np.floor((phi + np.pi) / (2 * np.pi) * n_bins).astype(int) % n_bins).ravel()
    r_int = rr.astype(int).ravel()
    rmax = int(rr.max())
    M = 5  # report harmonics m = 0..4

    def _amp_zmid(l: int) -> np.ndarray:
        pos = world[l]
        re = (project_carrier_re(pos[:, 0:3] - ref[:, 0:3]) - off_re[l]).reshape(nx, ny, nz)
        im = (pos[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[l]).reshape(nx, ny, nz)
        return np.sqrt(re[:, :, zc] ** 2 + im[:, :, zc] ** 2)  # (nx, ny)

    asym = np.zeros(n_T)
    centroid = np.zeros(n_T)
    pm_frac = np.zeros((n_T, M))
    r_peaks = np.zeros(n_T)
    A_phi_0 = np.zeros(n_bins)

    for l in range(n_T):
        amp = _amp_zmid(l)
        flat = amp.ravel()
        # radial profile -> donut radius
        rsum = np.bincount(r_int, weights=flat, minlength=rmax + 1)
        rcnt = np.maximum(np.bincount(r_int, minlength=rmax + 1), 1)
        prof = rsum / rcnt
        r_peak = float(np.argmax(prof))
        r_peaks[l] = r_peak
        band = max(2.0, 0.3 * max(r_peak, 1.0))
        sel = np.abs(rr.ravel() - r_peak) <= band
        # azimuthal average A(phi) over the radial band
        bcnt = np.maximum(np.bincount(bin_idx[sel], minlength=n_bins), 1)
        bsum = np.bincount(bin_idx[sel], weights=flat[sel], minlength=n_bins)
        A = bsum / bcnt
        F = np.fft.rfft(A)
        P = np.abs(F) ** 2
        Ptot = float(P.sum()) + 1e-40
        pm_frac[l, : min(M, len(P))] = P[: min(M, len(P))] / Ptot
        asym[l] = float(P[1:].sum() / Ptot)
        # energy-weighted core centroid over the band
        w = (amp ** 2) * (np.abs(rr - r_peak) <= band)
        wsum = float(w.sum()) + 1e-40
        cxw = float((w * xs).sum() / wsum)
        cyw = float((w * ys).sum() / wsum)
        centroid[l] = float(np.hypot(cxw - cx, cyw - cy))
        if l == 0:
            A_phi_0 = A.copy()

    # CSV
    csv_path = out_dir / "ring_azimuth.csv"
    rows = np.column_stack([times, asym, centroid, r_peaks, pm_frac])
    np.savetxt(
        str(csv_path), rows, delimiter=",",
        header="time,asymmetry_m>=1,centroid_offset,r_peak,"
               + ",".join(f"P{m}_frac" for m in range(M)),
        comments="",
    )

    # PNG
    _apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("D9 — Ring azimuthal asymmetry  (Fourier decomposition of |Psi| around the donut)",
                 fontweight="bold")

    amp0 = _amp_zmid(0)
    ax = axes[0, 0]
    im = ax.imshow(amp0.T, origin="lower", extent=[0, nx, 0, ny], aspect="equal", cmap="inferno")
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(cx + r_peaks[0] * np.cos(th), cy + r_peaks[0] * np.sin(th), "c--", lw=1.0,
            label=f"ring r={r_peaks[0]:.1f}")
    ax.plot([cx], [cy], "c+", ms=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("|Psi| z-midplane (slice 0) + sampling ring"); ax.legend(fontsize=8)

    ax = axes[0, 1]
    phis = np.linspace(-np.pi, np.pi, n_bins, endpoint=False)
    ax.plot(phis, A_phi_0, color="tab:purple")
    ax.set_xlabel("azimuth phi [rad]"); ax.set_ylabel("|Psi|(phi) on ring (slice 0)")
    ax.set_title(f"Azimuthal amplitude profile  (asym={asym[0]:.3f})")

    ax = axes[1, 0]
    ax.bar(range(M), pm_frac[0], color="tab:purple")
    ax.set_xlabel("azimuthal harmonic m"); ax.set_ylabel("power fraction P_m / ΣP (slice 0)")
    ax.set_title("Harmonic spectrum (m=1 = dipole / uneven ring)")

    ax = axes[1, 1]
    ax.plot(times, asym, label="asymmetry (Σ_{m≥1} P_m / ΣP)", color="tab:red")
    ax.plot(times, pm_frac[:, 1], label="m=1 (dipole) fraction", color="tab:orange", ls="--")
    ax.plot(times, centroid, label="core centroid offset [a]", color="tab:green", ls=":")
    ax.set_xlabel("time"); ax.legend(fontsize=8); ax.set_title("Asymmetry & core drift over the loop")

    _savefig(fig, out_dir / "ring_azimuth.png")

    return {
        "asymmetry_mean": float(np.mean(asym)),
        "asymmetry_max": float(np.max(asym)),
        "m1_dipole_frac_mean": float(np.mean(pm_frac[:, 1])),
        "centroid_offset_mean": float(np.mean(centroid)),
        "centroid_offset_max": float(np.max(centroid)),
        "r_peak_mean": float(np.mean(r_peaks)),
        "note": (
            "Azimuthal Fourier decomposition of |Psi| on the donut ring (z-midplane). "
            "Symmetric Y_1^1 seed -> asymmetry ~ 0 (all power in m=0); m>=1 power and "
            "a nonzero core-centroid offset quantify ring deformation (the 'one side "
            "brighter' effect) under relaxation. Progress tracker, no pass/fail."
        ),
        "csv": str(csv_path),
        "png": str(out_dir / "ring_azimuth.png"),
    }


# ---------------------------------------------------------------------------
# Device dispatch (shared by the real-state pass and the vacuum calibration)
# ---------------------------------------------------------------------------

_DEVICE_FNS = {
    "energy": device_energy,
    "confinement": device_confinement,
    "winding": device_winding,
    "berry": device_berry,
    "em_fields": device_em_fields,
    "color_channels": device_color_channels,
    "spectra": device_spectra,
    "ring_azimuth": device_ring_azimuth,
}


def _call_device(
    name: str,
    world: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    out_dir: Path,
    *,
    m_expected: int | None,
    n_t: int | None,
) -> dict[str, Any]:
    """Run one device, threading the two config-derived kwargs where needed.

    Used identically by the real-state measurement pass and the vacuum
    calibration pass, so a device is never exercised by two different code paths.
    """
    if name == "binding_probe":
        return device_binding_probe(world, ref, lattice, params, out_dir, n_t=n_t)
    if name == "winding":
        # Calibration passes m_expected=None (it only inspects the raw nulls);
        # the real pass passes the config m so the falsifier targets winding_z≈m.
        return device_winding(world, ref, lattice, params, out_dir, m_expected=m_expected)
    return _DEVICE_FNS[name](world, ref, lattice, params, out_dir)


def _run_calibration(
    world_vac: np.ndarray,
    ref: np.ndarray,
    lattice: SpacelikeLattice,
    params: ActionParams,
    n_t: int | None,
    verbose: bool,
) -> dict[str, dict]:
    """Run each vacuum-fixture device on the A=0 world into a throwaway dir.

    Returns device name -> metric dict (or ``{"error": ...}``).  Artifacts
    (PNG/CSV) are written into a TemporaryDirectory and discarded — only the
    returned metric dicts are inspected by the ledger.
    """
    calib: dict[str, dict] = {}
    with tempfile.TemporaryDirectory(prefix="branesim_calib_") as tmp:
        tmp_dir = Path(tmp)
        for name, contract in CONTRACTS.items():
            if contract.calib_fixture != "vacuum":
                continue
            try:
                # m_expected=None for the calibration winding run (raw nulls only).
                calib[name] = _call_device(
                    name, world_vac, ref, lattice, params, tmp_dir,
                    m_expected=None, n_t=n_t,
                )
            except Exception as exc:  # a device that can't survive A=0 fails calibration
                if verbose:
                    print(f"    [calib] device {name} raised on vacuum: {exc}")
                calib[name] = {"error": str(exc)}
    return calib


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_measurements(
    run_dir: str | Path,
    *,
    devices: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, str]:
    """Run all diagnostic devices on a finished worldvolume run folder.

    Parameters
    ----------
    run_dir : path-like
        Path to the run folder containing config.json and seed_world.npz.
    devices : list of str or None
        Subset of device names to run (default: all).
        Options: "energy", "confinement", "winding", "berry",
                 "em_fields", "color_channels", "spectra".
    verbose : bool
        Print progress.

    Returns
    -------
    dict mapping device name -> output path(s) summary string.
    """
    run_dir = Path(run_dir).resolve()
    diag_dir = run_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"[run_measurements] run_dir = {run_dir}")
        print(f"[run_measurements] diag_dir = {diag_dir}")

    data = _load_run(run_dir)
    world = data["world"]
    ref = data["ref"]
    lattice = data["lattice"]
    params = data["params"]
    config = data["config"]

    if verbose:
        print(f"  world shape: {world.shape}  dtype: {world.dtype}")
        print(f"  alpha={data['alpha']}  grid={data['grid_shape']}")

    # Extract n_t from config for device_binding_probe.
    # ActionParams is frozen=True — do NOT mutate it.  Instead pass n_t as a
    # kwarg directly to device_binding_probe via the special-case dispatch below.
    vp = config.get("vortex_params", {})
    config_n_t: int | None = None
    if "n_t" in vp:
        try:
            config_n_t = int(vp["n_t"])
        except (TypeError, ValueError):
            pass
    config_m: int | None = None
    if "m" in vp:
        try:
            config_m = int(vp["m"])
        except (TypeError, ValueError):
            pass

    all_devices = ["energy", "confinement", "winding", "berry",
                   "em_fields", "color_channels", "spectra", "ring_azimuth",
                   "binding_probe"]
    if devices is None:
        devices = all_devices

    results: dict[str, dict] = {}
    paths: dict[str, str] = {}

    for name in devices:
        if name not in _DEVICE_FNS and name != "binding_probe":
            print(f"  [WARNING] Unknown device {name!r}, skipping.")
            continue
        if verbose:
            print(f"  Running device {name} ...")
        try:
            res = _call_device(
                name, world, ref, lattice, params, diag_dir,
                m_expected=config_m, n_t=config_n_t,
            )
            results[name] = res
            paths[name] = str(diag_dir)
            if verbose:
                print(f"    -> OK  ({', '.join(str(v) for k,v in res.items() if k in ('csv','png'))})")
        except Exception as exc:
            print(f"    [ERROR] Device {name} failed: {exc}")
            import traceback
            traceback.print_exc()
            results[name] = {"error": str(exc)}

    # --- Trust ledger: calibrate every vacuum-fixture device on the A=0 state in
    #     this same run, then assemble the layer-ordered ledger (contracts.py).
    if verbose:
        print("  Running vacuum calibration for the trust ledger ...")
    world_vac = build_vacuum_world(ref, data["n_slices"], params.r_t)
    calib_results = _run_calibration(world_vac, ref, lattice, params, config_n_t, verbose)

    sr = config.get("solver_report") or {}
    converged = sr.get("converged")
    state = config.get("iter_label") or config.get("iter_kind") or "unknown state"
    if converged is not None:
        state = f"{state} (converged={converged})"
    ctx = {"converged": bool(converged), "m": config_m, "n_t": config_n_t,
           "state": state, "iter_kind": config.get("iter_kind")}
    ledger_rows, first_fail = build_ledger(results, calib_results, ctx)
    ledger_md = render_ledger_md(ledger_rows, first_fail, ctx)

    # Write report
    if verbose:
        print("  Writing report.md ...")
    report_path = _write_report(run_dir, diag_dir, config, results, ledger_md=ledger_md)
    paths["report"] = str(report_path)
    if verbose:
        print(f"  report.md -> {report_path}")

    if verbose:
        print(f"\n[run_measurements] DONE.  diagnostics/ = {diag_dir}")

    return paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m branesim.diagnostics.run_measurements <run_dir>")
        sys.exit(1)
    run_dir_arg = Path(sys.argv[1])
    if not run_dir_arg.exists():
        print(f"Error: {run_dir_arg} does not exist")
        sys.exit(1)
    run_measurements(run_dir_arg, verbose=True)
