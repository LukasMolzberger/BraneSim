"""Experiment: inject U(1) vortex-ring seed and render it.

TASK: "inject the object and see it" — validates the injection map, the
phase single-valuedness, and the render pipeline.  NO eigen-solve.

Physics spec: EXPERIMENT.md
Layer: E (experiments) — orchestrates init → visualize only.
Principles: §2 (no back-reaction), §3.2 (no clamps), §7.2 (no circular deps).

Geometry: single vortex RING (smoke ring).  The ring core is a circle of
radius R_ring in the z=z_c midplane; a torus of energy surrounds it.  The
U(1) phase winds m=1 times around the tube cross-section (meridionally).
Net winding through every periodic plane = 0 (contractible ring).

Usage (from repo root)::

    python -m branesim.experiments.vortex_seed_render

Output::

    runs/vortex_seed_<YYYY-MM-DD>_<HHMMSS>/
        config.json
        seed_world.npz          # the injected worldvolume (compact, float32)
        winding_closure.json    # net winding per periodic plane
        seed_snapshot_t0.png    # mid-plane slice at t=0 (phase->RGB + amplitude)
        renders/
            volume_phase.mp4    # 3D voxel, hue=U(1) phase, opacity=|u|
            slice_xy_phase.mp4  # XY midplane, phase->RGB
            slice_xz_phase.mp4  # XZ midplane, phase->RGB
            slice_yz_phase.mp4  # YZ midplane, phase->RGB
            volume_amp.mp4      # 3D voxel, amplitude only (inferno)
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.vortex_worldtube import (
    CARRIER_IM,
    CARRIER_RE,
    VortexParams,
    inject_vortex_worldtube,
    measure_winding_closure,
)
from branesim.visualization.volume_render import (
    create_slice_animation,
    create_volume_animation,
    phase_to_rgb,
)


# ---------------------------------------------------------------------------
# Experiment parameters
# ---------------------------------------------------------------------------

# Box: 48^3 with ring radius 6a, tube width 2.5a => torus fits in [~11a, ~37a]
# giving ~11a vacuum margin on all sides.
GRID_SHAPE = (48, 48, 48)
N_SLICES = 32                       # time slices; ~16 carrier periods at omega*dt=0.5
SPACING = 1.0                       # lattice unit

ALPHA = 0.7
DT = 0.25
BETA = 1.0
R_T = ALPHA * BETA * DT            # = 0.175

VORTEX_PARAMS = VortexParams(
    A0=0.3,                         # peak strain ~0.3
    w=2.5,                          # tube width (minor radius) ~2.5a
    R_ring=6.0,                     # ring major radius ~6a
    m_wind=1,                       # meridional winding number
    omega=2.0,                      # omega*dt = 0.5 rad/step
    geometry="vortex_ring",
)

# Render settings
FPS = 15
DPI = 100
DENSITY_THRESHOLD = 0.04
ALPHA_SCALE = 0.85
GAMMA = 0.65


# ---------------------------------------------------------------------------
# Helper: build displacement amplitude + phase fields for rendering
# ---------------------------------------------------------------------------


def _extract_amp_phase(
    world: np.ndarray,
    ref: np.ndarray,
    grid_shape: tuple[int, int, int],
) -> tuple[list[np.ndarray], list[np.ndarray], list[float]]:
    """Extract per-slice amplitude and phase of the carrier field.

    amplitude = sqrt(Re(u)^2 + Im(u)^2)
    phase     = atan2(Im(u), Re(u))

    Returns
    -------
    amps, phases : list of (nx,ny,nz) arrays
    times : list of floats (physical time = slice * DT)
    """
    nx, ny, nz = grid_shape
    amps = []
    phases = []
    times = []

    for l in range(world.shape[0]):
        pos = world[l]
        re_disp = pos[:, CARRIER_RE] - ref[:, CARRIER_RE]
        im_disp = pos[:, CARRIER_IM] - ref[:, CARRIER_IM]

        amp = np.sqrt(re_disp ** 2 + im_disp ** 2).reshape(nx, ny, nz)
        ph = np.arctan2(im_disp, re_disp + 1e-300).reshape(nx, ny, nz)

        amps.append(amp)
        phases.append(ph)
        times.append(float(l) * DT)

    return amps, phases, times


# ---------------------------------------------------------------------------
# Mid-plane snapshot (static PNG)
# ---------------------------------------------------------------------------


def _save_snapshot(
    amp: np.ndarray,
    phase: np.ndarray,
    grid_shape: tuple[int, int, int],
    spacing: float,
    R_ring: float,
    w: float,
    out_path: str,
) -> None:
    """Save a three-panel static snapshot of the t=0 seed.

    Panels:
      Left:   XZ midplane — phase->RGB, opacity = amplitude  (shows ring cross-section)
      Centre: XY midplane — phase->RGB, opacity = amplitude  (shows ring plan view)
      Right:  XZ midplane — amplitude only (inferno)

    The XZ midplane at y=y_c cuts through the torus cross-section and should
    show TWO bright spots (where the ring enters and exits the midplane) — the
    signature of a single ring rather than two separate tubes.
    """
    nx, ny, nz = grid_shape
    amp_max = float(np.max(amp))
    if amp_max < 1e-30:
        amp_max = 1.0

    # ---- Extract slices ----
    # XZ midplane (y = ny//2): shows ring cross-section
    sl_xz_amp = amp[:, ny // 2, :]          # (nx, nz)
    sl_xz_ph = phase[:, ny // 2, :]

    # XY midplane (z = nz//2): shows ring plan view (donut shape)
    sl_xy_amp = amp[:, :, nz // 2]          # (nx, ny)
    sl_xy_ph = phase[:, :, nz // 2]

    def _to_rgba(sl_amp, sl_ph, amax):
        rgb = phase_to_rgb(sl_ph)
        alpha = np.clip(sl_amp / amax, 0, 1)[..., None]
        return np.concatenate([rgb, alpha], axis=-1)

    rgba_xz = _to_rgba(sl_xz_amp, sl_xz_ph, amp_max)   # (nx, nz, 4)
    rgba_xy = _to_rgba(sl_xy_amp, sl_xy_ph, amp_max)   # (nx, ny, 4)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(
        f"Vortex ring seed  t=0 | R_ring={R_ring}a  w={w}a  m=1\n"
        "Left: XZ slice (ring cross-section — expect 2 spots)  "
        "Centre: XY slice (plan view — expect donut ring)  "
        "Right: XZ amplitude",
        fontsize=9,
    )

    ext_xz = [0, nx * spacing, 0, nz * spacing]
    ext_xy = [0, nx * spacing, 0, ny * spacing]

    ax = axes[0]
    ax.imshow(rgba_xz.swapaxes(0, 1), origin="lower", extent=ext_xz, aspect="equal")
    ax.set_xlabel("x"); ax.set_ylabel("z")
    ax.set_title("XZ midplane — phase→RGB  (ring cross-section)")

    ax = axes[1]
    ax.imshow(rgba_xy.swapaxes(0, 1), origin="lower", extent=ext_xy, aspect="equal")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title("XY midplane — phase→RGB  (plan view)")

    ax = axes[2]
    im = ax.imshow(
        sl_xz_amp.T, origin="lower", extent=ext_xz, aspect="equal",
        cmap="inferno", vmin=0, vmax=amp_max,
    )
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x"); ax.set_ylabel("z")
    ax.set_title("XZ midplane — amplitude")

    plt.tight_layout()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=DPI)
    plt.close(fig)
    print(f"  Snapshot saved: {out}")


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------


def run() -> None:
    # --- Run folder ---
    ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    run_dir = Path(__file__).parents[2] / "runs" / f"vortex_seed_{ts}"
    renders_dir = run_dir / "renders"
    run_dir.mkdir(parents=True, exist_ok=True)
    renders_dir.mkdir(parents=True, exist_ok=True)

    print(f"Run folder: {run_dir}")

    # --- Lattice and action params ---
    lattice_params = LatticeParams(
        grid_shape=GRID_SHAPE,
        spacing=SPACING,
        periodic_axes=(True, True, True),
    )
    action_params = ActionParams(
        k_s=1.0,
        alpha=ALPHA,
        rho=1.0,
        dt=DT,
        n_slices=N_SLICES,
        m_ambient=4,
        r_t=R_T,
        beta=BETA,
    )
    lattice = SpacelikeLattice(lattice_params)

    # --- Config ---
    config = {
        "experiment": "vortex_seed_render",
        "geometry": "vortex_ring",
        "timestamp": ts,
        "grid_shape": list(GRID_SHAPE),
        "n_slices": N_SLICES,
        "spacing": SPACING,
        "alpha": ALPHA,
        "dt": DT,
        "beta": BETA,
        "r_t": R_T,
        "vortex_params": VORTEX_PARAMS._asdict(),
        "carrier_re_component": CARRIER_RE,
        "carrier_im_component": CARRIER_IM,
        "render": {
            "fps": FPS,
            "dpi": DPI,
            "density_threshold": DENSITY_THRESHOLD,
            "alpha_scale": ALPHA_SCALE,
            "gamma": GAMMA,
        },
        "notes": (
            "Single vortex ring (smoke ring) geometry. Ring core is a circle of "
            f"radius R_ring={VORTEX_PARAMS.R_ring}a in the z=z_c midplane. "
            "Torus of energy (tube width w) surrounds the core. "
            "U(1) phase winds m=1 times meridionally around the tube. "
            "Net winding through any periodic plane = 0 (contractible ring). "
            "XZ midplane shows 2 bright spots (where ring pierces the plane). "
            "XY midplane shows the donut/ring plan view."
        ),
    }
    (run_dir / "config.json").write_text(json.dumps(config, indent=2))
    print("Config written.")

    # --- Injection ---
    print("Injecting vortex-ring seed...")
    t0 = time.perf_counter()
    world, seed_meta = inject_vortex_worldtube(
        lattice, action_params, VORTEX_PARAMS, N_SLICES
    )
    dt_inject = time.perf_counter() - t0
    print(f"  Injection done in {dt_inject:.2f} s  "
          f"world.shape={world.shape}  dtype={world.dtype}")

    # Save seed worldvolume (float32 for speed; diagnostic only)
    np.savez_compressed(str(run_dir / "seed_world.npz"), world=world.astype(np.float32))
    print("  seed_world.npz saved.")

    # --- Winding closure verification ---
    print("Measuring winding closure...")
    ref = lattice.reference_positions(4)
    winding = measure_winding_closure(world, lattice, slice_index=0)
    winding_report = {
        "geometry": VORTEX_PARAMS.geometry,
        "slice_index": 0,
        "windings": winding,
        "expected": 0.0,
        "max_abs_winding": max(abs(v) for v in winding.values()),
        "closure_ok": all(abs(v) < 0.1 for v in winding.values()),
    }
    (run_dir / "winding_closure.json").write_text(json.dumps(winding_report, indent=2))
    print(f"  Winding: {winding}")
    print(f"  Closure OK: {winding_report['closure_ok']}")

    # --- Field extraction ---
    print("Extracting amplitude + phase fields for rendering...")
    amps, phases, times_list = _extract_amp_phase(world, ref, GRID_SHAPE)

    peak_amp = max(float(np.max(a)) for a in amps)
    print(f"  Peak amplitude: {peak_amp:.4f}  (target ~{VORTEX_PARAMS.A0 * 0.607:.3f})")

    # --- Static mid-plane snapshot ---
    print("Saving mid-plane snapshot (t=0)...")
    _save_snapshot(
        amps[0], phases[0], GRID_SHAPE, SPACING,
        VORTEX_PARAMS.R_ring, VORTEX_PARAMS.w,
        str(run_dir / "seed_snapshot_t0.png"),
    )

    # --- Renders ---
    print("Rendering volume (phase->RGB)...")
    t0 = time.perf_counter()
    create_volume_animation(
        frames_amplitude=amps,
        frames_phase=phases,
        grid_shape=GRID_SHAPE,
        spacing=SPACING,
        output_path=str(renders_dir / "volume_phase.mp4"),
        times=times_list,
        fps=FPS,
        dpi=DPI,
        density_threshold=DENSITY_THRESHOLD,
        alpha_scale=ALPHA_SCALE,
        gamma=GAMMA,
        title_prefix="U(1) vortex ring",
    )
    print(f"  volume_phase.mp4 done ({time.perf_counter()-t0:.1f}s)")

    print("Rendering volume (amplitude, inferno)...")
    t0 = time.perf_counter()
    create_volume_animation(
        frames_amplitude=amps,
        frames_phase=None,
        grid_shape=GRID_SHAPE,
        spacing=SPACING,
        output_path=str(renders_dir / "volume_amp.mp4"),
        times=times_list,
        cmap_name="inferno",
        fps=FPS,
        dpi=DPI,
        density_threshold=DENSITY_THRESHOLD,
        alpha_scale=ALPHA_SCALE,
        gamma=GAMMA,
        title_prefix="U(1) vortex ring | amplitude",
    )
    print(f"  volume_amp.mp4 done ({time.perf_counter()-t0:.1f}s)")

    for plane in ("xy", "xz", "yz"):
        print(f"Rendering slice {plane} (phase->RGB)...")
        t0 = time.perf_counter()
        create_slice_animation(
            frames_field=amps,
            frames_phase=phases,
            grid_shape=GRID_SHAPE,
            spacing=SPACING,
            plane=plane,
            output_path=str(renders_dir / f"slice_{plane}_phase.mp4"),
            times=times_list,
            fps=FPS,
            dpi=DPI,
            title_prefix="U(1) vortex ring",
        )
        print(f"  slice_{plane}_phase.mp4 done ({time.perf_counter()-t0:.1f}s)")

    # --- Summary ---
    print("\n=== DONE ===")
    print(f"Run folder : {run_dir}")
    print(f"Geometry   : {VORTEX_PARAMS.geometry}  R_ring={VORTEX_PARAMS.R_ring}a  "
          f"w={VORTEX_PARAMS.w}a  m={VORTEX_PARAMS.m_wind}")
    print(f"Winding closure: {winding}")
    print(f"Closure OK : {winding_report['closure_ok']}")
    print(f"Peak amplitude : {peak_amp:.4f}")
    print(f"Snapshot   : {run_dir}/seed_snapshot_t0.png")
    renders = sorted(renders_dir.glob("*.mp4"))
    for r in renders:
        size_kb = r.stat().st_size / 1024
        print(f"  {r.name}  ({size_kb:.0f} kB)")


if __name__ == "__main__":
    run()
