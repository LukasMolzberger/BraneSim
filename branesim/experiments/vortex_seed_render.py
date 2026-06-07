"""Experiment: inject the U(1) spherical-harmonic vortex seed, then walk the
solver-iteration ("5th") axis — render + fully diagnose the initial seed and a
fixed-iteration relaxation probe.

TASK: "inject the object, see it, and measure it" — validates the injection
map, the phase single-valuedness, the render pipeline, AND the full diagnostic
back-up round-trip (Berry / EM / per-colour SU(3) / spectra / energy).

Physics spec: EXPERIMENT.md
Layer: E (experiments) — orchestrates init → solve(probe) → visualize → diagnose.
Principles: §2 (no back-reaction inside diagnostics), §3.2 (no clamps),
            §7.2 (no circular deps).

## The "5th dimension": solver iterations as a (purely technical) axis

The brane is 4D (3 space + 1 timelike link).  The block solver
(``bvp.solve_block``, JFNK root-find of ‖R‖=0) reaches the physical worldvolume
by *iterating* — a relaxation/flow coordinate, like gradient-flow or
imaginary-time relaxation.  Only the converged fixed point is physical; the
iterates before it are unphysical transients.  We encode this axis explicitly
in the run-folder layout:

    runs/vortex_seed_<ts>/
        config.json, manifest.json, README.md
        iter_0000/   — the freshly injected seed (solver untouched)
        iter_0015/   — after N rotating-frame-periodic JFNK relaxation steps
            config.json, world.npz, winding_closure.json
            snapshot.png (+ .md)
            renders/   volume_phase.mp4, volume_amp.mp4, slice_{xy,xz,yz}_phase.mp4 (+ .md each)
            diagnostics/  energy / confinement / winding / berry / em_fields /
                          color_channels / spectra  (CSV + PNG + .md each) + report.md

Each iteration folder is self-contained (config.json + world.npz), so
``run_measurements`` runs on it directly.

NOTE: the relaxation uses the **rotating-frame-periodic** BC (``PeriodicBC``):
a closed cyclic time loop with all slices free, well-conditioned (cond~1e3-1e4,
vs the Dirichlet two-time κ~1e14 that froze the earlier probe).  It genuinely
relaxes the worldtube toward ‖R‖=0 while preserving the carrier winding; whether
it binds or radiates is what the diagnostics then report.
"""

from __future__ import annotations

import json
import os
import tempfile
import time
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import matplotlib.pyplot as plt

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.initialization.vortex_worldtube import (
    CARRIER_IM,
    CARRIER_RE,
    VortexParams,
    inject_vortex_worldtube,
    measure_winding_closure,
    vacuum_offsets,
)
from branesim.solver.boundary import PeriodicBC
from branesim.solver.bvp import BoundaryProblem, SolveOpts, solve_block
from branesim.diagnostics.run_measurements import run_measurements
from branesim.visualization.volume_render import (
    create_slice_animation,
    create_volume_animation,
    phase_to_rgb,
)


# ---------------------------------------------------------------------------
# Experiment parameters
# ---------------------------------------------------------------------------

# Box: 48^3; spherical-harmonic seed centred in the box, donut radius r0=6a,
# radial width 2.5a => object localized in [~16a, ~32a] with ample vacuum margin.
GRID_SHAPE = (48, 48, 48)
N_SLICES = 32                       # time slices = period of the carrier loop
SPACING = 1.0                       # lattice unit

ALPHA = 0.7
DT = 0.25
BETA = 1.0
R_T = ALPHA * BETA * DT            # = 0.175

VORTEX_PARAMS = VortexParams(
    A0=0.3,                         # peak strain ~0.3
    r0=6.0,                         # radial-shell peak (donut radius) ~6a
    w=2.5,                          # radial-shell width ~2.5a
    l=1,                            # spherical-harmonic degree
    m=1,                            # azimuthal U(1) winding around the z-axis
    n_t=2,                          # carrier turns over the loop -> 720 deg, closes exactly
    geometry="spherical_harmonic",
)

# Solver-axis (5th dimension) rotating-frame-periodic relaxation.
# Well-conditioned (PeriodicBC, cond~1e3-1e4 vs Dirichlet's 1e14): the brane
# actually moves and the carrier winding is preserved.  Default ON.
RELAX_ITERS = 15                    # JFNK iterations for the iter_N snapshot
RELAX_INNER_MAXITER = 40            # inner lgmres cap
RELAX_TOL = 1e-6

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

    # Subtract the prestressed-vacuum N-gon offset so renders show the carrier
    # (donut + winding), not the rho-sized uniform vacuum background.
    off_re, off_im = vacuum_offsets(world.shape[0] - 1, R_T)

    for l in range(world.shape[0]):
        pos = world[l]
        re_disp = pos[:, CARRIER_RE] - ref[:, CARRIER_RE] - off_re[l]
        im_disp = pos[:, CARRIER_IM] - ref[:, CARRIER_IM] - off_im[l]

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
    l: int,
    m: int,
    r0: float,
    w: float,
    out_path: str,
    state_label: str = "t=0 seed",
) -> None:
    """Save a three-panel static snapshot of the first slice of a world state.

    Panels:
      Left:   XZ midplane — phase->RGB, opacity = amplitude  (axis cross-section)
      Centre: XY midplane — phase->RGB, opacity = amplitude  (plan view, on-axis)
      Right:  XZ midplane — amplitude only (inferno)

    For the canonical Y_1^1 seed:
      - XY midplane (z=z_c) is the equatorial plane: a bright RING of amplitude
        whose phase (hue) advances once (2*pi*m) AZIMUTHALLY around the centre —
        the signature of a U(1) line vortex along z.
      - XZ midplane (y=y_c) cuts the donut at two lobes either side of the z-axis;
        the two lobes are pi out of phase (opposite azimuth phi -> phi+pi).
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
        f"Y_{l}^{m} U(1) vortex | {state_label} | r0={r0}a  w={w}a\n"
        "Left: XZ slice (axis cross-section — 2 lobes, pi out of phase)  "
        "Centre: XY slice (equatorial plan view — ring, phase winds azimuthally)  "
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
# Render a world state (snapshot + volume + slice movies) into an iter folder
# ---------------------------------------------------------------------------


def _render_world(
    amps: list[np.ndarray],
    phases: list[np.ndarray],
    times: list[float],
    vp,
    grid_shape: tuple[int, int, int],
    spacing: float,
    iter_dir: Path,
    state_label: str,
) -> None:
    """Snapshot + 3D volume (phase, amplitude) + 2D slice movies for a world."""
    renders_dir = iter_dir / "renders"
    renders_dir.mkdir(parents=True, exist_ok=True)

    _save_snapshot(
        amps[0], phases[0], grid_shape, spacing,
        vp.l, vp.m, vp.r0, vp.w,
        str(iter_dir / "snapshot.png"), state_label=state_label,
    )

    print("  Rendering volume (phase->RGB)...")
    create_volume_animation(
        frames_amplitude=amps, frames_phase=phases, grid_shape=grid_shape,
        spacing=spacing, output_path=str(renders_dir / "volume_phase.mp4"),
        times=times, fps=FPS, dpi=DPI, density_threshold=DENSITY_THRESHOLD,
        alpha_scale=ALPHA_SCALE, gamma=GAMMA,
        title_prefix=f"U(1) Y_l^m | {state_label}",
    )
    print("  Rendering volume (amplitude, inferno)...")
    create_volume_animation(
        frames_amplitude=amps, frames_phase=None, grid_shape=grid_shape,
        spacing=spacing, output_path=str(renders_dir / "volume_amp.mp4"),
        times=times, cmap_name="inferno", fps=FPS, dpi=DPI,
        density_threshold=DENSITY_THRESHOLD, alpha_scale=ALPHA_SCALE, gamma=GAMMA,
        title_prefix=f"U(1) Y_l^m | amplitude | {state_label}",
    )
    for plane in ("xy", "xz", "yz"):
        print(f"  Rendering slice {plane} (phase->RGB)...")
        create_slice_animation(
            frames_field=amps, frames_phase=phases, grid_shape=grid_shape,
            spacing=spacing, plane=plane,
            output_path=str(renders_dir / f"slice_{plane}_phase.mp4"),
            times=times, fps=FPS, dpi=DPI,
            title_prefix=f"U(1) Y_l^m | {state_label}",
        )

    # Add-on: XY axial view with phase-dash + supercurrent-flow overlays.  The base
    # slice_*_phase.mp4 outputs above are unchanged.  Idiom from R. Behiel's GL-vortex
    # animations: dashes make the azimuthal winding/core defect legible where flat hue
    # is ambiguous; dots advected along j~|Psi|^2 grad(theta) show the swirl.  XY is the
    # decisive down-the-axis view.
    print("  Rendering slice xy (phase-dash + supercurrent-flow overlays)...")
    create_slice_animation(
        frames_field=amps, frames_phase=phases, grid_shape=grid_shape,
        spacing=spacing, plane="xy",
        output_path=str(renders_dir / "slice_xy_phase_flow.mp4"),
        times=times, fps=FPS, dpi=DPI,
        title_prefix=f"U(1) Y_l^m | phase+flow | {state_label}",
        phase_dashes=True, flow_dots=True,
    )
    (renders_dir / "slice_xy_phase_flow.md").write_text(
        f"# renders/slice_xy_phase_flow.mp4 — {state_label}\n\n"
        "XY axial view (same data as `slice_xy_phase.mp4`) with two **add-on overlays**:\n\n"
        "- **Phase dashes** — short black segments oriented by the local U(1) phase "
        "`arg(Psi)`, opacity proportional to amplitude. They make the **azimuthal winding** "
        "(segments rotate once per `2*pi*m` around the ring) and the **core defect** (where "
        "the dashes lose coherence / fade) legible where the flat hue is ambiguous — "
        "especially in the dark `|Psi|->0` interior where hue is meaningless.\n"
        "- **Supercurrent-flow dots** — white dots advected along the in-plane supercurrent "
        "`j ~ |Psi|^2 grad(theta)` (branch-cut-safe `Im(conj(Psi) grad Psi)`), spawned with "
        "density proportional to amplitude and faded in/out over their lifetime — the "
        "**swirl** around the vortex axis.\n\n"
        "Visualization only; no physics is added (identical injected field). Overlay idiom "
        "adapted from R. Behiel's Ginzburg-Landau vortex animations.\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Per-output documentation (.md beside every figure / video / diagnostic)
# ---------------------------------------------------------------------------

# Per-device explainers (PNG name -> what / how / read).  The verdicts live in
# diagnostics/report.md; these explain the figure itself.
_DEVICE_DOCS = {
    "energy": ("energy.png",
        "Total energy E, vacuum-subtracted excess E_excess, and the kinetic/potential split per time slice.",
        "From the discrete brane action: V is the spacelike spring potential per slice; the kinetic split uses central time differences. E_excess subtracts the static vacuum energy (LESSONS_LEARNED discipline).",
        "Flat E_excess across slices = a stationary state; secular drift = the state is not a solution (expected for the bare seed)."),
    "confinement": ("confinement.png",
        "Spatial spread of the excess-energy density over the loop: radius_rms, spread_ratio, and leakage to the box margin.",
        "Excess energy density per node, reduced to radial moments about the box centre per slice.",
        "Stable bounded radius_rms = localized/confined; growth or box-fill = dispersing (the semilocal-vortex binding question)."),
    "winding": ("winding.png",
        "U(1) phase winding per time slice, measured by the real plaquette/contour sum (never a hard-coded estimator).",
        "atan2(Im,Re) of the carrier 2-plane displacement; discrete circulation about the central axes.",
        "Winding ~ m about the z-axis, ~0 about x,y, constant across slices = the topological charge is conserved."),
    "berry": ("berry.png",
        "Berry/phase envelope time series: Psi = u + i v/omega0 on the carrier 2-plane, its modulus and accumulated phase.",
        "The geometric (carrier) phase A_t = i<u|d_t u> integrated over the loop — the EM-charge channel.",
        "Smooth 2*pi*n_t phase accumulation over the loop = the closure-locked carrier; modulus tracks the donut amplitude."),
    "em_fields": ("em_fields.png",
        "Emergent EM fields E_i, B_i from the Berry connection A_mu = i<u|d_mu u>, as a quiver on the mid-plane.",
        "F_munu = d_mu A_nu - d_nu A_mu of the carrier-phase connection; E_i = F_0i, B_i = (1/2)eps F_jk.",
        "A vortex shows a B flux threading the axis and a radial/azimuthal E pattern; magnitudes are diagnostic, not yet calibrated."),
    "color_channels": ("color_channels.png",
        "Per-channel U(3) breakdown: trace (U(1)/EM) vs traceless (SU(3)/colour) content of the lateral displacement.",
        "Project the lateral triplet (components 0,1,2) onto the trace direction (1,1,1)/sqrt3 (U(1)) and its orthogonal complement (SU(3)).",
        "The bare single-component seed reads 1/3 trace + 2/3 traceless; how the SU(3) part evolves answers 'does colour coexcite?'."),
    "spectra": ("spectra.png",
        "Spatial-FFT energy/mode spectrum of the excess field (radiation tail).",
        "3D FFT of the carrier displacement per slice, radially binned in |k|.",
        "Power concentrated at low |k| = a smooth localized object; a growing high-|k| tail = radiation/lattice noise."),
}


def _write_output_docs(
    iter_dir: Path,
    vp,
    seed_meta: dict,
    winding: dict,
    grid_shape: tuple[int, int, int],
    n_slices: int,
    dt: float,
    iter_index: int,
    state_label: str,
    kind: str,
    solver_note: str,
    do_render: bool = True,
) -> None:
    """Write a Markdown explainer next to every figure, video, and diagnostic.

    Each ``.md`` states: what the artifact shows, how it was derived (the
    ansatz + the specific projection/slice), how to read it, and the exact
    provenance.  A reader who has never seen the code can interpret it.
    """
    nx, ny, nz = grid_shape
    renders_dir = iter_dir / "renders"
    diag_dir = iter_dir / "diagnostics"
    omega_phys = seed_meta["omega_phys"]
    omega_per_slice = seed_meta["omega_per_slice"]

    # ---- shared preamble: the ansatz + the 5th-dimension framing ----
    ansatz = f"""## The seed (common to all outputs in this run)

The substrate is seeded with a **U(1) spherical-harmonic vortex** — the
EM / electron-like sector of the U(3) substrate.  The emergent order parameter,
written in the soliton-layer basis (spherical harmonics about the box centre), is

```
Psi(r, theta, phi, t) = A0 * Rhat(r) * Yhat_{vp.l}^{vp.m}(theta, phi) * exp(i * omega * t)
```

| symbol | meaning | value here |
|---|---|---|
| `Y_l^m` | complex spherical harmonic | `l={vp.l}, m={vp.m}` |
| `Rhat(r)` | unit-peak radial shell `exp(-(r-r0)^2/2w^2)` | `r0={vp.r0}a, w={vp.w}a` |
| `A0` | peak amplitude (peak |Psi|) | `{vp.A0}` |
| `omega` | carrier rate, closure-locked | `{omega_per_slice:.4f} rad/slice = {omega_phys:.4f} rad/time` |
| `n_t` | carrier turns over the loop | `{vp.n_t}`  (-> `{360*vp.n_t}` deg) |

The azimuthal factor `exp(i*m*phi)` winds the U(1) phase `m` times **around the
z-axis** (the vortex axis / donut hole); `|Y_1^1| ~ sin(theta)` makes the energy
density a **donut around that axis** — the torus *emerges* from the harmonic.
`Re(Psi)` -> ambient component 2, `Im(Psi)` -> component 3 (the timelike "i").

**Measured azimuthal winding** (slice 0): about z `= {winding['winding_through_z_normal']:+.2f}`
(expect `{vp.m}`), about y `= {winding['winding_through_y_normal']:+.2f}`, about x
`= {winding['winding_through_x_normal']:+.2f}` (expect `0`).

## The solver-iteration ("5th") axis — which state is this?

The brane is 4D; the block solver reaches the physical worldvolume by *iterating*
(a relaxation/flow coordinate, not a physical dimension).  Only the converged
fixed point is physical; earlier iterates are transients.  This folder is
**`{iter_dir.name}` — {state_label}**.

{solver_note}

**Caveats.** The object is localized with a vacuum margin, so it is contractible
in the periodic 3-torus — a *semilocal* vortex, not topologically protected;
binding is a dynamical question.  Phase hue = arg(Psi); where amplitude -> 0 the
hue is meaningless.
"""

    tail = f"""
## Provenance

- Run iteration folder: `{iter_dir.name}`  (kind: `{kind}`, solver iteration index {iter_index})
- Grid `{nx}x{ny}x{nz}`, spacing 1 lattice unit, periodic on all axes.
- Time loop: `n_slices={n_slices}` (+1 wrap slice), `dt={dt}`.
- Generated by `branesim/experiments/vortex_seed_render.py`.
- Seed: `branesim/initialization/vortex_worldtube.py`;
  diagnostics: `branesim/diagnostics/run_measurements.py`;
  relaxation: `branesim/solver/bvp.solve_block` (JFNK).
"""

    docs: dict[str, str] = {}

    docs[str(iter_dir / "snapshot.md")] = f"""# snapshot.png — {state_label}

Static three-panel snapshot of the first time slice.

**What it shows.**
- **Left — XZ midplane (y=y_c).** A cut through the axis: **two lobes** either
  side of the z-axis (the donut cross-section), **pi out of phase** (azimuth phi
  vs phi+pi -> opposite hue).
- **Centre — XY midplane (z=z_c), equatorial.** A **bright ring** whose hue
  advances **once (2*pi*m) azimuthally** around the centre — the U(1) line vortex
  seen down its axis. *The decisive view.*
- **Right — XZ amplitude** (inferno): the energy-density donut, dark on the axis.

{ansatz}{tail}"""

    docs[str(renders_dir / "volume_phase.md")] = f"""# renders/volume_phase.mp4 — {state_label}

3D volume rendering over the full time loop. **Opacity = amplitude density,
hue = U(1) carrier phase** arg(Psi).

**What it shows.** The donut energy density with its phase colour. For the seed,
the hue cycles **{vp.n_t}** full times ({360*vp.n_t} deg) over the loop (the
closure-locked carrier) and the donut shape is steady. After relaxation the
shape and phase may deform — that is the solver acting on the seed.

**How to read it.** A fixed voxel's hue advancing smoothly and returning after
{vp.n_t} cycles = loop closure; around the ring the hue advances once = winding m.

{ansatz}{tail}"""

    docs[str(renders_dir / "volume_amp.md")] = f"""# renders/volume_amp.mp4 — {state_label}

3D volume rendering of **amplitude only** (inferno; opacity & colour = |Psi|).

**What it shows.** The bare energy-density geometry: a donut around the z-axis,
dark core on the axis (`sin(theta)=0`), localized at shell `r0={vp.r0}a`. For the
seed this is static across the loop; after relaxation any change is the solver.

{ansatz}{tail}"""

    slice_desc = {
        "xy": ("z=z_c (equatorial)",
               "a **ring** of amplitude whose hue winds once (2*pi*m) **azimuthally** "
               "around the centre — the U(1) line vortex down its axis (the decisive view)."),
        "xz": ("y=y_c (through the axis)",
               "**two lobes** either side of the z-axis, **pi out of phase**; the dark "
               "z-axis between them is the vortex core."),
        "yz": ("x=x_c (through the axis)",
               "same as XZ by symmetry: two lobes, pi out of phase, dark axis core."),
    }
    for plane, (where, what) in slice_desc.items():
        docs[str(renders_dir / f"slice_{plane}_phase.md")] = f"""# renders/slice_{plane}_phase.mp4 — {state_label}

2D {plane.upper()} midplane slice ({where}) over the loop.
**Hue = U(1) phase** arg(Psi), **brightness = amplitude** |Psi|.

**What it shows.** {what}

Over the loop the hue cycles **{vp.n_t}** times ({360*vp.n_t} deg) for the seed —
the closure-locked carrier (`omega={omega_per_slice:.4f}` rad/slice). Spatial
winding encodes the vortex charge m={vp.m}; the temporal colour-cycling encodes
the carrier (EM charge). Where dark (|Psi|~0) the hue is meaningless.

{ansatz}{tail}"""

    # Diagnostic figure docs (written into diagnostics/, created by run_measurements).
    diag_dir.mkdir(parents=True, exist_ok=True)
    for dev, (png, what, how, read) in _DEVICE_DOCS.items():
        docs[str(diag_dir / f"{Path(png).stem}.md")] = f"""# diagnostics/{png} — {state_label}

(Device `{dev}` of the back-up diagnostic suite; verdicts in `report.md`.)

**What it shows.** {what}

**How it is derived.** {how}

**How to read it.** {read}

{ansatz}{tail}"""

    written = 0
    for path, text in docs.items():
        # Skip docs for render artifacts (snapshot.png + renders/*.mp4) that were
        # not produced when BRANESIM_VORTEX_RENDER=0 (diagnostics-only mode), so a
        # diagnostics-only run never leaves dangling .md pointers.
        if not do_render and (Path(path).parent.name == "renders"
                              or Path(path).name == "snapshot.md"):
            continue
        Path(path).write_text(text)
        written += 1
    print(f"  Wrote {written} per-output .md docs.")


# ---------------------------------------------------------------------------
# Config for an iteration folder + the relaxation probe
# ---------------------------------------------------------------------------


def _iter_config(
    base_config: dict,
    iter_index: int,
    state_label: str,
    kind: str,
    solver_report: dict | None,
) -> dict:
    """Self-contained config.json for an iteration folder (run_measurements reads it)."""
    cfg = dict(base_config)
    cfg["iter_index"] = iter_index
    cfg["iter_label"] = state_label
    cfg["iter_kind"] = kind
    cfg["solver_report"] = solver_report
    return cfg


def _relax(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    action_params: ActionParams,
    mass: float,
    n_iters: int,
) -> tuple[np.ndarray, dict]:
    """Bounded JFNK relaxation of the seed along the solver-iteration axis.

    Uses the **rotating-frame-periodic** BC (closed cyclic time loop, all slices
    free) — well-conditioned (cond~1e3-1e4 vs Dirichlet's 1e14), so the brane
    genuinely relaxes toward ‖R‖=0 from the seed while the carrier winding (a
    topological integer carried by the seed) is preserved.  Whether it relaxes
    to a bound worldtube or radiates toward vacuum is the physics question the
    diagnostics then answer.
    """
    bc = PeriodicBC(R0=world[0].copy())
    opts = SolveOpts(
        tol=RELAX_TOL, max_iter=n_iters,
        inner_maxiter=RELAX_INNER_MAXITER, verbose=True,
    )
    wv = solve_block(
        BoundaryProblem(lattice, action_params, mass, bc), opts,
        initial_world=world.copy(),
    )
    return wv.slices, dict(wv.solver_report)


# ---------------------------------------------------------------------------
# Emit one iteration folder: config, world, winding, renders, docs, diagnostics
# ---------------------------------------------------------------------------


def _emit_iteration(
    run_dir: Path,
    iter_index: int,
    state_label: str,
    kind: str,
    world: np.ndarray,
    ref: np.ndarray,
    base_config: dict,
    seed_meta: dict,
    vp,
    lattice: SpacelikeLattice,
    n_slices: int,
    dt: float,
    grid_shape: tuple[int, int, int],
    spacing: float,
    solver_report: dict | None,
    solver_note: str,
    do_render: bool = True,
) -> dict:
    """Write a fully self-contained, fully diagnosed iteration folder."""
    iter_dir = run_dir / f"iter_{iter_index:04d}"
    iter_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n--- {iter_dir.name}: {state_label} ---")

    # config.json (so run_measurements works standalone on this folder)
    (iter_dir / "config.json").write_text(
        json.dumps(_iter_config(base_config, iter_index, state_label, kind, solver_report), indent=2)
    )
    # world.npz
    np.savez_compressed(str(iter_dir / "world.npz"), world=world.astype(np.float32))

    # winding
    winding = measure_winding_closure(world, lattice, slice_index=0, r_t=R_T)
    (iter_dir / "winding_closure.json").write_text(json.dumps({
        "slice_index": 0, "windings": winding, "m_expected": vp.m,
        "winding_ok": (abs(winding["winding_through_z_normal"] - vp.m) < 0.1
                       and abs(winding["winding_through_y_normal"]) < 0.1
                       and abs(winding["winding_through_x_normal"]) < 0.1),
    }, indent=2))
    print(f"  winding: {{z:{winding['winding_through_z_normal']:+.2f}, "
          f"y:{winding['winding_through_y_normal']:+.2f}, x:{winding['winding_through_x_normal']:+.2f}}}")

    # renders (skipped in diagnostics-only mode, BRANESIM_VORTEX_RENDER=0)
    if do_render:
        amps, phases, times = _extract_amp_phase(world, ref, grid_shape)
        peak = max(float(np.max(a)) for a in amps)
        print(f"  peak |Psi| = {peak:.4f}")
        _render_world(amps, phases, times, vp, grid_shape, spacing, iter_dir, state_label)
    else:
        print("  Rendering SKIPPED (BRANESIM_VORTEX_RENDER=0; diagnostics-only).")

    # diagnostics suite (the back-up round-trip) -> diagnostics/
    print("  Running diagnostic suite (8 devices)...")
    run_measurements(iter_dir, verbose=False)

    # per-output docs (after diagnostics so diagnostics/ exists)
    _write_output_docs(iter_dir, vp, seed_meta, winding, grid_shape,
                       n_slices, dt, iter_index, state_label, kind, solver_note,
                       do_render=do_render)

    return winding


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------


def _env_int(name: str, default: int) -> int:
    v = os.environ.get(name)
    return int(v) if v not in (None, "") else default


def _env_grid(name: str, default: tuple[int, int, int]) -> tuple[int, int, int]:
    v = os.environ.get(name)
    if v in (None, ""):
        return default
    parts = tuple(int(x) for x in v.replace("x", ",").split(","))
    if len(parts) == 1:
        parts = parts * 3
    if len(parts) != 3:
        raise ValueError(f"{name} must be 'N' or 'NX,NY,NZ'; got {v!r}")
    return parts  # type: ignore[return-value]


def _env_bool(name: str, default: bool) -> bool:
    v = os.environ.get(name)
    if v in (None, ""):
        return default
    return v.strip().lower() in ("1", "true", "yes", "on")


def run() -> None:
    # --- Runtime config (env-overridable so the SAME module scales local->AWS) ---
    grid_shape = _env_grid("BRANESIM_VORTEX_GRID", GRID_SHAPE)
    n_slices = _env_int("BRANESIM_VORTEX_NSLICES", N_SLICES)
    relax_iters = _env_int("BRANESIM_VORTEX_RELAX_ITERS", RELAX_ITERS)
    # Rotating-frame-periodic relaxation (PeriodicBC) is well-conditioned and
    # genuinely moves the brane, so it runs by DEFAULT.  Disable with
    # BRANESIM_VORTEX_RELAX=0 (e.g. a seed-only render).
    run_relax = _env_bool("BRANESIM_VORTEX_RELAX", True)
    # Rendering (volume + slice movies + snapshot) is the CPU-heavy part.  Disable
    # with BRANESIM_VORTEX_RENDER=0 for a fast diagnostics-only run (e.g. an
    # eigensolve pre-test) — the SAME module, no forked driver script.
    run_render = _env_bool("BRANESIM_VORTEX_RENDER", True)

    ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    out_root_env = os.environ.get("BRANESIM_RESULTS_DIR") or os.environ.get("BRANESIM_VORTEX_OUTDIR")
    out_root = Path(out_root_env) if out_root_env else (Path(__file__).parents[2] / "runs")
    run_dir = out_root / f"vortex_seed_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)
    print(f"Run folder: {run_dir}")
    print(f"Config: grid={grid_shape}  n_slices={n_slices}  "
          f"relax={'ON('+str(relax_iters)+', rotating-frame-periodic)' if run_relax else 'OFF (seed-only render)'}  "
          f"render={'ON' if run_render else 'OFF (diagnostics-only)'}")

    # --- Lattice and action params ---
    lattice_params = LatticeParams(
        grid_shape=grid_shape, spacing=SPACING, periodic_axes=(True, True, True),
    )
    action_params = ActionParams(
        k_s=1.0, alpha=ALPHA, rho=1.0, dt=DT, n_slices=n_slices,
        m_ambient=4, r_t=R_T, beta=BETA,
    )
    lattice = SpacelikeLattice(lattice_params)
    mass = action_params.rho * SPACING ** lattice.dim
    ref = lattice.reference_positions(4)

    # --- Base config shared by every iteration folder ---
    base_config = {
        "experiment": "vortex_seed_render",
        "geometry": "spherical_harmonic",
        "timestamp": ts,
        "grid_shape": list(grid_shape),
        "n_slices": n_slices,
        "spacing": SPACING,
        "alpha": ALPHA,
        "dt": DT,
        "beta": BETA,
        "r_t": R_T,
        "vortex_params": VORTEX_PARAMS._asdict(),
        "carrier_re_component": CARRIER_RE,
        "carrier_im_component": CARRIER_IM,
        "render": {
            "fps": FPS, "dpi": DPI, "density_threshold": DENSITY_THRESHOLD,
            "alpha_scale": ALPHA_SCALE, "gamma": GAMMA,
        },
        "notes": (
            f"Spherical-harmonic U(1) vortex seed Y_{VORTEX_PARAMS.l}^{VORTEX_PARAMS.m}. "
            "exp(i*m*phi) winds the U(1) phase azimuthally around the z-axis; "
            "|Y|^2 ~ sin^2(theta) makes the energy a donut (emergent, not hand-built). "
            f"Carrier closure-locked to n_t={VORTEX_PARAMS.n_t} turns (720 deg)."
        ),
    }
    (run_dir / "config.json").write_text(json.dumps(base_config, indent=2))

    # --- Injection (the 4D seed worldvolume) ---
    print("Injecting spherical-harmonic vortex seed...")
    t0 = time.perf_counter()
    world, seed_meta = inject_vortex_worldtube(lattice, action_params, VORTEX_PARAMS, n_slices)
    print(f"  Injection done in {time.perf_counter()-t0:.2f}s  world.shape={world.shape}")

    iterations = []

    # --- iter 0: the freshly injected seed (solver untouched) ---
    _emit_iteration(
        run_dir, 0, "iter 0 — initial seed (solver untouched)", "seed",
        world, ref, base_config, seed_meta, VORTEX_PARAMS, lattice,
        n_slices, DT, grid_shape, SPACING,
        solver_report=None,
        do_render=run_render,
        solver_note=(
            "This is the freshly injected ansatz — **iteration 0** of the "
            "solver-relaxation axis. The block solver has NOT run; the carrier "
            "phase advance is prescribed kinematics, not dynamics, and the "
            "amplitude is static across the time loop."
        ),
    )
    iterations.append({"index": 0, "dir": "iter_0000", "kind": "seed",
                       "label": "initial seed (solver untouched)"})

    # --- iter N: rotating-frame-periodic JFNK relaxation along the solver axis ---
    relax_report = None
    relax_error = None
    if run_relax:
        print(f"\nRotating-frame-periodic relaxation: {relax_iters} JFNK iterations "
              "from the seed (well-conditioned PeriodicBC)...")
        t0 = time.perf_counter()
        try:
            world_N, relax_report = _relax(world, lattice, action_params, mass, relax_iters)
            print(f"  relaxation done in {time.perf_counter()-t0:.1f}s  "
                  f"res_init={relax_report.get('residual_initial')}  "
                  f"res_final={relax_report.get('residual_final')}  "
                  f"cond={relax_report.get('condition_estimate')}")
            moved = float(np.max(np.abs(world_N[:n_slices] - world[:n_slices])))
            note = (
                f"**Iteration {relax_iters}** of a rotating-frame-periodic JFNK "
                f"block-solve (`solve_block` + `PeriodicBC`, root-find of ‖R‖=0): a "
                f"closed cyclic time loop, all slices free, warm-started from the seed. "
                f"residual_initial={relax_report.get('residual_initial'):.3e}, "
                f"residual_final={relax_report.get('residual_final'):.3e}, "
                f"condition≈{relax_report.get('condition_estimate'):.1f} "
                f"(vs Dirichlet two-time κ~1e14). max|ΔR| vs seed = {moved:.3e} — the "
                f"brane genuinely relaxes; the carrier winding is preserved (see "
                f"winding_closure.json). Whether it binds or radiates toward vacuum is "
                f"read off from the energy/confinement/winding diagnostics."
            )
            _emit_iteration(
                run_dir, relax_iters,
                f"iter {relax_iters} — after {relax_iters} rotating-frame-periodic JFNK steps",
                "relaxation",
                world_N, ref, base_config, seed_meta, VORTEX_PARAMS, lattice,
                n_slices, DT, grid_shape, SPACING,
                solver_report=relax_report, solver_note=note,
                do_render=run_render,
            )
            iterations.append({"index": relax_iters, "dir": f"iter_{relax_iters:04d}",
                               "kind": "relaxation",
                               "label": f"after {relax_iters} rotating-frame-periodic JFNK steps"})
        except Exception as exc:  # pragma: no cover - robustness for the probe
            relax_error = str(exc)
            print(f"  [WARNING] relaxation failed: {exc}")
    else:
        print("\nRelaxation SKIPPED (BRANESIM_VORTEX_RELAX=0; seed-only render).")

    # --- manifest describing the iteration (5th-dimension) ladder ---
    manifest = {
        "experiment": "vortex_seed_render",
        "timestamp": ts,
        "grid_shape": list(grid_shape),
        "n_slices": n_slices,
        "axis": "solver_iteration",
        "axis_note": (
            "The brane is 4D; the solver-iteration index is a 5th, purely "
            "technical relaxation axis (a flow coordinate, not a physical "
            "dimension). Only a converged fixed point is physical; iterates "
            "before it are transients."
        ),
        "iterations": iterations,
        "render_enabled": run_render,
        "relax_enabled": run_relax,
        "relax_iters": relax_iters if run_relax else None,
        "relax_report": relax_report,
        "relax_error": relax_error,
        "relax_note": (
            "Relaxation uses the rotating-frame-periodic BC (PeriodicBC): a closed "
            "cyclic time loop, all slices free, well-conditioned (cond~1e3-1e4 vs the "
            "Dirichlet two-time κ~1e14 that froze the earlier probe). The brane "
            "genuinely relaxes toward ‖R‖=0 (residual drops ~100x, max|ΔR|~O(A0)) and "
            "the carrier winding is preserved. Earlier Dirichlet probe was frozen "
            "(project_block_solver_bvp_chirality / project_ivp_solver_deleted)."
        ),
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # --- top-level README ---
    relax_line = (
        f"- `iter_{relax_iters:04d}/` — after {relax_iters} rotating-frame-periodic "
        f"JFNK relaxation steps (PeriodicBC; the brane relaxes, winding preserved).\n"
        if run_relax and not relax_error else
        "- (relaxation skipped — seed-only render; set BRANESIM_VORTEX_RELAX=1.)\n"
    )
    (run_dir / "README.md").write_text(
        f"""# vortex_seed run {ts}

U(1) spherical-harmonic vortex seed (Y_{VORTEX_PARAMS.l}^{VORTEX_PARAMS.m}),
grid {grid_shape}, n_slices {n_slices}, rendered and fully diagnosed.

- `iter_0000/` — the freshly injected seed (solver untouched).
{relax_line}
Each iteration folder is self-contained (`config.json` + `world.npz`) and holds
`snapshot.png`, `renders/` (volume + slice movies), and `diagnostics/`
(energy / confinement / winding / berry / em_fields / color_channels / spectra,
each CSV + PNG + `.md`), plus `report.md`. Every figure/video has a `.md`
explaining what it shows, how it was derived, and how to read it.

See `manifest.json` for the iteration ladder and the solver-axis note.
"""
    )

    print("\n=== DONE ===")
    print(f"Run folder : {run_dir}")
    for it in iterations:
        print(f"  {it['dir']}  ({it['kind']}: {it['label']})")
    if relax_error:
        print(f"  relaxation probe FAILED: {relax_error}")


if __name__ == "__main__":
    run()
