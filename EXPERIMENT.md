# Current Experiment — the U(1) carrier-phase vortex (EM / electron-like object)

The single, fully-instrumented experiment we are building. Working mode: one
experiment at a time, one parameter set, **all** ingredients and **all** measurement
devices in place, pre-tested locally then run at scale on owner-operated AWS; no
parameter sweeps without explicit authorization. (See `LESSONS_LEARNED.md`.)

## Target object
A **protected `U(1)` carrier-phase vortex** — the electromagnetic / electron-like
sector of the `U(3)` substrate.
- **3D slice:** a donut (torus) of energy around a central **axis line** = the
  topological defect (the donut hole). The `U(1)` phase winds `2πm` azimuthally
  around the axis; the amplitude → 0 *on* the axis (regularizing the singular phase)
  and peaks off it → the donut.
- **4D:** a `U(1)` vortex is **codimension-2**, so the defect is a **2D worldsheet**
  (the axis swept over the time loop). The 3D time-slice of that worldsheet is the
  donut axis. ("Extend vortex to 4D" = a worldsheet, not a line.)
- **Two distinct holonomies:**
  - *spatial tumble* — the axis direction `n̂(t)` reorienting over the loop → a
    **twisted worldsheet**; a `2π` tumble is an `SO(3)` rotation → `ℤ₂`
    Finkelstein–Rubinstein holonomy → **spin-½**.
  - *temporal carrier* — the `U(1)` phase advancing `ωt` → the **geometric (Berry)
    phase / EM charge**.
- **Deterministic extended structure** — no probabilistic "point electron"; its size
  and shape are physical.

## Why this object, and the U(1)/SU(3) split (BACKBONE #25)
The carrier-complexified displacement is a `U(3)` field, `U(3)=U(1)×SU(3)`, split by
the rest length `α`. Topology differs sharply: `π₁(U(1))=ℤ` → line vortices (this
object); `π₁(SU(3))=0` (no vortices), `π₃(SU(3))=ℤ` → textures (the color/baryon — a
**different** object, the `π₃` Skyrmion). **Protection ≠ participation:** the winding
is held by `U(1)`, but the `SU(3)` sector is dynamically coexcited and may be
load-bearing for stabilization even at *net-zero color* (a color-singlet electron can
still be `SU(3)`-stabilized). So we do **not** freeze `SU(3)` — the full `U(3)` field
is dynamical, `U(1)` is the topological seed, `SU(3)` relaxes freely and is measured.

## Injection ansatz (soliton layer → substrate; the inverse of the diagnostic layer)
Order parameter `Ψ(x,t) = A(ρ)·exp(i[ m·χ + ω·t ])`:
- `m=1` winding, `ω` carrier rate, `χ` = angle around the axis, `ρ` = distance from it.
- donut profile `A(ρ) = A₀ (ρ/w) e^{−ρ²/2w²}` (0 on axis, peak `~w`).
- axis `n̂(t)` tumbles over the loop (the twisted worldsheet → spin).

Translate down to the substrate (`R = R_vacuum + u`) with the carrier in the
**trace direction** — `Ψ` populates the **symmetric `(1,1,1)` lateral direction**
(all three lateral components equally, sharing the winding phase), so the object is a
pure **EM/`U(1)`** vortex. (Construction check 2026-06-06: writing the phase into a
*single* lateral component is **not** pure EM — it reads `⅓` trace + `⅔` traceless
(color) by the `U(3)` projection. The bare seed populates only the trace; the full
`U(3)` field is then **free to relax**, so `SU(3)` can coexcite — which the per-color
diagnostic measures.)

**Phase must close two ways:** spatially (`2π`-single-valued winding, integer `m`) and
around the time loop (`ωT + tumble = 2π·integer`). So **`T` is the existence /
quantization condition**, not a free knob; spin-½ implies the `4π` subtlety.
*Open detail to pin before building:* the concrete periodic-consistent worldsheet
(axis wrapping a periodic spatial direction + the tumble providing closure) so the
phase is single-valued — else fall back to a contractible config.

## The one parameter set
`α=0.7`, `dt=0.25`, `β=1 → r_t=0.175` (prestressed substrate); `m=1`; `ω≈2`
(`ω·dt≈0.5 rad/step`, set jointly with the `T`-closure condition); `w≈3a` (donut
width); axis/donut geometry `≈6–8a`; `A₀ ~ O(0.3)` peak strain (tuned in the
pre-test, not swept). **Periodic** spatial BC (clamps the prestress);
**rotating-frame-periodic** temporal BC (closed worldtube loop).

## Solve — eigenstate, not seed-and-watch
Initialize the **full 4D worldvolume** with the ansatz and root-find `𝓡=0` with
`solve_block` over the periodic loop. The forward chiral march is kept only as a cheap
"what does it do dynamically" probe.

## Measurement suite (all devices in place — verify the base model)
Per-run folder `runs/<experiment>_<YYYY-MM-DD>_<HHMMSS>/` with `config.json`,
`manifest.json`, `diagnostics/` (CSV **and** paper-ready PNG per device), `renders/`
(multi-color volume movie + 2D slice movies xy/xz/yz), `report.md` (verdict +
figures). Devices — recover/adapt the legacy implementations from commit `c0f1aaa7`
(`volume_render.py`, `diagnostics/berry.py` incl. phase→RGB videos, `dispersion.py`,
`christoffel_6nn.py`, the viz driver):
- **Energy & consistency:** total `E`, vacuum-subtracted `E_excess`, kinetic/potential
  split; conservation over slices, `‖𝓡‖/DOF`, vacuum stability (`A=0 → E_excess≡0`).
- **Dispersion & spectra:** `ω(k)`, `c_T/c_L=√(1−α)`, anisotropy vs `christoffel_6nn`;
  spatial-FFT energy/mode spectrum (radiation tail).
- **Gauge layer:** Berry connection (`Ψ=u+iv/ω₀`) + phase-RGB videos; EM
  `A_μ=i⟨u|∂_μu⟩`, `F_μν → E,B` quiver/streamlines; **per-color-channel `SU(3)`/QCD
  breakdown** (does color coexcite/stabilize?). **Screening length / photon-mass check
  (OPEN_PROBLEMS D5):** fit the `B`/supercurrent radial tail — a power-law (`λ→∞`,
  unscreened) tail confirms a *massless* long-range EM photon; exponential screening
  (finite penetration depth `λ`) would mean the substrate Meissner-screened its own
  `U(1)` (a `W/Z`-like, not EM, object) and **falsifies the EM identification**.
- **Confinement & topology:** `spread_ratio`/`radius_rms`/leakage; winding via the
  **real signed-solid-angle degree** (never the hard-coded `F_inner=π` estimator).
- **Volume render:** opacity = excess-energy density, hue = `U(1)` phase (multi-color).

## Scale ladder (same parameters)
Local pre-test `48³×64` (validate injection + solve + measurement pipeline) → AWS
production `96³×128`.

## What this experiment answers
1. Does the substrate **bind** the `U(1)` vortex worldtube — which, with `SU(3)`
   unfrozen, is a **semilocal vortex** (the full vacuum is simply connected, so it is
   *not* topologically protected)? Binding is therefore a **dynamical** condition (the
   `β<1` / hardening-bound-mode regime, BACKBONE #25), on a clean periodic-clamped,
   `r_t>0`, eigenstate footing.
2. Does the **`SU(3)` sector coexcite and stabilize** the color-neutral electron — the
   non-Abelian/semilocal-vortex core moduli (BACKBONE #25; OPEN_PROBLEMS D4)?
3. The **geometric (carrier) phase** and the **spin-½ `ℤ₂`** holonomy, measured directly.

## Build status (2026-06-06)
- **Injection layer** (`branesim/initialization/vortex_worldtube.py`): single vortex
  **ring** (smoke-ring core, donut/torus energy), periodic-clamp-consistent
  (contractible → net winding machine-zero). Renders as one torus.
- **Renderers** (`branesim/visualization/volume_render.py`, recovered from `c0f1aaa7`):
  multi-color volume movie (opacity = energy, hue = `U(1)` phase) + 2D slice movies.
- **Measurement suite** (`branesim/diagnostics/run_measurements.py`): 7 devices
  (energy/consistency, confinement, winding, Berry/phase, EM `A_μ`/`F_μν`,
  per-color-channel `SU(3)`/QCD, spectra) → CSV + paper-ready PNG + `report.md`.
- **Run folder**: `runs/<exp>_<date>_<time>/` with config, worldvolume, diagnostics, renders.
- **Next**: (i) re-inject in the **trace (EM) direction**; (ii) add the **tumble**
  (spin) + set the **`ωT` quantization** (BVP eigenvalue); (iii) the eigen-solve
  (`solve_block`, periodic + rotating-frame-periodic BC); (iv) local `48³×64` pre-test;
  (v) AWS production (new generic bucket + owner-run scripts).
