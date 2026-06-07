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

## Extension — the U(1)↔SU(3) binding probe (OPEN_PROBLEMS D6)
Why a proton does not split into a color core + a free positive charge. EM/charge =
trace `(1,1,1)/√3` phase, color = internal orientation, of the **one** `Ψ∈ℂ³` lump.
The probe measures whether the trace (charge) content stays co-located with the
traceless (color) content. Derivation + falsifiers:
`paper/derivations/{u1_su3_binding,time_link_binding}.md`.

**BUILT (2026-06-07): device D8 `binding_probe`** (`branesim/diagnostics/binding_probe.py`,
registered in `run_measurements.py`; audit-clean; CSV+PNG+report stanza). Read-only,
single-worldvolume. Five probes (P1 sector centroids+separation `d(l)`; P2 longitudinal
stretch `Δu_∥(ρ)` sign+profile; P3 per-slice carrier rate `ω(l)` vs closure-locked; P4
trace-sector loop holonomy `γ_Γ`; P5 antisymmetric Kähler ratio `R_kähler=|Im𝒢|/|Re𝒢|`).
Pre-tested on the seed (`runs/vortex_seed_2026-06-06_190350/`): **all signals correctly
null** — `sep_d~1e-14`, `Δu_∥~0`, `R_kähler~1e-18`, `γ_Γ=−12.565≈−2π·n_t`, interior
`ω≈ω_ref`. The seed is static; the binding signals (`Δu_∥<0`, coexcited SU(3), nonzero
`R_kähler`) appear **only on a converged PeriodicBC worldtube** — so the device is ready
and the next step is the eigen-solve, then re-run D8.

- **Sub-probe 1 — co-location baseline (single `α=0.7` set).** On the relaxed
  eigenstate, compute the energy-weighted **sector centroids** `c_tr`, `c_⊥` and their
  separation `d(t)` over the worldtube. Stays co-located ⇒ binding; drifts apart ⇒ no.
- **Sub-probe 2 — forced-separation force law (Channel B; needs authorized `α`-ladder
  `{0.2,0.5,0.7,0.8}`).** Rigidly displace the two sector profiles by `d` in the seed,
  re-relax with the shift constrained (diagnostic energy scan, no back-reaction), read
  `E(d)`, finite-difference `F(d)=−dE/dd`. **Predict `F(d)∝α`** (ratios `α/0.5`),
  linear-in-`d` for `d≪w`, peak near `d≈√2w`. ⚠ Watch for the **cubic** vortex-core
  term (`F∝αΔu_∥`) dominating the quartic.
- **Sub-probe 3 — charge-strip (Channel C.2 GW, single `α`).** Add a pure-gauge trace
  phase gradient, re-relax. `Q_{U(1)}−Q_matter` snaps back to a quantized `c₁B` ⇒
  **topological** binding (a confining EM–color flux tube on forced separation);
  continuous drift ⇒ **energetic-only**. (Gated on a stable `B≠0` texture, C2.)
- **Sub-probe 4 — triality lock (Channel C.1, dynamics-free).** Three seeds — `B=3`
  color-singlet, single color-triplet (`t=1`), octet (`t=0`) — measure `Q_{U(1)} mod 1`
  via the trace holonomy. Predict singlet/octet → `0`, triplet → `±⅓`. Deviation of
  singlet/octet from integer **falsifies that the field is genuinely `U(3)`**.

**Read-only (built in D8) vs multi-solve (needs authorization):** D8 computes, per
worldvolume, the *per-run* quantities — sector separation `d`, `Δu_∥(ρ)`, `ω(l)`, `γ_Γ`,
`R_kähler`. The discriminating *falsifiers* are multi-solve campaigns (each effectively a
sweep, so owner-authorized): the force-law `F(d)∝α` (constrained re-relaxation over `d`
and the α-ladder), `Δu_∥∝ω²` (closure-index `n` ladder), `γ_Γ∝B¹` at fixed amplitude
(texture-winding ladder), and the triality `Q mod 1` over three seeds. D8 emits exactly
the scalars these campaigns collect; the campaign driver is not built.

## Build status (2026-06-06)
- **Injection layer** (`branesim/initialization/vortex_worldtube.py`): single vortex
  **ring** (smoke-ring core, donut/torus energy), periodic-clamp-consistent
  (contractible → net winding machine-zero). Renders as one torus.
- **Renderers** (`branesim/visualization/volume_render.py`, recovered from `c0f1aaa7`):
  multi-color volume movie (opacity = energy, hue = `U(1)` phase) + 2D slice movies.
- **Measurement suite** (`branesim/diagnostics/run_measurements.py`): 8 devices
  (energy/consistency, confinement, winding, Berry/phase, EM `A_μ`/`F_μν`,
  per-color-channel `SU(3)`/QCD, spectra, **D8 binding_probe**) → CSV + paper-ready PNG
  + `report.md`. Shared plot helpers in `diagnostics/_plot_helpers.py`.
- **Run folder**: `runs/<exp>_<date>_<time>/` with config, worldvolume, diagnostics, renders.
- **DONE (2026-06-07)**: (i) re-inject in the **trace (EM) direction** — the bare
  seed is now pure U(1)/EM (`u1_fraction → 1`; E1 resolved, all diagnostics route
  through `project_carrier_re`). Seed object resized ×2 (`r0 6→12a`, `w 2.5→4a`,
  the donut shell now ~10 cells across) on the unchanged 48³ brane; the time
  direction (`n_t`, `n_slices`, `dt`) is **not** scaled.
- **Next**: (ii) add the **tumble** (spin) + set the **`ωT` quantization** (BVP
  eigenvalue); (iii) the eigen-solve (`solve_block`, periodic +
  rotating-frame-periodic BC); (iv) local `48³×64` pre-test; (v) AWS production.

## Single experimental setup (no forked drivers)
There is **one** experiment module — `branesim/experiments/vortex_seed_render.py` —
parameterized entirely by env vars, so the **same** code runs locally and on AWS and
collects every result under `runs/`. Knobs: `BRANESIM_VORTEX_GRID`, `_NSLICES`,
`_RELAX` (`0`=seed-only), `_RELAX_ITERS`, `_RENDER` (`0`=fast diagnostics-only, e.g. an
eigensolve pre-test — replaces the old ad-hoc `_binding_eigensolve_pretest.py` fork),
`BRANESIM_RESULTS_DIR`. AWS recipe: `orchestration/aws/RUNBOOK.md` (the single runbook);
launcher default remote-command runs this experiment; results sync to S3 bucket
`branesim-runs-493652700851` and are fetched back into `runs/`. Origin is marked by the
run-dir suffix (`_aws` / `_local`). The retracted breather sweep lives in `archive/`
(experiment + runbook); only `branesim/solver/breather.py` is kept as infrastructure.
