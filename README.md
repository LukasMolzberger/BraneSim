# BraneSim

A numerical laboratory for a **substrate hypothesis**: that relativity, gauge
fields, gravity, and particles are all *emergent* features of one classical
object — a 4D elastic brane lattice embedded in 4D Euclidean ambient space. There
are no fundamental fields beyond the brane; everything else (the metric, the EM
four-potential, color, mass) is a diagnostic read off the substrate's geometry.

This repository holds the theory (a LaTeX paper + a curated set of design
documents), a Python solver package (`branesim/`), and the validation/experiment
machinery.

> **Status:** active research code. The physics is a working hypothesis under
> test, not an established result. See `OPEN_PROBLEMS.md` for what is *not* yet
> derived and the paper's §2 "Non-claims" for the honest scope.

---

## The idea in one paragraph

The fundamental object is a 4D world-volume that is a stationary point of an
elastic **brane action** `S[R]`. One of the four lattice directions plays the
role of time because the action carries Lorentzian sign structure (the kinetic
term enters `+`, the spring potential `−`); the ambient itself is fully
symmetric. Long-wavelength waves give emergent relativistic kinematics; the
complex carrier envelope of a band-isolated wavepacket carries a `U(3) = U(1) ×
SU(3)` gauge structure (EM + color); and localized, non-radiating solitons are
the candidate particles. The single prestress parameter `α := rest_length /
spacing` (default `0.2`) is the one physically-meaningful dial.

The canonical statement of all of this is **`BACKBONE.md`** — read that first.

---

## Repository map

### Documents (the theory and its guardrails)
| File | Role |
|---|---|
| `BACKBONE.md` | **Canonical, non-negotiable backbone of the theory.** Start here. |
| `PRINCIPLES.md` | Non-negotiable engineering/physics rules (substrate-only, no back-reaction, no hand clamps, layer separation). |
| `ARCHITECTURE.md` | Block-solver-centric code blueprint and design decisions. |
| `OPEN_PROBLEMS.md` | Index of unsolved derivations; full entries now live in each paper's `derivations/<bridge>/status.md`. |
| `LESSONS_LEARNED.md` | Mistakes not to repeat + results we trust. |
| `EXPERIMENT.md` | Spec of the current single instrumented experiment (the U(1) carrier-phase vortex). |
| `archive/VALIDATION_ROADMAP.md` | Sprint-organized validation subtasks (linear → gauge → Lorentz → solitons → gravity). |
| `archive/BARYON_SIMULATION_ROADMAP.md` | Soliton/baryon search program and ansatz menu. |
| `DEPLOYMENT.md` | Running large block solves on AWS (memory sizing, cost-safe scaffolding). |
| `paper/` | The LaTeX manuscript (`paper.tex` master) + `paper/derivations/` (math notes). |

### Code (`branesim/` package)
```
branesim/
  core/            # dimension-agnostic physics primitives (no I/O)
    lattice.py       #   4D lattice topology (6-neighbor axial spacelike + temporal)
    action.py        #   energies V, T, the action S, and the spacelike spring force
    residual.py      #   𝓡 = m·∂_τ²R − F  (the shared primitive; matrix-free; routes on r_t)
    conventions.py   #   α, units, light-cone helpers
  solver/
    ivp.py           #   forward Störmer–Verlet march (the r_t=0 / Cauchy special case)
    bvp.py           #   block root-find of 𝓡=0 over a 4D world-volume (JFNK; never minimizes S)
    boundary.py      #   chiral / two-time boundary conditions
    breather.py      #   time-periodic eigen-BVP (soliton search vehicle)
  initialization/    # boundary-data / seed generators (seeds.py, vortex_worldtube.py)
  diagnostics/       # read-only measurements (energy, confinement, winding, Berry, EM, color, spectra)
  visualization/     # volume + slice movie renderers
  io/                # versioned file contracts (the inter-component API)
  run_experiment.py  # config-driven entry point (writes worldvolume.zip + summary.json)
orchestration/       # pipeline driver + AWS launch/watch scaffolding + JSON configs
tests/               # pytest suite
```

**Layer separation is a hard rule** (`PRINCIPLES.md`): `core/`, `solver/`, and
`initialization/` never import `diagnostics/`, `visualization/`, or
`experiments/`; diagnostics are read-only and never feed forces back into the
solver.

---

## Install

Requires Python ≥ 3.10.

```bash
pip install -e .            # core (numpy + scipy)
pip install -e ".[viz,dev]" # + matplotlib renderers + pytest
```

## Run

The entry point is config-driven:

```bash
python -m branesim.run_experiment \
  --config orchestration/configs/branesim_ivp_smoke.json \
  --output-dir out/

# or, after install, the console script:
branesim-run --config orchestration/configs/branesim_bvp_dirichlet.json --output-dir out/
```

Outputs per run: `worldvolume.zip` (solved slices + `manifest.json` with the
solver report) and `summary.json`. For large block solves on AWS, see
`DEPLOYMENT.md`.

## Test

```bash
pytest            # 54 tests (core physics + BVP solver)
```

---

## Conventions (must match the code and derivations)

- **Prestress:** `α := rest_length / spacing`; `α = 1` is no prestress, `α = 0`
  is maximum prestress; **default `α = 0.2`**. A *single* `α` governs all four
  lattice directions (no spatial subscript).
- **Substrate:** one 4D-isotropic central-force spring lattice, parameterized by
  the temporal rest length `r_t`. `r_t = 0` is the linear/Verlet limit (the safe
  default); `r_t = α·β·Δt` is the prestressed canonical substrate.
- **Lattice:** 6-neighbor axial-only cubic (no diagonal shells); lab-frame cubic
  anisotropy is real and load-bearing, not retuned away.
- **Units:** dimensionless `k_s = a = ρ = 1`; light-cone `c_L² = k_s a²/m`,
  `c_T² = (1−α) k_s a²/m`.
- **The action is a saddle** (Lorentzian, unbounded below): the foundational
  solver **root-finds `∇S = 0`**, it does **not** minimize `S`.

---

## Where to start reading

1. `BACKBONE.md` — the theory, as 25 numbered non-negotiables.
2. `PRINCIPLES.md` — the rules any code change must respect.
3. `ARCHITECTURE.md` — how the solver is structured and why.
4. `paper/` — the full written argument.
5. `OPEN_PROBLEMS.md` — what is still missing.