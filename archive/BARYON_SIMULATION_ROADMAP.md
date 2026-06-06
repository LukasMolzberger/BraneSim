# BARYON_SIMULATION_ROADMAP.md

> **ARCHIVED / historical** — see `archive/README.md`. Superseded by the
> post-retraction program (`OPEN_PROBLEMS.md` C2/D4, `EXPERIMENT.md`). Kept for
> the VSH ansatz menu (Candidates 1–5), still cited by `BACKBONE.md` #20.

## Goal

Adapt the current BraneSim code so that it can search for stable baryon-like localized modes on the cubic substrate,
with proton/neutron configurations treated as the first decisive particle target.

> **Solver vehicle (2026-06-04).** The bound-particle search uses the
> TIME-SYMMETRIC time-periodic eigen-BVP `branesim/solver/breather.py::solve_breather(mode="topological")`
> (cyclic `R^P≡R^0` Skyrme common-carrier breather, root-finding `‖ℛ‖=0`), **not**
> the forward Verlet IVP march. IVP reintroduces a causal time-direction (violating
> the time-symmetry stance, OPEN_PROBLEMS A) and cannot find a bound eigenstate, so a
> non-eigenstate seed radiates regardless of parameters. Decisive verdicts:
> `converged` AND `floquet_multipliers` spectral radius ≤ 1 (stable) AND
> `harmonic_resonance_check` radiationless. Harness:
> `test-runs/sprint4b_skyrme_corrected/run_breather_sweep.py`. See
> `[[project_baryon_search_vehicle_breather]]`. IVP/Verlet remains the right tool
> for wave-propagation / dispersion, not bound states.

## What the current code already gets right

- The core force model already contains the **essential geometric nonlinearity**.
  `SpringForceComputer` computes forces from full 4D Euclidean link lengths, so the nonlinear feedback is already
  present through geometry.
- The solver is already dimension-agnostic enough for controlled 3D experiments.
- There is already infrastructure for Berry analysis and for localized electron-like initial conditions.

## Main limitations of the current codebase

1. **Electron-specific physical calibration is hardcoded**
   - `physical_constants.py` hardcodes the electron mass and electron Compton wavelength.
   - `analysis.py` hardcodes `omega0 = 2π c / lambda_C` using the electron scale.
   - Several initializers assume electron/Compton naming and interpretation.

2. **No baryon-specific initializer exists yet**
   - Existing localized initializers are electron-centered or generic photon packets.
   - There is no seed that explicitly targets a three-component axis-aligned confined mode.

3. **Diagnostics are not yet aligned with baryon questions**
   - There is no standard measurement for triplet confinement, color mixing, U(3) holonomy, charge-trace cancellation,
     or anisotropy of the localized mode itself.

## Required software changes

### Phase 1 — Make the physical scale species-generic

Create a new species/scale object, for example `ParticleScale` or `CarrierScale`, containing:
- name
- rest mass
- target rest energy
- characteristic carrier wavelength
- characteristic angular frequency
- optional target radius

Then:
- replace direct electron assumptions in `physical_constants.py` with a generic species API,
- make `analysis.py` accept `omega0` from the experiment config instead of hardcoding the electron Compton value,
- rename electron-specific helper names where they are really species-generic.

### Phase 2 — Add a partial-wave-organized baryon initializer

Create a new module such as `baryon_initialization.py` providing seed families
organized by their (J, L) vector-spherical-harmonic content
(see backbone #20 and paper §6.3, `sec:soliton-labels`).
Each seed family corresponds to a representation of the emergent SO(3) at
soliton scale; on the cubic substrate these descend to O_h irreps with small
splittings at sub-leading order in `a / R_seed`.

Every initializer must declare its (J, L) content explicitly in its
docstring or metadata so that diagnostics can be projected onto the correct
channel.

#### Candidate 1 — Hedgehog (J=0, L=1)

```
ξⁱ(x) = f(r) · x̂ⁱ        with f(0) = A, f(∞) = 0
```

- Color index `i` locked to spatial angular direction `x̂ⁱ`.
- Entire structure lives in the SU(3) traceless sector.
- U(1) trace averages to zero angularly → trace-neutral far field.
- No topological protection (winding 0) — expected to radiate unless
  geometric quartic alone suffices.
- Sweep parameters: amplitude `A`, profile width `w`, prestress `α`,
  profile shape (Gaussian / sech / power-tail).

#### Candidate 2 — Skyrme-twisted hedgehog (J=0, winding B=1)

```
ξⁱ(x)  = A · x̂ⁱ · sin F(r)
δX⁴(x) = A · cos F(r)
F(0) = π,  F(∞) = 0
```

- Topologically protected by π₃ winding number `B=1` on the
  (lateral triplet + X⁴) → S³ map.
- Uses the gravity channel `X⁴` (backbone #19) as the fourth direction the
  field wraps; topological winding and gravity channel share a structure.
- This is the most likely actually-stable baryon candidate per the
  topological-stability requirement of backbone #17.
- Sweep parameters: amplitude `A`, characteristic width `w` (defined by
  the radius at which `F(r) = π/2`), prestress `α`, profile shape.
- **Corrected seed target (2026-06-04).** The Derrick balance with the
  lattice-exact quartic (`W_4 ∝ k_s α/a`) gives `R_h/a = κ(A/a)√(α/(1−α))`,
  `κ = O(1)` (see `paper/derivations/vsh_channel_decomposition.md` §2.5a and
  `geometric_nonlinearity_alpha_scaling.md`). Radius is **linear in amplitude `A`**
  and **increasing in α** — the opposite of the old §2.5 `(ℓ₀/u₀)^{2/3}` target.
  Seed at **α = 0.7, A/a ≈ 10** (⇒ `R_h/a ≈ 8`); bracket **α ∈ {0.5, 0.7, 0.8}** at
  fixed `A/a ≈ 10` to trace the falsifiable curve `R_h(α)/R_h(0.5) = √(α/(1−α))`
  (0.65 / 1.00 / 1.53 / 2.00 at α = 0.3 / 0.5 / 0.7 / 0.8). NB the May-2026 sprint4
  sweep used tiny `u₀ ∼ 0.003–0.025` at α=0.2 — the collapse corner on *both* axes,
  which is why it dispersed. Use a box ≫ seed (periodicity masked dispersion).

#### Candidate 3 — Hedgehog with U(1) trace admixture (J=0, L=1 + L=0)

```
ξⁱ(x) = f₁(r) · x̂ⁱ  +  (1/√3) · (1,1,1)ⁱ · f₀(r)
```

- Adds a scalar trace component on top of the hedgehog.
- `f₀` gives nonzero U(1) far field → proton-like.
- Can be combined with the Skyrme twist of Candidate 2 to give a
  topologically stable proton-class seed.
- Sweep parameters: trace fraction `f₀(0) / f₁(0)`, amplitudes, widths,
  prestress.

#### Candidate 4 — Axis-triplet (negative control, J=1, L=0)

```
ξⁱ(x) = δⁱ_k · A_k · g(r)   for k ∈ {1,2,3} with independent weights A_k
```

- Three independent axis-polarized scalars.
- **No color-spatial locking; not a baryon ansatz.**
- Expected to radiate. Included only as a negative control to confirm
  that the search distinguishes partial-wave-locked configurations from
  unlocked ones. If candidate 4 confines but candidates 1-3 do not, the
  framework is wrong.

#### Candidate 5 — Hopfion / linked vortex ring (Hopf charge Q_H=1)

```
ξⁱ(x) = A · g(r) · n̂ⁱ(x/w)
```

with `n̂ ∈ S²` the charge-1 Hopf texture, written via the stereographic coordinate
`w_s = (n̂₁ + i n̂₂)/(1 − n̂₃)` of the lateral-triplet *direction*:

```
w_s(x) = 2 w (x + i y) / (2 w z + i (r² − w²)),    r² = x²+y²+z²
```

so `n̂ → ẑ` at infinity; `g(r)` localizes the amplitude (Gaussian / sech).

- **Topology: `π₃(S²) = ℤ`.** Preimages of any two target directions are **linked
  circles** — the "donut + central-axis line" structure. The Hopf charge is the
  linking number.
- **Distinct from Candidate 2.** The Skyrme-twisted hedgehog is `π₃(S³)` and wraps
  into the `X⁴` gravity channel; the Hopfion lives **entirely in the lateral-triplet
  direction field** (no `X⁴` escape required). If the amplitude `g` stays finite the
  texture is smooth (no core singularity); if `g → 0` on the central ring it
  degrades to a genuine vortex-ring `π₁` core.
- Stabilization against Derrick collapse needs a Skyrme-type quartic, which the
  geometric link-norm nonlinearity supplies (`Q = 8 k_s α / a² > 0` hardens —
  `[[project_pythagorean_breather_go]]`). Faddeev–Niemi Hopfions are the canonical
  stable knotted 3D solitons, so this may confine where the C2 Skyrme hedgehog
  disperses (`[[project_c2_skyrme_no_confinement]]`).
- Spin-statistics: under a `2π` rotation the Hopfion phase is tied to `Q_H` and to
  its trace-`U(1)` charge (spin-from-isospin, `OPEN_PROBLEMS.md` D1).
- Sweep parameters: amplitude `A`, ring scale `w`, profile shape, prestress `α`,
  optional internal twist (`Q_H > 1`), trace admixture (charged vs neutral).

This is a working simulation hypothesis, not a derived claim.

#### Working proton/neutron hypothesis

- **proton**: Candidate 3 (or 2+3 combined) — hedgehog + nonzero U(1)
  trace admixture, giving a nonzero far-field trace sector.
- **neutron**: Candidate 1 or Candidate 2 — pure hedgehog or
  Skyrme-twisted hedgehog with trace-cancelled far field.

This is a working simulation hypothesis, not yet a derived claim.

### Phase 3 — Add baryon diagnostics

Add diagnostics for:
- total energy stability,
- effective radius (RMS and percentile-based),
- triplet-component energy fractions,
- inter-axis mixing strength,
- U(1) trace holonomy versus traceless U(3) holonomy,
- radiation leakage rate,
- long-time confinement score,
- anisotropy score of the localized mode itself.

### Phase 4 — Create a dedicated 3D baryon experiment harness

Add a script such as `baryon_3d_experiment.py` that:
- builds a 3D grid with periodic or large-box boundaries,
- initializes a baryon seed,
- evolves it for long times,
- stores sparse checkpoints,
- computes diagnostics automatically,
- ranks parameter sets by confinement quality.

### Phase 5 — Parameter search

Sweep at least:
- prestretch `alpha`,
- neighbor-shell couplings `(k1, k2, k3)`,
- seed radius,
- seed amplitude,
- triplet phase offsets,
- degree of axis mixing,
- boundary conditions,
- damping-free relaxation / annealing preparation methods.

## Recommended order of work

1. species-generic scaling refactor
2. baryon initializer
3. baryon diagnostics
4. dedicated experiment harness
5. coarse parameter sweep
6. refined sweep near stable pockets
7. only then: proton/neutron interpretation sharpening

## Important modeling guardrails

- Do **not** add ad hoc confinement forces.
- Do **not** add artificial nonlinear saturation just to keep the packet together.
- Do **not** hide anisotropy by hand.
- Keep the confinement mechanism geometric and substrate-native.
- Let failure be visible if the substrate cannot support the mode.

## Practical recommendation

Do not begin with a fully calibrated SI-accurate proton run.
Begin with a **dimensionless baryon-mode existence search** on the cubic substrate.
Once a stable triplet mode exists in dimensionless form, perform the physical back-calibration.
