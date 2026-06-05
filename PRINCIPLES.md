# BraneSim Principles (Canonical)

**This file is the single source of truth for project principles.**
Read it before making **any** physics, mapping, measurement, or experiment changes.

---

## Quick Reference Card (Non-Negotiables)

1. **Substrate-Only Evolution**
   Evolve only the microscopic brane DOFs: node positions `R_p(t)` and velocities.
   Forces must come only from the brane's mechanical energy: `F_p = -∂U/∂R_p`.

2. **No Back-Reaction (Diagnostic-Only Emergence)**
   Anything "emergent" (metric, EM fields, potentials, Berry connection, curvature, etc.)
   is **measured** from the substrate state and **must never** be fed back into the solver
   as extra forces, constraints, damping, or special rules.

3. **Pure Geometric Coupling (Default)**
   Amplitude↔lateral coupling must come from geometry (e.g., 4D Euclidean distances /
   induced metric effects), not from an added coupling field or hand-crafted interaction.

4. **No Artificial Cutoffs / Clamps**
   No hard amplitude clamps, no piecewise "if energy > … then …", no imposed thresholds.
   If thresholds exist, they must **emerge** from smooth dynamics and the chosen energy.

5. **Ontology & Interpretation**
   - Fundamental entity: a **4D brane lattice embedded in 4D Euclidean ambient**
     (codimension 0). The lattice IS the substrate. The "amplitude" direction
     of older 3D-in-4D formulations is absorbed into the time direction.
   - **The ambient is symmetric; asymmetry lives in the brane.** One of the
     four lattice directions is picked out as **time** by the brane action's
     Lorentzian sign structure (`T − V` Lagrangian), not by anything intrinsic
     to the ambient.
   - **Causality is a property of solitons, not of the substrate.** The brane
     action is time-symmetric. The arrow of time is selected at the soliton
     level by chirality (matter = forward worldtube, antimatter = backward
     worldtube).
   - **The gauge / gravity split is observer-relative.** The 4-component
     displacement field decomposes into 3 spatial (U(3) gauge: EM + colour)
     + 1 timelike (gravity, via the contraction mechanism) only after the
     inside observer picks a timelike direction. Cosmological boundary
     conditions (Big Bang vertex) make this split globally consistent.
   - Relativity, EM, gravity, particles, quantization are **emergent
     descriptions** seen by the inside observer.
   - The current Verlet pipeline produces Cauchy slices through the 4D
     world-volume; it remains a valid forward-evolution diagnostic, not the
     foundational solver.

6. **Dimensionality-Agnostic Implementation**
   All core implementations (solver, forces, initialization, diagnostics) must work for 1D, 2D, and 3D
   without hard-coded dimension checks or separate code paths. Exceptions: visualization, geometry
   setup, and experiment orchestration may be dimension-specific.

**Quick check before a PR/change:**
- [ ] Only changing substrate evolution or clearly separated diagnostics/initialization?
- [ ] No emergent back-reaction introduced (even indirectly)?
- [ ] No clamps/cutoffs/hand-imposed thresholds?
- [ ] Coupling remains geometric (distance/metric), not an added "field"?
- [ ] Layer boundaries respected (core vs mapping vs measurements vs experiments)?
- [ ] Implementation works for 1D/2D/3D (or is in visualization/geometry/experiments)?

---

## 1. Core Model Commitments

### 1.1 Substrate definition (conceptual)

**Full ontology (foundational):**
- Material coordinates: `x = (x¹, x², x³, x⁴) ∈ Ω ⊂ ℝ⁴`. All four coordinates
  are on equal geometric footing in the ambient. The fourth is the
  **timelike** direction selected by the brane action's Lorentzian sign
  structure (per backbone #21).
- Embedding map: `X(x) ∈ ℝ⁴` — codimension 0. The brane IS the ambient;
  there is no extra perpendicular direction. The amplitude direction of
  older 3D-in-4D formulations is absorbed into the timelike direction.
- Discretization: nodes `p` indexed by 4D lattice coordinate `(i,j,k,l)`,
  positions `R_p ∈ ℝ⁴`.
- The full 4D world-volume is a stationary point of the brane action
  `S[R]`; the foundational solver finds it as a 4D boundary-value problem.

**Working pipeline (Cauchy slices, what the code currently implements):**
- The current implementation treats `x⁴` as the timelike direction
  explicitly and discretizes it as a sequence of 3D slices indexed by an
  integer time step.
- At this slicing, the working state is `{R_p(t), V_p(t)}` with
  `R_p(t) ∈ ℝ⁴` (the 4D ambient position of a node `p` at slice `t`) and
  `V_p(t) ∈ ℝ⁴` (its velocity, which is the discrete time-direction
  difference of `R_p`).
- The Verlet integrator (Substrate-Only Evolution rule) produces a
  forward-evolution diagnostic — a Cauchy slice through the 4D world-volume
  rather than the full stationary configuration.
- The reduction from "node has 4D position `R_p ∈ ℝ⁴`" to
  "node has 3D in-brane position `(x¹, x², x³)` plus a 1D
  amplitude `u`" used in older derivations remains a valid description
  **from the inside observer's frame**, where the chosen timelike direction
  is identified with the amplitude axis (per backbone #22).

The solver evolves the substrate state `{R_p(t), V_p(t)}` only. The 4D
foundational reading is what gives this evolution its physical meaning;
the practical computation is unchanged.

### 1.1a Structural facts about the lattice (must be respected)

These are mathematical consequences of the central-force pair-spring energy
`U_link = ½ k_δ (|R_{p+δ} − R_p| − α a |δ|)²`. They are not policy — they are
identities the project relies on.

**Scope note.** The facts below describe the **spacelike-slice geometry** of
the substrate — the 3D-cubic lattice seen at a fixed time slice by the
inside observer. They remain the correct description of what the current
Verlet pipeline simulates. The full 4D-in-4D ontology (backbone #1, #21,
#22) adds a timelike fourth lattice direction; the connectivity and
elastic structure along that fourth direction (e.g. whether links along
time are explicit lattice bonds with rest-length spacing, or are encoded
implicitly via the kinetic term `½ m v²` in the current Verlet pipeline)
is a separate design question being deferred. None of the facts below
change.

- **Canonical stencil: 6-neighbor axial-only (spacelike slice).** Each node
  connects to its six nearest axial neighbours `±êᵢ` in the three spacelike
  directions. Diagonal-shell bonds (face-diagonal, body-diagonal) are
  intentionally absent. The dynamical matrix `D(k)` is then **diagonal in
  the Cartesian basis** at every `k` and every `α`, with k-independent
  eigenvectors `(ê_x, ê_y, ê_z)`. See
  `paper/derivations/lattice_to_continuum.md` for the algebra and
  `test-runs/sprint2_subtask9_d_of_k_diagonal/` for the structural certificate.
- **α convention** (matches code and derivations): `α := rest_length / spacing`,
  so `α = 1` ↔ no prestress (rest length = held distance) and `α = 0` ↔
  maximum prestress (rest length zero). Default operating point: `α = 0.2`.
- **Cubic anisotropy is a structural feature, not a defect.** On the
  6-neighbor axial-only lattice the long-wavelength dispersion is direction-
  dependent: along `[100]` the speeds are `c_L = 1, c_T = √(1−α)`; along
  `[111]` all three lateral eigenvalues are exactly degenerate at every `α`.
  This anisotropy is acknowledged at the lab-observer level and is the
  structural source of the U(3) gauge sector under the dual-observer
  framework (`paper/backbone.md` #8, #15, #16, #19).
- **Prestress α runs the U(1)³ → U(3) crossover via eigenvalue degeneracy
  (Mechanism ii).** What runs with `α` is the eigenvalue spread of `D(k)`,
  *not* the eigenframe (which is k-independent). At `α = 1` the transverse
  modes are at zero frequency and the three lateral channels are dynamically
  decoupled → closest to `U(1)³`. At `α = 0` the lattice energy is purely
  geometric `(1−α)|Δu|²` and `D(k) ∝ I`: the three lateral channels are
  fully degenerate and the carrier-triplet gauge group is the full
  `U(3) = U(1) × SU(3)`. The default `α = 0.2` sits close to the degenerate
  end (≈ 10.6 % L–T gap along `[100]`, 0 % along `[111]`).
- **Cubic anisotropy provides the *basis* for the gauge sector, not the
  gauge group.** `U(3)` on a 3-dim complex envelope is generic; what the
  cubic anisotropy provides is the preferred lateral subspace `(ξ¹, ξ², ξ³)`,
  the preferred axis-aligned basis within it (giving operational meaning to
  the eight SU(3) traceless generators as "colour" labels), and the natural
  identification of the U(1) trace with the EM-like sector. See backbone #19.
- **Geometric quartic = Skyrme-class stabilization.** The induced-metric
  contribution `∂_i u ∂_j u` gives a `(|∇u|²)²` term, which under Derrick
  scaling `u(x) → u(λx)` scales as `λ^{+1}` in 3D and resists collapse.
  Combined with the `λ^{−1}` quadratic-gradient term, this is the standard
  Skyrme balance and *defeats Derrick's no-go at the continuum level*.
  Lattice spacing `a` provides an independent UV cutoff.
  **Coefficient (lattice-exact, ∝ α).** The quartic prefactor is `∝ k_s α/a`
  (continuum `W_4 = (k_s α/8a)|∂u_⊥|⁴`), derived from the norm term
  `−k_s αa|ΔR|` of the central-force spring — it vanishes at α→0. A naive StVK
  identification `(μ/4ℓ₀⁴)(|∇u|²)²` locks the prefactor to `μ ∝ (1−α)` and has
  the α-scaling **inverted**; do not use it for soliton sizing. The Derrick
  `λ`-scaling above is α-independent and unaffected. See
  `paper/derivations/geometric_nonlinearity_alpha_scaling.md`.
- **Two soliton regimes.** Width `~ a` (lattice-stabilized) → strong
  Peierls–Nabarro pinning, soliton mass tied to the UV scale, no continuum
  limit. Width `≫ a` (Skyrme-stabilized) → exponentially suppressed PN,
  Lorentz-respecting motion, derivable physical scale. **The program targets
  width ≫ a**; experiments must report the width-to-spacing ratio explicitly.

### 1.2 Dynamics (implementation truth)

**Foundational view.** The substrate is defined by a brane action `S[R]`
over the full 4D world-volume. The physical configuration is a stationary
point: `δS/δR = 0`. The action has Lorentzian sign structure (backbone
#21) — time-direction terms enter with opposite sign from spacelike terms,
recovering standard `T − V` Lagrangian mechanics in the inside-observer's
frame.

The explicit discrete form is the 4D-cubic link sum with `6 + 2` links per
node — six spacelike central-force springs entering `−V`, and two temporal
links (`±` in time) entering `+T` — derived in full in
`paper/derivations/discrete_4d_brane_action.md`. Two facts from that
derivation are load-bearing for how the solver must be understood:

- **The Euler–Lagrange stencil IS Störmer–Verlet.** Stationarity at an
  interior node, `δS/δR_p^l = 0`, is term-for-term the leapfrog update
  `m (R^{l+1} − 2R^l + R^{l−1})/Δt² = F^l = −∂V^l/∂R^l`. Verlet is the
  *discrete variational (symplectic) integrator* of this action, not an
  approximation to it. The local stencil is therefore identical whether the
  model is solved as a forward initial-value problem (the working pipeline)
  or as a 4D block boundary-value problem (the foundational reading); only
  the global solution philosophy and boundary conditions differ.
- **The action is a saddle, not a minimum** (Lorentzian ⇒ unbounded
  below). A foundational block solver must **root-find `∇S = 0`**
  (Newton–Krylov), never gradient-descend `S`; relaxation/minimization
  would silently solve the *Euclidean* heat-equation problem instead.
  The two-time block BVP is also not unconditionally well-posed
  (resonances → non-uniqueness). Both are tracked in `OPEN_PROBLEMS.md`
  §A1/§A2, with soliton chirality (§1.5) as the conjectured selection
  principle. Time-symmetric solvers that find the full 4D stationary
  configuration under past+future boundary data are a separate development
  track.

**Working pipeline (what the code currently implements).** The simulation
is driven by a mechanical energy `U` on each 3D Cauchy slice plus the
kinetic term `½ m v²` along the timelike direction. This is the discrete
Lagrangian whose stationarity is the action-stationarity condition above,
specialised to forward time-marching from an initial slice.

**Typical structure (conceptual, not prescriptive):**
- Stretching/edge energy: depends on deviations of neighbor distances from rest lengths
- Optional bending/regularization energy: penalizes sharp angular changes / curvature proxies
- Total energy: `U = U_str + U_bend (+ optional smooth regularizers)`

**Forces:**
- `F_p = -∂U/∂R_p` (computed from the substrate configuration only)

**Integration:**
- Use a stable time integrator (often symplectic / Verlet-family) and track energy.

> If a change introduces a force that is not a gradient of the chosen substrate energy,
> it violates Substrate-Only Evolution.

### 1.3 Emergence hierarchy (how to reason about dependencies)
This project treats "physics we recognize" as **read-outs** / effective descriptions.
A practical dependency order used in debugging and testing:

1. **Linear wave regime & isotropy** (dispersion tests, group speed, stability)
2. **Emergent kinematics / relativity tests** (effective cone structure, Lorentz-like behavior)
3. **Berry / holonomy diagnostics** (phase transport on an isolated eigen-branch)
4. **EM-like diagnostics** (connections/curvatures compared to Maxwell-like structure)
5. **Solitons / particle analogs** (stability, topology, chirality)
6. **Gravity-like diagnostics** (lateral contraction / induced metric interpretations)

Don't "skip levels" by hard-coding behavior at a lower level.

### 1.4 Lab vs inside observer (dual-observer framework)

The substrate has two natural observers. They see different things, and
several project commitments depend on keeping them distinct.

**Lab observer (the one writing the simulation).**
- Sees the 4D Euclidean ambient and the brane lattice within it.
- Sees the brane action's Lorentzian sign structure as a property of the
  action, not of the ambient.
- Sees the lattice's 4D-cubic anisotropy (which, on a spacelike Cauchy
  slice, reduces to the 3D-cubic anisotropy of §1.1a) as a real, measurable
  feature.
- Has access to the full 4D state and can read off any global property.

**Inside observer (an excitation propagating in the brane).**
- Built from substrate excitations; their rods, clocks, and signal cones
  are renormalized coherently with the local wave speed.
- Picks one of the four lattice directions as their **timelike direction**.
  Globally, cosmological boundary conditions (Big Bang vertex) make every
  inside observer pick the same direction; locally, different inside
  observers are related by emergent Lorentz transformations on the
  perpendicular spacelike 3-subspace.
- Sees an **effective Lorentzian metric**, modulated by the contraction
  field sourced by displacements in the timelike direction (the
  gravity-channel mechanism).
- Is conjectured to be blind to the lab-frame cubic anisotropy on their
  spacelike slice (backbone #8). This is the load-bearing assumption for
  emergent Lorentz kinematics; it must be tested, not assumed.
- Sees the displacement 4-vector decompose into **3 spatial channels
  (carrying U(3) gauge: U(1) EM + SU(3) colour) plus 1 timelike channel
  (carrying gravity, via contraction)** — observer-relative but globally
  consistent (backbone #22).

**Information flow.** All "emergent" descriptions — Lorentzian metric,
Maxwell-like fields, Berry/Wilczek–Zee connections, soliton dynamics —
live at the inside-observer level. The lab observer measures them as
diagnostics but does not feed them back into the solver (the
No-Back-Reaction rule, §0 quick-reference card #2).

### 1.5 Bell's theorem and the retrocausal worldtube interpretation

BraneSim is manifestly local and deterministic. Bell's theorem rules out
theories that are simultaneously local, deterministic, and
measurement-independent. The project commits to the **retrocausal
loophole**: the future measurement context propagates back along the
particle's worldtube. Concretely:

- The brane action is time-symmetric (its Lorentzian signature picks out a
  *kind* of direction, not a *direction*). The arrow of time is selected
  at the soliton level by the **chirality of soliton solutions**: matter
  solitons are forward-propagating worldtubes; antimatter solitons are
  backward-propagating worldtubes (Feynman–Stueckelberg-consistent).
- **Entanglement is V-branching of one worldtube.** Spacelike correlations
  between entangled pairs are continuity of one extended 4D object whose
  worldtube branches at a V-vertex in the shared past.
- This places the model in the same interpretive family as Aharonov's
  TSVF, Cramer's transactional interpretation, Price & Wharton's
  retrocausal models, and Sutherland's time-symmetric Bohmian model.

**Open derivations** (not yet established; flagged as future work):
Tsirelson's bound (`2√2`), the no-signalling theorem, and the
baryon-to-photon ratio. These are constraints on what the project *can*
eventually claim, not established results. The live, detailed list lives in
`OPEN_PROBLEMS.md` §B (the central tracker) — keep it there, not in the
paper or the theory-structure diagram.

---

## 2. Strict Separation of Concerns (Layers)

### 2.1 Layers (what goes where)

**A) Core substrate solver (dimensionless / numerical core)**
- State: positions/velocities, neighbor topology, energies, forces, integrator
- Must remain free of "physics interpretations" (no EM/gravity forces)

**B) Mapping / calibration (SI ↔ simulation units)**
- Converts between physical scales and internal dimensionless parameters
- Purpose: numerical conditioning + transparent calibration
- Must not smuggle in emergent forces

**C) Initialization (allowed to be "bidirectional")**
- Builds initial conditions `{R_p(0), V_p(0)}` from higher-level descriptions
  (e.g., a wavepacket, a tubular mode, an EM-like field used as a generator).
- Once `t > 0`, substrate evolution is one-way only.

**D) Measurements / diagnostics (read-only)**
- Computes derived fields from `{R_p(t), V_p(t)}`: induced metric proxies,
  coarse-grained potentials, Berry phases, holonomies, dispersion curves, etc.
- Must be read-only w.r.t. the solver

**E) Experiments**
- Orchestrate: setup → run → measure → export → visualize
- Experiments should not implement "new physics" inside the solver

**F) Visualization**
- Reusable plotting / rendering utilities live here
- No physics logic beyond formatting/aggregation

### 2.2 Information flow rule
**During simulation:**
`substrate state → measurements/plots` only (never the other direction)

**During initialization:**
High-level description → substrate initial conditions is allowed

---

## 3. Nonlinearity Policy

### 3.1 What nonlinearity is allowed
- Geometric nonlinearity from the embedding and distance/metric dependence
- Smooth constitutive nonlinearity inside the substrate energy is allowed
  (if clearly part of the mechanical model and not a clamp)

### 3.2 What is forbidden
- Any piecewise clamp ("cap X⁴", "limit strain to …", "if amplitude > … then …")
- Any manually injected "collapse rule" in the integrator
- Any EM/gravity back-reaction forces computed from diagnostics

If a threshold/quantization/collapse-like behavior exists, it must arise as an emergent,
dynamical consequence of smooth substrate mechanics and coarse-graining.

---

## 4. Electromagnetism & Berry-Phase Program (Current Direction)

### 4.1 What we are *not* doing
- We do not identify "the FFT phase" with a gauge potential.
- We do not hard-code Maxwell forces into the solver.
- We do not assume a single global phase without isolating a mode/eigen-branch.

### 4.2 What we are doing (diagnostic path)
We treat EM-like structure as potentially emerging from **phase transport** on an
isolated, quasi-linear eigen-branch of the brane's small-amplitude dynamics.

Operationally:
1. Create/track a wavepacket that stays mostly within a narrow spectral band
2. Build a **projector / filter** that isolates that band (the "eigen-branch")
3. From the filtered complex mode field `u`, compute Berry objects (diagnostic only):
   - Berry connection (example form): `a_μ = i⟨u | ∂_μ u⟩`
   - Berry curvature: `F = d a` (or discrete analogue)
4. Compare patterns and invariants to EM expectations in vacuum-like regimes

**Critical requirement:** the "band isolation" step is not optional. If we mix modes,
the phase is not gauge-meaningful and Berry diagnostics become unreliable.

### 4.3 Narrowband is local, not global

The complex-envelope description (`ξⁱ ≈ Re[Ψⁱ e^{−iω₀ t}]`) is required to promote
the real lateral triplet `(ξ¹, ξ², ξ³)` to the complex `Ψ ∈ ℂ³` on which `U(3)`
acts. This **does not** require a globally coherent narrowband sector. The
description is per-wavepacket: each band-isolated excitation carries its own
local complex envelope and its own Wilczek–Zee bundle. A global carrier
postulate is not part of this project's commitments.

### 4.4 SU(3) emergence chain (record what is doing the work)

The candidate path from microphysics to color-like internal structure is:

1. Cubic central-force lattice → three real lateral displacement channels.
   The current `1/|δ|²` shell weights do **not** give an exactly isotropic
   leading-order IR tensor; an `SO(3)`-like triplet is only justified after
   shell retuning or after quantitatively showing that residual anisotropy is
   below the tolerance of the regime being used.
2. Prestress `α` changes the acoustic/static tensors and the degree of
   branch mixing through the same substrate energy. Its proposed role in the
   `U(1)³ → U(3)` crossover must be measured through band isolation and
   Wilczek--Zee diagnostics, not assumed from the Cauchy relation alone.
3. Per-wavepacket narrowband ansatz promotes real `(ξ¹, ξ², ξ³)` to complex
   `Ψ ∈ ℂ³`, on which `U(3) = U(1) × SU(3)` acts.
4. The `U(1)` trace sector is the EM-like channel; the `SU(3)` traceless
   sector is the candidate color channel.

What this chain does NOT provide:
- **Color confinement** (asymptotic freedom, area-law Wilson loops, hadron
  spectrum). That is a separate dynamical question; the SU(3) here is
  kinematic gauge structure, not QCD's full nonperturbative content.
- **Discrete cubic group → continuous SU(3) representation theory.** SU(3)
  emerges from the *near-degenerate continuum triplet*, not from `O_h`.
- **Half-integer spin.** This requires spinorial holonomy of a 2-level
  polarization subspace inside the triplet; it does not come for free.

### 4.5 Berry phase outputs expected from experiments
- Local Berry phase along a 1D path (per position / per segment)
- Holonomy on closed loops (2π vs 4π behaviors where relevant)
- Robustness under gauge choices (re-phasing of `u`)

All of this remains **measurement**.

---

## 5. Particle / Soliton Program (W&vdM lineage)

- "Particles" are modeled as candidate **stable localized excitations** (soliton-like),
  potentially with topology (e.g., tubular/toroidal structures).
- Chirality and internal phase structure may encode charge-like behavior.
- Any EM-to-brane mapping used here is an **initialization generator** only.

Again: soliton stability must come from substrate dynamics, not added constraints.

---

## 6. Testing & Validation Expectations

### 6.1 Always export the basics
Each experiment should, at minimum, export:
- energy vs time (K, U, total)
- CFL / stability parameters used
- dispersion/isotropy results when relevant
- the diagnostic quantity under test (e.g., Berry phase profile)

### 6.2 Falsifiability mindset
Prefer tests that could fail clearly:
- isotropy error bounds in linear regime
- Berry phase consistency under gauge changes
- stability / decay rates of localized excitations
- sensitivity to resolution and timestep

---

## 7. Project Hygiene Rules (for Claude Code / contributors)

### 7.1 File placement
- Reusable plotting/rendering → `visualization/`
- Diagnostics / derived field calculators → `measurements/` (or similar)
- Experiment scripts / orchestrators → `experiments/`
- Solver + energies + integrators → core module (keep minimal and clean)
- Mapping/calibration utilities → `mapping/`

(Use the repository's existing structure if it differs; keep the layer boundaries.)

### 7.2 No circular dependencies
- Core must not import experiments/visualization
- Measurements must not modify the solver state
- Experiments may depend on core + measurements + visualization

### 7.3 When changing physics-core code
A change is "physics-core" if it touches:
- energies, forces, neighbor topology, integration, boundary conditions

For any such change:
- document which non-negotiables are relevant
- add/extend at least one test or experiment that would detect regressions

### 7.4 Clean refactoring (no backwards compatibility mess)
When refactoring code:
- Do the refactoring completely or don't do it at all
- No backwards compatibility if-statements or deprecated parameters with warnings
- No conditional logic to support old and new APIs simultaneously
- Make a clean break with the old API

Technical debt from mixed old/new APIs is confusing and error-prone. Clean code is easier to maintain.

### 7.5 Delete obsolete code immediately after refactoring
When creating new implementations that replace old code:
- **Delete the old implementation files immediately** after the refactoring is complete
- Do not leave multiple versions of the same functionality in the codebase
- Ensure experiments and tests use the new implementation, then remove old files
- If backward compatibility wrappers are needed, they should be thin re-exports, not duplicated implementations

**Rationale:** Dead code wastes developer time. Fixing bugs in obsolete code that's no longer used is extremely frustrating. Keep one clear implementation path.

**Checklist after refactoring:**
- [ ] Verified all experiments/tests use new implementation
- [ ] Deleted old implementation files (not just deprecated them)
- [ ] Updated imports in any remaining backward compatibility layers
- [ ] No duplicate implementations of the same functionality exist

### 7.6 Dimensionality-agnostic implementation (detailed)

**Core principle:** The brane can be 1D, 2D, or 3D in its intrinsic dimensionality. All core
infrastructure must handle any of these cases automatically, without dimension-specific code paths.

**What MUST be dimension-agnostic:**
- **Solver core:** forces, energy computation, integration steps
- **Initialization:** carrier compilation, velocity initialization, polarization selection
- **Diagnostics:** Berry phases, holonomy, spectrum analysis, induced metric
- **Grid operations:** neighbor lookups, boundary conditions, gradients
- **State management:** positions, velocities, accelerations updates

**Implementation guidelines:**
- Use `len(grid_shape)` or `state.dimension.value` to determine dimensionality at runtime
- Loop over spatial dimensions: `for axis in range(ndim):`
- Use dynamic slicing: `slice_list = [slice(None)] * ndim; slice_list[axis] = ...`
- Tensor operations should broadcast naturally across dimensions
- Grid reshaping: `field.view(*grid_shape)` adapts to any dimension

**Example (correct):**
```python
# Dimension-agnostic gradient computation
ndim = len(grid_shape)
for axis in range(ndim):
    center = [slice(None)] * ndim
    center[axis] = slice(1, -1)
    grad[tuple(center)] = (field[tuple(plus)] - field[tuple(minus)]) / (2*h)
```

**Example (incorrect):**
```python
# Hard-coded for 3D only - WRONG
grad_x = (field[1:-1, :, :] - field[:-2, :, :]) / h
grad_y = (field[:, 1:-1, :] - field[:, :-2, :]) / h
grad_z = (field[:, :, 1:-1] - field[:, :, :-2]) / h
```

**What MAY be dimension-specific:**
- **Visualization:** plotting 1D vs 2D vs 3D requires different matplotlib calls
- **Geometry setup:** torus knots are inherently 3D, 1D chains are inherently 1D
- **Experiment orchestration:** specific experiments may target a particular dimension
- **Test fixtures:** individual tests can fix dimensionality for clarity

**When dimension-specific code is unavoidable:**
Use clear conditional structure and document why:
```python
if ndim == 1:
    # 1D-specific visualization
    plt.plot(x, field)
elif ndim == 2:
    # 2D-specific visualization
    plt.imshow(field.reshape(nx, ny))
else:  # 3D
    # 3D-specific visualization (slice or 3D plot)
    plt.imshow(field.reshape(nx, ny, nz)[:, :, nz//2])
```

**Validation:**
Before merging code that operates on grids or fields:
- [ ] Test manually with 1D, 2D, and 3D inputs (or add dimension-parametrized tests)
- [ ] No hard-coded axis indices (0, 1, 2) except in dimension-specific sections
- [ ] No assumptions like "grid_shape has exactly 3 elements"
- [ ] Dynamic slicing used for multi-dimensional array operations

---

## 8. Summary

**One substrate. Mechanical evolution only.**
Everything else is either:
- initialization convenience, or
- diagnostic measurement / coarse-grained interpretation.

If an effect cannot be obtained without adding an emergent back-reaction force,
it does not count as "emergent" in this project.

---
