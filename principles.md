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
   - Fundamental entity: a **3D brane embedded in 4D** (the 4th coordinate is "amplitude").
   - Time `t` is an **external evolution parameter** (not a geometric coordinate in the substrate).
   - Relativity, EM, gravity, particles, quantization are **emergent descriptions**.

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
- Material coordinates: `x = (x¹, x², x³) ∈ Ω ⊂ ℝ³`
- Embedding map: `X(x,t) ∈ ℝ⁴`
- Discretization: nodes `p` with positions `R_p(t) ∈ ℝ⁴`

The solver evolves the substrate state `{R_p(t), V_p(t)}` only.

### 1.1a Structural facts about the lattice (must be respected)

These are mathematical consequences of the central-force pair-spring energy
`U_link = ½ k_δ (|R_{p+δ} − R_p| − α a |δ|)²`. They are not policy — they are
identities the project relies on:

- **Cauchy relations vs cubic isotropy (do not confuse).** Any cubic Bravais
  lattice with purely central pair-wise interactions, expanded around its
  stress-free reference, satisfies the Cauchy relation `C_{1122} = C_{1212}`
  (in Voigt notation: `C_{12} = C_{44}`). This reduces the number of
  independent cubic elastic constants from 3 to 2 and is automatic.
  However, **cubic isotropy** — the condition `C_{1111} − C_{1122} = 2 C_{1212}`
  — is a *separate* requirement and is NOT implied by the Cauchy relation.
  Cubic isotropy depends on the specific shell weights and is NOT automatic
  with the `1/|δ|²` weights used in `components/simulation/grid.py`.
  See `paper-v4/derivations/lattice_to_continuum.md` for the closed-form
  computation. With the current weights and prestress `α`, the residual
  cubic-anisotropy index is `η_cub = −7α / [2(39 − 22α)]`, which evaluates to
  roughly `−21%` at `α = 1` and `−2%` at `α = 0.2` (the current default).
  This anisotropy is leading-order in `ka`, not an `O((ka)⁴)` correction.
- **Implications of leading-order anisotropy.** Until shell weights are
  retuned (or the project commits to a finite static anisotropy at this
  level), claims about "operational rotational invariance" or "emergent
  Lorentz" must be quantitative and acknowledge this baseline. Subtask 2
  in `paper-v4/validation_roadmap.md` is the empirical check; the dispersion
  experiment must measure direction-dependent `c_L` and `c_T` and compare
  against the predicted `η_cub`.
- **Prestress α controls inter-axis coupling.** The cross-term in the link
  energy proportional to `(1−α)·(Δξ)²/(2a)` is the only linear coupling
  between the three lateral channels `ξ¹, ξ², ξ³` from the diagonal-shell
  springs. At `α = 1` the channels decouple at the linear level (gauge group
  reduces to `U(1)³`); at `α < 1` the coupling activates and the gauge group
  is the full `U(3) = U(1) × SU(3)` of the complex narrowband triplet.
  Therefore `α` is the dial that runs the U(1)³ → U(3) crossover.
- **Geometric quartic = Skyrme-class stabilization.** The induced-metric
  contribution `∂_i u ∂_j u` enters the StVK energy as `(μ/4ℓ₀⁴)(|∇u|²)²`,
  which under Derrick scaling `u(x) → u(λx)` scales as `λ^{+1}` in 3D and
  resists collapse. Combined with the `λ^{−1}` quadratic-gradient term, this
  is the standard Skyrme balance and *defeats Derrick's no-go at the
  continuum level*. Lattice spacing `a` provides an independent UV cutoff.
- **Two soliton regimes.** Width `~ a` (lattice-stabilized) → strong
  Peierls–Nabarro pinning, soliton mass tied to the UV scale, no continuum
  limit. Width `≫ a` (Skyrme-stabilized) → exponentially suppressed PN,
  Lorentz-respecting motion, derivable physical scale. **The program targets
  width ≫ a**; experiments must report the width-to-spacing ratio explicitly.

### 1.2 Dynamics (implementation truth)
The simulation is defined by a mechanical energy (or discrete Lagrangian) of the substrate.

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
