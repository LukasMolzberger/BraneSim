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

**Quick check before a PR/change:**
- [ ] Only changing substrate evolution or clearly separated diagnostics/initialization?
- [ ] No emergent back-reaction introduced (even indirectly)?
- [ ] No clamps/cutoffs/hand-imposed thresholds?
- [ ] Coupling remains geometric (distance/metric), not an added "field"?
- [ ] Layer boundaries respected (core vs mapping vs measurements vs experiments)?

---

## 1. Core Model Commitments

### 1.1 Substrate definition (conceptual)
- Material coordinates: `x = (x¹, x², x³) ∈ Ω ⊂ ℝ³`
- Embedding map: `X(x,t) ∈ ℝ⁴`
- Discretization: nodes `p` with positions `R_p(t) ∈ ℝ⁴`

The solver evolves the substrate state `{R_p(t), V_p(t)}` only.

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

### 4.3 Berry phase outputs expected from experiments
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

---

## 8. Summary

**One substrate. Mechanical evolution only.**
Everything else is either:
- initialization convenience, or
- diagnostic measurement / coarse-grained interpretation.

If an effect cannot be obtained without adding an emergent back-reaction force,
it does not count as "emergent" in this project.

---