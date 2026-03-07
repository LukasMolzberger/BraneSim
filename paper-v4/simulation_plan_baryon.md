# Baryon-scale simulation plan for the brane substrate model

This plan is designed to preserve the core theory bridges rather than remove them:

1. microscopic cubic lattice dynamics,
2. continuum/hyperelastic interpretation,
3. axis-aligned triplet mixing,
4. Berry/Wilczek--Zee diagnostics,
5. localized nonlinear bound states.

The goal is to make those bridges calculable in software.

## 1. Why proton/neutron-like states are the better first target

At the current model stage, proton- and neutron-like states are better first targets than an electron-like state because the lattice hypothesis already contains a natural **three-channel internal sector** (`x`, `y`, `z` axis-aligned carrier branches). That means the first decisive numerical question is not “can we make any localized bump?” but:

> Can the substrate support a stable, femtometer-scale, three-channel localized mode with robust long-range diagnostics?

That is much closer to a baryon-style problem than to a single-channel electron problem.

## 2. What the current code already supports

The current uploaded code already gives a useful base:

- `state.py` supports full **3D intrinsic grids with 4D embedding coordinates**.
- `grid.py` already supports **3D neighborhoods** and **periodic axes**.
- `solver.py` gives a clean **Velocity Verlet** real-time integrator.
- `dimensions.py` and `dimensional_mapping.py` already separate physical and simulation units cleanly.

So the project is *not* starting from zero.

## 3. What the current code cannot yet do

The current code is still too close to a generic linear spring-wave simulator for a proton/neutron search.
The main missing pieces are:

### 3.1 No anisotropic shell physics

`forces.py` uses one global spring constant and one global rest length for all links.
That is too limited for the present theory, because the baryon program needs at least:

- nearest-neighbor axis links,
- face-diagonal links,
- body-diagonal links,
- possibly separate constitutive laws per shell.

Without this, the model cannot cleanly express how the cubic lattice generates axis-triplet structure while still tuning toward approximate isotropy at long wavelength.

### 3.2 No nonlinear confinement mechanism

The present `SpringForceComputer` is linear Hooke-like.
That is not enough for a realistic localized bound-state search.
A proton/neutron-style mode needs a **smooth nonlinear energy law** so that self-guidance and confinement can emerge from the substrate energy itself.

Important: this must remain smooth and variational.
No artificial thresholds, no hard clamps, no piecewise “capture” rule.

### 3.3 No stationary-state solver

`solver.py` supports time evolution, but not:

- damped energy minimization,
- constrained relaxation,
- continuation in parameters,
- stability analysis around a candidate bound state.

Trying to discover a baryon purely by random real-time initial conditions would be extremely inefficient.

### 3.4 No 3D localized triplet initializers

`initial_conditions.py` is still wave-packet oriented and mostly 1D.
For baryon-scale work, you need seeded 3D configurations such as:

- spherically localized triplet seeds,
- toroidal/ring-like seeds,
- parity-odd/even triplet combinations,
- proton-like vs neutron-like partner seeds with different far-field phase structure.

### 3.5 No spectral / Berry diagnostic pipeline

A working theory requires direct calculation of the bridges between:

- lattice eigenmodes,
- triplet subspace,
- mode mixing,
- non-Abelian holonomy,
- localized-mode labels.

The current code has no dedicated modules yet for dispersion calculation, Brillouin-zone transport, or Berry/Wilczek--Zee observables.

## 4. Recommended software architecture changes

## 4.1 Replace the force model, not the whole solver

Keep `VelocityVerletSolver` as the main real-time integrator.
But replace `SpringForceComputer` with a more expressive physics layer.

Suggested new structure:

- `link_potentials.py`
  - abstract smooth potential interface
  - `energy(strain)`
  - `force_prefactor(strain)`

- `forces.py`
  - `AnisotropicNonlinearForceComputer`
  - shell-specific parameters `(k1, k2, k3)`
  - shell-specific rest-length scaling
  - optional shell-specific nonlinear laws

- `grid.py`
  - expose neighbor shell IDs explicitly
  - precompute masks for axis / face-diagonal / body-diagonal links

This keeps the clean current separation:

- grid topology,
- state tensors,
- force computation,
- time stepping.

## 4.2 Add a stationary-state search pipeline

Add new modules instead of forcing this into the Verlet class:

- `relaxation.py`
  - damped relaxation / gradient-flow-like search
  - optional projection to preserve chosen invariants

- `continuation.py`
  - continue a converged state while varying `alpha`, shell couplings, or target energy

- `stability.py`
  - linearize around a converged state
  - compute unstable eigenvalues / Floquet multipliers if periodically breathing

This is essential.
A publishable working theory will stand or fall on whether it can **find and hold** candidate localized states.

## 4.3 Add a spectral branch-analysis layer

Suggested modules:

- `spectrum.py`
  - linearized dynamical matrix on periodic lattices
  - eigenvalues/eigenvectors for selected `k`
  - branch tracking across the Brillouin zone

- `triplet_analysis.py`
  - identify the axis-aligned triplet subspace
  - quantify branch mixing as a function of `alpha`
  - measure effective `g_mix` directly from overlaps/eigenvectors

This is the missing bridge between the paper’s continuum narrative and the actual microscopic lattice model.

## 4.4 Add Berry / Wilczek--Zee diagnostics explicitly

Suggested modules:

- `berry.py`
  - link variables on the discretized Brillouin zone
  - plaquette Berry curvature
  - non-Abelian Wilson loops / holonomies

- `observables.py`
  - charge-like far-field labels
  - holonomy around localized cores
  - channel occupation and mixing fractions

This is where the “color-like” idea becomes either measurable or empty.

## 4.5 Add localized triplet seed generators

Suggested module:

- `localized_initializers.py`

Seed families:

1. **Spherical triplet core**
   - same radial envelope in all three channels
   - different relative phases

2. **Toroidal triplet seed**
   - ring-shaped energy density
   - one or two winding numbers

3. **Parity partner seeds**
   - to search for proton-like / neutron-like neighboring states

4. **Breather seeds**
   - static core plus oscillatory carrier content

The purpose is not to impose the answer, but to search the solution space intelligently.

## 5. Parameter-scaling guidance

For baryon-scale simulations, the dimensional mapper should be re-centered on femtometer physics.

Suggested first pass:

- choose `h_phys` in the range of roughly `0.05 fm` to `0.1 fm`,
- choose the computational box large enough to separate core and far field,
- keep `c_light_sim = 1` through the current mapper design,
- use target rest energies near proton/neutron scale only as **post-search validation targets**, not as hard constraints inserted into the dynamics.

The crucial point is:

- **size** and **energy** should emerge from the nonlinear substrate configuration,
- but the search algorithm may still use them as convergence targets or ranking metrics.

## 6. Suggested milestone sequence

## Milestone A - Linear lattice spectroscopy

Deliverables:

- dispersion surfaces,
- branch polarization content,
- identification of the triplet subspace,
- direct measurement of how `alpha` changes triplet mixing.

This is the first place where the theory can become quantitatively sharp.

## Milestone B - Berry geometry of the triplet sector

Deliverables:

- Abelian Berry phases in the decoupled regime,
- non-Abelian holonomies in the mixed regime,
- plots of curvature vs `alpha` and shell couplings.

This makes the gauge bridge concrete.

## Milestone C - Static or breathing localized triplet states

Deliverables:

- converged localized solutions,
- conserved-energy report,
- core-radius measurement,
- channel occupation report,
- robustness under perturbation.

This is the make-or-break milestone.

## Milestone D - Proton-like / neutron-like partner search

Deliverables:

- at least two nearby stable localized families,
- similar spatial scale,
- close but distinct total energy,
- different long-range holonomy / phase diagnostics.

Even if these are only “proton-like” and “neutron-like” analogues at first, this would be a major result.

## 7. Minimal concrete coding order

If development time is limited, I would implement in exactly this order:

1. shell-aware neighbor classification in `grid.py`
2. anisotropic nonlinear force model in `forces.py`
3. periodic spectral analysis module `spectrum.py`
4. triplet mixing diagnostics `triplet_analysis.py`
5. Berry / Wilson-loop module `berry.py`
6. constrained relaxation solver `relaxation.py`
7. 3D localized seeds `localized_initializers.py`
8. stability analysis `stability.py`
9. proton-like / neutron-like ranking metrics in `observables.py`

## 8. What would count as success

A meaningful first success is not “we simulated the proton exactly.”
It is this:

- the cubic prestretched substrate admits a real three-channel localized mode,
- the mode is stable or long-lived,
- its scale is femtometer-like after dimensional mapping,
- its energy is in the baryon range,
- its long-range diagnostics are not trivial,
- and nearby partner states exist.

That would show that the theory is not merely a philosophical reinterpretation but a calculable dynamical model.
