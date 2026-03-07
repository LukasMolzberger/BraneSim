# baryon_simulation_roadmap.md

## Goal

Adapt the current BraneSim code so that it can search for stable baryon-like localized modes on the cubic substrate,
with proton/neutron configurations treated as the first decisive particle target.

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

### Phase 2 — Add a baryon initializer

Create a new module such as `baryon_initialization.py` with at least two seed families:

1. **Axis-triplet standing-wave seed**
   - three coupled components associated primarily with x, y, z carrier channels,
   - common spherical or slightly cubic envelope,
   - controllable relative phases and amplitudes,
   - option for locked or weakly mixed triplet states.

2. **Confined triplet tube / shell seed**
   - a more structured seed with internal circulation or standing-wave nodes,
   - explicit control of parity-like and charge-trace-like combinations,
   - intended for proton vs neutron comparisons.

Working proton/neutron hypothesis for the first experiments:
- **proton**: confined triplet with nonzero far-field U(1) trace sector,
- **neutron**: confined triplet with near-cancelled far-field U(1) trace sector,
- both share a confined internal triplet/color structure.

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
