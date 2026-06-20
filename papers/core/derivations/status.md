# Substrate Bridge — Status

## HAVE

- Discrete 4D brane action with codimension-zero embedding; Lorentzian sign as a property of the action (A5), not the ambient.
- Geometric nonlinearity derivation: the entire anharmonic sector is the norm term −k_s αa|ΔR|, exactly ∝ α; quadratic geometric stiffness ∝ (1−α); quartic coefficient +k_s α/(8a²) > 0 (hardening).
- Lattice-to-continuum reduction: acoustic tensor A_abcd = (k_s/a)[α Q_abcd + (1−α) δ_ac δ_bd]; wave speeds c_L² = k_s a²/m, c_T² = (1−α) k_s a²/m.
- Spring constitutive continuum limit: central-force spring ½k_s(|ΔR|−αa)² established as the constitutive law; StVK demoted to quadratic-order proxy (inverted quartic α-scaling).
- Cross-check of discrete and continuum wave speeds via conventions.py closed-form dispersion.
- Two-time (block) BVP well-posedness **resolved in the linear regime**: a two-past-slice chiral boundary condition gives an N-independent, bounded condition number (OPEN_PROBLEMS A2, `in-progress`, linear regime resolved 2026-05-30; stated in §substrate-model of the paper).
- See: discrete_4d_brane_action.md, geometric_nonlinearity_alpha_scaling.md, lattice_to_continuum.md, spring_constitutive_continuum_limit.md.

## MISSING

- Saddle-point solver for the full Lorentzian action (OPEN_PROBLEMS A1, `open`); the linear-regime well-posedness above does not yet extend to the nonlinear two-time problem with the time-link quartic restored.
- Proof that the 6-neighbor stencil with retuned shell weights produces an isotropic acoustic tensor; or explicit acceptance of the leading-order cubic anisotropy as a feature with quantified magnitude. (The paper takes the latter stance: the cubic anisotropy is the structural source of the gauge sector, not a defect to tune away.)
- Verification of α_t = α consistency for the canonical prestressed vacuum (OPEN_PROBLEMS A4a, `adopted; verification open`).
- Marsden–West reference still uncited: §substrate-model attributes the Störmer–Verlet = variational-integrator result to Marsden–West in prose, but references.bib has no matching entry.