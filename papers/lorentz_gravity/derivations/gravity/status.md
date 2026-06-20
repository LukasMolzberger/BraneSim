# Gravity Bridge — Status

## HAVE

- **Poisson equation (derived).** ∇²Φ_grav = −4πG_eff ρ_eff, with ρ_eff = w_⊥/c_T² (local transverse
  energy density), sourced by the geometric eigenstrain ∂ᵢu ∂ⱼu already in gᵢⱼ propagated through W₂
  equilibrium. See G_eff_derivation.md §3a.
- **Probe bending → Newton (derived).** Eikonal ray equation in the ξ-background reduces to ẍ = −∇Φ_grav.
  See G_eff_derivation.md §3b.
- **Equivalence principle (derived).** θ = ∇·ξ enters every probe branch through the common metric
  trace — all matter sectors fall the same way. See G_eff_derivation.md §3b.
- **G_eff explicit (derived).** G_eff = (κ/16π)(1−α)k_s/(a³ρ_m²), with κ = 2/3 (O(1) metric-trace
  factor). See G_eff_derivation.md §3c.
- **Falsifiable scaling (sharp).** Observable contraction amplitude ∝ α(1−α); slope
  ∂ ln Φ_grav^obs / ∂ ln α = (1−2α)/(1−α) → zero crossing at α = 0.5, slope +0.75 at α = 0.2.
  Monotonic measurement falsifies the mechanism. See G_eff_derivation.md §3d.
- **Architecture.** Time-link = Lorentz sign + perpendicular axis only. Quartic = spatial self-confinement
  ∝ α. Gravitational time dilation = kinematic: slow ξ modifies probe dispersion relation; Lorentz sign
  converts to frequency shift (backbone #22, not yet derived).
- **Single-dial co-variation (A4a).** G_eff and κ_bind both carry exactly one factor of α from the
  geometric eigenstrain — not independently tunable.

## MISSING

- Derivation of gravitational time dilation (kinematic consequence of Lorentz structure + ξ-modified
  dispersion relation; mechanism identified, computation deferred).
- Coupling to stress-energy tensor: T_μν source equation (relativistic generalization).
- Numerical verification: measure Φ_grav^obs vs α ∈ {0.2, 0.35, 0.5, 0.65, 0.8} with fixed-shape
  transverse bump u; confirm zero crossing near α = 0.5. Requires `diagnostics/contraction_field.py`.
- Exact O(1) factor κ (needs full tensor solve beyond monopole coarse-graining).
- Constitutive-law independence: current G_eff specific to 6-neighbor central-force (λ=0, ν=0);
  check which parts survive a neo-Hookean or Ogden law.