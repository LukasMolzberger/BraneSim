# Gravity Bridge — Status

> The "derived" items below are **conditional** on the assumptions listed under RESTS ON; they
> are not closed, unconditional results. The paper section (§gravity-open / §gravity-open-reduction)
> correspondingly presents the scalar reduction and G_eff as a priority open problem. See
> G_eff_derivation.md §5 for the full DERIVED / ANSATZ / LOAD-BEARING-AXIOM ledger.

## HAVE (conditional)

- **Poisson structure (derived under Assumptions 1–6 of G_eff_derivation.md).** ∇²Φ_grav = −4πG_eff ρ_eff,
  sourced by the geometric eigenstrain ∂ᵢu ∂ⱼu already in gᵢⱼ propagated through W₂ equilibrium. The
  Poisson *structure* is derived; the identification ρ_eff = w_⊥/c_T² (local transverse energy density)
  *as a gravitational mass density* is an ansatz. See G_eff_derivation.md §3a, §5.
- **Probe bending → Newton (derived, eikonal regime).** Eikonal ray equation in the ξ-background reduces
  to ẍ = −∇Φ_grav. See G_eff_derivation.md §3b.
- **Equivalence principle (derived, same regime).** θ = ∇·ξ enters every probe branch through the common
  metric trace — all matter sectors fall the same way. See G_eff_derivation.md §3b.
- **G_eff explicit (conditional).** G_eff = (κ/16π)(1−α)k_s/(a³ρ_m²); the prefactor uses κ = 2/3 from a
  monopole metric-trace **ansatz** (exact value needs the full tensor solve with the traceless/tidal
  source restored). See G_eff_derivation.md §3c, §5.
- **Falsifiable scaling (sharp).** Observable contraction amplitude ∝ α(1−α); slope
  ∂ ln Φ_grav^obs / ∂ ln α = (1−2α)/(1−α) → zero crossing at α = 0.5, slope +0.75 at α = 0.2.
  Monotonic measurement falsifies the mechanism. See G_eff_derivation.md §3d.
- **Architecture.** Time-link = Lorentz sign + perpendicular axis only. Quartic = spatial self-confinement
  ∝ α. Gravitational time dilation = kinematic: slow ξ modifies probe dispersion relation; Lorentz sign
  converts to frequency shift (backbone #22, not yet derived).
- **Single-dial co-variation (A4a).** G_eff and κ_bind both carry exactly one factor of α from the
  geometric eigenstrain — not independently tunable.

## RESTS ON (assumptions not themselves derived)

- **Load-bearing axiom: inside-observer isotropy (backbone #8).** Without it the "Newtonian 1/r"
  exterior is only "elastic near-field along each lattice axis." Not proven; to be checked numerically.
- **Monopole coarse-graining ansatz:** the traceless/tidal part of the eigenstrain source is dropped;
  κ = 2/3 and ρ_eff = w_⊥/c_T² are ansätze, not exact results.
- **Quasi-static (Newtonian) limit** ∂ₜₜξ ≪ c_T²∇²ξ and small slope |∇u|, |∇ξ| ≪ 1.
- **6-neighbor central-force stencil (λ = 0, ν = 0):** the α(1−α) exponent and the α = 0.5 zero crossing
  are stencil-specific; only the existence of the Poisson structure, the inward sign under prestress, and
  the ∝ α source coupling are claimed law-independent.

## MISSING

- Discharge of the isotropy axiom: numerical confirmation that the cubic far field is isotropic to the
  inside observer (otherwise the 1/r exterior does not follow).
- Full tensor solve replacing the monopole ansatz (exact κ; traceless/tidal source restored).
- Derivation of gravitational time dilation (kinematic consequence of Lorentz structure + ξ-modified
  dispersion relation; mechanism identified, computation deferred).
- Coupling to stress-energy tensor: T_μν source equation (relativistic generalization).
- Numerical verification: measure Φ_grav^obs vs α ∈ {0.2, 0.35, 0.5, 0.65, 0.8} with fixed-shape
  transverse bump u; confirm zero crossing near α = 0.5. Requires `diagnostics/contraction_field.py`.
- Exact O(1) factor κ (needs full tensor solve beyond monopole coarse-graining).
- Constitutive-law independence: current G_eff specific to 6-neighbor central-force (λ=0, ν=0);
  check which parts survive a neo-Hookean or Ogden law.