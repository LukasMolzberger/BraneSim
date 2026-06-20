# Derivation of G_eff from the spatial contraction channel

## 1. Layer interface

**Layer A (substrate / continuum elasticity):** state variable = embedding X(x,t) ∈ ℝ⁴, decomposed near
identity as X = (h★ xⁱ + ξⁱ, u) where ξⁱ(x,t) is the in-brane displacement triplet and u = X⁴ the
transverse channel. Governing energy: central-force spring W = W₂ + W₄.

**Layer B (inside-observer Newtonian gravity):** state variable = scalar contraction potential Φ_grav(x,t)
felt by a probe wavepacket; probe obeys ẍ = −∇Φ_grav.

**Coarse-graining map:** Φ_grav ∝ θ := ⅓ ∇·ξ (trace of the slow in-brane strain), which is the dilation
field sourced by the transverse-gradient eigenstrain ∂ᵢu ∂ⱼu via quasi-static elastic equilibrium.

---

## 2. Assumptions (numbered)

1. **Slow/static separation.** ξⁱ relaxes fast compared to the timescale on which the source u moves:
   ∂ₜₜξ ≪ c_T² ∇²ξ. (Newtonian limit; to be checked numerically.)
2. **Small slope.** |∇u| ≪ 1 and |∇ξ| ≪ 1: W truncated at W₂ for ξ; eigenstrain source kept to leading
   order ∂ᵢu ∂ⱼu. Leading neglected term O(|∇u|⁴) = soliton self-confinement W₄.
3. **Prestress α ∈ (0,1) required.** At α → 0 the geometric anharmonic sector vanishes; all anharmonicity ∝ α.
4. **6-neighbor axial-only stencil ⇒ Lamé λ = 0, μ = (1−α)k_s/a, Poisson ratio ν = 0.**
   Results that survive any isotropic hyperelastic law: existence of the Poisson-type equation, inward sign
   under prestress, eikonal bending, anharmonic source coupling ∝ α. The exact α(1−α) exponent and the
   zero-crossing at α = 0.5 are 6-neighbor specific.
5. **Isotropic far field (inside-observer, backbone #8 conjecture).** Cubic anisotropy invisible to the
   inside observer. To be checked numerically.
6. **Localized, finite-energy source.** u has compact support, so a 1/r exterior solution exists.
7. **Probe = band-isolated linear wavepacket** in eikonal regime.

---

## 3. Derivation

### (a) Poisson equation for the contraction field

**Step 1 — Quasi-static elastic equilibrium with the geometric eigenstrain.**

The total in-brane strain splits into a displacement part and a transverse-gradient eigenstrain:

    εᵢⱼ = ½(∂ᵢξⱼ + ∂ⱼξᵢ)  +  ½ ∂ᵢu ∂ⱼu
          ε^ξ_ij                ε*_ij

The eigenstrain ε* is the ∂ᵢu ∂ⱼu term already in gᵢⱼ; it enters W₂ as a strain offset.
With λ = 0 (Assumption 4), the stress is:

    σᵢⱼ = 2μ(ε^ξ_ij − ε*_ij)

**Step 2 — Force balance ∂ᵢσᵢⱼ = 0.**

    μ(∇²ξⱼ + ∂ⱼ(∇·ξ)) = μ ∂ᵢ(∂ᵢu ∂ⱼu)

**Step 3 — Take the divergence; define θ := ∇·ξ.**

    2∇²θ = ∂ᵢ∂ⱼ(∂ᵢu ∂ⱼu)

**Step 4 — Identify the source.** Define the local transverse energy density:

    w_⊥(x) := ½ k_⊥ |∇u|²,    k_⊥ = 2μ = 2(1−α)k_s/a

So |∇u|² = w_⊥/μ. Setting Φ_grav := −½ c_T² θ (contraction θ < 0 ⟹ attractive well Φ < 0):

    ∇²Φ_grav = −4πG_eff ρ_eff

with

    ρ_eff := w_⊥ / c_T²

No force was added; the source is the ∂ᵢu ∂ⱼu term already in gᵢⱼ, propagated through W₂.

### (b) Probe bending (eikonal → Newton)

**Step 5 — Probe dispersion on the ξ-background.** The background dilation θ shifts the local metric trace;
the probe's local speed renormalizes:

    c²_loc(x) = c_T²(1 − κ θ(x)),    κ = 2/3   (from the metric trace)

**Step 6 — Eikonal.** Hamiltonian H = c_loc(x)|k|, ray equation:

    ẋ = ∂H/∂k,    k̇ = −∂H/∂x = −|k| ∇c_loc

**Step 7 — Reduce to Newton.** For a non-relativistic probe (|v| ≪ c_T):

    ẍ ≈ −½ c_T² κ ∇θ = −∇Φ_grav,    Φ_grav = ½ c_T² κ θ

**Universality (equivalence principle):** θ enters every probe branch through the common metric trace,
independent of polarization — all matter sectors fall the same way.

### (c) G_eff in substrate parameters

Matching prefactors in Steps 3–4 with μ = (1−α)k_s/a, c_T² = (1−α)k_s a²/m, ρ_m = m/a³:

    c_T² / (μ a²) = a/m

    G_eff = (κ / 16π) · (c_T² a / m) = (κ / 16π) · (1−α) k_s a³ / m²

    G_eff = (κ / 16π) · (1−α) k_s / (a³ ρ_m²)

with κ = 2/3 from the metric trace (O(1) geometric factor; exact value needs full tensor solve).

**Note on observable field strength.** G_eff itself ∝ (1−α) (response modulus). But the *source coupling*
— the eigenstrain ε* that drives contraction — carries one factor of α from the geometric nonlinearity
(all anharmonicity ∝ α). So the **observable contraction amplitude** for a fixed soliton excitation scales as:

    Φ_grav^obs ∝ α(1−α)

### (d) Falsifiable scaling prediction

    ∂ ln Φ_grav^obs / ∂ ln α = (1 − 2α)/(1 − α)

| α   | predicted slope |
|-----|-----------------|
| 0.2 | +0.75           |
| 0.5 |  0.00 (zero crossing) |
| 0.8 | −3.00           |

**Failure threshold:** if the measured Φ^obs(α) is monotonic across [0.2, 0.8] with no extremum near α = 0.5,
the contraction-channel mechanism is falsified.

**Measurement procedure:** build a fixed-shape transverse bump u (width ≫ a, periodic-clamped vacuum),
measure θ = ∇·ξ via the induced-metric trace diagnostic, at α ∈ {0.2, 0.35, 0.5, 0.65, 0.8}. Fit
ln Φ^obs vs ln α.

**Co-variation with binding (A4a).** The source coupling of both G_eff and κ_bind carry exactly one factor
of α (from the geometric eigenstrain). Both the contraction-field source and the binding compression
Δu_∥ ∝ α come from the same geometric nonlinearity. The single-dial claim: binding strength and
gravitational coupling are not independently tunable.

---

## 4. Summary

    ∇²Φ_grav = −4πG_eff ρ_eff

    ρ_eff = w_⊥ / c_T²

    G_eff = (κ/16π) · (1−α) k_s / (a³ ρ_m²)

    ẍ = −∇Φ_grav     (Newton)

    G_eff^obs ∝ α(1−α),   zero slope at α = 0.5

---

## 5. What remains open

- **DERIVED:** Poisson structure from W₂ equilibrium; eikonal → Newton reduction; (1−α) response modulus;
  α source coupling; α(1−α) observable scaling; equivalence-principle argument.
- **ANSATZ:** monopole coarse-graining (traceless/tidal source dropped); κ = 2/3 (exact value needs full
  tensor solve); ρ_eff = w_⊥/c_T² identified as gravitational mass density.
- **LOAD-BEARING AXIOM:** inside-observer isotropy (backbone #8). Without it "Newtonian 1/r" is only
  "elastic near-field along each lattice axis."
- **NOT addressed:** gravitational time dilation (kinematic, follows once (b) holds — slow strain ξ
  modifies probe dispersion relation; Lorentz sign of time-axis converts this into a frequency shift);
  relativistic/T_μν-source generalization; whether G_eff is the same constant probed by binding vs free-fall.

---

## 6. Relevant files

- `paper/03_continuum_substrate_model.tex` — §gravity-channel, §constitutive-law, cone-speeds
- `paper/04_emergent_relativity.tex` — analogue-metric interval (Step 5)
- `paper/derivations/substrate/geometric_nonlinearity_alpha_scaling.md` — α source-coupling scaling
- `paper/derivations/mass/time_link_binding.md` — Part 3, binding↔gravity single-dial co-variation
- `branesim/diagnostics/binding_probe.py` — extend or add `diagnostics/contraction_field.py`