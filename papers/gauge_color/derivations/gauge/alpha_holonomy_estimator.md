# α from the EM-vs-colour coupling ratio via fibre-internal holonomy

**Question.** Can the empirical EM-vs-strong coupling ratio fix the prestress
factor `α := rest_length/spacing`?

**Bottom line.** **No — α is normalization-undetermined at linear order.** The
linear spectrum fixes the *α-scaling* of the trace-`U(1)` / `SU(3)`-WZ holonomy
ratio, but the *overall coupling hierarchy* lives in the connection normalization,
which is set by the nonlinear fibre metric (the geometric quartic), absent from the
linear theory. A normalization-independent falsifiable prediction survives.

Owner: physics-derivation. Tracked: `status.md` D2 (gauge bridge).

---

## 1. Layer interface

L3/L4 continuum complex-envelope (linear quasi-particle) ↔ L4-gauge effective
gauge-field diagnostics (the `U(3)` connection a band-isolated wavepacket sees).

- Lower state: real lateral triplet `ξⁱ(x,t) ∈ ℝ³` on the 6-neighbor cubic slice.
- Upper state: per-wavepacket complex envelope `Ψ(x,t) ∈ ℂ³`, via
  `ξⁱ = Re[Ψⁱ e^{i(k₀·x − ω₀t)}]` (backbone #18).
- Coarse-graining: band-isolation projector `Π_{k₀}` onto a narrow shell about
  `k₀`, then analytic-signal lift to `Ψ`. The gauge object is the **fibre-internal**
  connection on the `ℂ³` bundle over the *physical* base `(x,t)`, split by
  `P_U(1) = (1/3)·𝟙𝟙ᵀ` (trace, EM) and `P_SU(3) = I − P_U(1)` (traceless, colour).

## 2. Key assumptions

1. Linear regime — drop the geometric quartic (StVK Skyrme term, backbone #17).
   *This is exactly the assumption that makes the answer "undetermined."*
2. Scale separation `ε = |∇Ψ|/(k₀|Ψ|) ≪ 1`, width `W ≫ a`, narrowband.
3. Diagonal `D(k)` — exact on the 6-neighbor axial stencil (backbone #16).
4. **k-space connection ≡ 0 ∀α** (real-symmetric `D` ⇒ trivial real eigenbundle;
   proven in BACKBONE #16 / `status.md` D2). This is *why* the
   only nontrivial connection is fibre-internal in `(x,t)`.
5. Carrier coherence of the triplet requires `k₀∥[111]` (exact degeneracy) or
   `α→0`, else the three branches dephase at rate `|ω_a−ω_b|`. (Tension — see §6.)
6. Adiabatic transport; no inter-band transitions.
7. No back-reaction (PRINCIPLES #2): the connection is measured, never re-injected.
8. Constitutive law = central-force pair-spring (Hookean link-norm). The kinematic
   split `P_U(1)/P_SU(3)` survives any isotropic hyperelastic law; the *number*
   `√3 α/(3−2α)` does not.

## 3. Derivation

**Fibre-internal connection (the only nonzero one).** On the `ℂ³` bundle over
physical base `xᵘ=(t,x,y,z)`,

    𝒜_μ = i⟨Ψ|∂_μ|Ψ⟩ ∈ 𝔲(3),
    ℱ_μν = ∂_μ𝒜_ν − ∂_ν𝒜_μ − i[𝒜_μ,𝒜_ν].

This is the object PRINCIPLES §4.2 sanctions — *not* the k-space connection
(which is zero, assumption 4); it lives over `(x,t)`.

**Sector split (gauge-invariant).**
`𝒜_μ = a_μ·(1/√3)𝟙 + Σ_{A=1}^{8} 𝒜^A_μ T_A` (Gell-Mann `T_A`, `tr T_A T_B = ½δ_AB`).
The U(1) part `a_μ = (1/√3) tr 𝒜_μ` is the EM connection (abelian — shifts by an
exact `∂_μ(phase)` under a fibre gauge); the traceless part is the colour
connection (covariant). Gauge-invariant scalars: the abelian flux
`Φ_U(1) = ∮_γ a_μ dxᵘ` and the `SU(3)` Wilson-loop holonomy angle `Φ_SU(3)`.

**Response ratio.** For a fixed small `(x,t)` loop `γ` of area `𝒮`,

    R(α) ≡ |∂Φ_U(1)/∂𝒮| / |∂Φ_SU(3)/∂𝒮|.

**Where the spectrum enters — and where it does NOT.** The eigenvalue split gives
the *spectral susceptibility* of each sector:

    λ_a − λ̄ = (2k_s/ρ) α (h_a − H/3)        (traceless ∝ α)
    λ̄       = (2k_s/ρ) H (1 − 2α/3)          (trace)
    ⇒ ρ_SU3(k̂,α) = g(k̂)·√3 α/(3 − 2α)        (★)   [locked, alpha_separability.py]

But `(★)` is a ratio of eigenvalue *splittings* (dispersion), not of couplings.
The coupling lives in the **connection normalization**: the canonical gauge kinetic
term `−(1/4g_A²)(ℱ^A)²` has `1/g_A² ∝ tr(P_A 𝒢)` with fibre metric
`𝒢_ij = ⟨∂_iΨ|∂_jΨ⟩`. **At linear order `𝒢 ∝ 𝟙`** (envelope norm is flat ℓ² on
`ℂ³`): rotating `Ψ` by `U(1)` or by an `SU(3)` generator costs the same
envelope-norm. The only sector difference is **dimension counting**
`tr P_U(1)=1` vs `tr P_SU(3)=2` — *not* anything α-dependent.

## 4. Result

    ┌─────────────────────────────────────────────────────────────────┐
    │  R(α) = (C / g(k̂)) · (3 − 2α)/(√3 α),                           │
    │  C = N_U(1)/N_SU(3) = tr(P_U(1)𝒢)/tr(P_SU(3)𝒢)                  │
    │      is normalization-undetermined at linear order (𝒢 ∝ 𝟙).     │
    └─────────────────────────────────────────────────────────────────┘

- `g(k̂) = √(Σ_a(h_a − H/3)²)/H`; `g([111]) = 0`.
- The **α-dependence** is carried entirely by the spectral factor (a spectral
  ratio, not a coupling ratio); the **overall scale** `C` is a number the linear
  theory cannot weight by physics.

**Identification + inversion.** The defensible empirical map is `R(α) ↔ g_EM/g_s`
(coupling *amplitudes* — a connection-response scales linearly with `g`, whereas
`α_fine = g²/4π` is the squared object). Inverting:

    α = 3 / (2 + √3 g(k̂) (g_EM/g_s)/C).

With `g_EM/g_s ≈ √((1/137)/0.12) ≈ 0.25`, `g([100]) = √6/3 ≈ 0.816`, and the bare
`C = 1/√2`, this gives `α ≈ 1.2 > 1` — **out of range**, itself a signal the linear
identification is invalid. Hence: **α is normalization-undetermined at linear
order.**

## 5. Regime of validity

`ε ≪ 1`, `W ≫ a`, adiabatic, `k₀∥[111]` (coherence). Leading neglected term: the
geometric quartic, whose continuum coefficient is **`∝ k_s α/a`** (not `μ/ℓ₀⁴`; see
`geometric_nonlinearity_alpha_scaling.md` — the entire anharmonic sector is the
norm term `−k_s αa|ΔR|`, exactly ∝ α). This is *precisely the term that supplies
the missing `C`* by making `𝒢` sector-dependent (the trace/dilatational direction
couples to the `X⁴` gravity channel; the traceless/shear directions do not).
Because the quartic ∝ α, the nonlinear correction enters as `C = C₀ + O(α·
amplitude²)`, so the nonlinear computation *can* determine α. `R`'s α-scaling is
reliable to `O(ε²)`; its absolute value is not defined until the quartic is
restored.

## 6. Falsifiable prediction (decisive, normalization-independent)

    R(α₂)/R(α₁) = [(3 − 2α₂)/α₂] / [(3 − 2α₁)/α₁]   (independent of C, g(k̂)).

For `(α₁,α₂) = (0.2, 0.5)`:  **R(0.5)/R(0.2) = 4.0/13.0 = 0.3077.**

Procedure (extend `branesim/diagnostics/berry_holonomy.py`):
1. Off [111] (e.g. [100], `k₀a = π/4`) build band-isolated `ℂ³` envelopes at
   `α = 0.2` and `α = 0.5`.
2. Transport each around a fixed small loop `γ` in `(x,t)` (NOT in `k`); compute
   `𝒜_μ = i⟨Ψ|∂_μΨ⟩`, split by `P_U(1), P_SU(3)`.
3. Form `R(α)` from curvature per unit loop area.
4. Compare `R(0.5)/R(0.2)` to `0.3077`.

**Failure threshold:** deviation > 10% (at `W/a ≥ 8`, `𝒮 → 0` extrapolated)
falsifies the spectral-susceptibility factorization — i.e. the linear split `(★)`
does not control the fibre-internal curvature, and the linear-layer U(1)/SU(3)
sector identification fails.

## 7. Open / ansatz / missing ingredient

- **Derived:** the α-scaling `(3−2α)/(√3 α g(k̂))`; `k`-space connection ≡ 0;
  gauge-invariance of `Φ_U(1), Φ_SU(3)`; dimension-counting `C` at linear order.
- **Ansatz (testable, not derived):** (a) that the curvature ratio factorizes as
  spectral-susceptibility × normalization (the §6 claim); (b) the identification
  `R ↔ g_EM/g_s`.
- **Sharpest open gap:** `k₀∥[111]` (coherence) **conflicts** with `g([111])=0`
  (colour content vanishes ⇒ `R→∞`); off-[111] the triplet dephases. **No linear
  direction is both coherent and colour-active** — colour curvature is
  intrinsically a soliton-layer (nonlinear) object.
- **Missing ingredient that would determine α:** the sector-dependent fibre metric
  `𝒢`, supplied only by the geometric quartic (backbone #17). A nonlinear
  computation; no phenomenological coupling was introduced to fill `C` (forbidden,
  PRINCIPLES #3).
- **Constitutive sensitivity:** kinematic split and the *existence* of `C` survive
  any isotropic hyperelastic law; the specific α-scaling is StVK/Hookean-specific.

---

**Files:** BACKBONE #16 / `status.md` D2 (gauge bridge; k-space connection
= 0); `branesim/diagnostics/alpha_separability.py` (closed-form `(★)`, `g(k̂)`);
`branesim/diagnostics/berry_holonomy.py` (transport harness to extend for §6);
`BACKBONE.md` #16/#17/#18/#19; `PRINCIPLES.md` §4.2–4.4.