# Color Bridge — Status

## HAVE

- Kinematic confinement argument: coherence and color-activity are mutually exclusive on the cubic lattice (k₀ ∥ [111] gives g([111]) = 0, killing color content; off-[111] the triplet dephases). No free colored asymptotic states.
- U(3) decomposition into U(1) trace direction ê_s = (1,1,1)/√3 (EM) and SU(3) traceless complement (color).
- ℤ₃ triality lock established (group-theoretic, dynamics-free): U(3) = (U(1)×SU(3))/ℤ₃ forces trace charge q ≡ −t (mod 3), binding charge fractionality to color representation.
- Holonomy ratio prediction R(0.5)/R(0.2) = 0.3077; normalization-independent, constitutive-law-independent.
- VSH channel decomposition: hedgehog (J=0, L=1) as canonical baryon-seed angular structure; cubic-anisotropy multiplet splitting Δω/ω ~ α(a/w)² at soliton scale.
- Soft winding-lock: Q = κ(A)·B with κ(A) ∝ A_tr A_⊥(ω_tr − ω_⊥), amplitude-graded, soliton-only.
- See: vsh_channel_decomposition.md, u1_su3_binding.md.

## MISSING

- Stiefel-Whitney topology clarification: whether the real substrate bundle has nontrivial w₁ or w₂.
- Derivation that U(3) triplet transport maps to genuine SU(3) color dynamics (not just group-theoretic availability).
- Color dynamics beyond kinematic confinement: asymptotic freedom, running coupling, hadron spectrum.
- Quark representations: explicit identification of the three axis-aligned directions with the SU(3) fundamental.

## Open derivations

*Relocated from the former central open-problems tracker (group D, gauge sector —
the `U(1)`↔`SU(3)` binding items D4 and D6). Original IDs retained. The
EM/`U(1)`-only items D2, D3, D5 live in the Gauge bridge.*

### D4. EM = `U(1)` vortex / color = `SU(3)` texture; do they bind, and does color stabilize the electron? — `open` (active experiment, 2026-06-06)
**Statement.** The `U(3)=U(1)×SU(3)` field carries two topologically distinct particle
sectors (BACKBONE #25): the EM/electron-like object is a `π₁(U(1))=ℤ` carrier-phase
**vortex** (codim-2 worldsheet, donut cross-section, spin-½ from the `2π` axis tumble,
charge from the temporal carrier); the color/baryon is a `π₃(SU(3))=ℤ` **texture**
(`π₁(SU(3))=0` — no vortex). They are one field, so **protection ≠ participation**:
`SU(3)` is dynamically coexcited around a `U(1)` vortex and may be load-bearing even at
net-zero color (a color-singlet electron `SU(3)`-stabilized).
**Open questions.** (i) Does the substrate **bind** the `U(1)` vortex worldtube — which,
with `SU(3)` unfrozen, is a **semilocal vortex** (simply-connected full vacuum, *not*
topologically protected), so binding is a **dynamical** condition (the `β<1` /
hardening-bound-mode regime; Vachaspati–Achucarro), on a clean periodic-clamped `r_t>0`
eigenstate (the periodic-BC eigenstate route, not seed-and-watch)? (ii) Does the `SU(3)` sector
coexcite and **stabilize** the color-neutral electron (its non-Abelian/semilocal core
moduli)? (iii) **Scale hierarchy:** does the size ratio of the `U(1)`
vortex vs the `SU(3)` texture track `f(α)` — the candidate origin of the QCD-vs-EM
scale split? *Caveat:* the canonical `U(1)/SU(3)` split is dynamical (running
couplings); a static `α` most naturally gives a coupling *ratio*, so a true *scale*
hierarchy is not guaranteed (no geometric-`α` precedent — see the prior-art note).
(iv) the worldsheet's periodic-closure / `(ω,T)` quantization condition. (v) EM =
the **trace `(1,1,1)`** lateral direction — a single-component seed is only `⅓` EM.
**Paper (2026-06-07).** §6.8 now presents this as one of **two complementary numerical
targets** (the `U(1)` carrier-phase vortex alongside the baryon/`π₃` texture), with the
vortex framed as a **semilocal** (non-topologically-protected) bound state stable only in
the type-I `β=(m_scalar/m_vector)²<1` regime (below the Bogomol'nyi point;
Vachaspati–Achúcarro). The framing — not the resolution — moved into the manuscript; the
binding question stays open.
**Status.** This is the active single experiment — full spec in `EXPERIMENT.md`
(object, injection ansatz, one parameter set, measurement suite, scale ladder).
Owner: physics-derivation (ansatz/closure) + the build.

### D6. U(1)↔SU(3) binding: what stops the proton splitting into a color core + a free charge? — `in-progress` (2026-06-07)
**Statement.** EM/charge is the trace (`(1,1,1)/√3`) phase and color is the internal
orientation of **one** `U(3)` triplet envelope `Ψ∈ℂ³`. If the sectors decoupled, a
proton would split into a localized `SU(3)` core and a free positive `U(1)` charge
that wanders off. They must interact; the *nature* of that interaction is the
question. **At quadratic (linear) order they are exactly decoupled** (block-diagonal
gradient energy; D2), so any binding is **nonlinear and must vanish as `α→0`**.
**Why it matters.** A bound charged baryon is required for the whole particle picture
(BACKBONE #25) and is the physical glue between the EM (D3/D4) and color (C2/#24)
sectors. Without it the `U(1)`-vortex / `SU(3)`-texture split (D4) would predict
unbound fragments.
**Three channels (derivation: `paper/derivations/u1_su3_binding.md`, 2026-06-07).**
- **B — energetic cross-vertex (DERIVED; net REPULSIVE, 2026-06-07).** The norm term
  generates a quartic vertex `V=(k_s α/4a²)ξ_s²|ξ_⊥|²` (envelope `g|Ψ_tr|²|Ψ_⊥|²`,
  `g=Θk_sα/a²`), exactly `∝α` — **the splitting parameter is also the coupling
  parameter.** The cubic-vs-quartic sign is now resolved: the true cubic is the
  *separable* `+(k_sα/2a²)Δu_∥(ξ_s²+|ξ_⊥|²)` (the earlier `Δu_∥ξ_s|ξ_⊥|²` was
  dimensionally wrong); slaving `Δu_∥*=−(α/2a²)|Δu_⊥|²<0` gives an attractive but
  `O(αε²)`-**subdominant** correction, so `g_net=(k_sα/4a²)[1−αε²]>0`: **the
  geometric-only channel is overlap-repulsive.** Binding survives *only* if `Δu_∥` is
  **externally sourced (negative) by the timelike-link prestress `r_t`** (see Next).
  Falsifier: sign+α-power of `F(d)` (`+α¹` repulsive = quartic; `+α²` restoring =
  cubic-dominated binding).
- **C.1 — ℤ₃ triality lock (ESTABLISHED, dynamics-free).** `U(3)=(U(1)×SU(3))/ℤ₃`
  forces `q≡−t (mod 3)`: a color triplet carries fractional (`⅓`) trace charge, a
  color singlet integer charge. Binds charge *fractionality mod 1* to color
  independent of dynamics — *given* the field is truly `U(3)`. This is the firmest
  binding result and doubles as a test of the `U(3)` declaration itself.
- **C.2 — Goldstone–Wilczek current (CONDITIONAL-NO, 2026-06-07).** `J^μ_{U(1)}=c₁B^μ`
  would lock *integer* charge to baryon winding (forced separation → confining
  EM–color flux tube). Verdict: the "i from time" **genuinely escapes** BACKBONE #16
  (nonzero `(x,t)` U(1) curvature, the required `ε^{μνρσ}` factorizes time×spatial-
  winding), but the **modulus-only norm vertex sources the P-even color *energy*
  `|Ψ_⊥|²`, not the P-odd *winding* `Im tr(LLL)`** — so **no GW term; charge tracks
  color energy continuously (`∝α`), not `B`; binding is energetic-only.** The single
  flip-to-YES condition: a nonzero **antisymmetric (Kähler) part of the quartic fibre
  metric `𝒢_ij`** (= the open D2 gap), sourced by the time-quadrature. Falsifier
  (P-GW-2): amplitude-rescale `Ψ_⊥→λΨ_⊥` at fixed `B` — no-GW gives `γ_Γ∝λ²`, GW gives
  flat/quantized.
**Falsifiers (the binding probe, `EXPERIMENT.md`).** (i) sector-centroid co-location
baseline; (ii) forced-separation `F(d)` — **sign + α-power** (`+α¹` repulsive =
quartic dominance, the derived expectation; `+α²` restoring = cubic-dominated binding
from an externally-sourced `Δu_∥`); needs an authorized 4-point `α`-ladder
`{0.2,0.5,0.7,0.8}`; (iii) amplitude-rescale holonomy `γ_Γ∝λ²` ⇒ energetic vs
topological (C.2); (iv) `Q mod 1` of singlet/triplet/octet seeds (C.1). Reuses
`confinement.py` (sector centroids via `P_tr=ê_sê_s†`) and `berry_holonomy.py`
(`Q_{U(1)}`, winding `B`), plus a `Δu_∥(ρ)` longitudinal-stretch readout.
- **TIME-LINK BINDING (DERIVED 2026-06-07, VERDICT CONDITIONAL-YES;
  `paper/derivations/time_link_binding.md`).** The two negative spatial channels are
  **both flipped by the single `r_t`**: (Part 1, energetic, DERIVED) with `ω` fixed
  externally by worldtube closure, the time-link norm term sources a definite-negative
  compression `Δu_∥^ext=−(mα/4k_sb²)(ωA)²<0` co-located with the lump → cubic vertex
  flips → `g_net=(k_sα/4a²)[1−χ]`, **binding iff `χ=2α(ω/ω_*)²(A_0/a)²>1`** (reachable
  at α≈0.5–0.8); (Part 2, topological, RESOLVED) the contraction lands on the
  **winding** `tr(L^3)`, not the energy `|Ψ_⊥|²` (P-odd `f` can only wedge into
  `f∧tr(L^3)`, the WZW descent) → `Q=κ(A)·B` — a **soft winding-lock**: charge tracks
  baryon number but with `κ(A)∝A_tr A_⊥` *amplitude-graded*, not a quantized integer,
  and only in the soliton sector (zero in the linear spectrum via `g([111])=0`).
  **Unification:** the same `r_t` sources gravity (#22), so **binding ↔ gravity co-vary
  in α, not independently tunable** — a new test of A4a (see core/Gravity bridges).
**Status.** in-progress — **time-link binding is CONDITIONAL-YES** (energetic leg
derived & reachable; topological leg conditional-maybe). The verdict rests on a
**closure-locked carrier**: `ω` is not a dynamical degree of freedom but is fixed
externally by worldtube periodicity, `ω̄=2π n_t/(P dt)` (period `P`, step `dt`, integer
closure index `n_t`); this is what makes the time-link contraction `Δu_∥^ext`
definite-negative and co-located with the lump. Two refinements to the premise: (a) only
the *period-averaged* `ω` is pinned, so per-slice `ω(l)` must be read on the converged
worldtube; (b) since `ω` is imposed, the *amplitude self-selects* to it, so the decisive
test is `Δu_∥∝ω²` across the *converged family* of closure indices `n` (dividing out the
measured `A`), not at fixed `A_0`. C.1 triality remains the dynamics-free positive (binds
fractionality, not co-location). **Part 2 contraction is CLOSED:** charge contracts onto
the winding `tr(L^3)` → soft winding-lock `Q=κ(A)·B`, `κ(A)∝A_tr A_⊥` (amplitude-graded,
not quantized), soliton-sector-only. Decisive falsifier: `γ_Γ∝B¹` at fixed amplitude
(energy-only is flat in `B`). **Remaining open:** (i) the soliton coherence-pinning (that
the hedgehog fixes `Δφ₀(x)` so `∮Im𝒢dt` doesn't self-cancel) — asserted from the texture
form, must be computed on a converged state; (ii) a stable `B≠0` texture to compute it on
(C2/D4 — the hard gate; the entire topological leg is vacuous without it). Owners:
physics-derivation (binding/contraction reductions) + soliton-hunter (stable texture).
See `[[project_u3_topology_scale]]`, `[[project_temporal_link_4d_spring]]`,
`paper/derivations/{u1_su3_binding,time_link_binding}.md`.

> *Apparatus note (not theory).* The read-only binding-probe diagnostic (device D8,
> `branesim/diagnostics/binding_probe.py`) emits the scalars these falsifiers need
> (`Δu_∥(ρ)`, `ω(l)`, force law, `γ_Γ` vs `B`/`λ`, `Im𝒢`); the multi-solve driver is
> owner-authorized and not yet built. See `EXPERIMENT.md`.