# How the U(1) (EM) and SU(3) (color) sectors interact — the binding channels

**Motivating question.** The carrier-complexified lateral displacement is one
`U(3)` triplet envelope `Ψ=(ψ_x,ψ_y,ψ_z)∈ℂ³`. EM/charge is its overall phase (the
trace / `(1,1,1)/√3` direction); color is its internal orientation (the traceless
`SU(3)` directions). If the two sectors decoupled, a proton would split into a
localized color core and a free positive EM charge that wanders off. They must
therefore **interact**. This note derives the interaction in three channels and
states, for each, what is established and what remains open.

Date: 2026-06-07. Owners: physics-derivation (channels) + the active build
(diagnostics). Primary model: the lattice central-force spring (NOT continuum
StVK). Builds on `geometric_nonlinearity_alpha_scaling.md` (the exact `∝α`
anharmonicity) and `paper/05b_effective_field_theory.tex` §5b (the
`U(3)=U(1)⊕su(3)` decomposition, EM = trace direction `ê_s=(1,1,1)/√3`).

**Reframing (load-bearing).** There are not two fields near each other; there is
one ℂ³ lump. The question is not "what force holds two objects together" but
"**what restoring coupling stops the single lump from deforming so its trace
content (charge) spreads to long range while its traceless content (color) stays
localized.**" That fixes what to measure: the co-location of the trace- vs
traceless-sector centroids and the restoring force against pulling them apart.

---

## Channel A — kinematic (shared core), quadratic order: DECOUPLED

At quadratic (linear) order the gradient energy is block-diagonal,

    Ψ†(−∇²)Ψ = |∇Ψ_tr|² + |∇Ψ_⊥|² ,

so trace and traceless **decouple** — the established "α-undetermined / sectors
decoupled at linear order" result (`../gauge/status.md` D2; `alpha_holonomy_estimator.md`).
Linear theory therefore predicts the charge *does* drift off, beyond the soft cost
of tearing one amplitude lump into two. **Any binding must be nonlinear**, and (by
the exact `∝α` structure of the anharmonicity) must vanish as `α→0`. This is the
consistency check every candidate below must pass.

---

## Channel B — energetic cross-vertex (∝ α): DERIVED

The only anharmonicity is the norm term `−k_s α a |ΔR|`. For a transverse link,
split the real lateral displacement `ξ = ξ_s ê_s + ξ_⊥` (trace + traceless,
orthogonal), so `|Δu_⊥|² = ξ_s² + |ξ_⊥|²` and
`|ΔR| = √(a² + ξ_s² + |ξ_⊥|²)` (longitudinal stretch set to its mean, see caveat
below). Expanding `√(1+x)`:

    −k_s α a |ΔR| = −k_s α a²[ 1 + (ξ_s²+|ξ_⊥|²)/2a² − (ξ_s²+|ξ_⊥|²)²/8a⁴ + … ].

- The **quadratic** piece `−(k_s α/2)(ξ_s²+|ξ_⊥|²)` is block-diagonal → confirms
  Channel A decoupling.
- The **quartic** piece `+(k_s α/8a²)(ξ_s²+|ξ_⊥|²)²` contains the genuine
  U(1)↔SU(3) vertex via the binomial cross product:

      ┌─────────────────────────────────────────────┐
      │   V_{U(1)×SU(3)} = (k_s α / 4a²) · ξ_s² |ξ_⊥|²   │   (microscopic, per link)
      └─────────────────────────────────────────────┘

Cycle-averaging the carrier (`ξ=Re[Ψ e^{−iω₀t}]`) promotes this to an envelope
vertex

    V̄ = g |Ψ_tr|² |Ψ_⊥|²,   g = Θ · k_s α / a²,   Θ ∈ [1/16, 3/32],

where `Θ` is an O(1) carrier-phase-correlation factor (undetermined here: `1/16`
if the two sectors carry independent carrier phases, `3/32` if they share one).
Continuum density: `(k_s α/4a) |∂u_tr|² |∂u_⊥|² ∝ k_s α/a`, matching the overall
geometric-quartic scaling.

**Established:** the vertex exists, is degree-4 in the field (quadratic in each
sector), and is **exactly `∝α`** — it vanishes at `α=0`, recovering exact
decoupling (consistency check PASS). The splitting parameter `α` is therefore also
the coupling parameter: one dial both splits the scales and couples the sectors.

**Falsifiable prediction (force law).** For Gaussian sector envelopes of common
width `w`, rigidly displacing the trace centroid from the traceless centroid by `d`
gives a cross energy `E_x(d) = g N₁N₂(πw²)^{3/2} e^{−d²/4w²}` and a force

    F(d) = (Θ k_s/a²) · α · C_prof · (d/2w²) · e^{−d²/4w²},   C_prof = N₁N₂(πw²)^{3/2}.

Three signatures: (1) **`F(d) ∝ α`** — predicted ratios `F_α/F_{0.5}=α/0.5`
(`0.4,1.0,1.4,1.6` at `α=0.2,0.5,0.7,0.8`), extrapolating to 0 as `α→0`;
(2) linear-in-`d` for `d≪w` (stiffness `κ=Θ α k_s C_prof/(2a²w²)`), peaking near
`d≈√2 w`, then exponentially decaying; (3) `F ∝ A_tr² A_⊥²` (quartic in amplitude).
**Failure threshold:** if `F(d)` does not scale linearly with `α`
(`|F_{0.2}/F_{0.5} − 0.4| > 0.1` after dividing out profile/amplitude factors), the
geometric-only origin of the coupling is falsified.

**Binding sign — RESOLVED (2026-06-07 derivation): geometric-only channel is REPULSIVE.**
- **The cubic vertex, exactly.** Keeping `Δu_∥≠0`, the √-expansion gives a *cubic*
  cross term `V_3 = +(k_s α/2a²) Δu_∥ (ξ_s² + |ξ_⊥|²)` — **separable**: both sectors
  couple to the longitudinal stretch `Δu_∥` with the *same* coefficient `+k_sα/2a²`,
  `∝α`, →0 at α→0. (This *corrects* the earlier flagged form `(k_sα/2a²)Δu_∥ξ_s|ξ_⊥|²`,
  which was dimensionally wrong — degree-4 and not produced by the norm.)
- **Core stretch sign.** Relaxing `Δu_∥` against its harmonic cost gives
  `Δu_∥* = −(α/2a²)|Δu_⊥|² < 0`: a transversely-displaced link **shortens** axially
  ("guitar-string pull-in"). In the donut, `Δu_∥*(ρ) ∝ −A(ρ)²`, most compressed on the
  ridge `ρ≈w`, zero on the axis hole.
- **Net sign.** Substituting the slaved `Δu_∥*` yields an *attractive* mediated cross
  term `−(k_sα²/4a⁴)ξ_s²|ξ_⊥|²`, but it is `O(αε_⊥²)≈0.08` (with `ε_⊥=|Δu_⊥|/a`,
  `A₀~0.3`) — an order of magnitude **smaller** than the repulsive bare quartic. So
  `g_net = (k_s α/4a²)[1 − αε_⊥²] > 0`: **net overlap-repulsive, the cubic does NOT
  flip the sign.** The earlier worry that the cubic dominates is retracted.
- **The one route to binding: an externally-sourced `Δu_∥`.** The cubic survives at
  `O(α)` (and *can* dominate the quartic, flipping to attractive) **only if `Δu_∥` is
  not slaved to `|Δu_⊥|²` but sourced independently with a definite negative sign — by
  the timelike-link/worldtube prestress `r_t`** (the worldtube channel, A4a/#22). This
  is the open knob that decides the binding sign; it is set by the time link, not by
  the spatial geometry. See Synthesis below.
- **Falsifier (sharpened).** Measure the sign **and α-power** of `F(d)`: `F∝α¹` and
  **positive (repulsive)** ⇒ quartic dominance (this derivation); `F∝α²` and
  **negative (restoring)** ⇒ cubic-dominated binding (would mean `Δu_∥` is externally
  sourced, falsifying the slaving assumption). Also: directly measured `Δu_∥(ρ)` must
  be **negative, peaking at `ρ≈w`, `∝αA₀²`** — a positive core stretch falsifies the analysis.

---

## Channel C — topological (Goldstone–Wilczek + ℤ₃ triality)

### C.1 ℤ₃ triality lock — ESTABLISHED (group-theoretic, dynamics-free)

The model declares the field to be `U(3) = (U(1)×SU(3))/ℤ₃` (the center `ℤ₃⊂SU(3)`
identified with cube-roots of unity in `U(1)`) — *not* an independent
`U(1)×SU(3)`. Single-valuedness then forces the trace charge `q` and the `SU(3)`
triality `t∈{0,1,2}` (`𝟙→0, 𝟑→1, 𝟑̄→2, 𝟖→0`) to satisfy

      ┌──────────────────────┐
      │   q ≡ −t   (mod 3)    │
      └──────────────────────┘

So a color **triplet** (single "quark", `t=1`) is *forced* to carry fractional
(`⅓`-integer) trace charge; a color **singlet** (`t=0`) carries integer charge.
This holds exactly, with no small parameter and independent of all dynamics, *given*
the field is genuinely `U(3)`. It binds charge **fractionality mod 1** to color.

**Falsifiable (P-triality):** build a `B=3` color-singlet, a single color-triplet
sub-lump (`t=1`), and an octet (`t=0`); measure `Q_{U(1)} mod 1` via the trace
holonomy. Predicted: singlet & octet → `0`; triplet → `±⅓`. Any deviation of the
singlet/octet from integer by `>0.1` **falsifies that the substrate field is truly
`U(3)`** (it would be `U(1)×SU(3)` with no shared center).

### C.2 Goldstone–Wilczek current — PLAUSIBLE, NOT DERIVED

A term `L_GW = c₁ a_μ B^μ` in `S_eff` (with `B^μ = (1/24π²) ε^{μνρσ} tr(L_νL_ρL_σ)`
the topologically conserved baryon current, `L_ν=U⁻¹∂_νU`) yields, on varying `a_μ`,

      ┌────────────────────────────────────────────┐
      │   J^μ_{U(1)} = c₁ B^μ + J^μ_matter           │
      │   ⇒  Q_{U(1)} = c₁ B + Q_matter              │
      └────────────────────────────────────────────┘

i.e. the Skyrme/GW relation `Q = c₁B + (non-topological)` (`c₁=½` in QCD). If
present, charge cannot leave without the winding `B`: binding becomes **topological**,
and forced separation draws a confining EM–color flux tube rather than two free
pieces. The integer/topological part of the charge complements the `mod 1`
fractionality lock of C.1 (the two are independent, not redundant).

**Substrate-native candidate mechanism.** The same geometric quartic of Channel B
couples the trace/dilatational mode to the traceless/shear modes; where the `SU(3)`
texture density `|u_⊥^{tl}|²` is large (nonzero winding), the trace mode is
energetically cheapest to *also* wind there — the spatial signature of GW. Treating
the cross term as the source gives a candidate scaling `c₁ ~ O(α)` (ansatz, not
derived) — making the topological binding vanish at `α→0` like everything else.

**The obstruction — RESOLVED (2026-06-07 derivation): VERDICT CONDITIONAL-NO.**
The cleanest GW origin (Wess–Zumino / anomaly inflow) needs a parity-odd carrier
response; BACKBONE #16 says the k-space curvature is identically zero ∀α (`D(k)`
real-symmetric) — the **BZ-anomaly route is blocked.** The `(x,t)`-fibre escape
([[project_complex_u1_from_time]]) is **genuine but insufficient:**
- The "i from time" (carrier time-quadrature `J: ξ↦ξ̇/ω₀`, `J²=−1`) supplies an
  orientation, and the required `ε^{μνρσ}` **factorizes** as `(time leg, from J) ×
  (spatial winding, from the SU(3) texture)`. So GW is *not* forbidden by reality —
  and this is exactly why a nonzero `(x,t)` U(1) curvature exists (good for EM/D3),
  escaping #16. The BZ-zero and the `(x,t)`-nonzero are **different base spaces**, no
  contradiction.
- **But the norm-only vertex sources the wrong partner.** `|ΔR|=√(a²+|Δu|²)` depends
  on the *modulus* of the real displacement, so the carrier feels
  `δω₀ ∝ (k_sα/a²)|Ψ_⊥|²` — the **P-even color energy density**, NOT the **P-odd
  winding density** `Im tr(LLL)` a baryon current `B^μ` requires. Hence `Q` tracks
  color *energy* (continuously, `∝α`), not winding `B`: **no GW term, binding is
  energetic-only, charge NOT topologically locked.**
- **The single condition that flips it to YES:** if the *quartic* fibre metric
  `𝒢_ij=⟨∂_iΨ|∂_jΨ⟩` develops a nonzero **antisymmetric (Kähler / J-twisted) part**
  sourced by the time-quadrature, a P-odd source `∝Im tr(LLL)` appears and a linking
  term `c₁ a_μ B^μ` (`c₁` = 4D linking `Lk(Σ₂,ℓ₁)`, topological) is generated. This is
  exactly the **open D2 fibre-metric gap** — the one unchecked algebraic step.
- **Falsifier (P-GW-2, amplitude-rescale at fixed winding):** rescale `Ψ_⊥→λΨ_⊥`
  holding `B`. No-GW (this verdict): trace holonomy `γ_Γ ∝ λ²` (tracks energy, slope
  `∝α`, ratio 2.5 across `α=0.2→0.5`). GW: `γ_Γ` flat/quantized (ratio 1.0). A
  continuous `λ²` scaling **confirms** energetic-only; flatness falsifies it.

**Falsifiable (P-GW, charge-stripping):** relax a `B≠0` texture; measure `B` (lattice
winding) and `Q_{U(1)}` (trace holonomy around the worldtube). Add a pure-gauge trace
phase gradient and re-relax. Topological binding → `Q_{U(1)}−Q_matter` snaps back to
the same quantized `c₁B`; energetic-only → `Q_{U(1)}` drifts continuously. With the
`c₁~α` ansatz, `c₁(0.5)/c₁(0.2)=2.5`. Failure threshold: post-relaxation residue
differing from its first value by `>10%` of one quantum falsifies topological GW
binding; a `c₁` ratio off `2.5` by `>20%` falsifies only the `∝α` ansatz.

---

## Summary table

| Channel | Status | Mechanism | α-scaling | Falsifier |
|---|---|---|---|---|
| A kinematic | decoupled (linear) | shared envelope only | — (vanishes) | n/a — baseline |
| B energetic | **derived — net REPULSIVE** | quartic `+` dominates the slaved cubic; `g_net=(k_sα/4a²)[1−αε²]>0` | `∝α` | sign+power of `F(d)`: `+α¹` repulsive |
| C.1 triality | **established** | `U(3)=(U(1)×SU(3))/ℤ₃`, `q≡−t mod 3` | dynamics-free | `Q mod 1` of singlet/triplet/octet |
| C.2 GW current | **soft winding-lock** (resolved via time link) | `f∧tr(L^3)`: charge tracks winding `Q=κ(A)·B`, but `κ(A)∝A_tr A_⊥` amplitude-graded, not quantized; soliton-only | `∝α` | `γ_Γ∝B¹` at fixed amplitude (vs flat=energy) |

**Bottom line (2026-06-07).** The sectors are *not* independent, but the two
*spatial-substrate* binding channels both come back **negative**: the energetic
quartic is net **repulsive** (the slaved cubic does not flip it), and the topological
GW current is **absent** from the modulus-only norm vertex (CONDITIONAL-NO). The only
firmly *positive* binding is the **ℤ₃ triality lock** — but that binds charge
*fractionality* to color representation, not the spatial *co-location* of the charge
and color lumps. So the question "what stops the proton fragmenting spatially" is
**not yet answered by the spatial geometry alone.**

**Synthesis — both open knobs route through the timelike link / worldtube (`r_t`).**
Remarkably, the single condition that could flip *each* negative to positive is the
same channel:
- Channel B binds only if `Δu_∥` is **externally sourced with a definite negative
  sign** — by the timelike-link prestress `r_t` (the worldtube), not by the slaved
  spatial response.
- Channel C.2 generates a GW term only if the **antisymmetric (Kähler) part of the
  quartic fibre metric `𝒢`** is nonzero — and that antisymmetry is sourced by the
  **time-quadrature `J`** (the "i from time"), i.e. again the timelike/worldtube
  structure.

Both point to the **time link as the binding agent**, tying directly to A4a (single
4D prestress `α_t=α`) and the unified-contraction picture (#22).

**Time-link derivation — DONE (2026-06-07): VERDICT CONDITIONAL-YES**
(`paper/derivations/time_link_binding.md`). The single `r_t` turns **both** knobs:
- **Part 1 (energetic, DERIVED).** With `ω` fixed externally by worldtube closure
  (`ωT=2πn`), the time-link norm term sources a definite-**negative** compression
  `Δu_∥^ext=−(mα/4k_sb²)(ωA(ρ))²<0` co-located with the lump — *externally sourced*,
  not slaved. It flips the cubic vertex: `g_net=(k_sα/4a²)[1−χ]`,
  `χ=2α(ω/ω_*)²(A_0/a)²`. **Binding when `χ>1`, reachable at α≈0.5–0.8.** This upgrades
  the spatial-only *repulsive* verdict above (which was the `r_t→0` limit). Linchpin
  axiom: `ω` set by closure, not relaxation — verify first in the solver.
- **Part 2 (topological, RESOLVED → "soft winding-lock").** The time-quadrature `J` + the
  α-split `ω_tr≠ω_⊥` give a nonzero Kähler `Im𝒢_sa∝(ω_tr−ω_⊥)∝α`. The contraction
  question is now closed: the P-odd `f=da` **cannot** wedge into the P-even scalar
  `|Ψ_⊥|²`; the only available P-odd 4-form is `f∧tr(L^3)` (WZW descent), so **charge
  contracts onto the winding** `Q=κ(A)·B`. **But** `κ(A)∝A_tr A_⊥` is amplitude-weighted,
  not a quantized integer — a **soft topological lock** (winding-tracked, amplitude-graded),
  nonzero only in the soliton sector (zero in the linear spectrum via `g([111])=0`), gated
  on a stable `B≠0` texture. So C.2: CONDITIONAL-NO → **soft winding-lock**.
- **Unification (sharp, falsifiable).** The same `r_t` sources gravitational time
  dilation (#22), so **binding strength and gravity strength co-vary in the single
  dial α** (`∂ln κ_bind/∂ln α = ∂ln c_grav/∂ln α + 1`) — not independently tunable; a
  new test of A4a.

**Net picture.** The proton binds spatially *because of the worldtube/time prestress*,
a single mechanism unifying EM–color binding with gravity — conditional on the closure
axiom (Part 1) and the winding-contraction (Part 2). Clean negative fork: if `Δu_∥` is
ω-independent and `Im𝒢` small, sectors separate spatially and only the ℤ₃ triality lock
survives. Tracked in `status.md` **D6** (color bridge); experiment in `EXPERIMENT.md` (binding probe;
decisive test = `Δu_∥(ρ)∝ω²` via the closure index).