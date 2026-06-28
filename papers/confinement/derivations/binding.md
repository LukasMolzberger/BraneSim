> **Archived 2026-06-28:** superseded by `derivations/audit_2026-06-28.md`; retained for provenance only.

# Confining the two sectors together: the binding question

**Status:** `conditional`. Spatial-only binding is a **closed negative** (three independent routes);
the only derived binder is the time link, conditional on the worldtube-closure axiom.
**Sources:** `gauge_color/derivations/color/u1_su3_binding.md`,
`matter_mass/derivations/mass/time_link_binding.md`, `[[project_u1_su3_binding]]`,
`[[project_temporal_link_4d_spring]]`. Spatial routes derived with the `physics-derivation` agent.

---

## 0. The reframing

There are not two fields near each other; there is **one `ℂ³` lump**. Binding is not "what force
holds two objects together" but "**what restoring coupling stops the single lump from deforming so
its trace content (charge) spreads to long range while its traceless content (colour) stays
localized.**" The measurable is the co-location of the trace- vs traceless-sector centroids and the
restoring force against pulling them apart.

## 1. Spatial-only binding is a closed NEGATIVE

All three spatial routes return "does not bind," for three structurally different reasons. The deep
common cause: **every anharmonicity is the P-even norm term `−k_sαa|ΔR|`** (`nonlinearity.md` §4),
which depends only on displacement *moduli*. A spatial binder needs a P-odd or phase-coherent cross
term; the modulus structure forbids it at every order.

### Route 1 — VSH coherence-pinning: **does not bind** (derived)
Could the hedgehog's geometric phase-locking make the cross term add coherently (negatively)? No.
Cycle-averaging the vertex `ξ_s²|ξ_⊥|²` gives `¼A_s²A_⊥²`, **phase-independent** (the phase-coherent
piece survives only at the resonance `ω_tr=ω_⊥`, killed by the α-split `ω_tr≠ω_⊥`). The cross energy
is non-negative pointwise:

    E_x = ¼ g ∫ A_s²(x) A_⊥²(x) d³x ≥ 0,   ∀ Δφ₀(x).

There is no sign for `Δφ₀` to flip. The only term where `Δφ₀` *could* enter is a trace–traceless
**bilinear**, but (a) the `L=1` hedgehog gives `∮ x̂ dΩ = 0` (no trace overlap on the sphere), and
(b) any such bilinear is quadratic-order, where the sectors are exactly block-diagonal (Channel A).
Closed at every order the norm term produces.

### Route 2 — semilocal type-I moduli: **does not bind** (mechanism unavailable)
Type-I (`β<1`) semilocal flux attraction needs a **gauged** `U(1)`: `m_v²=e²v²` from `|D_μΨ|²`. But
the substrate `U(1)` is **ungauged** — the connection `a_μ=i⟨u|∂_μu⟩` is diagnostic-only
(No-Back-Reaction). No dynamical gauge field ⇒ no `m_v` ⇒ `β` undefined. The substrate is a *global*
`U(3)` σ-model with a gradient-quartic; the would-be co-localizing term reduces to the
already-**repulsive** Channel-B quartic. (We decline to introduce `m_v` by hand.)

### Route 3 — kinematic shared-core tearing: **does not bind** (derived)
With the cross-vertex set to zero, the quadratic gradient energy is block-diagonal and **separable**:
`E = E_tr[Ψ_tr] + E_⊥[Ψ_⊥]`. Rigidly translating the trace centroid by `d` leaves `E_tr` invariant
(translation symmetry):

    dE/dd |_{g=0} = 0    for d ≳ w.

The neck cost `E₂ ~ k_s(1−α)A²w` is a one-time, finite barrier to *initiate* separation, not a
confining potential — it does not grow with `d`, and the free trace charge is energetically favored
to *delocalize* (spreading lowers gradient energy). This is the original "charge wanders off"
problem.

**Summary (closed):**

| Route | Decisive relation | Verdict |
|---|---|---|
| 1 coherence-pinning | P-even `E_x=¼g∫A_s²A_⊥²≥0`; `∮x̂dΩ=0` | does not bind |
| 2 semilocal type-I | `U(1)` ungauged ⇒ `β` undefined; → repulsive Channel B | does not bind |
| 3 gradient tearing | `dE/dd\|_{g=0}=0` for `d>w`; neck cost finite | does not bind |

Spatial geometry alone does **not** confine charge to colour.

## 2. The only derived binder: the time link / worldtube closure — `conditional`

The single timelike-link prestress `r_t` turns the two knobs the spatial routes could not.

**Part 1 — energetic (derived).** With the carrier frequency `ω` fixed **externally** by worldtube
closure `ωT = 2πn` (not by spatial relaxation), the time-link norm term sources a definite-**negative**
longitudinal compression co-located with the lump:

    Δu_∥^ext(ρ) = −(m α / 4 k_s b²)(ω A(ρ))² < 0,   b = β Δt.

This is exactly the externally-sourced negative `Δu_∥` that Route 3/Channel B could not supply. It
flips the cubic vertex:

    g_net = (k_s α / 4a²)[ 1 − χ ],   χ = 2α (ω/ω_*)² (A₀/a)²,   ω_* = c_L/b.

    ┌──────────────────────────────────────────────┐
    │  Binding iff  χ > 1  (reachable at α ≈ 0.5–0.8) │
    └──────────────────────────────────────────────┘

The linchpin axiom is the **closure-locked carrier**: `ω` is imposed by worldtube periodicity, not a
relaxed degree of freedom. (Refinement: only the period-averaged `ω` is pinned; the decisive test is
`Δu_∥ ∝ ω²` across the converged family of closure indices `n`.)

**Part 2 — topological soft winding-lock (resolved).** The time-quadrature `J` ("i from time") plus
the α-split `ω_tr≠ω_⊥` give a nonzero Kähler part `Im𝒢 ∝ (ω_tr−ω_⊥) ∝ α`. The P-odd `f=da` cannot
wedge into the P-even `|Ψ_⊥|²`; the only available P-odd 4-form is `f∧tr(L³)` (Wess–Zumino–Witten
descent), so charge contracts onto the **winding**:

    Q = κ(A)·B + Q_matter,   κ(A) ∝ A_tr A_⊥ (ω_tr−ω_⊥) ∝ α.

This is a **soft** lock: charge tracks baryon number `B`, but `κ(A)` is amplitude-graded (not a
quantized integer), and nonzero only in the soliton sector (zero in the linear spectrum via
`g([111])=0`).

**Unification with gravity.** The same `r_t` sources gravitational time dilation (`#22`), so
**binding strength and gravity strength co-vary in the single dial `α`** —
`∂ln κ_bind/∂ln α = ∂ln c_grav/∂ln α + 1` — not independently tunable. A new falsifiable test of the
`α_t=α` consistency (A4a): if a simulation tunes binding without tuning gravity, the single-prestress
hypothesis fails.

## 3. The surviving closed statement: ℤ₃ triality

The one firmly *positive*, dynamics-free confinement statement is the `ℤ₃` triality lock
`q ≡ −t (mod 3)` (`su3_texture.md` §6): it binds charge **fractionality** to colour
**representation**, forcing quarks into singlets. It is silent on spatial co-location — that is what
the time link must supply.

## 4. Ledger

| Item | Status |
|---|---|
| Spatial-only binding (3 routes) | `closed` NEGATIVE |
| Time-link energetic binding `χ>1` | `conditional` (closure axiom) |
| Soft winding-lock `Q=κ(A)B`, `κ∝α` | `conditional` (gated on stable `B≠0` texture) |
| Binding↔gravity α-co-variation | `conditional` (A4a test) |
| `ℤ₃` triality `q≡−t (mod 3)` | `closed` (representation, not position) |
| `χ>1` reachable in a converged eigenstate | `open` (numerical) |
| Coherence-pinning of `Δφ₀` on a real texture | `open` (asserted from form) |
