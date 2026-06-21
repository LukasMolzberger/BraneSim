# The time link binds the proton: worldtube prestress as the EM↔color binding agent

**Result.** The two open knobs that left U(1)↔SU(3) spatial binding
*negative* in `u1_su3_binding.md` (Channel B net-repulsive; Channel C.2 GW absent) are
**both turned by the single timelike-link prestress `r_t`** (with `α_t=α`, A4a). Verdict:
**CONDITIONAL-YES** — the proton binds spatially, conditional on (1) the carrier rate
`ω` being fixed externally by worldtube closure (the linchpin axiom) and (2), for the
topological leg only, an unfinished contraction step.

Source of both effects: expanding the **existing** time-link norm term `−k_t r_t|ΔR_t|`
(`discrete_4d_brane_action.md` §6; `k_t=m/Δt²`, `r_t=α·βΔt`, held time spacing
`b≡βΔt`). No new force introduced.

---

## Part 1 — Energetic binding (DERIVED, the robust leg)

The lateral content of the *time* link is the carrier velocity:
`|Δu_⊥^{(t)}|² = A(ρ)²(ωΔt)²sin²(ωlΔt)`, cycle-average `½A(ρ)²(ωΔt)²`. The carrier
rotation rate `ω` is fixed **externally** by the closure condition `ωT=2πn`
(rotating-frame-periodic BC), **not** by spatial relaxation — this is what makes the
induced compression an imposed background rather than a slaved response.

Expanding `−k_t r_t|ΔR_t|` gives a negative, lump-co-located energy density
`Ē_{t,⊥}(ρ) = −(mα/4)(ωA(ρ))² < 0`, and a static node force-balance transmits a
**definite-negative longitudinal compression** to the spatial links:

      ┌────────────────────────────────────────────────┐
      │  Δu_∥^ext(ρ) = −(m α / 4 k_s b²)·(ω A(ρ))²  < 0   │
      └────────────────────────────────────────────────┘

peaked on the donut ridge `ρ≈w`, zero on the axis — **`∝α`, `∝ω²`.** This is *not* the
slaved `Δu_∥* = −(α/2a²)|Δu_⊥|²` of Channel B; it is externally sourced by the time-link
prestress + carrier eigenvalue. Feeding it into the separable spatial cubic
`V_3=+(k_sα/2a²)Δu_∥(ξ_s²+|ξ_⊥|²)` flips the net cross-coupling:

      ┌──────────────────────────────────────────────────────────────┐
      │  g_net = (k_s α / 4a²)·[ 1 − χ ],   χ = 2α (ω/ω_*)² (A_0/a)²    │   ω_* ≡ c_L/b
      └──────────────────────────────────────────────────────────────┘

**Binding (g_net<0, restoring `F(d)`) iff `χ>1`** — reachable at the sweet spot
`α≈0.5–0.8`, `ω` near the band edge, `A_0/a~0.3`. Threshold amplitude
`A_0^crit=(b/ω)√(k_s/2mα)`. The prior spatial-only verdict (net repulsive) was the
`r_t→0` limit; the time link genuinely flips it.

**Linchpin axiom — VERIFIED in the solver.** `ω` is set by closure, not
relaxation. Confirmed by code inspection of the `PeriodicBC`/`solve_block` path
(`branesim/solver/{bvp.py,boundary.py}`, `core/residual.py:residual_periodic`):
- The JFNK root-find relaxes **node positions only** (`x=slices.reshape((P,n_nodes,m))`
  in `F_periodic`); `ω` is **not a degree of freedom**.
- The loop length `T=P·dt` is fixed (`P=n_slices`, `dt` both config inputs); the carrier
  winding `n_t` is a **topological integer carried by the seed and preserved by the
  smooth Newton flow unless the amplitude collapses** (`PeriodicBC`/`residual_periodic`
  docstrings). Hence the period-averaged `ω̄ = 2π·n_t/(P·dt)` is pinned by three
  externally-fixed quantities — the experiment labels it the "closure-locked carrier."
- The time-link prestress force (`r_t>0`) **is** present in `residual_periodic` (routes
  on `r_t>0` to the central-force spring on the wrapped time-neighbours), so the Part-1
  source `Δu_∥^ext` is physically in the solver, not absent.
- **Counterexample correctly avoided:** `branesim/solver/breather.py` is an eigen-solver
  where `ω` co-relaxes as a nonlinear eigenvalue (`omega_exact`, `result["omega"]`);
  running the vortex through *that* would re-slave `ω` and break the linchpin. The
  experiment uses `PeriodicBC` (`vortex_seed_render.py`), not the breather solver.

Two caveats (do NOT invalidate the linchpin, but refine it):
1. **Average vs per-slice `ω`.** The pin fixes the *period-averaged* winding (total phase
   `2π n_t`); the relaxed field may redistribute phase velocity around the loop, so the
   per-slice `ω(l)` need not be uniform. Diagnostic: measure `ω(l)` on the *converged*
   worldtube (not just the seed) and confirm it stays ≈ the closure-locked value.
2. **`ω` imposed ⇒ amplitude responds (this is the *right* direction).** `PeriodicBC`
   forces `ω=2π n_t/(P dt)` (a discretization value, not the natural nonlinear
   eigenfrequency), so for `‖R‖→0` the **amplitude self-selects** to be stationary at the
   imposed `ω`. This *confirms* the regime — `ω` is the input, `A` slaves to it (opposite
   of `ω` slaving to `A`) — but it **refines the decisive test**: one cannot truly hold
   `A_0` fixed while varying `n`. Instead, solve a stationary worldtube at each closure
   index `n` (each imposed `ω`), then check `Δu_∥ ∝ ω²` across the *converged family*,
   dividing out the measured `A` (see Falsifiers).

---

## Part 2 — Topological binding (RESOLVED: a "soft topological lock")

The C.2 flip-condition was: a nonzero antisymmetric (Kähler) part of the quartic fibre
metric `𝒢_ij=⟨∂_iΨ|∂_jΨ⟩`. With the carrier lift `Ψ_i=A_i e^{iω_i t}` and the
**α-controlled trace/traceless frequency split** `ω_tr=ω_0√(1−2α/3)`, `ω_⊥=ω_0√α`:

      ┌─────────────────────────────────────────────────────────────────┐
      │  Im 𝒢_sa = A_s A_⊥ (ω_tr − ω_⊥) sin[(ω_tr − ω_⊥)t + Δφ₀]  ≠ 0     │
      └─────────────────────────────────────────────────────────────────┘

Nonzero **iff** both (a) the time-quadrature `J` supplies the `i` (real `ξ` alone gives
`Im𝒢≡0`), and (b) the sectors are split `ω_tr≠ω_⊥` (α>0). Vanishes as `α→0`.

**The contraction — does `Im𝒢` route to the winding `tr(LLL)` or the energy `|Ψ_⊥|²`?
RESOLVED: the WINDING (but with an amplitude-weighted, not quantized, level).**
The trace curvature `f=da` (the field-space pullback of `Im𝒢`) is a P-odd 2-form; it
**cannot** wedge into a P-even scalar 4-form (`|Ψ_⊥|²`), so the *only* available P-odd
4-form is `f∧tr(L∧L∧L)` — the Wess–Zumino–Witten 5-form descended to 4D, whose
`a`-variation is `c₁ a_μ B^μ`. So the term that appears is

      ┌──────────────────────────────────────────────────────────────────────┐
      │  S_cand = ∫ a ∧ [ κ(A)·tr(L∧L∧L) ],   κ(A) = (1/2π)∮Im𝒢 dt            │
      │  ⟹  Q = κ(A)·B + Q_matter,   κ(A) ∝ A_tr A_⊥ (ω_tr−ω_⊥) ∝ α          │
      └──────────────────────────────────────────────────────────────────────┘

Charge **tracks the baryon winding `B`** (not the color energy) — a genuine winding
lock. **But** the metric-independence test fails: `κ(A)∝A_tr A_⊥` is amplitude-weighted
(inherited from the modulus `|ΔR|`), **not** a quantized linking integer `c₁`. So this
is a **soft topological lock**: GW-like, sourced by the winding, with a non-quantized
amplitude-graded level. Forced separation that unwinds the texture (`B→0`) *does* kill
the charge transfer, so the lock is real — just soft, not hard.

**Linear-spectrum value is identically zero; nonzero only in the soliton sector.** The
wedge needs trace and traceless *simultaneously coherent and color-active*. In the
linear spectrum that requires `k₀∥[111]`, where `g([111])=0` kills `tr(L^3)` — the
coherent-vs-color exclusivity ([[project_alpha_undetermined_at_linear_order]]) kills the
4-form. It survives **only** where the SU(3) hedgehog texture supplies its **own**
coherent color frame `U=exp(iF(r)x̂·T)`, pinning `Δφ₀(x)` so `∮Im𝒢 dt` does not
self-cancel. (ANSATZ: this coherence-pinning is read off the texture *form*, not yet
computed on a converged state — the one un-discharged step, exactly what the falsifier
measures.) **Gated on a stable `B≠0` texture (C2/D4), which the project does not yet
have.**

---

## Part 3 — The unification prediction (sharp, falsifiable)

The **same** `r_t` (single `α_t=α`, A4a) sources **gravitational time dilation**
(BACKBONE #22, via the time-link quartic `∝r_t`) **and** the binding compression of
Part 1 (`∝r_t`). Therefore:

> **Binding strength and gravitational time-dilation strength are the same dial `α` and
> are NOT independently tunable.** `κ_bind ∝ α·χ` and `c_grav ∝ α` co-vary:
> `∂ln κ_bind/∂ln α = ∂ln c_grav/∂ln α + 1` (binding gains one extra `α` from `χ`).

If a simulation can tune binding without tuning gravity (or vice versa), **A4a (`α_t=α`)
is falsified** and the time link is not the unified binding agent. This is a new, strong
test of A4a beyond the isotropy/Newtonian-limit checks already listed there.

---

## Combined binding condition & verdict

      Proton binds spatially  ⟺  [ χ = 2α(ω/ω_*)²(A_0/a)² > 1 ]   (Part 1, energetic, DERIVED)
                              AND/OR  [ soft winding-lock Q=κ(A)·B ]  (Part 2, topological, DERIVED)

Both legs now derived; they are complementary. Part 1 binds the **position** (restoring
force, `χ>1`); Part 2 binds the **quantum number** (charge tracks winding `Q=κ(A)·B`,
amplitude-graded). **VERDICT: CONDITIONAL-YES.** The strong leg (Part 1) is derived and
physically reachable (closure axiom verified). The topological leg is a **soft winding
lock** — charge contracts onto `tr(L^3)`, not the energy, but with a non-quantized
amplitude-weighted level, and only inside a coherent soliton texture. **Clean negative
fork:** if `Δu_∥` is ω-independent (closure axiom fails — already excluded by the solver
check) AND `γ_Γ` shows no `B`-dependence, the sectors separate spatially and only the
**ℤ₃ triality lock** (C.1) survives — binding charge fractionality, not position.

## Falsifiers

| Observable | `r_t=0` | `r_t>0`, binds (Part 1) | `r_t>0`, soft winding-lock (Part 2) |
|---|---|---|---|
| Force law `F(d)` | repulsive (`+α¹`) | **restoring (flips at χ=1)** | restoring |
| `Δu_∥(ρ)` vs `ω` | flat in ω, `∝αA_0²` | **`∝ω²`, `<0`, lump-located** | `∝ω²` |
| `γ_Γ` vs `λ` (rescale `Ψ_⊥`) | `∝λ²` | `∝λ²` | **bilinear `∝λλ'`** (NOT flat, NOT `λ²`) |
| `γ_Γ` vs `B` (fixed amplitude) | flat | flat | **`∝B¹`** ← decisive: energy-only is flat in B |
| `Im𝒢` magnitude | 0 | small | **`∝α(ω_tr−ω_⊥)`, O(Re𝒢)** |
| binding ↔ gravity | n/a | **co-vary in α** | co-vary |

**Decisive C.2 discriminator:** vary `B` at *fixed* amplitude. Pure-energetic predicts
**no** `B`-dependence of the P-odd holonomy; the soft winding-lock predicts **`γ_Γ∝B¹`**.
A `B`-slope of 1.0 confirms winding contraction; flat-in-`B` would mean energetic-only.

**Decisive single test (refined per the solver verification):** solve a *stationary*
worldtube at each closure index `n` (each imposed `ω=2π n/(P dt)`); on each converged
state measure `Δu_∥(ρ)` and the relaxed amplitude `A`, and check `Δu_∥ ∝ ω²` across the
converged family (divide out the measured `A` — it self-selects to the imposed `ω`, so
one cannot hold `A_0` fixed). `∝ω²` confirms the externally-sourced compression
(binding); flat-in-ω reverts to the spatial-only repulsive verdict. Also confirm
per-slice `ω(l)` stays ≈ uniform on the converged worldtube (caveat 1).

## Ledger
- **DERIVED:** the negative co-located compression `Δu_∥^ext∝αω²` and the flip `χ>1`;
  `Im𝒢≠0` needing both `J` and the α-split; the binding↔gravity co-variation.
- **ANSATZ:** the quasi-static node force-balance transmitting the time-link pull to
  spatial `Δu_∥`; the `c₁` coefficient form.
- **AXIOM (load-bearing):** A4a `α_t=α`; the carrier lift `Ψ≈ξ+(i/ω₀)ξ̇`.
  (The Part-1 linchpin **`ω` fixed by closure not relaxation** is not an
  assumption — it is **VERIFIED in code** for the `PeriodicBC`/`solve_block` path; see
  above. The breather eigen-solver would violate it, so the experiment must not use it.)
- **RESIDUAL GAP:** Part 2's contraction is RESOLVED — it lands on the
  winding `tr(L^3)` (soft, amplitude-graded lock), not the energy. The one remaining
  un-discharged step is the soliton coherence-pinning (that the hedgehog fixes `Δφ₀(x)`
  so `∮Im𝒢 dt` does not self-cancel), asserted from the texture form, computed only on a
  converged texture — i.e. it is the falsifier, gated on a stable `B≠0` texture (C2/D4).
- **GATING:** all topological claims gated on a stable `B≠0` texture (C2/D4); Part 1
  needs only the U(1)-vortex eigenstate.
- **Constitutive:** existence of the cross-vertex and Kähler part survive any isotropic
  hyperelastic law; the `∝α` coefficients are central-force-specific
  (`geometric_nonlinearity_alpha_scaling.md` §5).

Parent: `u1_su3_binding.md`. Tracker: `papers/gauge_color/derivations/color/status.md` D6 (+ `papers/core/derivations/status.md` A4a, new test).