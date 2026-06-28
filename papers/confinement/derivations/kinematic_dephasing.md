> **Archived 2026-06-28:** the dephasing estimate survives only as linear wave physics, not confinement.

# Kinematic colour confinement, quantitative: the dephasing length and the conjugacy bound

**Status:** `closed` for the dephasing length and the divergence at `[111]`; `conditional` for the
figure-of-merit bound (rests on the spectral-susceptibility proportionality, itself an ansatz in
`alpha_holonomy_estimator.md`). Derived with the `physics-derivation` agent.
**Sources:** `gauge_color/sections/03_color.tex` (`subsec:color-confinement`, `eq:holonomy-ratio`),
`gauge_color/derivations/gauge/alpha_holonomy_estimator.md`, `vsh_channel_decomposition.md` §3.

> Upgrades kinematic colour confinement from a qualitative either/or (coherence only on `[111]`,
> colour only off it) to a **closed-form trade-off with a number**: how far colour survives, and a
> direction-/`α`-independent bound on how much colour can be carried over one coherence length.

---

## 1. Branch-frequency split (the spectral driver)

Long-wavelength dispersion (diagonal in Cartesian, `D(k)` real-symmetric):
`ω_a(k₀) = c k₀ √((1−α) + α k̂_a²)`. Expanding about the mean stiffness `s = 1 − 2α/3`:

    δω/ω = √3 α g(k̂) / (2(3−2α)),     g(k̂) = √(Σ_a(k̂_a²−⅓)²) at small k.

- `[100]`: `g=√6/3≈0.8165` ⇒ `δω/ω = α/(√2(3−2α))` (= 0.054 at α=0.2);
- `[111]`: `g=0` exactly ⇒ `δω/ω=0` (degenerate triplet);
- near `[111]` (`k̂=[111]/√3+ε`): `g≈(2/√3)|ε|` ⇒ `δω/ω→0` **linearly** in the angular deviation.

`∝α` and vanishes at `α=0` (where `D∝I`) — passes the single-nonlinearity check.

## 2. The dephasing length

The relative phase among the three components **is** the operationally resolvable colour content.
Components advance at their own branch frequency; the relative-phase spread after distance `L` is
`Δφ_rel(L) = L k₀ (δω/ω)`. Setting `Δφ_rel(L_deph)=1`:

    ┌──────────────────────────────────────────────────────────┐
    │  L_deph(α,k̂,k₀) = 2(3−2α) / ( √3 α g(k̂) k₀ )  ~ w·(c/δc)  │
    └──────────────────────────────────────────────────────────┘

For a packet of width `w` (`Δk~1/w`), inter-branch dephasing beats intra-branch dispersion
(`L_deph ≪ L_disp`) whenever `(k₀w)² · √3αg/(2(3−2α)) ≫ 1` — easily met for `w≫a` off `[111]`. So
the colour dies **as colour** (relative phase slip), not merely by smearing.

## 3. The confinement statement, quantitative

`L_deph` is **finite for every direction off `[111]`** and **diverges as `k̂→[111]`**
(`L_deph → (3−2α)/(α k₀|ε|) → ∞`). Coherence is bought only by sending the colour content to zero.

**Figure of merit (the conjugacy bound).** The resolvable `SU(3)` holonomy rate per unit length is
set by the same branch split, `κ_color ∝ k₀(δω/ω)`. The colour accumulated over one dephasing length:

    ℳ ≡ κ_color · L_deph = ( k₀ δω/ω ) · ( k₀ δω/ω )⁻¹ = O(1),

with `α`, `g(k̂)`, `k₀` **cancelling identically**. Colour rotation and phase coherence are
**conjugate**: their product per propagation length is pinned to unity by the same gap `δω` that
defines both. A free linear excitation can have arbitrarily long coherence *or* a well-defined
colour rate, **never an `O(1)` amount of both over a coherence length** — the quantitative form of
"no free colored asymptotic state."

## 4. Adiabaticity ties the same knot

The WZ connection is well-defined only adiabatically, `Ω_mix ≪ Δω_gap`, with
`Δω_gap ∝ α g(k̂)` — the **same** spectral quantity. Adiabaticity wants a large gap (off `[111]`,
large `α`); coherence wants a small gap (toward `[111]`). The adiabaticity-breakdown loop size
`ℓ_ad = v_g/Δω_gap ~ L_deph`: the largest loop on which WZ transport is valid and the largest length
over which colour survives are the same scale. No clean window — the conjugacy again.

## 5. Falsifiers

- **Direction scan:** at fixed `α,k₀`, `L_deph(k̂)·g(k̂) = const`; e.g. `L_deph[110]/L_deph[100] = 2.0`
  exactly. Deviation `>20%` (at `k₀a≤0.3`) falsifies the `δω∝αg` driver.
- **The bound (decisive):** `ℳ = κ_color·L_deph = O(1)`, independent of `α,k̂,k₀`. If `ℳ` at
  `([100],α=0.2)` and `([110],α=0.5)` differ by more than ~2×, the conjugacy cancellation — the heart
  of the claim — fails.
- Predicted `L_deph[100]·k₀ ≈ 18.4` (≈`23a` at `k₀a=π/4`, order-of-magnitude at that `k`).

## 6. Caveats / ledger

| Item | Status |
|---|---|
| `δω/ω = √3αg/(2(3−2α))`; `→0` at `[111]` | `closed` (from `Λ_a`) |
| `L_deph = 2(3−2α)/(√3αg k₀)`; finite off `[111]`, diverges on it | `closed` (diagonal `D(k)`) |
| Figure of merit `ℳ = κ_color L_deph = O(1)` | `conditional` (κ_color∝k₀δω/ω is the ansatz) |
| `Δω_gap∝αg`; `ℓ_ad~L_deph` | `closed` |
| Absolute numbers at `k₀a=π/4` | order-of-magnitude (long-`λ` expansion); scalings survive |
| `κ_color` proportionality (the `★` susteptibility) | `open` (test via `ℳ` direction/α scan) |

This is a no-go on asymptotic states with a closed coherence length — **not** a string tension or
area law (the QCD non-claim stands).
