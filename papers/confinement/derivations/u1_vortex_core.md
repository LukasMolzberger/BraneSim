> **Archived 2026-06-28:** the localized zero-vacuum ansatz has no protected `U(1)` winding and is not a static solution of the stated energy.

# The U(1) centerline: a lattice-regulated phase defect (no Higgs vortex)

**Status:** the no-`|Ψ|⁴`-potential / no-VEV / no-Meissner facts are `closed`; the core-energy
result is a **derived scaling within a localized-envelope ansatz** (not a stable-eigenstate
existence claim); stability is `conditional` (semilocal, unit-winding unstable at the Derrick
scale); existence of a stable unit-winding particle is `open` (numerical). *Per the external
critique: do not speak of the single-charge object as an established particle.*
**Sources:** `core/sections/03_substrate_model.tex` (`eq:cone-speeds`, `eq:W4`),
`gauge_color/sections/03_color.tex` (`eq:em-trace`, `subsec:topology`, Target-1 run),
`gauge_color/derivations/color/vsh_channel_decomposition.md` (§2.5 Derrick radius).
Derived with the `physics-derivation` agent.

---

## 1. Object and ansatz

`π₁(U(1))=ℤ` line vortex of the trace phase `arg Ψ_tr`: a centerline (codim-2) in a 3-D slice, a
**2-D worldsheet** in 4-D. Static, straight, axisymmetric ansatz (cylindrical `ρ,φ` about the
centerline, winding `n∈ℤ`):

    Ψ_tr(ρ,φ) = f(ρ) e^{i n φ},   f(0)=0  (single-valuedness).

This is the `m=1`/`Y₁¹`-type angular sector: the axial symmetry reduces the 2-D transverse problem
to a **1-D radial problem** for `f(ρ)`. This is the spherical-harmonic reduction in the `U(1)`
sector (contrast the `SU(3)` texture, `su3_texture.md`, where SH is not the right basis). The price
of the reduction is the curvilinear connection: the `n²/ρ²` centrifugal barrier in `|∇Ψ_tr|²` below
is the Levi-Civita/Christoffel term of the cylindrical metric `g_{φφ}=ρ²`
(`connections_holonomy.md` §2) — **moving-frame covariantization** (a connection term from the
rotating basis), not an added potential and not a new constitutive nonlinearity.

## 2. The honest point: there is NO `|Ψ|⁴` potential

In Ginzburg–Landau / Nielsen–Olesen the vortex core is regulated by a potential
`V = (β/4)(|Ψ|²−v²)²`, which (i) selects a VEV `v` and (ii) heals `f→v` over `ξ_GL=1/√β v`.

**The central-force spring produces no such term.** The only anharmonicity is the *gradient*-quartic
`W₄ ∝ |∇Ψ|⁴` (`nonlinearity.md` §2); there is no algebraic `V(|Ψ|)`. Consequences (`closed`,
derived from `eq:spring-potential-exact`):

- vacuum is `Ψ=0`, **flat** — no spontaneously selected amplitude `v`;
- **no mass term** `m²|Ψ|²` ⇒ the trace `U(1)` is massless;
- **no Meissner screening** of the trace photon (resolves the screening-tension worry,
  `[[project_massive_photon_screening_tension]]`): the substrate is *not* in a broken/condensate
  phase for the trace sector — consistent with the "carrier-phase is rotation along the timelike
  axis, not a thermodynamic condensate" picture (`[[project_complex_u1_from_time]]`).

We do **not** add a phenomenological `|Ψ|⁴` potential (project rule). We derive what the
gradient-theory gives.

## 3. Closed core energy with the lattice as the physical cutoff

Coarse-grained energy density for the trace envelope:

    W = ½ K |∇Ψ_tr|² + (1/8) Q |∇Ψ_tr|⁴,
    K = ρ_m c_T² = (1−α) k_s/a   (transverse stiffness),
    Q = k_s α/a                  (geometric quartic).

With the ansatz, `|∇Ψ_tr|² = f'² + n²f²/ρ² ≡ G`. In the London limit (`f≈f₀` away from core),
`G = n²f₀²/ρ²` and

    E/L = ∫ 2πρ dρ [ ½ K n²f₀²/ρ² + (1/8) Q n⁴f₀⁴/ρ⁴ ]
        = π K n²f₀² ∫ dρ/ρ  +  (π/4) Q n⁴f₀⁴ ∫ dρ/ρ³ .

The first integral is the standard logarithm. The second integral *would* diverge as `ρ→0` in the
continuum — but the lattice caps it: `ρ ≥ a` and `|ΔR|=√(a²+|Δu|²)` is analytic there
(`nonlinearity.md` §5), so the inner cutoff is `ρ_min = a` **physically, not by fiat**:

    ┌─────────────────────────────────────────────────────────────────────┐
    │  E/L = π(1−α)(k_s/a) n² f₀² ln(R/a)  +  E_core,                       │
    │  E_core(α) = (π/16)(k_s α / a³) n⁴ f₀⁴   ∝ α.                         │
    └─────────────────────────────────────────────────────────────────────┘

**No true singularity:** the term that would diverge is cut at `ρ=a`. `E_core ∝ α` vanishes as
`α→0` (passes the single-nonlinearity check) and is finite for all `α>0`. The log prefactor is the
transverse stiffness `∝(1−α)`; the core energy is the geometric quartic `∝α` — the two ends of the
`α` dial.

## 4. Radial ODE and healing length

Full Euler–Lagrange (quasilinear, field-dependent stiffness `K_eff = K + ½QG`):

    d/dρ [ ρ(K + ½QG) f' ] = ρ(K + ½QG) n²f/ρ²,   G = f'² + n²f²/ρ² .

Every term carries a derivative — there is **no algebraic restoring term**, confirming §2. At large
`ρ` (`G→0`) this reduces to `f'' + f'/ρ − n²f/ρ² = 0` with solutions `f ~ ρ^{±n}` — there is **no
`f→f_∞=const` plateau** (a gradient theory has no length to heal to). The finite-norm physical
object is therefore a vortex on a **localized envelope** (`f→0` at infinity), i.e. a lump carrying a
phase winding — not a condensate vortex. This matches the simulation reality (the Target-1 seed is a
localized shell, `gauge_color/sections/03_color.tex` §6).

The only length is where the quartic-gradient balances the quadratic-gradient at the core:

    ξ = (|n| f₀ / 2) √(α/(1−α))   (units of a).

This is the **same** `√(α/(1−α))` combination as the Derrick hedgehog radius (`eq:derrick-radius`,
`vsh_channel_decomposition.md` §2.5) — an internal consistency check across the two defect sectors.
At the operating point `α=0.2, f₀≈0.3, n=1`: `ξ≈0.075a < a` — the healing length is **sub-lattice**,
so the core pins at the lattice cutoff `ρ_min=a`. To get `ξ≳a` needs `|n|f₀√(α/(1−α))≳2`
(e.g. `α=0.7, f₀≈2`).

## 5. Semilocal stability and the β window (`conditional`)

With `SU(3)` unfrozen the full vacuum `U(3)` is simply connected, so the trace vortex is
**semilocal** (Vachaspati–Achúcarro): topologically protected in the `U(1)` factor, dynamically
(un)stable against spreading into the simply-connected `SU(3)` directions. Standard criterion
`β=(m_s/m_v)²<1` (type-I, stable).

**Caveat (load-bearing):** the substrate `U(1)` is **ungauged** — the Berry connection
`a_μ=i⟨u|∂_μu⟩` is a *diagnostic* read out of the state, never fed back as a force
(No-Back-Reaction, PRINCIPLES §2). So there is **no gauge mass `m_v`** and the GL `β` cannot be
imported literally; the substrate is a *global* `U(3)` σ-model with a gradient-quartic, whose
"vortex" is a global semilocal/`CP²` lump with power-law (`1/r`) tails, not a screened Nielsen–Olesen
profile. Building an effective `β` by dimensional matching of the two stiffnesses (amplitude channel
`Q`, winding channel `K`) at the core scale `ξ`:

- at the single-scale Derrick balance `ξ` cancels and `β_eff = 4/n²`: **unit winding `n=1` gives
  `β_eff=4>1` (unstable to spreading)** — exactly the orientational drift seen in the Target-1 run
  (core drifted ~19 units, winding flipped). `|n|≥2` gives `β_eff≤1` (stable).
- for a lattice-pinned core (`ξ<a`, the actual `α=0.2` situation), `β_eff^pinned = (α/(1−α))(f₀/a)²`,
  which **is** `<1` for small amplitude across essentially the whole `α` window (`α≲0.9` at
  `f₀≈0.3a`); it only crosses to unstable when `f₀≳a`.

So the small-amplitude / lattice-pinned electron vortex is stable; the Derrick-scaled unit-winding
vortex is not. The single-charge `n=1` object sits at the boundary and needs branch-selection
machinery — this is the genuine `open` gap (and the precise reason the existing run did not
converge). **Note:** the `β` construction is a dimensional model, not a derived fluctuation
spectrum; the load-bearing open step is the linearized scalar-vs-vector mode analysis about the
relaxed core.

## 6. Falsifiers

- **Core-energy scaling:** `E_core(α) ∝ α` ⇒ `E_core(0.5)/E_core(0.2)=2.5`, `E_core(0.7)/E_core(0.2)=3.5`
  at fixed `n,f₀`. Flat/decreasing in `α`, or deviation `>15%`, falsifies the geometric-quartic
  regulator (would indicate a hidden `|Ψ|⁴` potential).
- **Tail test (massless EM):** isolated trace charge far field `∝1/r` (power-law), **not**
  `∝e^{−m r}` (screened). An exponential beating power-law (ΔAIC>10) would mean self-screening —
  falsifying the EM identification (flags a W/Z-class object).
- **β verdict:** small-amplitude (`A₀≈0.3`) `n=1` core stays on-axis; raising `A₀≳1` destabilizes
  it (`β_eff→4`). If the small-amplitude `n=1` core also drifts, then `β_eff>1` always at `n=1`,
  forcing `|n|≥2` for the electron sector — a substantive physical consequence.

## 7. Ledger

| Item | Status |
|---|---|
| No `\|Ψ\|⁴` potential / no VEV / no mass / no Meissner | `closed` (derived) |
| `E/L` scaling with lattice `a` as physical cutoff, no singularity | `closed-in-ansatz` (localized envelope) |
| `E_core(α)=(π/16)(k_sα/a³)n⁴f₀⁴ ∝ α` | `closed-in-ansatz` |
| Stable unit-winding particle exists | `open` (numerical) |
| Healing length `ξ=(\|n\|f₀/2)√(α/(1−α))`; sub-lattice at op. point | `closed` |
| Vortex lives on a localized envelope (`f→0`), not a condensate | `closed` |
| `β_eff=4/n²` (Derrick scale) / `(α/(1−α))(f₀/a)²` (pinned) | `conditional` (dimensional model) |
| Converged unit-winding eigenstate / branch-selection | `open` (numerical) |
| Linearized scalar-vs-vector fluctuation spectrum (proper β) | `open` |
