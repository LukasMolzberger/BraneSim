> **Archived 2026-06-28:** valid spring algebra is restated without the unsupported gauge interpretation in `derivations/surviving_results.md`.

# Where the nonlinearity lives, and how it couples the two gauges

**Status:** `closed` (assembly of established results, with one new framing — the lattice as
physical UV regulator).
**Sources:** `core/sections/03_substrate_model.tex` (`eq:spring-potential-exact`, `eq:W4`),
`core/derivations/geometric_nonlinearity_alpha_scaling.md`,
`gauge_color/derivations/color/u1_su3_binding.md` (Channel B vertex).

---

## 1. The single anharmonicity

Each link carries `E_link = ½ k_s (|ΔR| − αa)²`. Split into squared-norm and Euclidean-norm:

    E_link = ½ k_s |ΔR|²  −  k_s α a |ΔR|  +  ½ k_s α² a² .

Because `|ΔR|² = a² + 2a Δu_∥ + |Δu|²` is **exactly quadratic** in node coordinates, the
squared-norm term `½k_s|ΔR|²` is a degree-2 polynomial: it contains no cubic or higher term and no
`α`. The constant `½k_s α²a²` is inert. Therefore:

> **Single-nonlinearity theorem (`closed`).** Every cubic and higher term in the substrate energy
> descends from the one term `−k_s α a |ΔR|`, and is therefore **exactly `∝α`**. At `α=0` the theory
> is exactly linear (Hookean), not approximately so.

This is the consistency check every confinement mechanism below must pass: *it must vanish as
`α→0`.*

## 2. The nonlinearity is intrinsically transverse

Write `Δu = Δu_∥ ê + Δu_⊥` (longitudinal/transverse to the link). Then

    |ΔR| = √(a² + 2aΔu_∥ + |Δu_⊥|²) .

Pure-longitudinal motion (`Δu_⊥=0`) gives `|ΔR| = a + Δu_∥` exactly — **linear**, no square root
bites. The square root (hence all anharmonicity) acts only through the **transverse** displacement
`Δu_⊥`. This is why the lateral triplet `Ψ∈ℂ³` (the transverse channels) is where all gauge
nonlinearity lives, and why the trace↔traceless coupling exists at all.

Expanding on the pure-transverse slice (`Δu_∥` set to its mean):

    −k_s α a |ΔR| = −k_s α a²[ 1 + |Δu_⊥|²/2a² − |Δu_⊥|⁴/8a⁴ + … ] .

- Quadratic piece `−(k_s α/2)|Δu_⊥|²`: feeds the transverse stiffness; **block-diagonal** in
  trace/traceless → sectors decouple at linear order (the established D2 result).
- Quartic piece `+(k_s α/8a²)|Δu_⊥|⁴`: the **geometric quartic**
  `W₄ = (k_s α/8a)|∂u_⊥|⁴` (continuum, `eq:W4`), hardening, `∝α`.

## 3. The trace↔traceless cross-vertex

Split the transverse displacement into trace + traceless, `ξ = ξ_s ê_s + ξ_⊥`,
`|Δu_⊥|² = ξ_s² + |ξ_⊥|²`. The quartic `(ξ_s²+|ξ_⊥|²)²` contains the cross product

    V_{U(1)×SU(3)} = (k_s α / 4a²) · ξ_s² |ξ_⊥|²   (per link),

which cycle-averages to the envelope vertex

    V̄ = g |Ψ_tr|² |Ψ_⊥|²,   g = Θ k_s α / a²,   Θ ∈ [1/16, 3/32].

**The splitting parameter is the coupling parameter:** the single dial `α` both splits the two
sectors' spectral scales and couples them. `g > 0` — the bare vertex is **repulsive** (see
`binding.md`). Geometrically, `g∝α` is the off-diagonal block that makes the internal (Wilczek–Zee)
connection irreducible: it is the dial of the `U(1)³↔U(3)` reducibility transition
(`connections_holonomy.md` §4), and the *same* dial that sets the scale-flow of the anharmonicity
(`scale_transition.md` §2).

## 4. Parity: the nonlinearity is P-even (modulus-based)

`|ΔR| = √(a²+|Δu|²)` depends only on the **modulus** of the displacement. Hence every term it
generates is P-even. Two structural consequences, used repeatedly downstream:

- **(a)** The cross-vertex is `|Ψ_tr|²|Ψ_⊥|²` — a product of moduli-squared, non-negative
  pointwise. No phase-coherent (P-odd) bilinear is produced at any order. *(This kills the
  coherence-pinning binding route — see `binding.md` Route 1.)*
- **(b)** The carrier feels `δω₀ ∝ (k_s α/a²)|Ψ_⊥|²` — the colour **energy** density, not the
  colour **winding** density `Im tr(L³)`. *(This is the Goldstone–Wilczek obstruction — the norm
  vertex sources the wrong, P-even partner; see `binding.md` Channel C.2.)*

## 5. The lattice is the physical UV regulator (no singularities) — `closed`, new framing

This is the analytic content of the owner's "limited resolution prevents infinitely tiny
singularities." The precise statement is about a **spatial phase-gradient / angular wave-number**,
not a temporal frequency.

### 5.1 It is the angular wave number that would diverge, not the carrier frequency

In a continuum curvilinear chart the angular variation of a defect corresponds to a **local spatial
wave number**

    k_ang(r) ~ √(ℓ(ℓ+1)) / r   (spherical harmonic ℓ),     k_φ(ρ) ~ n/ρ   (vortex phase e^{inφ}).

Mapping back to Cartesian and letting `r→0` (or `ρ→0`) makes this angular gradient appear to
diverge. Read naively through a dispersion relation `ω~ck`, that *looks* like an infinite-frequency,
infinite-energy core. **This reading is a category error.** The temporal carrier frequency `ω` is a
separate quantity — fixed externally by worldtube closure (`ωT=2πn`, see `binding.md` §2), not free
to run to infinity. What actually blows up in the continuum is the **spatial** phase gradient, i.e.
the angular wave number `k_ang/k_φ` — a property of the *chart*, not of the dynamics.

### 5.2 The lattice caps the spatial gradient

The norm `|ΔR| = √(a² + |Δu|²)` is **analytic and bounded for every configuration**. A field
gradient on the lattice is a finite difference of node displacements over one link of length `a`:

    |∇Ψ| ≤ |ΔΨ|_max / a ~ f₀/a   (saturates at separation ~a).

The angular circumference `2πr` eventually falls below what the lattice can resolve: below `r~a`
there is no continuum ring of independent angular samples, so `k_ang` cannot exceed `~1/a`.
Concretely, the discrete circulation around the smallest plaquette advances the trace phase by `2πn`
over a perimeter `~8a`, so the maximal per-link phase gradient is `≤ 2πn/8a` — finite. A topological
defect therefore has a **finite** core energy (the would-be continuum `∫dρ/ρ³` divergence is cut at
`ρ=a` *physically*, not by a hand-chosen regulator; `u1_vortex_core.md` §3) and **no** true
singularity at the centerline / texture point.

### 5.3 The core claim

> A particle is **not** a continuum singularity at the center of a spherical-harmonic expansion. It
> is a **lattice-regulated topological defect**: the place where the continuum angular description
> demands more resolution than the substrate possesses. The apparent `1/r` or `1/ρ` divergence is
> therefore not a physical infinity but the **failure of the spherical/cylindrical chart at the
> core**. The lattice replaces the continuum point by a finite core of radius `O(a)`, while the
> topological winding remains well-defined outside that core.

This also replaces, on the UV side, the continuum Derrick collapse: a soliton cannot shrink below
one cell, because energy concentrated at a single node is immediately dispersed by the springs
(Paper IV `subsec:derrick`). The lattice spacing `a` is the hard inner length of every confinement
formula in this paper.

## 6. Ledger

| Item | Status |
|---|---|
| All anharmonicity `∝α` from the one norm term | `closed` |
| Anharmonicity is transverse-only | `closed` |
| Cross-vertex `V̄=g|Ψ_tr|²|Ψ_⊥|²`, `g>0` repulsive | `closed` |
| P-even ⇒ no phase-coherent bilinear; sources energy not winding | `closed` |
| Lattice `a` = physical inner cutoff, no singularities | `closed` |
| `Θ` (carrier-phase correlation `1/16`–`3/32`) | `open` (O(1), undetermined) |
