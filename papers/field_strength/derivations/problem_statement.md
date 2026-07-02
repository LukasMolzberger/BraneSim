# Field-Strength Tensors from the Spring Stencil (Paper VII) — Problem Statement

A drill-down of Paper III (`gauge_color`). Paper III *identifies* the electromagnetic field tensor
with the Berry curvature and the colour field tensor with the Wilczek–Zee curvature, and treats this
identification as formal. This paper does the one job that identification leaves implicit: **construct
the carrier on which those curvatures live directly from the substrate's spring stencil, then build the
Faraday tensor `F_μν` and the QCD field-strength tensor `G^a_μν` as the continuum limit of its
plaquette holonomies**, and show the construction reproduces (and grounds) the Berry/Wilczek–Zee
identification of Paper III.

This is a *kinematic* paper. It builds the carrier, the two antisymmetric rank-2 curvature objects, and
their algebraic properties (antisymmetry, Bianchi identity, sector split, non-Abelian self-coupling).


## What the second attempt fixes

The first attempt built the plaquette dictionary correctly but left two things hand-wavy. This attempt
makes them load-bearing:

1. **The complex carrier is earned, not assumed.** The substrate state is real
   (`R(n), Ṙ(n) ∈ ℝ⁴`). The complex carrier is the *phase-space amplitude* of the oscillator,
   `ψ^a = δR^a + i δṘ^a/ω` — position and velocity, 90° out of phase, packaged into one complex number.
   The `i` is then literally the carrier's rotation along the timelike worldtube axis
   (`[[project_complex_u1_from_time]]`): velocity is the time-derivative, so complexifying *is* encoding
   the time-link rotation. See `carrier_construction.md`. This answers the reviewer's question "where
   did `i` come from?" without invoking magic.

2. **`α` activates the colour sector — one parameter, one direction.** The rest-length parameter
   `α := ℓ₀/h⋆ ∈ (0,1]` controls *both* which gauge group the carrier carries *and* how strongly the
   colour curvature appears, and both increase **together** with `α`. At `α=0` (the linear limit) the
   three lateral carriers decouple: the connection is Abelian `U(1)³` and the traceless holonomy,
   weighted `∝α`, vanishes. As `α→1` they lock into the degenerate triplet with full `U(3)` frame
   freedom and the traceless `SU(3)` holonomy turns on, weighted `∝α` against the trace's `∝(3−2α)`.
   These are two descriptions of one monotonic activation, not competing effects: there is no preferred
   interior value of `α`. The clean theoretical statement is the two limits, not a fitted point. See
   `field_tensors.md` §4.

## The substrate and its stencil (interface from Paper I)

The substrate is a 4D hypercubic brane lattice, spacing `a`. Each node `n` carries a real embedding
coordinate `R(n) ∈ ℝ⁴` and is connected by **central-force springs along 8 axial links**:

- 6 spacelike links `±ê_x, ±ê_y, ±ê_z` (the canonical 6-neighbour axial stencil, `PRINCIPLES.md`);
- 2 temporal links `±ê_t` (the timelike link is a 4D-isotropic spring like the spatial ones,
  parameterized by `r_t`; `[[project_temporal_link_4d_spring]]`).

Together: **8 axial links per node**, one ordered link family per signed lattice direction `±μ`,
`μ∈{t,x,y,z}`. Each spring has a stress-free **rest length `ℓ₀`**; with the held ground-state scale
`h⋆` this defines the **prestress (rest-length) parameter** `α := ℓ₀/h⋆ ∈ (0,1]` — the single physical
knob. The only anharmonicity is the Euclidean norm term `−k_s α a |ΔR|`, exactly `∝α`
(`[[project_geometric_nonlinearity_alpha_scaling]]`); `α=0` is the exactly-linear limit.

## The central object: a field strength is a plaquette curvature

A field-strength tensor `F_μν` is, structurally, the curvature 2-form of a connection: the holonomy of
parallel transport around an infinitesimal oriented loop in the `(μ,ν)` plane. On a lattice this loop
is a **plaquette** — the ordered product of four link variables around an elementary square. The program
is the standard lattice-gauge dictionary, run *backwards* from a physical substrate:

```
spring links → phase-space carrier ψ → carrier link variables U_μ(n) → plaquette □_μν(n)
             → log □_μν = ia²F_μν + O(a⁴)
```

with `μ,ν` ranging over the 8-link directions. The two physical field tensors are the two sectors of
this single `U(3)` construction (the `U(3)=U(1)⊕SU(3)` carrier split of Paper III).

## The load-bearing subtlety (must be confronted, not hidden)

The substrate displacement field `R(n)` is **real**, and a plaquette built from the *displacement
geometry* carries trivial curvature: the dynamical matrix `D(k)` is real symmetric, so the
**Brillouin-zone Berry curvature is identically zero ∀α** (BACKBONE #16;
`[[project_spin_half_is_soliton_layer]]`). A naïve "plaquette of the spring lattice" gives `F≡0`. The
nonzero field strength cannot live in the displacement geometry or in the Brillouin zone.

It lives in the **carrier phase over physical `(x,t)`**. The phase-space carrier `ψ` carries a genuine
complex phase (the time-link rotation), and the link variable is the *ordered overlap of the carrier
polarization eigenframe between adjacent nodes*. The plaquette of these overlaps is the object whose
continuum log is `a²F_μν`. This is exactly why Paper III insists the connection's base is `(x,t)`, not
the BZ. The temporal link is not decoration: the `E`-field components `F_{0i}` are time–space
plaquettes, impossible without the 2 temporal links (`[[project_complex_u1_from_time]]`: two time
slices are minimal).

## Questions this paper must answer in closed form

1. **Carrier.** Construct the complex carrier `ψ` from the real `(R, Ṙ)` phase space; show the real
   degenerate triplet has only an `O(3)` rotation freedom and that complexification is what promotes it
   to `U(3)`; distinguish the single-vector `U(1)` phase from the transported three-frame `U(3)`
   connection (do **not** conflate a vector in `ℂ³` with a matrix in `u(3)`).
2. **Link variables.** From the carrier and the 8-link stencil, define `U_μ(n)` rigorously with its
   gauge (carrier-rephasing) law; show spacelike and timelike links enter on the same footing.
3. **Plaquette → curvature.** Show `□_μν=exp(ia²F_μν+O(a⁴))` with the leading term **antisymmetric** in
   `μν` — the structural origin of the tensor's antisymmetry, not an assumption.
4. **Faraday tensor (trace `U(1)`).** `F_μν=∂_μA_ν−∂_νA_μ`, `E_i=F_{0i}`, `B_i=½ε_{ijk}F_{jk}`;
   homogeneous Maxwell `∂_{[λ}F_{μν]}=0` as the exact lattice (closed-cube) Bianchi identity.
5. **QCD field-strength tensor (traceless `SU(3)`).** `G^a_μν=∂_μA^a_ν−∂_νA^a_μ+f^{abc}A^b_μA^c_ν`; the
   `f^{abc}` self-coupling as the non-commutativity of the ordered link product; non-Abelian Bianchi.
6. **The role of the rest length.** Show the single parameter `α` activates the colour sector
   monotonically: `α=0` (linear) ⇒ Abelian `U(1)³`, traceless weight `∝α → 0`; `α→1` ⇒ degenerate
   triplet, full `U(3)`, traceless `SU(3)` on, weight `(3−2α):α`. Degeneracy and curvature-weight move
   the same way — no tension, no interior sweet spot. State the linear-order caveat (D2: `α` fixes the
   *relation/ratio*, not `α_EM/α_s`).
7. **Reconciliation with Paper III.** Prove the continuum limit of the `(x,t)` plaquette holonomy equals
   the Berry curvature `f_μν` / Wilczek–Zee curvature `𝓕_μν`; restate why the BZ curvature vanishes.

## Scope boundary

- **In scope:** the carrier construction; `F_μν` and `G^a_μν` as lattice-plaquette curvatures; their
  antisymmetry, sector split, Bianchi identities; the `α`-activation of the colour sector; equivalence
  to Berry/WZ.
- **Out of scope (deferred, not claimed):** field *dynamics* — `−¼F²`/`−¼G²` action, inhomogeneous
  equations `d⋆F=J`, masslessness, propagation speed `c` (D3, open); the coupling-constant hierarchy /
  value of `α` (D2, undetermined at linear order, `[[project_alpha_undetermined_at_linear_order]]`);
  confinement (archived 2026-06-28).
- This is **not** a derivation of QCD: no running coupling, asymptotic freedom, string tension, or
  hadron spectrum. Only the field-strength *tensor* — the curvature 2-form — is constructed here.