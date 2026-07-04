# Field-Strength Tensors from the Spring Stencil (Paper VII) — Problem Statement

A drill-down of Paper III (`gauge_color`). Paper III *identifies* the electromagnetic field tensor
with a Berry curvature and the colour field tensor with a Wilczek–Zee curvature, and treats that
identification as formal. This paper does the job the identification leaves implicit: **construct the
carrier on which those curvatures live directly from the substrate's spring stencil, build the Faraday
tensor `F_μν` and the QCD field-strength tensor `G^a_μν` as the continuum limit of its plaquette
holonomies, and show the construction reproduces — and thereby grounds — the Berry/Wilczek–Zee
identification of Paper III.**

This is a *kinematic* paper: it builds the carrier, the two antisymmetric rank-2 curvature objects, and
their algebraic properties (antisymmetry, Bianchi identity, sector split, non-Abelian self-coupling). It
does not build their dynamics.

## Where this sits in the layered structure

The emergence program is stratified. Each object must be assigned to the layer where it actually lives,
and the field strength is a *linear-layer* object:

- **Substrate layer.** The 4D brane lattice: real node embeddings `R(n) ∈ ℝ⁴`, central-force springs,
  the discrete action. Ontic (Paper I).
- **Linear / carrier layer.** Bloch modes, polarization branches, the phase-space carrier `ψ`, its Berry/
  Wilczek–Zee connection over `(x,t)`. **`F_μν` and `G^a_μν` are built here.** This is a description of
  small excitations on the substrate, not a new structural commitment of the substrate.
- **Soliton layer.** Bound, non-radiating localized modes: vortex cores (charge, `π₁(U(1))=ℤ`), colour
  textures (`π₃(SU(3))`), and the spin-½ `ℤ₂` framing holonomy (`π₁(SO(3))=ℤ₂`). **Everything topological
  and everything about spin-½ is here, and is out of scope for this paper** (`matter_mass`, Paper IV;
  `[[project_spin_half_is_soliton_layer]]`, `[[project_soliton_layer_description_language]]`).

The discipline this buys: the field strength is the curvature of a *linear-layer* connection; the
`ℤ₂` objects that look adjacent to it (Zak phase of the spin-1 envelope; spin-½ of a bound soliton) are
not field strengths and belong to other layers. Keeping the layers apart is what keeps the paper honest.

## The substrate and its stencil (interface from Paper I)

The substrate is a 4D hypercubic brane lattice, spacing `a`. Each node `n` carries a real embedding
coordinate `R(n) ∈ ℝ⁴` and is connected by **central-force springs along 8 axial links**:

- 6 spacelike links `±ê_x, ±ê_y, ±ê_z` (the canonical 6-neighbour axial stencil, `PRINCIPLES.md`);
- 2 temporal links `±ê_t` (the timelike link is a 4D-isotropic spring like the spatial ones,
  parameterized by `r_t`; `[[project_temporal_link_4d_spring]]`).

Together: **8 axial links per node**, one ordered link family per signed lattice direction `±μ`,
`μ ∈ {t,x,y,z}`. The Lorentzian structure enters through the *intrinsic lattice metric*
`η_{μν} = diag(−1,+1,+1,+1)` on the stencil directions, i.e. through the sign `s_μ` in the link coupling
— **not** through making the ambient `ℝ⁴` itself Lorentzian. Each spring has a stress-free **rest length
`ℓ₀`**; with the held ground-state scale `h⋆` this defines the **prestress (rest-length) parameter**
`α := ℓ₀/h⋆ ∈ (0,1]` — the single physical knob. The only anharmonicity is the Euclidean norm term
`−k_s α a |ΔR|`, exactly `∝α` (`[[project_geometric_nonlinearity_alpha_scaling]]`); `α=0` is the
exactly-linear limit.

## The central object: a field strength is a plaquette curvature

A field-strength tensor `F_μν` is, structurally, the curvature 2-form of a connection: the holonomy of
parallel transport around an infinitesimal oriented loop in the `(μ,ν)` plane. On a lattice this loop
is a **plaquette** — the ordered product of four link variables around an elementary square. The program
is the standard lattice-gauge dictionary, run *backwards* from a physical substrate:

```
spring links → Bloch carrier ψ → carrier link variables U_μ(n) → plaquette □_μν(n)
             → log □_μν = −ia²F_μν + O(a⁴)
```

with `μ,ν` ranging over the 8-link directions. The two physical field tensors are the two sectors of
this single `U(3)` construction (the `U(3) = U(1) ⊕ SU(3)` carrier split of Paper III).

## The chain from the substrate, stated carefully

The construction has a subtle centre, and it is worth stating the whole chain the way it actually goes,
because two natural shortcuts are wrong.

1. **Bloch modes carry one `U(1)` phase, not four.** A linearized perturbation is a Bloch mode
   `δR_p^A(k) = u^A(k)\,e^{ik_μ p^μ}`, `k_μ = (−ω, k_1, k_2, k_3)`. The exponential `e^{ik_μ p^μ}` is a
   *single* complex phase `θ` living in one circle `θ ∼ θ + 2π`; it has four directional gradients
   `∂_μ θ = k_μ`, but it is **one `U(1)`, not `U(1)⁴`.** (See `carrier_construction.md` §1.)

2. **The Berry connection comes from the eigenvector's phase freedom, not the raw plane-wave phase.** The
   polarization eigenvector `u(k)` is defined only up to a phase `u(k) → e^{iχ(k)}u(k)`; *that* phase
   freedom is the gauge `U(1)`, and the Berry connection is `A_μ(k) = i⟨u(k)|∂_{k_μ}u(k)⟩`. A pure scalar
   Bloch wave with constant `u(k)` has trivial curvature. Curvature needs `u(k)` to genuinely twist as
   `k` moves. (See `carrier_construction.md` §1.)

3. **The rest length is the knob that makes `u(k)` twist — necessary, not sufficient.** The rest length
   `ℓ₀` sets the *transverse* bond stiffness through the pre-tension factor `(1 − ℓ_μ/L_μ)` in the bond
   matrix `K_μ`. Non-commuting `K_μ` cannot be diagonalized in one basis, so `u(k)` becomes a genuinely
   `k`-dependent mixture of the ambient displacement components — the ingredient a Berry connection needs.
   **But this is necessary, not sufficient:** `D(k) = Σ_μ f_μ(k) K_μ` is real symmetric, so `u(k)` can be
   chosen real at every `k`, and a real (even if `k`-dependent) eigenvector gives
   `A_μ = i⟨u|∂u⟩ = 0`. The Brillouin-zone `U(1)` curvature is identically zero `∀α`
   (BACKBONE #16; `[[project_spin_half_is_soliton_layer]]`). (See `carrier_construction.md` §1a.)

4. **The genuine curvature lives in the complex carrier over physical `(x,t)`.** The complex structure
   the gauge fields need is not in `R` and not in the BZ. It is the phase-space amplitude
   `ψ^a = δR^a + i\,δṘ^a/ω`: position and velocity, 90° out of phase, packaged into one complex number.
   The `i` is literally the carrier's rotation along the timelike worldtube axis — velocity is the
   time-derivative, so complexifying *is* encoding the time-link rotation (`[[project_complex_u1_from_time]]`).
   This is why the connection's base is `(x,t)`, not the BZ, and why the `E`-field components `F_{0i}`
   are time–space plaquettes, impossible without the 2 temporal links. (See `carrier_construction.md`
   §§2–6, `berry_reconciliation.md` §3.)

The one-line summary of the subtlety: **the rest-length/stiffness geometry supplies the `k`-mixing that a
Berry connection needs, but only the complex time-carrier supplies a nonzero curvature.** Both halves are
in the paper, honestly joined.

## Questions this paper must answer in closed form

1. **Carrier.** From the real `(R, Ṙ)` phase space, construct the complex carrier `ψ`; show the real
   degenerate triplet has only `O(3)` rotation freedom and that complexification is what promotes it to
   `U(3)`; distinguish the single-vector `U(1)` phase from the transported three-frame `U(3)` connection
   (do **not** conflate a vector in `ℂ³` with a matrix in `u(3)`).
2. **Link variables.** From the carrier and the 8-link stencil, define `U_μ(n)` rigorously with its gauge
   (carrier-rephasing) law; show spacelike and timelike links enter on the same footing.
3. **Plaquette → curvature.** Show `□_μν = exp(−ia²F_μν + O(a⁴))` with the leading term **antisymmetric**
   in `μν` — the structural origin of the tensor's antisymmetry, not an assumption.
4. **Faraday tensor (trace `U(1)`).** `F_μν = ∂_μA_ν − ∂_νA_μ`, `E_i = F_{0i}`, `B_i = ½ε_{ijk}F_{jk}`;
   homogeneous Maxwell `∂_{[λ}F_{μν]} = 0` as the exact lattice (closed-cube) Bianchi identity.
5. **QCD field-strength tensor (traceless `SU(3)`).** `G^a_μν = ∂_μA^a_ν − ∂_νA^a_μ + f^{abc}A^b_μA^c_ν`;
   the `f^{abc}` self-coupling as the non-commutativity of the ordered link product; non-Abelian Bianchi.
6. **The role of the rest length.** From the triplet spectrum `λ_A = κ[(1−α)|k|² + α k_A²]`, show the
   trace susceptibility `∝ (3−2α)` and the traceless splitting `∝ α g(k̂)`, so the colour curvature weight
   grows `∝ α` off `[111]` and vanishes at `α=0`. Note the *near-degeneracy tension*: exact degeneracy
   (and `U(3)` frame freedom) is the `α=0` limit, and the splitting that turns colour on is what breaks it
   — degeneracy and colour move in opposite directions (`[111]` coherence-vs-colour, kinematic
   confinement). Derive the normalization-independent ratio `R(0.5)/R(0.2) = 0.3077`, and state the D2
   caveat (`α` fixes the *relation/ratio*, not `α_EM/α_s`).
7. **Reconciliation with Paper III.** Prove the continuum limit of the `(x,t)` plaquette holonomy equals
   the Berry curvature `f_μν` / Wilczek–Zee curvature `𝓕_μν`; restate why the BZ curvature vanishes.
8. **Diagnostics (candidate section).** The same plaquette/Wilson-loop objects, computed on the existing
   substrate graph, are a reusable diagnostic toolbox borrowed from lattice gauge theory — *not* a second
   lattice. What is directly reusable, partially reusable, and not (see `lgt_diagnostics.md`).

## Scope boundary

- **In scope:** the carrier construction; `F_μν` and `G^a_μν` as lattice-plaquette curvatures; their
  antisymmetry, sector split, Bianchi identities; the `α`-activation of the colour sector; equivalence to
  Berry/WZ; the LGT-diagnostic correspondence on the same lattice.
- **Out of scope (deferred, not claimed):** field *dynamics* — `−¼F²`/`−¼G²` action, inhomogeneous
  equations `d⋆F = J`, masslessness, propagation speed `c` (D3, open); the coupling-constant hierarchy /
  value of `α` (D2, undetermined at linear order, `[[project_alpha_undetermined_at_linear_order]]`);
  confinement (archived); **soliton-layer objects — vortex cores, colour textures, and the spin-½
  `ℤ₂` framing holonomy — which belong to Paper IV** (`[[project_spin_half_is_soliton_layer]]`).
- This is **not** a derivation of QCD: no running coupling, asymptotic freedom, string tension, or hadron
  spectrum. Only the field-strength *tensor* — the curvature 2-form — is constructed here.