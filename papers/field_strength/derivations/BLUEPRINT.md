# Field-Strength Paper (Paper VII) — Reconstruction Blueprint

The minimal set of pieces sufficient to build the paper. Each piece names the claim, its key equation, and
which derivation note owns it. Read top to bottom: it is the logical spine in order.

## 0. Thesis (one sentence)

Paper III *identifies* the EM field tensor with a Berry curvature and the colour field tensor with a
Wilczek–Zee curvature; this paper **constructs the carrier those curvatures live on directly from the
spring stencil, builds `F_μν` and `G^a_μν` as the continuum limit of its plaquette holonomies, and proves
the construction equals the Berry/WZ object** — a *kinematic* result (tensors + algebra), not dynamics.

## 1. Layer discipline (frame the whole paper)

The field strength is a **linear / carrier-layer** object. Substrate layer = the lattice (ontic).
Soliton layer = vortex cores, colour textures, and the spin-½ `ℤ₂` framing holonomy (Paper IV). Keeping
the layers apart is what lets the paper build `F`, `G` without dragging in charge quantization or spin-½.
(`[[project_spin_half_is_soliton_layer]]`, `[[project_soliton_layer_description_language]]`.)

## 2. Substrate interface (given, from Papers I, III)

- 4D hypercubic brane lattice, spacing `a`. Node carries real `R(n) ∈ ℝ⁴`, velocity `Ṙ(n) ∈ ℝ⁴`.
- **8 axial central-force springs/node**: 6 spacelike `±ê_{x,y,z}` + 2 temporal `±ê_t` (timelike link is a
  4D-isotropic spring, `r_t`; `[[project_temporal_link_4d_spring]]`). Lorentzian sign `s_μ=(−1,+,+,+)` on
  the *intrinsic* stencil metric, not on ambient `ℝ⁴`.
- One knob: rest-length / prestress `α := ℓ₀/h⋆ ∈ (0,1]`. Sole anharmonicity `−k_sαa|ΔR| ∝ α`; `α=0` is
  the exactly-linear limit (`[[project_geometric_nonlinearity_alpha_scaling]]`).

## 3. The Bloch → Berry scaffold and the obstruction (owns: `carrier_construction.md` §§1–1b)

- **Bloch mode** `δR_p^A(k)=u^A(k)e^{ik_μp^μ}`, `k_μ=(−ω,k_1,k_2,k_3)`. The exponential is **one `U(1)`**
  (four gradients `∂_μθ=k_μ`, one circle) — *not* `U(1)⁴`.
- **Berry connection is the eigenvector's phase freedom** `u(k)→e^{iχ(k)}u(k)`, `A_μ=i⟨u|∂_{k_μ}u⟩` — not
  the raw plane-wave phase. Needs `u(k)` to twist.
- **Rest length makes `u(k)` twist (necessary).** Bond stiffness
  `K_μ^{AB}=κ_μ[n̂^An̂^B+(1−ℓ_μ/L_μ)(δ^{AB}−n̂^An̂^B)]`; transverse weight `(1−α)`. Non-commuting `K_μ`
  ⇒ `k`-dependent polarization mixing (`D(k)=Σ_μ f_μ(k)K_μ`).
- **But not sufficient.** `R` real ⇒ `D(k)` real symmetric ⇒ eigenvectors real ⇒ **BZ Berry/WZ curvature
  ≡ 0 ∀α** (BACKBONE #16). A real, even `k`-dependent, `u(k)` gives `A_μ=0`. Real `k`-bundles carry at
  most a `ℤ₂` (Zak) holonomy — the **spin-1 envelope's**, not spin-½. **So the field strength cannot live
  in the displacement geometry or the BZ — it must live in the carrier phase over physical `(x,t)`.**

## 4. Earn the complex carrier from phase space — the load-bearing piece (owns: `carrier_construction.md` §§2–6)

1. **Phase-space amplitude:** `ψ^a_n = δR^a_n + i δṘ^a_n/ω`. The `i` is not assumed — `δṘ` is the time-
   derivative, so complexifying *is* the carrier's rotation along the timelike worldtube axis
   (`[[project_complex_u1_from_time]]`). Free mode `ψ~e^{−iωt}`. Two time slices minimal.
2. **Carrier space:** lateral triplet ⇒ `ψ_n=(ψ^x,ψ^y,ψ^z)∈ℂ³`. Linear envelope is spin-1, so `ℂ³` is
   natural.
3. **`O(3)→U(3)`:** a *real* degenerate triplet has only `O(3)` frame freedom; complexifying to `ℂ³`
   enlarges it to `U(3)` (`ψ→Uψ`, `U†U=𝟙`). Complexification *is* what upgrades `O(3)` to `U(3)`.
4. **Vector ≠ connection:** a single carrier vector's overall phase → `U(1)`; a transported three-frame's
   overlap `M_{ij}=⟨u_i|u_j⟩` → `U(3)` connection (`u(3)`-valued). Carrier stays a vector.
5. **Dynamic phase ≠ curvature:** the clock `e^{−iωt}` gives plaquette holonomy `=1`, `F=0`. Curvature
   needs the *relative* frame between neighbours path-dependent — which it is because the complex link is
   the time-rotation.

## 5. Link → plaquette → curvature (owns: `link_holonomy.md`)

- **Link** on `n→n+ê_μ` (all μ, space- and timelike on the same footing): `U_μ(n)=⟨u(n)|u(n+ê_μ)⟩/|·|`
  (rank-1) or polar-unitary part of `M_{ij}` (rank-`d`, ∈`U(d)`). Reverse `U_{-μ}(n+ê_μ)=U_μ(n)^†`.
- **Gauge law** (carrier rephasing `u→Vu`): `U_μ(n)→V(n)U_μ(n)V(n+ê_μ)^†` — gauge symmetry *is* the local
  eigenframe redundancy.
- **Plaquette** `□_μν=U_μ(n)U_ν(n+ê_μ)U_μ(n+ê_ν)^†U_ν(n)^†`. **Continuum:** `□_μν=exp(−ia²F_μν+O(a⁴))`,
  `F_μν=∂_μA_ν−∂_νA_μ−i[A_μ,A_ν]`.
- **Antisymmetry structural:** `□_νμ=□_μν^†` ⇒ `F_νμ=−F_μν`. A 2D oriented loop can only make a rank-2
  antisymmetric tensor — *why* both field tensors are antisymmetric.

## 6. The two sectors from the one `U(3)` plaquette (owns: `field_tensors.md`)

Split `A_μ=A^0_μ(𝟙/√6)+A^a_μT^a`, `T^a=½λ^a`.

- **Faraday (trace `U(1)`).** `F_μν=∂_μA_ν−∂_νA_μ`; `E_i=F_{0i}` (needs a temporal link), `B_i=½ε_{ijk}F_{jk}`.
  **Homogeneous Maxwell = exact lattice Bianchi:** six face-plaquettes of an elementary 3-cube multiply to
  `1` ⇒ `∂_{[λ}F_{μν]}=0` exactly (`∇·B=0`, `∇×E+∂_tB=0`).
- **Yang–Mills (traceless `SU(3)`).** `G^a_μν=∂_μA^a_ν−∂_νA^a_μ+f^{abc}A^b_μA^c_ν`. The `f^{abc}` self-
  coupling **is** `−i[A_μ,A_ν]` projected to the adjoint — the gluon vertex = non-commutativity of the
  ordered link product. Only structural difference from EM: `d=3` vs `d=1`. Non-Abelian Bianchi
  `D_{[λ}G_{μν]}=0` by the same closed-cube argument with path-ordering.

## 7. The `α` spectrum — colour and degeneracy move OPPOSITE ways (owns: `field_tensors.md` §4)

Triplet spectrum `λ_A(k) = κ[(1−α)|k|² + α k_A²]`; trace `λ̄ ∝ (3−2α)`, traceless `λ_A−λ̄ ∝ α g(k̂)`.

```
 α = 0 : EXACT degeneracy → U(3) frame freedom, but traceless splitting = 0 ⇒ colour curvature FLAT (EM on, colour OFF)
 α → 1 : MAXIMAL splitting → non-degenerate ⇒ traceless SU(3) curvature MAX (∝ α) (colour ON)
```

- **Degeneracy is maximal at `α=0` and the splitting/colour grows with `α` — opposite directions.** The
  `U(3)` frame freedom needs degeneracy; the colour curvature needs the splitting that breaks it. This is a
  genuine **near-degeneracy tension** (`[111]` coherence-vs-colour, `g([111])=0` — kinematic confinement),
  not a smooth co-activation. WZ is the standard *near-degenerate* construction.
- The one monotonic, unambiguous thing is the **colour curvature weight `∝ α`** off `[111]`, vanishing at
  `α=0`. No interior sweet spot, no fitted `α`.
- **Weighting derived here** from our own `K_μ`: `R(α) ∝ (3−2α)/α`, normalization-independent, giving
  `R(0.5)/R(0.2) = 0.3077` — a genuine result of this paper that reproduces the colour bridge.
- **Caveat (D2):** `R(α)` is a spectral (dispersion) ratio; identifying it with a physical coupling ratio
  `α_EM/α_s` needs the normalization `C` (undetermined at linear order). `α` fixes the *relation/ratio*,
  not the coupling hierarchy (`[[project_alpha_undetermined_at_linear_order]]`).

## 8. Equivalence theorem (owns: `berry_reconciliation.md`)

Expand the overlap link: `U_μ(n)=exp(−ia·i⟨u|∂_μu⟩+O(a²))=exp(−ia·a_μ+…)`, so `□_μν` is the Wilson loop
of the Berry connection and (Stokes / non-Abelian Stokes):
`log □_μν = −ia²f_μν+O(a⁴)` (Abelian), `−ia²𝓕_μν+O(a⁴)` (non-Abelian). **Therefore `F_μν=f_μν`,
`G_μν=𝓕_μν`** — plaquette field strength = Berry/WZ curvature; one object, integral vs differential.
BZ-vanishing restated from this side (real links ⇒ trivial plaquette; complex link = time-rotation).

## 9. Diagnostics on the same graph (owns: `lgt_diagnostics.md`) — candidate section

Effective links `U_i^{eff}(x)=M_i(x)(M_i†M_i)^{−1/2}∈U(N)` from the carrier eigenframe overlap on the
*existing* substrate links (not a second lattice). `N=1`→`U(1)`; `N=3`→`U(3)≃U(1)×SU(3)/ℤ₃`. Plaquette
`U_{ij}^{eff}` = Berry-curvature diagnostic; Wilson loops, topological charge, finite-size scaling.
Same lattice artifacts as LQCD (BZ, cutoff) but ontology, not regulator (`a→0` not mandatory). Reusable:
Wilson loops, smearing (diagnostic only), blocking/RG, finite-size scaling, topological charge, HPC
infra. Partial: HMC (init/BC ensembles only), gauge fixing (compute invariants), Wilson-action fit
(derived, not postulate). Not: Haar path integral, Euclidean gauge action as microscopic law, fundamental
quarks, mandatory `a→0`. **Section-vs-cross-reference placement is open.**

## 10. Scope

- **In:** carrier construction; `F_μν`, `G^a_μν` as plaquette curvatures; antisymmetry; sector split;
  Bianchi; the `α`-activation of colour; equivalence to Berry/WZ; LGT diagnostics on the same graph.
- **Out (deferred, not claimed):** dynamics (D3 — `−¼F²`/`−¼G²`, `d⋆F=J`, masslessness, `c`); coupling /
  `α` hierarchy (D2); confinement (archived); soliton-layer objects (vortex cores, colour textures,
  spin-½ — Paper IV). **Not** a derivation of QCD.

## 11. Open items (status.md "Open derivations")

1. Link-variable definedness across carrier band-crossings (band-isolation; for colour, essentially the
   soliton layer — `g([111])=0` has no defined holonomy).
2. `O(a⁴)` is a genuine higher-derivative term (clover/improvement check).
3. Non-Abelian ordering convergence (`f^{abc}` the unique non-commuting remnant).
4. Numerical equivalence: transport a `ℂ³` carrier around a fixed `(x,t)` plaquette, confirm
   `log □ = −ia²F` vs `branesim/diagnostics/berry_holonomy.py`. Entry point for the `lgt_diagnostics.md`
   §6 pipeline.
5. Identification of the spectral ratio `R(α) ∝ (3−2α)/α` as a physical coupling ratio (D2) — needs the
   connection normalization `C`, undetermined at linear order. The `0.3077` *ratio* is derived; the
   *coupling* reading is deferred.

## 12. Honesty headline

Closes the **kinematic construction** of `F_μν` and `G^a_μν` from the phase-space carrier on the 8-link
stencil, proves it is Paper III's Berry/WZ object, and shows the LGT diagnostic toolbox applies on the
same substrate graph. Does **not** close **dynamics** (no action/field equations/propagation) or
**coupling** (`α` undetermined; fixes relation, not `α_EM/α_s`).

## File map

`problem_statement.md` (scope+motivation, layer discipline) · `carrier_construction.md` (§§3–4) ·
`link_holonomy.md` (§5) · `field_tensors.md` (§§6–7) · `berry_reconciliation.md` (§8) ·
`lgt_diagnostics.md` (§9) · `status.md` (props, open items, risks).
