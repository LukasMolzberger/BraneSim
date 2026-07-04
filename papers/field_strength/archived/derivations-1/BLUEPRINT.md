# Field-Strength Paper (Paper VII) — Reconstruction Blueprint

The minimal set of pieces sufficient to rebuild the paper from scratch. Each piece names the claim, its
key equation, and which derivation note owns it. Read top to bottom: it is the logical spine in order.

## 0. Thesis (one sentence)

Paper III *identifies* the EM field tensor with a Berry curvature and the colour field tensor with a
Wilczek–Zee curvature; this paper **constructs the carrier those curvatures live on directly from the
spring stencil, builds `F_μν` and `G^a_μν` as the continuum limit of its plaquette holonomies, and
proves the construction equals the Berry/WZ object** — a *kinematic* result (tensors + algebra), not
dynamics.

## 1. Substrate interface (given, from Papers I, III)

- 4D hypercubic brane lattice, spacing `a`. Node carries real `R(n) ∈ ℝ⁴`, velocity `Ṙ(n) ∈ ℝ⁴`.
- **8 axial central-force springs/node**: 6 spacelike `±ê_{x,y,z}` + 2 temporal `±ê_t` (timelike link is
  a 4D-isotropic spring, `r_t`; `[[project_temporal_link_4d_spring]]`).
- One knob: rest-length / prestress `α := ℓ₀/h⋆ ∈ (0,1]`. Sole anharmonicity `−k_sαa|ΔR| ∝α`;
  `α=0` is the exactly-linear limit (`[[project_geometric_nonlinearity_alpha_scaling]]`).

## 2. The obstruction that forces the whole construction

`R(n)` is real ⇒ `D(k)` real symmetric ⇒ eigenvectors choosable real ⇒ **BZ Berry/WZ curvature ≡ 0
∀α** (BACKBONE #16; `[[project_spin_half_is_soliton_layer]]`). A naïve plaquette of the displacement
lattice gives `F≡0`. Note the subtlety that must be pre-empted: the rest length *does* set the transverse
bond stiffness `(1−α)` and can make the polarization eigenvector `u(k)` genuinely `k`-dependent, but that
`k`-mixing is **necessary, not sufficient** — real-symmetric `D(k)` keeps `u(k)` real, so `A_μ=i⟨u|∂u⟩=0`
still, and only a `ℤ₂`/Zak holonomy (not a `U(1)` curvature) can survive in `k`. **So the field strength
cannot live in the displacement geometry or the Brillouin zone — it must live in the carrier phase over
physical `(x,t)`.** (Owns: `berry_reconciliation.md` §3, `carrier_construction.md` §1a.)

## 3. Earn the complex carrier from phase space — THE load-bearing new piece

(Owns: `carrier_construction.md`.) The chain, in five moves:

1. **Phase-space amplitude.** Per embedding direction, package position+velocity (90° out of phase):
   ```
   ψ^a_n = δR^a_n + i δṘ^a_n / ω .
   ```
   The `i` is **not** assumed — `δṘ` is the time-derivative, so complexifying *is* encoding the carrier's
   rotation along the timelike worldtube axis (`[[project_complex_u1_from_time]]`). Free mode runs
   `ψ ~ e^{−iωt}`. Two time slices minimal.
2. **Carrier space.** Lateral triplet (Paper III's three transverse polarizations) ⇒
   `ψ_n = (ψ^x,ψ^y,ψ^z) ∈ ℂ³`. Linear envelope is spin-1 vector, so `ℂ³` is natural, not imposed.
3. **`O(3) → U(3)`.** A *real* degenerate triplet has only `O(3)` frame freedom (preserves real energy).
   Complexifying to `ψ ∈ ℂ³` enlarges the norm-preserving frame group to `U(3)` (`ψ → Uψ`, `U†U=𝟙`).
   **The complexification is what upgrades `O(3)` to `U(3)`** — without it there is no complex gauge
   structure.
4. **Vector ≠ connection (avoid the `ℂ³`/`u(3)` slip, critique #4).** A *single* carrier vector's
   overall phase → `U(1)` (trace/EM). A transported *three-frame* `{u_1,u_2,u_3}` → overlap matrix
   `M_{ij}=⟨u_i(n)|u_j(n+ê_μ)⟩` whose unitary part ∈ `U(3)`; the **connection** is `u(3)`-valued. The
   carrier stays a vector.
5. **Dynamic phase ≠ curvature.** The clock `e^{−iωt}` alone gives plaquette holonomy `=1`, `F=0`. A
   field strength needs the *relative* carrier frame between neighbours to be path-dependent (nonzero
   loop holonomy). The construction supplies the connection; nonzero curvature is the frame genuinely
   rotating across `(x,t)` — which it does because the complex link is the time-rotation.

## 4. Link → plaquette → curvature (shared spine)

(Owns: `link_holonomy.md`.)

- **Link variable** on `n→n+ê_μ` (all μ, spacelike and timelike on the same footing):
  ```
  U_μ(n) = ⟨u(n)|u(n+ê_μ)⟩ / |·|                      (rank-1, Abelian)
  U_μ(n) = unitary part of M_{ij}=⟨u_i(n)|u_j(n+ê_μ)⟩  (rank-d, ∈ U(d))
  ```
  Reverse link `U_{-μ}(n+ê_μ)=U_μ(n)^†`. **Gauge law** (carrier rephasing `u→Vu`, `V∈U(d)`):
  `U_μ(n) → V(n)U_μ(n)V(n+ê_μ)^†` — gauge symmetry *is* the local eigenframe redundancy.
- **Plaquette** `□_μν(n) = U_μ(n)U_ν(n+ê_μ)U_μ(n+ê_ν)^†U_ν(n)^†`.
- **Continuum limit** (`U_μ=exp(−iaA_μ)`): `□_μν = exp(−ia²F_μν + O(a⁴))`,
  `F_μν = ∂_μA_ν − ∂_νA_μ − i[A_μ,A_ν]`.
- **Antisymmetry is structural:** orientation reversal `□_νμ=□_μν^†` ⇒ `F_νμ=−F_μν`. A 2D oriented loop
  can only make a rank-2 antisymmetric tensor — this is *why* both field tensors are antisymmetric.

## 5. The two sectors from the one `U(3)` plaquette

(Owns: `field_tensors.md`.) Split `A_μ = A^0_μ(𝟙/√6) + A^a_μ T^a`, `T^a=½λ^a`.

- **Faraday (trace `U(1)`).** `F_μν=∂_μA_ν−∂_νA_μ`; `E_i=F_{0i}` (time–space plaquette — needs a
  temporal link), `B_i=½ε_{ijk}F_{jk}` (spatial plaquette). **Homogeneous Maxwell = exact lattice
  Bianchi:** the six face-plaquettes of an elementary 3-cube multiply to `1` ⇒ `∂_{[λ}F_{μν]}=0` exactly
  (no `a→0`); these are `∇·B=0` and `∇×E+∂_tB=0`. (`E` requires the 2 temporal links — structural.)
- **Yang–Mills (traceless `SU(3)`).** Non-commuting matrix links keep the commutator:
  `G^a_μν=∂_μA^a_ν−∂_νA^a_μ+f^{abc}A^b_μA^c_ν`. The `f^{abc}` self-coupling **is** `−i[A_μ,A_ν]`
  projected to the adjoint — the gluon vertex = non-commutativity of the ordered link product. The only
  structural difference from EM is `d=3` vs `d=1` (links commute or not). Non-Abelian Bianchi
  `D_{[λ}G_{μν]}=0` by the same closed-cube argument with path-ordering.

## 6. The α convention — MONOTONIC, ONE DIRECTION (do not break this)

(Owns: `field_tensors.md` §4; guardrail `[[feedback_alpha_direction_consistency]]`.)
Everything `α` does points the same way; degeneracy and curvature-weight move *together*:

```
 α = 0 (linear) : carriers decouple ⇒ Abelian U(1)³ ; traceless/colour OFF
 α → 1          : carriers lock, degenerate triplet ⇒ full U(3) ; traceless SU(3) MAX
```

- Holonomy weighting (INPUT from Paper III §3 — *not re-derived here*; provenance of `(3−2α)` is an open
  cite to chase): trace `∝(3−2α)`, traceless `∝α`, `R(α)=(3−2α)/α`. Falsifier `R(0.5)/R(0.2)=0.3077`.
- **Why consistent:** a traceless `SU(3)` needs the triplet degenerate (off-diagonal frame mixing), so
  "traceless present ⟺ degenerate"; with traceless `∝α`, degeneracy must also grow with `α`. The colour
  sector is exactly what the anharmonicity `∝α` switches on.
- **Forbidden** (past bug): no "two roles in tension," no interior sweet spot `α≈0.5–0.8`, no fitting `α`
  from experiments. State the two clean limits.
- **Caveat (D2):** `α` fixes the *relation/ratio* between sectors, not physical `α_EM/α_s`
  (`[[project_alpha_undetermined_at_linear_order]]`).

## 7. Equivalence theorem (the payoff)

(Owns: `berry_reconciliation.md`.) Expand the overlap link:
`U_μ(n)=exp(−ia·i⟨u|∂_μu⟩+O(a²))=exp(−ia·a_μ+...)`, so `□_μν` is the Wilson loop of the Berry connection
and (Stokes / non-Abelian Stokes):
```
log □_μν = −ia² f_μν + O(a⁴)  (Abelian) ,   −ia² 𝓕_μν + O(a⁴)  (non-Abelian).
```
**Therefore `F_μν=f_μν`, `G_μν=𝓕_μν`** — plaquette field strength = Berry/WZ curvature; one object,
integral vs differential. Grounds Paper III's identification in the stencil. BZ-vanishing restated from
this side (real links ⇒ trivial plaquette; complex link = time-rotation).

## 8. Scope

- **In:** carrier construction; `F_μν`, `G^a_μν` as plaquette curvatures; antisymmetry; sector split;
  Bianchi identities; the `α`-activation of colour; equivalence to Berry/WZ.
- **Out (deferred, not claimed):** dynamics (D3 — `−¼F²`/`−¼G²` action, `d⋆F=J`, masslessness,
  propagation `c`); coupling value / `α` hierarchy (D2); confinement (archived); soliton-layer objects
  (vortex cores, spin-½ — Paper IV). **Not** a derivation of QCD (no running coupling, asymptotic
  freedom, string tension, hadrons).

## 9. Open items (status.md "Open derivations")

1. Link-variable definedness across carrier band-crossings (band-isolation; for colour, essentially the
   soliton layer — `g([111])=0` has no defined holonomy).
2. `O(a⁴)` is a genuine higher-derivative term (clover/improvement check; no spurious symmetric/parity-odd
   piece).
3. Non-Abelian ordering convergence (`f^{abc}` is the unique non-commuting remnant).
4. Numerical equivalence: transport a `ℂ³` carrier around a fixed `(x,t)` plaquette, confirm
   `log □ = ia²F` vs `branesim/diagnostics/berry_holonomy.py` (same loop as the `R` falsifier).
5. **Provenance of `(3−2α)`** — trace it in Paper III / `gauge_color`; currently an inherited input.

## 10. Honesty headline

Closes the **kinematic construction** of `F_μν` and `G^a_μν` from the phase-space carrier on the 8-link
stencil and proves it is Paper III's Berry/WZ object. Does **not** close **dynamics** (no action/field
equations/propagation) or **coupling** (`α` undetermined; fixes relation, not `α_EM/α_s`).

## File map

`problem_statement.md` (scope+motivation) · `carrier_construction.md` (§3) · `link_holonomy.md` (§4) ·
`field_tensors.md` (§5–6) · `berry_reconciliation.md` (§7) · `status.md` (props, open items, risks).