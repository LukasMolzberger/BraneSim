# ARCHITECTURE.md — Layered structure of the theory (paper backbone)

This file describes the **conceptual stack** of the theory and how the paper is
organized around it. `BACKBONE.md` holds the non-negotiable *principles*; this
file holds the *derivational architecture* — the ordered chain of layers, each
built strictly on the one below, that the manuscript must present.

The organizing idea: a single ontic substrate at the bottom, and a sequence of
**bridges** that each recover one piece of known physics as an *emergent*
description living on top of that substrate. Nothing above Layer 0 is a new
fundamental postulate; each higher layer is a derived, effective object.

```
 L5  Particles / spin-1/2 solitons        (self-confined standing waves, Z2 holonomy)
 L4  Field strengths: Faraday + QCD        (Berry/WZ curvature = F_munu, G^a_munu)
 L3  Berry / Wilczek-Zee connections       (holonomy of the eigenbundle)
 L2  Bloch / phase layer                   (4D plane-wave analysis, band eigenbundles)
 L1  Linear branch structure               (dispersion, c_L vs c_T, α_s-controlled mixing)
 L0  Brane-lattice substrate               (8-link stencil, Lorentzian-sign action)
 L-  Gravity channel (transverse)          (induced-metric ∂u·∂u eigenstrain; runs parallel)
```

**Gauge-sector derivation chain** (the spine of L2→L4):

```
static 4D hypercubic brane
   ⇓
Lorentzian-sign stationary action   (η_μ = prestress sign pattern; α_s, α_t, γ_t magnitudes)
   ⇓
anisotropic stiffness tensor        (C_{nμ}[R̄]: Pythagorean Hessian about nonlinear R̄)
   ⇓
4D Bloch spectral problem           (D_{R̄}(k): C_{nμ} as matrix-valued hopping coeff.)
   ⇓
complex eigenbundle over the 4D Brillouin zone   (rank-3 transverse fluctuation bundle;
                                                   4th/time-amplitude mode split off by α_t)
   ⇓
Berry / Wilczek–Zee holonomy        (𝒜_μ^a_b ∈ u(3): μ=1..4 base, a,b=1..3 carrier)
   ⇓
U(1) × SU(3) gauge structure        (trace = EM Faraday, traceless = color)
```

The complex phase is intrinsic to the 4D spectral problem (the character of the
lattice translation group) — it is **not** reconstructed from motion through a
preferred time direction. No analytic-signal / worldtube construction enters.

---

## Layer 0 — The substrate (ontology, fully deterministic)

The fundamental and *only* fundamental object.

- A **static** 4D hypercubic (quartic) brane lattice embedded in a 4D
  **Euclidean** ambient space (codimension 0). Node positions `R_n^A ∈ ℝ^4` are
  real-valued; the ambient dot product is `δ_AB = diag(1,1,1,1)`. The four ambient
  directions are geometrically equivalent, and the four internal lattice
  directions are treated symmetrically at the stencil and Bloch level — the
  time-like direction is one of the four internal directions, **not** an external
  evolution parameter.
- Local **8-link stencil**: 6 spacelike axial neighbours (`±ê_x,±ê_y,±ê_z`) plus
  2 temporal neighbours (`±ê_t`).
- Single compact action (this is the whole theory at the bottom):

  ```
  S_4[R] = ½ Σ_n Σ_{μ=1..4} η_μ κ_μ ( L_{nμ} − r_μ )²,
  L_{nμ} = sqrt( Σ_A (R_{n+ê_μ}^A − R_n^A)² ),
  η_μ = (−1,−1,−1,+1),
  r_1=r_2=r_3 = α_s a,   r_4 = α_t a  (0 < α_t < 1, likely α_t = 1−ε_t),
  κ_1=κ_2=κ_3 = κ_s,     κ_4 = κ_t.
  ```

- **Prestress: sign vs magnitude are separate.** `η_μ` carries only the Lorentzian
  *sign*; the *magnitudes* need not match in space and time. The signed prestress is
  `ρ_μ ~ η_μ κ_μ a (1 − α_μ)`. The **minimal dimensionless control set** is
  `{α_s, α_t, γ_t}`, `γ_t = κ_t/κ_s`, with held spacing `a_s = a_t = a` (preserves
  the quartic substrate). Their jobs: `α_s` — spatial branch mixing, soliton shape,
  color-sector activation; `α_t` — temporal/time-amplitude coupling, gravity channel,
  and separation of the fourth carrier mode (Layer 3); `γ_t` — effective light-cone
  calibration. This **replaces the earlier single-`α`, `r_t ≈ 0` formulation**: the
  old kinetic limit made the time-like sign look like a kinetic-sign artifact rather
  than genuine prestress (the T12 tension — now **resolved**, `derivations/T12_*`: the
  `(1−α_t)` anisotropy of the temporal link's continuum kinetic term is the fingerprint
  of genuine prestress, and `η_4=+1` makes `S=T−V` → Verlet wave dynamics). `α = 1` is the *zero-mismatch* limit (no
  prestress, transverse stiffness `∝(1−α) → 0`), and `α_t = 1` is forbidden because
  it over-decouples the time/amplitude direction from spatial geometry.

- **The Lorentzian signature is the sign pattern of the directional prestress —
  not an independent postulate, not an ambient metric, not a special treatment of
  time.** The prestress state fixes a *directional* strength whose sign is `η_μ`,
  while its magnitude is set separately by `α_μ, κ_μ`:

  ```
  prestress state  →  ρ_μ = η_μ κ_μ a (1 − α_μ)  →  Lorentzian-sign brane action.
  ```

  The action sign is exactly `η_μ = (−1,−1,−1,+1)` (three spacelike directions one
  way, the time-like direction opposite). The ambient stays Euclidean; the relative
  sign of temporal vs spacelike links makes the stationary equation hyperbolic
  (wave-like) rather than elliptic. This removes one more independent-looking
  assumption: signature is *derived* from prestress, not imposed alongside it.
- **The only nonlinearity is the Pythagorean link length** (the square root).
  Even a Hookean scalar spring law is nonlinear in the node coordinates.
- **Spatial prestress `α_s ∈ (0,1]` is the master spatial control parameter.** It
  sets the mismatch between held spacing `a` and stress-free length `α_s a`, and
  thereby the transverse/longitudinal stiffness ratio — the knob that runs from
  *complete axis separation* (`α_s→1`) to *complete axis mixture* (`α_s→0`) in all
  higher layers. `α_t` and `γ_t` play the analogous role for the time-like direction
  and set the gap that isolates the rank-3 color carrier (Layer 3).
- Stationarity `∂S_4/∂R_n^A = 0` yields, depending on boundary conditions:
  - initial-value view → Störmer–Verlet time march;
  - 4D block view → a **saddle-point** root-find `∇S = 0` (S = T − V is indefinite,
    so *never minimize*).

## Layer 1 — Linear branch structure

Linearize `R = R_★ + u` about the tensioned vacuum.

- Per-axis tangent stiffness `M_i^{AB} = (1−α_s)δ^{AB} + α_s δ_i^A δ_i^B`
  (temporal: `M_4^{AB} = (1−α_t)δ^{AB} + α_t δ_4^A δ_4^B`).
- Branch speeds: longitudinal `c_L² = κ_s a²/m`; transverse
  `c_T² = (1−α_s) κ_s a²/m`.
- **α_s interpolation, made concrete here:**
  `α_s→0` ⇒ `c_T = c_L`, isotropic, axes fully **mixed/degenerate**;
  `α_s→1` ⇒ `c_T→0`, axes fully **separated**. This is the microscopic meaning of
  the "separation ↔ mixture" transition that the whole gauge story rides on.
- The tensioned vacuum is stable **only with periodic BC** (see LESSONS_LEARNED):
  open/free BC lets the prestressed lattice relax globally.
- **`M_i^{AB}` is the bare-vacuum tangent operator only — a reference case, not the
  arbiter of the gauge sector.** It is the Hessian about a straight, axis-aligned
  background and linearizes away the Pythagorean nonlinearity. The load-bearing
  object for L2–L4 is the **anisotropic stiffness tensor** — the fluctuation Hessian
  about a *finite-amplitude nonlinear* background `R̄`:

  ```
  C_{nμ}^{AB} = η_μ κ_μ [ (1 − r_μ/L̄_{nμ}) δ^{AB} + (r_μ/L̄_{nμ}) Q̂_{nμ}^A Q̂_{nμ}^B ],
  Q̂_{nμ} = (R̄_{n+ê_μ} − R̄_n) / L̄_{nμ}.
  ```

  (isotropic `δ` part + rank-one `Q̂Q̂ᵀ` part: a spring is stiffer along its own link
  than transverse to it.) In the trivial vacuum `Q̂_{nμ}^A = δ_μ^A` and this reduces
  to the diagonal `M`; with finite transverse deformation `Q̂Q̂ᵀ` gains off-diagonal
  components, so the tensor becomes a **background-dependent 4D frame tensor**.
  Because different links carry different `Q̂`, generally `[C_{nμ}, C_{mν}] ≠ 0` —
  this noncommutativity is the *mechanical seed* of a non-Abelian transported
  subspace. This `C_{nμ}^{AB}` is the concrete object that carries the Pythagorean
  nonlinearity forward: in Layer 2 it becomes the **matrix-valued Bloch hopping
  coefficient** of `D_{R̄}(k)` (a link crossing supercell displacement `d`
  contributes `C_ℓ^{AB} e^{ik·d}`). The leading extra coupling is cubic —
  `V^{(3)} ~ (α_s κ/2a)(Δ_μ u^μ) Σ_{A≠μ}(Δ_μ u^A)²`, longitudinal strain ↔
  transverse-amplitude² (BACKBONE #4) — and must **not** be discarded: a negative
  gauge-sector result that relies on dropping it is an approximation artifact, not a
  theorem (LESSONS_LEARNED #3).

## Layer 2 — Four-dimensional Bloch / phase layer

Apply Bloch analysis directly to the full 4D hypercubic brane lattice — in **all
four** lattice directions, treated symmetrically.

- Lattice index `n = (n¹,n²,n³,n⁴) ∈ ℤ⁴`; node positions stay real,
  `R_n^A ∈ ℝ⁴`.
- Small perturbations about the stationary substrate are expanded in 4D Bloch
  modes:

  ```
  u_n^A = ε^A(k) · exp( i Σ_{μ=1..4} k_μ n^μ ).
  ```

- The fourth lattice direction is **not** an external evolution parameter — it is
  one of the four internal brane directions. Its only distinguished role is the
  opposite sign `σ_4 = +1` in the brane action. The Lorentzian character therefore
  enters through the action, not through the ambient metric and not through a
  time-evolution postulate.
- The complex phase `e^{i k·n}` is **not** a new fundamental complex field: it is
  the character of the 4D lattice translation group, used to diagonalize
  translations. Physical configurations stay **real**, recovered from conjugate
  Bloch-mode pairs `u_n^A = ε^A(k) e^{ik·n} + ε^{A*}(k) e^{−ik·n}`. No
  analytic-signal / worldtube construction is needed.
- **Bloch analysis is applied to fluctuations about a stationary background `R̄`**
  (`∂S/∂R[R̄] = 0`), not only the trivial vacuum. The resulting operator is
  `D_{R̄}(k)`; for the gauge sector `R̄` is the **universal nonlinear periodic vacuum
  carrier** — a finite-amplitude periodic supercell state (the vacuum microtexture),
  *not* the bare straight lattice and *not* a soliton core, because **color is a
  universal vacuum gauge sector** (matter enters only at Layer 5). So `D_{R̄}(k)` is
  assembled from the anisotropic stiffness tensors `C_{nμ}[R̄]` (L1) as
  **matrix-valued hopping terms** `C_ℓ^{AB} e^{ik·d}`. In the trivial vacuum these hopping matrices commute and share
  fixed eigenvectors (trivial bundle); in a nonlinear background they do not commute,
  so the eigenvectors of `D_{R̄}(k)` acquire nontrivial `k`-dependence — exactly what
  Berry/WZ geometry needs. The trivial-vacuum operator `D_vac(k;α_s)` is only a
  reference limit.
- Output: polarization eigenspaces over the 4D Brillouin zone — the carrier
  **eigenbundles** `k ↦ ε_n(k)` on which Berry / Wilczek–Zee holonomy (L3) is
  defined.

## Layer 3 — Berry / Wilczek–Zee connections

Holonomy of the Layer-2 eigenbundles.

- Non-degenerate band → Abelian Berry connection `a_μ = i⟨ε|∂_{k_μ}ε⟩`.
- Degenerate triplet → non-Abelian **Wilczek–Zee** connection
  `[𝒜_μ]_{ab} = i⟨ε_a|∂_{k_μ}ε_b⟩`, `U(3)`-valued.
- **Index structure (two independent spaces).** `𝒜_μ^a_b` has an *external* Bloch/
  brane index `μ = 1..4` (the base) and an *internal* carrier index `a,b = 1..3`
  (the fiber). Four-dimensionality of the base does **not** enlarge the gauge group:
  `f_μν` and `𝒢_μν` are 4D antisymmetric tensors whose *values* are 3×3 internal
  matrices. Geometrically `𝒜 ∈ Ω¹(M₄, u(3))`.
- Algebraic split `U(3) ≃ (U(1)×SU(3))/ℤ_3`:
  `𝒜_μ = a_μ I_3 + B_μ`, `Tr B_μ = 0`. **Interpretation A (selected):** the trace
  `a_μ = ⅓ Tr 𝒜_μ` is the electromagnetic U(1) (common phase of the whole triplet),
  the traceless `B_μ` is the color su(3) (relative internal twist). Candidate
  identification until T3/T4 are derived.
- **The natural rank-3 carrier is the complexified transverse fluctuation bundle**
  of the nonlinear background: at each link the direction `Q̂ ∈ ℝ⁴` leaves a 3D
  orthogonal complement `T⊥ = {v : v·Q̂ = 0}` (dim 3 because the ambient is 4D),
  which Bloch-complexifies to `ℂ³`. This *moving* triplet — not the fixed spatial
  axes `span(e_x,e_y,e_z)` — is the color-carrier candidate; its projector
  `P_3(k) = Σ_a |ε_a(k)⟩⟨ε_a(k)|` can vary over `k`, so the WZ curvature
  `𝓕_{μν} = i P_3[∂_{k_μ}P_3, ∂_{k_ν}P_3]P_3` need not vanish.
- **Why U(3), not SU(4) — the fourth mode is spectrally separated, not deleted.**
  The substrate and Bloch base stay fully 4D; only the *transported carrier* is
  rank 3. `α_t` (with `γ_t`) lifts the fourth (time/amplitude) polarization mode out
  of the near-degenerate spatial triplet. The carrier condition is
  `λ_1 ≈ λ_2 ≈ λ_3` with `|λ_4 − λ_triplet| ≫ |λ_a − λ_b|` (large outside gap,
  small internal splitting). The exclusion is *primarily geometric* (the fourth
  direction is the local longitudinal/time-amplitude mode) and *secondarily
  spectral* (`α_t` makes the gap dynamically stable). Rank 3 → U(3); rank 4 would be
  U(4)/SU(4). Near-degeneracy with a large outside gap (`ε ≪ Δ_outside`) — not exact
  degeneracy — is the requirement (settled); the composite rank-3 subspace carries
  the U(3).
- **Settled — the real su(3) test (T1), now PROVED:** a real transverse frame alone
  yields only `so(3) ⊂ su(3)` (3 generators, confirmed by a control calculation).
  Full `su(3)` (8 generators) is obtained once the complex Bloch phases act on the
  noncommuting real anisotropic stiffness frames: they supply the 5
  symmetric-traceless + Cartan directions that complete `so(3)→su(3)`. The genuine
  criterion `Lie⟨𝒢_μν, [𝒢_μν, 𝒢_ρσ], …⟩ = su(3)` (rank 8, gauge-invariant, not a
  pure-gauge `u(3)` artifact) is met by an exact stationary periodic carrier — see
  `derivations/T1_derivation_and_proof.md`, `derivations/t1_su3_witness.py`. The
  bare-vacuum "no SU(3)" is a trivial-bundle special case (LESSONS_LEARNED #3), not
  a no-go; the color lives in `D_{R̄}(k)`, and `so(3)`-only is a measure-zero
  real-frame case.
- **`{α_s, α_t, γ_t}` control the sector split.** Tune the stiffness/gap hierarchy so
  the trace `U(1)` mode is soft/long-range and the traceless `SU(3)` modes are
  stiff/gapped/short-range, with weak residual coupling ε≪1.

## Layer 4 — Field strengths (the two demonstrations we owe)

- **Faraday tensor (U(1)):** trace curvature `f_{μν} = ∂_μ a_ν − ∂_ν a_μ`. **T3 now
  derived** (`derivations/T3_*`): Bianchi/homogeneous Maxwell automatic; U(1) gauge
  invariance + locality + power counting force Maxwell `−¼e⁻²f²` (massless photon;
  `1/e²` = positive/finite substrate abelian quantum metric); inhomogeneous
  `∂^μf_μν=e²J_ν` with `∂·J=0`. The pristine vacuum is **EM-flat** (abelian Berry
  curvature ≈0 by a PT-like symmetry — no background EM field), with EM switching on
  under deformation/matter; consistent with `su(3)≠0` (non-abelian off-diagonal
  sector). Charge = quantized **U(1) vortex winding** (`π₁(U(1))=ℤ`); the photon is
  continuous / non-quantized (see below).
- **The photon is a continuous, non-quantized mode** (BACKBONE #13): the
  propagating U(1) Berry excitation carries continuous amplitude and has no
  intrinsic quantum. Its apparent `hν` granularity is *inherited* — it only shows
  up at emission/absorption, gated by the discrete confinement-equilibrium
  transitions of the L5 spin-1/2 solitons it couples to. Quantization enters the
  radiation sector solely through this coupling, never as a property of `f_{μν}`
  itself.
- **QCD field strength (SU(3)):** traceless curvature
  `𝒢_{μν} = ∂_μ B_ν − ∂_ν B_μ − i[B_μ,B_ν]` (4D antisymmetric in `μ,ν`, valued in
  traceless 3×3 internal matrices). **T4 now derived** (`derivations/T4_*`): the
  emergent local `SU(3)` invariance (T1) + locality + power counting force the
  leading effective action to be uniquely `−¼g⁻² G^a_{μν} G^{aμν}` — the gluon
  mass term `Tr(BB)` is forbidden (gluon massless; pure-gauge `𝒢=0` verified), and
  `1/g²` is the substrate's positive, finite, strain-dependent integrated
  non-Abelian quantum metric (derives the coupling-vs-strain of Decision J). The
  `−i[B_μ,B_ν]` term is the algebraic source of the self-interaction, whose
  `su(3)` structure constants (`f^{acd}f^{bcd}=2N δ`, `N=3`) fix the cubic/quartic
  gluon vertices with the single coupling `g`. Universal field over empty space
  (Decision M); colored matter enters at Layer 5.
- **Color is a universal vacuum gauge sector; nonlinear strain sets its *coupling
  strength*, not its existence.** The U(3) carrier — hence the SU(3) structure —
  exists over the whole substrate, because the gauge vacuum is the universal
  nonlinear periodic carrier whose projector `P_3(k)` moves in `k` (universal
  baseline `𝒢_{μν} ≠ 0`). Matter (solitons) enters only at Layer 5, carrying color
  charge into an already-universal field. What *scales with local nonlinear strain*
  is the magnitude of the traceless curvature (Large_Energies discussion): in weak,
  long-wavelength regions `∂_k P_3` is small and the SU(3) sector is soft /
  perturbative; where the link frames `Q̂_{nμ}` vary strongly — high energy density,
  short length scales (`|Δu|/a ~ 1`), topological confinement — it becomes large and
  nonperturbative. So strong-coupling / confinement-like SU(3) *dynamics* is
  short-range (near matter) while color itself is universal, analogous to QCD
  asymptotic freedom. Any coupling-vs-strain scale must come from `α_s, α_t, γ_t, a`,
  never imposed by hand (BACKBONE #5).

## Layer 5 — Particles: spin-1/2 solitons

- Particles are extended, self-confined standing-wave / solitonic patterns with
  **spherical-harmonic internal structure and constraints** — never point-like;
  point behaviour is an interaction/localization effect.
- Rest mass = energy of a self-confined closed-loop / toroidal mode.
- Charge = U(1) vortex winding, `π_1(U(1)) = ℤ`.
- **Spin-1/2 = a separate ℤ_2 condition** (Spin_One_Half discussion): a 2π
  winding around the core is necessary but *not sufficient*; spin-1/2 requires the
  soliton texture to trace the nontrivial element of `π_1(SO(3)) = ℤ_2` under a 2π
  rotation *relative to the brane far-field*, trivial only at 4π (belt-trick).
- **Solitons are the sole source of intrinsic quantization** (BACKBONE #13):
  discreteness (mass, charge, spin, energy levels) is the confinement-equilibrium
  spectrum of the self-trapped mode, a property of the bound equilibrium — not of
  the substrate field. Every quantum appearing in the radiation sector (L4) is
  inherited from these discrete soliton transitions. This layer is therefore where
  all `hν`-like granularity originates.

## Layer − (parallel) — Gravity channel

- The induced metric `g_{ij} = ∂_i X · ∂_j X` contains the transverse term
  `∂_i u ∂_j u` — the same Pythagorean nonlinearity as Layer 0.
- This term acts as an **eigenstrain**: localized transverse excitation energy
  sources a slow in-brane contraction/dilation → candidate gravity-like channel.
- Runs alongside, not on top of, the gauge stack (it uses the *geometric* channel,
  they use the *phase/holonomy* channel).

---

## Dual-observer split (load-bearing, spans all layers)

- **Outside/lab observer** knows the lattice, sees hypercubic axis anisotropy and
  direction-dependent wave speeds. This anisotropy is *kept*, not tuned away — it
  is the structural seed of the SU(3) sector.
- **Inside/soliton observer** is built from substrate excitations and is
  *conjectured* blind to the anisotropy (rods, clocks, cones renormalize with the
  local wave speed) → measures an isotropic `c`. This is the load-bearing
  assumption for emergent Lorentz kinematics.
- **T2 status** (`derivations/T2_*`): the effective metric's **Lorentzian signature
  (3,1) is derived from the `η` sign pattern** (not imposed), the branch speeds
  `c_L,c_T,c_4(α_s,α_t,γ_t)` are calibrated in closed form, the long-wavelength cone
  is massless and Lorentz-invariant per sector (slower modes causal/subluminal), and
  the `e_4` amplitude mode is exactly isotropic. But the lab is genuinely
  anisotropic/birefringent for `α_s>0` (spread→0 only as `α_s→0`, which disables the
  gauge sector), so the *single universal* inside-observer cone is not achievable by
  tuning — it stays the conjecture above, now with a measured obstruction and
  concrete resolution routes (nonlinear-carrier gauge-mode dispersion; Layer-5
  soliton-scale isotropy; single-sector rod/clock renormalization).