# Field-Strength Tensors from the Spring Stencil (Paper VII) — Status

Scaffolding stage. Statuses marked `closed` / `in-ansatz` / `inherited` / `open` per project convention.
The spine: the complex carrier is *earned* from Bloch phase space; the rest length supplies the `k`-mixing
a Berry connection needs but cannot supply the curvature (necessary, not sufficient); the curvature is the
complex `(x,t)` time-carrier; the plaquette of that carrier is the field strength and is Paper III's
Berry/WZ object. The colour curvature weight grows `∝ α` (splitting), while exact degeneracy sits at
`α=0` — they move in opposite directions (near-degenerate WZ / confinement tension). The LGT toolbox
applies as diagnostics on the same substrate graph.

## Layer assignment (keeps the paper honest)

- **Linear / carrier layer (this paper):** Bloch modes, polarization branches, phase-space carrier `ψ`,
  its `U(3)` Berry/WZ connection over `(x,t)`, `F_μν`, `G^a_μν`, the Zak `ℤ₂` of the spin-1 envelope.
- **Soliton layer (Paper IV, out of scope):** vortex cores (`π₁(U(1))=ℤ`, charge), colour textures
  (`π₃`), and the spin-½ `ℤ₂` framing holonomy (`π₁(SO(3))=ℤ₂`, gear `χ=½(θ_env−θ_lat)`). Cross-referenced
  only. `[[project_spin_half_is_soliton_layer]]`, `[[project_soliton_layer_description_language]]`.

## HAVE (imported / established — close as restatement)

- **8-link 4D stencil.** 6 spacelike axial + 2 temporal axial links per node; central-force springs;
  single anharmonicity `−k_sαa|ΔR| ∝ α`. From Papers I, III. Lorentzian sign `s_μ=(−1,+,+,+)` on the
  intrinsic stencil metric, not on the ambient `ℝ⁴`.
- **BZ curvature vanishes `∀α`.** `D(k)` real symmetric ⇒ `k`-space Berry/WZ curvature `≡ 0`, even when
  the rest length makes `u(k)` `k`-dependent. Field strength lives over physical `(x,t)`. BACKBONE #16.
- **Carrier phase = time rotation.** The complex `U(1)` `i` is carrier rotation along the timelike link
  (`[[project_complex_u1_from_time]]`); two time slices minimal. *Why* the temporal links are mandatory
  for `F_{0i}` (the `E`-field).
- **Sector split.** `u(3) = u(1) ⊕ su(3)`: trace = EM, traceless = colour (Paper III).

## TARGET (this paper's content)

- **`carrier_construction.md`** — Bloch mode → one `U(1)` phase (not `U(1)⁴`) → Berry connection from the
  eigenvector's phase freedom (not the raw plane-wave phase); rest length → bond stiffness `K_μ` (pre-
  tension weight `(1−α)`) → `k`-dependent `u(k)` (necessary); real-symmetric `D(k)` ⇒ zero BZ curvature
  (not sufficient); `ψ^a = δR^a + iδṘ^a/ω` from phase space (`i` = time-link rotation, sufficient); real
  degeneracy `O(3)` ⇒ complex carrier `U(3)`; single-vector `U(1)` phase vs transported three-frame `U(3)`
  connection; dynamic phase = raw material, not yet curvature. **Load-bearing.**
- **`link_holonomy.md`** — `U_μ(n)` from the carrier eigenframe; gauge law; plaquette `□_μν`;
  `log □_μν = −ia²F_μν + O(a⁴)`; antisymmetry from orientation reversal.
- **`field_tensors.md`** — both sectors from the one `U(3)` plaquette: Faraday (trace `U(1)`, `E`/`B`,
  homogeneous Maxwell = exact Bianchi) and Yang–Mills (`G^a_μν`, `f^{abc}` from non-commuting links,
  non-Abelian Bianchi); the single-direction `α`-activation of the colour sector.
- **`berry_reconciliation.md`** — continuum limit of the `(x,t)` plaquette holonomy `=` Berry / WZ
  curvature of Paper III; the two constructions are one object; BZ-vanishing restated.
- **`lgt_diagnostics.md`** — the LGT toolbox (Wilson loops, plaquette curvature, topological charge,
  finite-size scaling) as diagnostics on the *same* substrate graph; effective links = polar-unitary part
  of the carrier eigenframe overlap; reusable / partially / not-reusable triage; ontology vs regulator.

## OUT OF SCOPE (deferred — do not claim)

- **Dynamics (D3).** `−¼F²`/`−¼G²` action, inhomogeneous equations `d⋆F = J`, masslessness, propagation
  speed `c`. Open in Paper III; this paper builds the tensor, not its action.
- **Coupling / α (D2).** `α` undetermined at linear order; no `1/137`, no `α_s`. `α` fixes the *relation*,
  not `α_EM/α_s`. `[[project_alpha_undetermined_at_linear_order]]`.
- **Confinement.** Archived; not referenced as a result.
- **Soliton-layer objects.** Vortex/texture cores, spin-½ holonomy — Paper IV / matter sector.

## Named results (to promote to proposition in the paper)

- **Prop 0 (carrier).** `ψ = δR + iδṘ/ω` is the phase-space carrier; the real degenerate triplet carries
  only `O(3)` frame freedom, and complexification promotes it to `U(3)`. The `i` is the time-link
  rotation. The Berry connection is the eigenvector's phase freedom, not the raw Bloch phase.
- **Prop 1 (link variable & gauge law).** `U_μ(n)` is well-defined off band-crossings and transforms as a
  lattice connection under carrier rephasing `u → Vu`, `V ∈ U(d)`.
- **Prop 2 (curvature & antisymmetry).** `□_μν = exp(−ia²F_μν + O(a⁴))`, `F_μν = −F_νμ`, antisymmetry
  inherited from plaquette orientation — the structural origin of the rank-2 antisymmetric tensor.
- **Prop 3 (Faraday + Bianchi).** Trace sector ⇒ Faraday tensor; homogeneous Maxwell = exact lattice
  Bianchi identity.
- **Prop 4 (Yang–Mills field strength).** Traceless sector ⇒ `G^a_μν` with `f^{abc}` self-coupling from
  ordered non-commuting links; non-Abelian Bianchi.
- **Prop 5 (`α` and the colour sector — near-degenerate).** From the triplet spectrum
  `λ_A = κ[(1−α)|k|² + α k_A²]`: trace `∝ (3−2α)`, traceless splitting `∝ α g(k̂)`. The colour curvature
  weight grows `∝ α` off `[111]` and vanishes at `α=0`. **Exact degeneracy (and `U(3)` frame freedom) is
  the `α=0` limit; the splitting that turns colour on is what breaks it — degeneracy and colour move in
  *opposite* directions** (`[111]` coherence-vs-colour, kinematic confinement). The holonomy ratio
  `R(α) ∝ (3−2α)/α` is normalization-independent, giving `R(0.5)/R(0.2) = 0.3077` — derived here from our
  own `K_μ`, reproducing the colour bridge.
- **Theorem (equivalence).** Continuum limit of the `(x,t)` plaquette holonomy equals the Berry / WZ
  curvature — the plaquette construction *is* Paper III's identification, grounded in the stencil.
- **Prop 6 (diagnostics, candidate).** The field strength is measurable with the LGT toolbox on the same
  substrate graph; effective links are the polar-unitary part of the carrier eigenframe overlap — no
  second lattice.

## Open derivations

*Per-paper bridge entries (`[[project_open_problems_tracker]]`). This bridge = `field-strength`.*

1. **Link-variable definedness across band-crossings.** `U_μ` needs band isolation; characterize where the
   carrier band is gapped enough that `U_μ` is single-valued (ties to Paper III adiabaticity).
2. **`O(a⁴)` lattice-artifact terms.** Confirm the subleading plaquette expansion is a genuine
   higher-derivative correction (no spurious symmetric or parity-odd piece) — clover/improvement check.
3. **Non-Abelian ordering convergence.** Verify the `f^{abc}` term is the unique non-commuting remnant and
   that the symmetric/path-ordering choice does not shift the continuum `G^a_μν`.
4. **Numerical equivalence check.** Transport a band-isolated `ℂ³` carrier around a fixed `(x,t)` plaquette
   and confirm `log □ = −ia²F` against `branesim/diagnostics/berry_holonomy.py` (the same loop the D2
   falsifier uses). Also the entry point for the `lgt_diagnostics.md` §6 pipeline.
5. **Coupling reading of `R(α)`.** The spectral ratio `R(α) ∝ (3−2α)/α` is derived here; identifying it
   with a physical coupling ratio `α_EM/α_s` (D2) needs the connection normalization `C`, undetermined at
   linear order. The `0.3077` *ratio* is closed; the *coupling* reading is deferred.

## Self-flagged risks (carry into critique when reviewed)

1. **Band-isolation is load-bearing.** `U_μ` is single-valued only off band-crossings; the colour `d=3`
   subspace is degenerate by construction and `g([111]) = 0` means the colour-coherent direction has no
   defined holonomy. Valid *where* the subspace is band-isolated — for colour, essentially the soliton
   layer (Paper IV), not a free linear direction.
2. **"Field-strength tensor" ≠ "the gauge theory."** Constructing `G^a_μν` is not deriving QCD. Keep the
   action, running coupling, and confinement (archived) explicitly out.
3. **The complex structure is earned.** The `i` is the timelike-link rotation, realized as the phase-space
   pairing `ψ = δR + iδṘ/ω`. Always trace it to the temporal spring; cite `[[project_complex_u1_from_time]]`.
4. **No `ℂ³`-as-`u(3)`-matrix slip.** The split is of the *connection* valued in `u(3)`, not of the carrier
   vector. Single vector ⇒ `U(1)` phase; transported three-frame ⇒ `U(3)` connection.
5. **`k`-mixing is necessary, not sufficient.** Do not let the rest-length/stiffness story imply a nonzero
   BZ curvature. Real-symmetric `D(k)` forces `A_μ = 0` in `k`; the curvature is the complex `(x,t)`
   carrier. Keep the three `ℤ₂`-adjacent objects apart (`U(1)` core winding; Zak `ℤ₂` of the spin-1
   envelope; spin-½ `ℤ₂` at the soliton layer).
6. **`O(a⁴)` honesty.** Lattice plaquettes carry higher-derivative artifacts; the claim is the *leading*
   term is `a²F_μν`, pending the clover/improvement check (Open derivation 2).
7. **Get the `α` direction right: degeneracy and colour move *opposite* ways.** Exact degeneracy (and
   `U(3)` frame freedom) is the `α=0` limit; the traceless splitting `∝ α g(k̂)` that turns colour on is
   what *breaks* the degeneracy. The colour curvature weight is monotonic `∝ α` (off `[111]`), vanishing at
   `α=0` — but do **not** say degeneracy grows with `α`, and do **not** deny the `[111]`
   coherence-vs-colour tension: that tension is the kinematic-confinement physics. No interior sweet-spot
   `α`, no fitting `α` from experiments.
8. **`(3−2α):α` is derived here, with a spectral-vs-coupling caveat.** The ratio and `R(0.5)/R(0.2)=0.3077`
   fall out of our own `K_μ` — present them as this paper's result. But `R(α)` is a ratio of eigenvalue
   *splittings* (dispersion); its identification as a physical coupling ratio `α_EM/α_s` needs the
   normalization `C` (D2, deferred). State that caveat whenever the number appears.
9. **LGT diagnostics are diagnostics, not dynamics.** The Wilson-action-fit is a *question* about an
   emergent effective action, never a replacement for the elastic action. Smoothing/gradient-flow applies
   to the diagnostic eigenbundle, never the physical substrate state.

## Honesty headline

This paper closes the **kinematic construction** of `F_μν` and `G^a_μν` from the phase-space carrier on
the 8-link spring stencil, proves it is the Berry/WZ object of Paper III, and shows the LGT diagnostic
toolbox applies on the same substrate graph. It does **not** close the **dynamics** (no action, no field
equations, no propagation) or the **coupling** (`α` undetermined — it fixes the *relation/ratio* between
sectors, not `α_EM/α_s`). Say so plainly in the scope section.
