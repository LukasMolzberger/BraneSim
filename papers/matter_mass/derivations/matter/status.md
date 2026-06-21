# Matter Bridge — Status

## HAVE

- Self-adjoint eigenproblem formulation for localized modes.
- VSH angular structure: hedgehog ansatz (J=0, L=1) as canonical baryon seed; Skyrme-twisted extension with F(r): [0,∞] → [π,0].
- Soliton labels: (J,P; SU(3) irrep; Q_U(1); B_winding).
- Kinematic confinement of color (coherence and color-activity mutually exclusive on cubic lattice).
- Derrick theorem analysis: anti-collapse (UV) blocked by lattice UV cutoff (spacing a; Derrick's λ→0 rescaling halts at one cell; single-node concentration immediately disperses). Anti-expansion (IR) *is* the confinement problem (color bridge #24/#25), not a separate problem; the geometric quartic ∝ k_s α/a (Pythagorean norm term) gives a Derrick λ⁻¹ balance R_h/a = κ (A/a) √(α/(1−α)) and is *part of* — but not shown sufficient for — the stabilization (C2 hardening sweep still disperses); SU(3) mechanics (kinematic confinement + π₃ texture) likely also load-bearing. Quartic and SU(3) confinement to be assessed together.
- Soliton sweet spot: α ≈ 0.5–0.8, operational seed target α=0.7, A/a ≈ 10 → R_h/a ≈ 8.
- Periodic BC requirement established: a prestressed lattice (rest length αa < a) is under tension; only a periodic box (fixed at N·a) holds that tension, giving a stable tensioned vacuum. Open/free boundaries let the network relax toward αa and collapse even the A=0 vacuum — so the tensioned vacuum, and any soliton on it, must be posed with periodic spatial BC.
- A bound particle is an eigenstate of the tensioned vacuum, not a marched/relaxed seed: stability and winding must be read from a self-consistent stationary state, not a seed-and-watch trajectory.

## MISSING

- Converged numerical solution for any localized mode with periodic-clamped tensioned vacuum.
- Rigorous Derrick's theorem write-up (cite theorem directly; state lattice UV cutoff argument cleanly; treat the IR/anti-expansion direction as the confinement problem — geometric quartic *and* SU(3) confinement together, not the quartic alone — and tie both to the effective field manifold).
- Proof of existence (not just Derrick balance — need topological stability argument independent of Derrick).
- Baryon-triplet simulation with correct BC and eigenstate approach (not seed-and-watch).

## Open derivations

*Constitutive law & soliton stability (group C), plus D1, the EM-charge ⇄ spin-½ bridge.*

### C1. StVK compression non-coercivity vs soliton stability — `open`
**Statement.** The St. Venant–Kirchhoff baseline is neither polyconvex nor
coercive in compression: its stored energy does not diverge as `det F → 0⁺`, so
it does not penalize local volume collapse / element inversion and can soften
under strong compression. This is **structural, not a continuum artifact**: the
microscopic central-force pair energy `½ k_δ(|R_{p+δ}−R_p| − αa)²` stays bounded
(→ `½ k_δ(αa)²`) as a link collapses `|R_{p+δ}−R_p| → 0`. The lattice spacing
`a` masks the instability (cannot compress below the grid); the continuum
exposes it.
**Why it matters.** The decisive numerical target is a stable baryon-like
soliton (backbone #12). A strongly localized mode could collapse through this
compression instability rather than being protected. The Skyrme/geometric
quartic (#17) guards only against **Derrick** (uniform-rescaling) collapse — a
distinct failure mode — so it does not by itself remove this.
**Candidate approach.** For soliton-scale runs either (i) monitor `det F` /
minimum link length and verify configurations stay non-inverting within the
StVK-valid regime, or (ii) substitute a polyconvex hyperelastic completion
(neo-Hookean / Mooney–Rivlin / Ogden) within the *same* induced-metric
framework, which leaves the geometric transverse–lateral coupling unchanged.
Stated in paper §3 (constitutive subsection). Owner: soliton-hunter.

### C1a. StVK quartic coefficient has the wrong α-scaling — `open`
**Statement.** The lattice spring `½k_s(|ΔR|−αa)²` is quadratic in *stretch*; StVK
is quadratic in Green–Lagrange *strain* `E`. They agree at quadratic order
(identical wave speeds) but disagree at quartic: the lattice-exact transverse
quartic is `+k_s α/(8a²) ∝ α` (continuum `W_4 ∝ k_s α/a`, vanishing at α→0; from the
norm term `−k_s αa|ΔR|`, see
`paper/derivations/geometric_nonlinearity_alpha_scaling.md` §4), whereas a naive
StVK identification locks the quartic to `μ ∝ (1−α)` — inverted α-scaling.
**Why it matters.** Any soliton-width / hedgehog-radius formula must use the
geometric quartic prefactor, not the StVK one. The hedgehog radius is
`R_h/a = κ(A/a)√(α/(1−α))` (`vsh_channel_decomposition.md` §2.5): radius grows
with α as `√(α/(1−α))` and is *linear* in amplitude `A`. The Derrick `λ`-scaling
(backbone #17) is α-independent and unaffected; the α≈0.5–0.8 sweet-spot story
carries no sign error.
**Status.** open. Owner: soliton-hunter. Feeds the sprint4 target
(α≈0.7, A/a≈10). See `[[project_geometric_nonlinearity_alpha_scaling]]`.

### C2. Rest baryon — `open`
**(BACKBONE #25, D4.)** The baryon is specifically the `SU(3)`
`π₃` **texture** (a localized Skyrmion/instanton, codim-3) — topologically a different
object from the `U(1)` carrier-phase **vortex** (the EM/electron-like object, codim-2,
`π₁`). The **active first experiment is the `U(1)` vortex** (`EXPERIMENT.md`/D4),
which is `π₁`-protected and the cleaner test; the `SU(3)` texture (this entry) follows.
**Statement.** The time-periodic, amplitude-breathing baryon ansatz (the
common-carrier Skyrme breather) is dynamically unstable, and prestress α does not
stabilize it: the breather *mechanism* is sound — Pythagorean hardening
`Q = 8k_sα/a² > 0` places a localized mode above the transverse band
([[project_pythagorean_breather_go]]) — but the breathing/scale degree of freedom
itself carries an instability (Floquet multiplier `>1`, α-independent). The
conclusion is structural: a pulsating-scale ansatz is the wrong vehicle, so the
rest baryon must be sought as a *static* topological soliton (see Replacement).
**Why it matters.** A stable baryon-like soliton is the decisive numerical target
(backbone #12). The breather ansatz ([[project_baryon_search_vehicle_breather]]) is
the wrong vehicle, so the program needs a static-soliton replacement.
**Replacement (proposed, not yet owner-locked).** The rest baryon is a **static**
topological soliton — an extremum of the *spatial* potential `V` at fixed winding
`B=1` — dressed by a rigid internal **iso-rotation** collective coordinate
(spin-from-isospin, D1), **not** a pulsating lump. Rationale: (i) it removes the
breathing DOF that carries the Floquet instability; (ii) it is the bound object
that the kinematic color-confinement conclusion (backbone #24) *requires* anyway —
color has no free linear state and lives only in a phase-locked soliton; (iii) it
is the same hedgehog that carries the spin-½ `ℤ₂` Finkelstein–Rubinstein holonomy
([[project_spin_half_is_soliton_layer]], D1); (iv) it is cheaper than the breather.
**Variational principle (consistency note).** Unlike the Lorentzian *action* `S`
(a saddle, unbounded below — must root-find `∇S=0`, principles §saddle), the
*static spatial energy* `V` on a single slice is bounded below, so the static
soliton is a genuine **constrained minimization of `V` at fixed `B=1`** (or
equivalently `∇V=0` + Hessian check). Derrick's theorem is the obstruction; the
geometric quartic (#17, coeff `∝ k_sα/a`) supplies the `λ^{+1}` balance.
**Regime (load-bearing).** This must live in the **lattice-scale `w~a`,
large-amplitude (O(1) strain)** corner — **not** the wide `w≫a` regime, which C1
closed negatively (no interior energy minimum; no hyperelastic rescue). A `w~a`
static soliton is Peierls–Nabarro grid-pinned, re-opening the PN / emergent-Lorentz
tension (#8, #17): the test **must report `w/a` and quantify the PN barrier**.
**Test spec.**
- (a) *Static extremum.* Constrained-minimize `V` (or root-find `∇V=0`) on a `B=1`
  hedgehog + `X⁴`-twist seed at α≈0.5–0.8, large `A`, **periodic** boundary (box at
  N·a, to clamp the prestress — see HAVE), box ≫ seed; verify `B(slice)=1` from the
  actual field and that the config stays localized (`spread_ratio ≪ 1`, not box-fill —
  use the metrics from [[project_c2_skyrme_no_confinement]], never
  `leakage_fraction`).
- (b) *Stability.* Spectrum of the static Hessian (second variation of `V`):
  require **no negative-eigenvalue directions** beyond the allowed zero modes
  (3 translations + the iso-rotation/orientation moduli). A Floquet test would
  apply only to a periodic orbit, which is not the vehicle here.
- (c) *Spin.* Quantize the `SO(3)`/iso collective coordinate; the rigid `2π`
  rotation must give the π (`ℤ₂`) FR holonomy via the **open-path** Fukui–Hatsugai
  product (closed loop → silent false-negative; see
  [[project_spin_half_is_soliton_layer]]). Already scoped:
  at L5 (D1; `[[project_spin_half_is_soliton_layer]]`).
**Status.** open. The rest baryon is to be sought as a static, periodic-BC
eigenstate (test spec a/b/c) rather than a periodic *orbit* (breather), per the
Statement and Replacement above. Owners: soliton-hunter (static extremum + Hessian) +
physics-derivation (iso-rotation quantization) + berry-validator (L5 spin holonomy).
See `archive/BARYON_SIMULATION_ROADMAP.md` Phase 2.

> *Apparatus note (not theory).* Soliton-scale runs must use periodic spatial BC to
> clamp the prestress and must read winding `B` from the actual field — two harness
> errors (open BC on a tensioned lattice; a hard-coded `F_inner=π` B-estimator)
> produce spurious "stable B=1" results. Detail lives in
> `LESSONS_LEARNED.md`; kept out of this derivation to avoid re-running the mistake.

### C3. How the geometric quartic and SU(3) confinement combine to block expansion — `(a) RESOLVED; (b) moot; (c) test derived`
**(Anti-expansion / IR direction of Derrick; the confinement problem, BACKBONE #24/#25.)**

**RESULT (scaling leg, physics-derivation — negative, division-of-labor).** There is
**no distinct Skyrme/Faddeev gradient-quartic on the central-force substrate, and no orientation
potential.** The norm term `−k_s αa|ΔR|` depends only on the *scalar* `|ΔR|²`, so its quartic is
the **symmetric** contraction `(∂ξ·∂ξ)²` — never the **antisymmetric** Faddeev–Skyrme
`(∂ξ∧∂ξ)²`, which would need a non-central (bending/torsional) interaction or a WZ term. Robust
across isotropic hyperelastic completions (neo-Hookean/Mooney–Rivlin/Ogden all give only
symmetric quartics); flips only under new non-central microphysics.
- **(a) channel structure — RESOLVED, option (i), in the strong form.** The SU(3) texture energy
  is *not even a separate additive channel* — it is literally the **traceless projection of the
  same `E₂+E₄`**: the σ-model 2-gradient piece `∝ A_⊥²|∂n̂|²` *is* `E₂^⊥`, and the only
  4-gradient operator is the geometric `E₄` (symmetric, positive & finite on the hedgehog,
  `E₄ = 2π k_s α a⁻¹ A_⊥⁴ I₄`). λ-powers are `(λ⁺¹, λ⁻¹)`, identical to the lateral sector. **The
  geometric quartic is the *only* `λ⁻¹` stabilizer; there is no independent IR length.**
- **(b) which dominates — moot.** Only one mechanism sets a length (the geometric quartic).
  SU(3)'s role is *not* a length-setting energy channel: it is **topological existence**
  (`π₃` winding) **+ kinematic confinement** (no free colored asymptotic state to disperse
  into) — these block dispersal but set *no* `R_h`. Clean division of labor: SU(3) ⇒ *that* it
  stays a localized bound object; geometric quartic ⇒ *how wide* it is.
- **(c) discriminating size law — DERIVED.** The single law `R_h/a = κ(A_⊥/a)√(α/(1−α))` stands
  (no competing texture law). Falsifiable on the C2 eigenstate: a genuine Faddeev–Skyrme
  stabilizer would absorb the amplitude into the unit field `n̂`, giving `R_h` **independent of
  `A_⊥`** and of `α/(1−α)`. So measure (1) `R_h` vs `α∈{0.3,0.5,0.7,0.8}` at fixed `A` — predict
  `R_h(α)/R_h(0.5) = √(α/(1−α)) = {0.65,1.00,1.53,2.00}` (flat ⇒ falsified); (2) `R_h` vs `A_⊥`
  — predict slope 1 on log–log (slope <0.5 ⇒ hidden texture scale ⇒ violates the all-anharmonicity-
  is-the-norm-term axiom, i.e. new substrate physics).

**One place the verdict could flip (carried-axiom):** kinematic confinement is asserted to set
*no* length (group-theoretic only). If a future derivation produced a confining *linear potential*
in the color sector (a `V(n̂)`, scaling `λ⁺³`), that would be option (ii) and would supersede the
quartic — but no microscopic source for it exists under the norm-term axiom. The numerical leg
(measuring `R_h(α,A)` on a converged C2 state) remains open and is the live test.

Full derivation returned by physics-derivation; key conventions in
`[[project_geometric_nonlinearity_alpha_scaling]]` and `vsh_channel_decomposition.md` §1.4–2.5.

---
**[original scoping, retained]**
**Statement.** The IR (anti-expansion) instability direction is the confinement problem,
and two distinct mechanisms are expected to contribute: (i) the **geometric quartic**
`W₄ ∝ k_s α/a` (Pythagorean norm term, Paper I), which supplies a Derrick `λ⁻¹` balance
against the gradient `λ⁺¹` and sets a finite equilibrium width `R_h/a = κ(A/a)√(α/(1−α))`;
and (ii) **SU(3) confinement** of the color sector (kinematic confinement — no free colored
asymptotic state, [[project_c2_skyrme_no_confinement]] / color bridge #24 — together with the
`π₃(SU(3))` topological winding). **Neither has been shown sufficient alone:** the C2 hardening
sweep with the quartic present still disperses to box-fill, and the kinematic-confinement
argument is so far group-theoretic/dynamics-free (it forbids free *linear* colored states but
has not been shown to set a *length scale* for the bound texture). The open question is how the
two combine **quantitatively**.
**Why it matters.** The decisive numerical target is a stable, finite-width baryon-like soliton
(backbone #12). If the answer is reported as "the geometric quartic stabilizes it," that is the
over-attribution this entry exists to correct: the quartic alone has not confined anything in
simulation. A wrong stabilization story predicts the wrong `R_h(α, A)` law and the wrong
sweet-spot regime, and mis-scopes which sector (geometric vs. color) must be active in the seed.
**Specific sub-questions.**
- (a) *Channel structure.* Are the quartic balance and the SU(3) confinement **additive
  λ-channels** in the Derrick functional (`E(λ) = λ E₂ − (corrections) + λ⁻¹ E₄ + E_SU(3)(λ)`),
  or does the color texture set an **independent IR length** that supersedes the quartic balance
  (so `R_h` is fixed by the texture, not by `√(α/(1−α))`)? Derive `E_SU(3)(λ)` — its Derrick
  `λ`-power for a `π₃` texture in 3D — and compare to the quartic's `λ⁻¹`.
- (b) *Which dominates, and where.* If both are `λ⁻¹` (or comparable), which prefactor wins in
  the `α ≈ 0.5–0.8`, `A/a ≈ 10`, `w ~ a` regime (C2's load-bearing corner)? Does the dominant
  mechanism switch with `α`?
- (c) *Falsifiable size law.* Quartic-only gives `R_h/a = κ(A/a)√(α/(1−α))` (linear in `A`).
  A texture-set scale would give a **different** `R_h(α, A)` — derive it and state the
  discriminating measurement on the converged eigenstate (the same periodic-clamped run as C2).
**Candidate approach.** Write the full Derrick scaling functional `E(λ)` for the hedgehog/Skyrme
seed including the color/texture energy explicitly (not just the transverse-quartic term);
identify each term's `λ`-power; minimize and read off `R_h`. Cross-check against the C2 converged
eigenstate once it exists (measure `R_h` vs `α`, `A` and test the two competing size laws).
This is gated on a converged localized solution (C2) for the numerical leg, but the scaling
derivation can proceed now. Owners: physics-derivation (Derrick functional incl. SU(3) term) +
soliton-hunter (the `R_h(α,A)` measurement). Feeds Paper II §matter (Derrick) and Paper IV
(confinement). See `[[project_c2_skyrme_no_confinement]]`, `[[project_geometric_nonlinearity_alpha_scaling]]`,
`[[project_pythagorean_breather_go]]`, and `papers/gauge_color/derivations/color/status.md`.

### D1. EM-charge ⇄ spin-½ bridge (spin-from-isospin) — `open`
**Statement.** Two emergent quantum numbers are currently scoped *independently*:
(i) the trace-`U(1)` holonomy is read as electromagnetism (backbone #19); (ii)
spin-½ is scoped as an L5 *soliton-rotation* effect — a `π₁(SO(3))=ℤ₂`
Finkelstein–Rubinstein phase under `2π` rigid rotation of a winding-odd hedgehog
(prediction P3; see `[[project_spin_half_is_soliton_layer]]`). Established
field theory (Jackiw–Rebbi 1976; Hasenfratz–'t Hooft 1976; Witten's dyon 1979)
makes these *not* independent: a bosonic soliton carrying one quantum of
`U(1)`/non-Abelian charge acquires **half-integer spin and Fermi statistics** when
the charge collective coordinate is quantized ("spin from isospin") — the same
mechanism that turns a charged Skyrmion/monopole into a fermion.
**Why it matters.** If it holds on the substrate, the `ℤ₂` spin phase is *sourced
by* the trace-`U(1)` charge rather than being an unrelated topological accident.
It yields a sharp, falsifiable correlation that strengthens P3: **only
trace-`U(1)`-charged solitons are fermions (spin-½); `U(1)`-neutral solitons stay
integer/spin-1.** With the working proton/neutron hypothesis
(`archive/BARYON_SIMULATION_ROADMAP.md` Phase 2) this predicts the charged (proton-class,
trace-admixed) and neutral (neutron-class, trace-cancelled) seeds carry the *same*
baryon winding, but the spin-statistics phase tracks the `U(1)` charge — testable
by the same `2π`-rotation holonomy experiment run on charged vs neutral seeds.
**Candidate approach.** Identify the collective coordinate on the lattice hedgehog
whose quantization carries the trace-`U(1)` charge; compute the induced
rotation/exchange phase (Witten-effect / `θ`-term analogue from the geometric
coupling); compare to the directly-measured `2π`-rotation holonomy of charged vs
neutral seeds. This is a *soliton-layer* (L5) statement — the linear spectrum
gives spin-1 at all α (proven, same derivation) — so it is gated on a confined
soliton. Owners: physics-derivation (the spin-from-isospin reduction) +
berry-validator (the L5 charged-vs-neutral rotation holonomy). See
`[[project_spin_half_is_soliton_layer]]`.