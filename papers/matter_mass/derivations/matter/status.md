# Matter Bridge — Status

## HAVE

- Self-adjoint eigenproblem formulation for localized modes.
- VSH angular structure: hedgehog ansatz (J=0, L=1) as canonical baryon seed; Skyrme-twisted extension with F(r): [0,∞] → [π,0].
- Soliton labels: (J,P; SU(3) irrep; Q_U(1); B_winding).
- Kinematic confinement of color (coherence and color-activity mutually exclusive on cubic lattice).
- Derrick theorem analysis: anti-collapse blocked by lattice UV cutoff (spacing a; Derrick's λ→0 rescaling halts at one cell; single-node concentration immediately disperses); anti-expansion blocked by geometric quartic ∝ k_s α/a (corrected from earlier inverted StVK scaling). Derrick balance (gradient vs. quartic) gives equilibrium size R_h/a = κ (A/a) √(α/(1−α)).
- Corrected soliton sweet spot: α ≈ 0.5–0.8, operational seed target α=0.7, A/a ≈ 10 → R_h/a ≈ 8.
- Periodic BC requirement established: open BC causes tension relaxation (prestressed lattice contracts to αa), collapsing even A=0 vacuum; periodic BC fixes box at N·a, holds tension, vacuum stable.
- Preliminary U(1) vortex relaxation attempted (non-converged; semilocal drift observed — expected behavior).
- ALL 2026-06-05 nonlinear runs RETRACTED (unclamped pretension + hard-coded B-estimator invalidated static/box-doubling/w-scan/worldtube runs).

## MISSING

- Converged numerical solution for any localized mode with periodic-clamped tensioned vacuum.
- Rigorous Derrick's theorem write-up (cite theorem directly; state lattice UV cutoff argument cleanly; tie quartic IR balance to effective field manifold).
- Proof of existence (not just Derrick balance — need topological stability argument independent of Derrick).
- Baryon-triplet simulation with correct BC and eigenstate approach (not seed-and-watch).

## Open derivations

*Relocated from the former central open-problems tracker (group C, constitutive law &
soliton stability; plus D1, the EM-charge ⇄ spin-½ bridge). Original IDs (C1, C1a, C2, D1) retained.*

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
**Why it matters.** Any soliton-width / hedgehog-radius formula built on the StVK
quartic prefactor must be re-derived. `vsh_channel_decomposition.md` §2.5's old
`R_h/a ~ (ℓ₀/u₀)^{2/3}` is superseded by `R_h/a = κ(A/a)√(α/(1−α))` (corrected
balance, 2026-06-04): radius grows with α as `√(α/(1−α))` and is *linear* in
amplitude `A` — the opposite of the old note on both axes. The Derrick `λ`-scaling
(backbone #17) is α-independent and unaffected; only the coefficient was wrong, so
no sign error in the α≈0.5–0.8 sweet-spot story.
**Status.** open. Owner: soliton-hunter. Feeds the sprint4 re-run target
(α≈0.7, A/a≈10). See `[[project_geometric_nonlinearity_alpha_scaling]]`.

### C2. Rest baryon — `open` · ⚠ ALL prior nonlinear runs RETRACTED 2026-06-05 (unclamped pretension + hard-coded B-estimator); restarting periodic-clamped + eigenstate
**Reframed 2026-06-06 (BACKBONE #25, D4):** the baryon is specifically the `SU(3)`
`π₃` **texture** (a localized Skyrmion/instanton, codim-3) — topologically a different
object from the `U(1)` carrier-phase **vortex** (the EM/electron-like object, codim-2,
`π₁`). The **active first experiment is the `U(1)` vortex** (`EXPERIMENT.md`/D4),
which is `π₁`-protected and the cleaner test; the `SU(3)` texture (this entry) follows.
**Statement.** The time-periodic, amplitude-breathing baryon ansatz (the
common-carrier Skyrme breather, `solve_breather(mode="topological")`) **exists
but is dynamically unstable, and prestress α does not stabilize it.** Converged
32³ orbits at α=0.80/0.85/0.90 give Floquet `ρ = 3.59 / 3.47 / 3.70` — robustly
`>1`, flat in α (ticks *up* at 0.9), and robust to `inner_maxiter` (physical, not
a solver artifact); the earlier apparent "α-stabilizing trend" (95→54→3.6) was an
artifact of non-converged orbits at α=0.5/0.7 (residual 8–27) and must be
excluded (sprint4b trend-sweep; data retired with `test-runs/`, result in `LESSONS_LEARNED.md`). The breather
*mechanism* is sound — Pythagorean hardening `Q = 8k_sα/a² > 0` places a localized
mode above the transverse band ([[project_pythagorean_breather_go]]) — but the
breathing/scale degree of freedom itself carries the instability.
**Why it matters.** A stable baryon-like soliton is the decisive numerical target
(backbone #12). The breather was the leading candidate vehicle
([[project_baryon_search_vehicle_breather]]); it is now closing negative, so the
program needs a replacement.
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
  hedgehog + `X⁴`-twist seed at α≈0.5–0.8, large `A`, **open** boundary, box ≫ seed;
  verify `B(slice)=1` preserved and the config stays localized
  (`spread_ratio ≪ 1`, not box-fill — use the corrected metrics from
  [[project_c2_skyrme_no_confinement]], never `leakage_fraction`).
- (b) *Stability.* Spectrum of the static Hessian (second variation of `V`):
  require **no negative-eigenvalue directions** beyond the allowed zero modes
  (3 translations + the iso-rotation/orientation moduli). This replaces the
  Floquet test, which only applies to the (now-abandoned) periodic orbit.
- (c) *Spin.* Quantize the `SO(3)`/iso collective coordinate; the rigid `2π`
  rotation must give the π (`ℤ₂`) FR holonomy via the **open-path** Fukui–Hatsugai
  product (closed loop → silent false-negative; see
  [[project_spin_half_is_soliton_layer]]). Already scoped:
  at L5 (D1; `[[project_spin_half_is_soliton_layer]]`).
**⚠ RETRACTED — 2026-06-05. ALL executed nonlinear-soliton runs in this thread are INVALID.**
Two independent harness errors invalidate every C2-thread numerical result — the
"corrected harness" `N_neg=0` minimum, the box-doubling "boundary-confined" verdict,
the w-scan, the saturation analysis, and the worldtube-tornado run:

1. **Pretension not clamped (dominant error).** The runs used **open spatial
   boundaries** on a *prestressed* lattice (bonds stretched to `a`, rest length
   `αa < a`, so the medium is under tension). An unpinned tensioned network relaxes
   its prestress — it contracts globally toward `αa`. The **`A=0` control proves it**:
   with *no soliton at all*, the open-BC march still collapses to `E_excess ≈ −V_vac`
   (every bond reaches rest length, releasing the entire tension). The "vacuum" was
   never held at the intended tensioned equilibrium, so the dynamics measured tension
   relaxation, not soliton physics. Test-spec (a)'s "**open** boundary" above is the
   error and is withdrawn.
2. **B was never measured.** `compute_B_analytic` hard-codes `F_inner = math.pi`, so
   it returns `B ≈ 1` whenever the boundary is at vacuum — *regardless of the interior
   field*. Every "B = 1 preserved" claim is this artifact. Under static FIRE the field
   in fact **unwound to the square-lattice vacuum**: the celebrated `N_neg=0, V*=273.8`
   minimum is exactly the `13³` vacuum prestress floor
   (`(k_s/2)·3·13²·12·(0.3)² = 273.78`) — the stable *vacuum*, mislabeled `B=1`.

**Corrected standard for the restart (non-negotiable):**
- **Clamp the pretension with periodic spatial BC.** The period fixes the box at
  `N·a`, holds the tension, and the tensioned lattice is then a *stable* vacuum —
  confirmed: `A=0` periodic → `E_excess ≡ 0` to machine precision.
- **A bound particle is an eigenstate**, not a marched/relaxed seed. Every
  seed-and-watch run here dispersed / unwound / contracted because the seed was not a
  stationary solution. Solve for the self-consistent (rotating worldtube) state.
- **Real diagnostics only:** vacuum-subtracted excess energy against the stable
  periodic vacuum, and a winding degree computed from the *actual* field. The
  hard-coded B-estimator is deleted.

**What survives (unaffected by the bug):** the linear dispersion/anisotropy layer
(machine-exact, periodic, small-amplitude), the `bvp_chiral` solver, the
geometric-nonlinearity-∝α decomposition, and kinematic color confinement (backbone
#24). The breather Floquet result (`ρ≈3.5`) is a distinct code path; held pending
re-validation under this standard, but its instability is not attributable to the
pretension bug.

**Status.** open — restarting from the periodic-clamped tensioned vacuum (step 1)
toward a periodic-BC **eigenstate** baryon (step 2). The static-min test spec (a/b/c)
above is superseded by the eigenstate approach.

**Competing live hypothesis (kept, ranked below).** Larger *breather* solitons
stabilize (48³/64³, AWS). Ranked below the static route because `ρ≈3.5` is flat in
α and shows no approach to 1 at the resolutions tried — no evidence stability is
hiding just past 32³.
**Owner.** soliton-hunter (static extremum + Hessian) + physics-derivation
(iso-rotation quantization) + berry-validator (L5 spin holonomy). See
`archive/BARYON_SIMULATION_ROADMAP.md` Phase 2.

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