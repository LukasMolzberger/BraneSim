# BraneSim — Open Problems & Open Derivations

**Purpose.** A single central place to track *unsolved* derivational debts and
mathematical gaps. These are **not** part of the theory and must **not** be
presented as structure in the theory-structure diagram or foregrounded in the
paper. If a problem remains unsolved when the paper is finalized, it may be
mentioned briefly there as future work — but the working record lives here, not
in the manuscript.

Conventions:
- **Status:** `open` · `in-progress` · `resolved` (link the result) · `parked`
- Keep each entry to: *statement → why it matters → candidate approach → status*.
- When an item is resolved, move its conclusion into `principles.md` /
  `paper/backbone.md` / the paper, and mark it `resolved` here with a pointer.

Related living documents (not duplicated here):
- `principles.md` — canonical non-negotiables
- `paper/backbone.md` — canonical project backbone
- `paper/validation_roadmap.md` — numerical validation sprints
- `paper/baryon_simulation_roadmap.md` — soliton search program
- Paper §2 "Non-claims" — honest scope limits (distinct from open *derivations*)

---

## A. Foundational solver — the 4D block-variational formulation

These arose from making the static 4D world-volume the foundational object
(brane action `S[R]`, time-symmetric, stationary point under boundary data)
rather than forward Verlet marching. The local Euler–Lagrange stencil is
unchanged (it *is* the Störmer–Verlet / variational-integrator stencil); what is
open is the **global solution strategy and its well-posedness**.

### A1. Saddle-point solver for the Lorentzian action — `open`
**Statement.** The discrete action `S = Σ Δt (T − V)` has Lorentzian sign
structure and is therefore **unbounded below** — a saddle, not a minimum.
Gradient descent / energy relaxation on `S` diverges along the kinetic
direction (and, if both link types were taken positive, would silently solve the
*Euclidean* heat-equation problem instead of Newtonian/wave dynamics). The
foundational solver must therefore **root-find `∇S = 0`** (e.g. Newton–Krylov on
the discrete d'Alembertian, or minimize `‖∇S‖²`), not minimize `S`.
**Why it matters.** Determines the entire architecture of the foundational
(non-Verlet) solver track. Forward Verlet is the well-posed IVP special case;
the block solve is the general case needed for the retrocausal/two-boundary
problem.
**Candidate approach.** Linear case → sparse linear solve of the discrete
hyperbolic operator with prescribed boundary data. Nonlinear (StVK/soliton)
case → nonlinear BVP via Newton–Krylov with the Verlet stencil as the interior
residual.

### A2. Well-posedness of the two-time (block) BVP — `in-progress` (linear regime resolved 2026-05-30)
**Statement.** Hyperbolic operators are naturally Cauchy/IVP problems. Posing
the brane equation as a two-time boundary-value problem (data on past **and**
future spacelike slices) is **not unconditionally well-posed**: when the
time-extent hits a normal-mode resonance of the spatial operator, the
homogeneous problem has nontrivial solutions and the BVP becomes singular /
non-unique. (The Euclidean analogue is always well-posed; the Lorentzian one is
not.)
**Why it matters.** This is the price of the block-universe picture. A naive
"prescribe `q⁰` and `qᴺ`, solve the interior" can be non-unique exactly in the
regimes of interest.
**Candidate approach / conjecture.** The physical selection principle that
removes the ambiguity may be the *same* one already in the ontology:
**soliton chirality / causal-direction selection** (matter = forward, antimatter
= backward worldtube). Test whether the chirality boundary condition and the
BVP well-posedness condition are the same statement in disguise. See
`[[retrocausal-worldtube-interpretation]]` (project memory).

**Resolution (linear regime, 2026-05-30 — derivation + numerical confirmation).**
The conjecture is **confirmed, as an exact algebraic identity.** Since the
6-neighbor `D(k)` is diagonal, the block BVP decouples per spatial mode into the
scalar recurrence `a^{l+1} − 2cosθ(k) a^l + a^{l−1} = 0`,
`θ(k)=arccos(1−Δt²ω²(k)/2)`. The Dirichlet two-time operator has determinant
`2i·sin(Nθ)` → singular at `Nθ=mπ` (resonances), generically ill-conditioned
(`κ ~ 1/min_k|sin(Nθ(k))|`; ~most time-extents fail at realistic mode counts —
verified: `|det|=2|sin(Nθ)|` to machine precision, 67/79 extents with `κ>10³` for
an evenly-spaced 32-mode spectrum). Splitting into characteristics
`a^l=a₊e^{−ilθ}+a₋e^{+ilθ}` and imposing **one characteristic per end** (fix
forward `a₊` from the past slice; no-incoming/Sommerfeld `a^N−e^{−iθ}a^{N−1}=0`
on the future, killing `a₋`) gives `κ=1` exactly ∀N — and this *is* the
matter=forward/antimatter=backward selection (exact change of basis). DC modes
(`θ→0`) route through Dirichlet (well-posed there). Implementation recipe in
`ARCHITECTURE.md` §2 (D2). **Still open:** uniqueness of the *nonlinear* block BVP
— the `κ=1` result is the conditioning of each linearized JFNK step (necessary,
not sufficient globally); the geometric quartic re-couples modes. See
`[[block-solver-bvp-chirality]]`.

**Refinement — verdict (a), 2026-05-31 (implemented, supersedes the
characteristic-future-condition recipe above).** Imposing "kill `a₋` per mode" is
correct 2×2 algebra but **wrong for a real field**: reality couples
`a₊(k)=conj(a₋(−k))`, so zeroing `a₋` per mode deletes both characteristics →
non-real garbage. The correct chiral BC is simply **two adjacent past slices
`(R⁰,R¹)` marched** (forward=matter, backward=antimatter); one slice can't encode
direction, two can. Reality automatic; condition bounded `O(10)`. Implemented in
`branesim/solver/{boundary,bvp}.py` (commit 4a2985c), tested to 8e-14 recovery at
a Dirichlet-resonant N. **Implication:** the well-posed chiral solve = Cauchy
march, so the **two-time Dirichlet BVP is the wrong vehicle for a bound
particle** — that requires a **time-periodic (`R^{l+P}=R^l`) + spatially-localized
nonlinear eigen-BVP** (the carrier+envelope worldtube, #18). New sub-problem to
track below.

### A3. `c²` tuning for emergent Lorentz invariance — `open`
**Statement.** With explicit temporal links, the ratio of temporal stiffness
(`m/Δt²`) to spatial stiffness (`k`) sets the lattice light-cone speed,
`c² ≈ k Δx² / m`. Emergent isotropic Lorentz invariance in the long-wavelength
limit is a **tuning condition** on that ratio and on `k·a` (small), not
automatic.
**Why it matters.** "Observe Lorentzian action but do not model it actively"
(principles §5 ontology) is only achievable in the tuned, long-wavelength
regime.
**Candidate approach.** Connect to the dispersion layer (`c_T²`, `c_L²`,
anisotropy vs `k·a`); quantify the tuning window. Owner: dispersion-analyst.

### A4. Temporal-link form — `resolved + implemented (2026-06-05): one 4D-isotropic spring`
**Resolved.** The timelike link is a central-force spring `½ k_t (|ΔR| − r_t)²` —
the same law as the three spacelike links — so the substrate is a fully symmetric
4D-cubic spring lattice. There is **one model parameterized by the temporal rest
length `r_t`**, not an "a vs b" fork:
- `r_t = α·β·dt` (single prestress, A4a): the canonical prestressed substrate.
- `r_t = 0`: the **linear/Verlet limit** — `½ (m/Δt²) |ΔR|²` kinetic, EL = plain
  Newton = the Verlet stencil. Used for small-amplitude wave/dispersion validation.

**Exact reduction.** At matched stiffness `k_t = m/Δt²` the spring force reduces
*term-for-term* to the Verlet stencil as `r_t → 0` (the energies differ by
`−k_t r_t |ΔR| + ½ k_t r_t² ∝ r_t`); verified machine-precision in the
implementation. See `discrete_4d_brane_action.md` §6.

**What `r_t` controls.** The distinguishing term carries the Euclidean *norm*
`|ΔR|` (not `|ΔR|²`), which generates a **geometric quartic** in the time link with
coefficient `∝ r_t` — the time-link analogue of the Skyrme/contraction term (#17)
and the natural home for backbone #22's unified contraction (gravitational time
dilation). The `r_t = 0` limit has no such quartic. *(The distinction is the
norm-nonlinearity, not transverse stiffness — `½ k_t|ΔR|²` already has isotropic
Hessian.)*

**Implemented (2026-06-05).** `ActionParams` exposes `r_t` (plus `k_t`, `beta`); the
residual routes on `r_t` — `r_t = 0` uses the linear Verlet stencil (the fast exact
path and the **default**), `r_t > 0` uses the central-force spring force (block
solve). The `r_t → 0` reduction is bit-identical to the prior linear path; principles
audit clean. The Lorentzian action is a saddle, so the foundational solve is
root-find `∇S = 0`, not minimization. **Consequences:** forward Verlet / IVP is the
`r_t = 0` limit, not the general dynamics — the block solver is mandatory for
`r_t > 0`; A2 well-posedness must be re-checked with the time-link quartic restored;
A4a opens. Folded into backbone #21/#22. Owner: contraction-channel (gravity
derivation).

### A4a. Single 4D prestress `α_t = α` — `adopted (2026-06-05); verification open`
**Decision.** A **single** prestress `α` governs all four lattice directions: every
link's rest length is `α ×` its held spacing (spatial `ℓ₀ = αa`, temporal `r_t = α·βΔt`).
This is the symmetric reading of the 4D-isotropic spring (A4) and ties soliton binding (spatial
geometric quartic `∝ α`) and gravitational time dilation (temporal quartic `∝ α`) to
**one** dial — the unified-contraction picture (#22) taken literally, and a strong
falsifiable structural prediction: gravity strength and soliton confinement are **not**
independently tunable.
**Open (verification, not the decision).** `α_t = α` is adopted as the working model but
is not yet proven consistent. Two checks must pass: (i) a single `α` keeps the emergent
long-wavelength light-cone **isotropic** (this is the `c²`-tuning of A3, now *constrained*
by `α` rather than free); (ii) it reproduces the correct **Newtonian limit** in the
contraction channel. If either fails, `α_t` must be freed from `α` again. Owner:
contraction-channel + physics-derivation.

---

## B. Quantum-foundations derivations (Bell stance)

Carried over from `principles.md` §1.5 and `paper/backbone.md` #23. These are
constraints the retrocausal worldtube interpretation places on the program —
honest debts, not results. **Do not claim as established in the paper.**

### B1. Tsirelson's bound `2√2` — `open`
Reproduce the exact CHSH maximum from V-branching soliton correlations on the
substrate. Retrocausal models notoriously under-/over-shoot; deriving the exact
bound (not merely nonzero correlation) is the most direct empirical handle on
the stance.

### B2. No-signalling theorem — `open`
Derive no-signalling as a *theorem* about branching-worldtube dynamics (no
one-branch observable depends on the other branch's measurement context),
rather than assuming it. This is what prevents the retrocausal channel from
sending a macroscopic bit backward in time.

### B3. Baryon-to-photon ratio `η ≈ 6×10⁻¹⁰` — `open`
Show whether the geometric-vertex picture (matter/antimatter expanding from the
Big Bang vertex in opposite time directions) can match `η` quantitatively,
rather than merely permitting opposite-direction expansion in principle.

---

## C. Constitutive law & soliton stability

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
excluded (`test-runs/sprint4b_skyrme_corrected/trend_sweep.csv`). The breather
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
  `test-runs/alpha_separability/L5_spin_half_target.md`.
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
`BAYRON_SIMULATION_ROADMAP.md` Phase 2.

---

## D. Gauge sector — emergent quantum numbers

### D1. EM-charge ⇄ spin-½ bridge (spin-from-isospin) — `open`
**Statement.** Two emergent quantum numbers are currently scoped *independently*:
(i) the trace-`U(1)` holonomy is read as electromagnetism (backbone #19); (ii)
spin-½ is scoped as an L5 *soliton-rotation* effect — a `π₁(SO(3))=ℤ₂`
Finkelstein–Rubinstein phase under `2π` rigid rotation of a winding-odd hedgehog
(prediction P3, `test-runs/alpha_separability/derivation_H_eff.md`). Established
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
(`BAYRON_SIMULATION_ROADMAP.md` Phase 2) this predicts the charged (proton-class,
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

### D2. Prestress α from the EM/colour coupling ratio — `resolved (linear: undetermined)` · `open (nonlinear)`
**Statement.** Can the empirical EM-vs-strong coupling ratio fix the prestress α?
The linear *spectrum* cannot: the traceless/`SU(3)` content `∝ α` and the trace/
`U(1)` content `∝ (1−2α/3)` stay comparable (≈13% apart at α=0.2) and both vanish
together as α→0, so the spectrum produces **no** EM/colour hierarchy — the famous
`α_EM/α_s` must live in the *connection/curvature normalization*, not eigenvalue
magnitudes. Compounding this, the **k-space Berry/WZ curvature is identically zero
∀α** (real-symmetric `D(k)`; `derivation_H_eff.md`), so any coupling ratio is a
*fibre-internal* `(x,t)` holonomy object.
**Why it matters.** A derived α would convert the project's one free linear
parameter into a prediction. A clean *negative* result (α undetermined at linear
order) is itself valuable: it localizes the EM/colour hierarchy in the nonlinear
normalization and rules out spectrum-ratio "derivations" of `1/137`.
**Lorentz-bound corollary (already usable).** Mutual photon/gluon Lorentz
invariance (no measured vacuum birefringence between EM and colour) does *not*
bound α small; it forces physical carriers onto the `[111]` isotropic locus where
all three lateral eigenvalues are exactly degenerate ∀α (`g([111])=0`). This
*frees* α to sit at O(0.1–0.5) and is the cleanest current "SM fact → lattice"
import.
**Resolution (linear order, 2026-06-04 — derivation; α NOT determined, by design).**
The fibre-internal holonomy ratio is `R(α) = (C/g(k̂))·(3−2α)/(√3 α)`, where
`g(k̂)=√(Σ_a(h_a−H/3)²)/H` is the directional anisotropy and `C = N_U(1)/N_SU(3)`
is the **relative gauge-kinetic normalization** = `tr(P_U(1)𝒢)/tr(P_SU(3)𝒢)` set by
the fibre metric `𝒢_ij=⟨∂_iΨ|∂_jΨ⟩`. At linear order `𝒢 ∝ 𝟙` (envelope norm is
flat ℓ² on ℂ³), so `C` reduces to pure dimension-counting (`1/√2`) with **no
α-dependence** — the EM/colour hierarchy is therefore *absent* from the linear
theory and the naive inversion gives `α≈1.2 > 1` (out of range, a self-flag that
the linear identification is invalid). **α is normalization-undetermined at linear
order.** The missing ingredient is the *sector-dependent* fibre metric supplied by
the geometric quartic (backbone #17), which couples the trace/dilatational
direction to the `X⁴` gravity channel differently from the traceless/shear
directions — a nonlinear computation. **Sharpest new gap surfaced:** the triplet is
coherent only for `k₀∥[111]`, but `g([111])=0` kills the colour content there
(`R→∞`); off-[111] the branches dephase. There is **no linear direction that is
both coherent and colour-active** — color curvature is intrinsically a
soliton-layer (nonlinear) object. **Falsifiable, normalization-independent
spin-off (decisive):** the α-scaling is `R(α₂)/R(α₁) = [(3−2α₂)/α₂]/[(3−2α₁)/α₁]`;
for `(α₁,α₂)=(0.2,0.5)` this is **0.3077** — measurable by extending
`branesim/diagnostics/berry_holonomy.py` to transport a band-isolated ℂ³ envelope
around a fixed `(x,t)` loop off-[111]; a >10% deviation falsifies the
spectral-susceptibility factorization. Full derivation:
`paper/derivations/alpha_holonomy_estimator.md`.

### D3. Dynamical Maxwell: effective action for the EM four-potential — `open`
**Statement.** Electromagnetism is the `U(1)`-trace Berry connection of the carrier — the
**electromagnetic four-potential** `A_μ = i⟨u|∂_μu⟩` — with **field tensor**
`F_μν = ∂_μA_ν − ∂_νA_μ` (`E_i = F_{0i}`, the time–space part; `B_i = ½ε_{ijk}F_{jk}`, the
space–space part). The **kinematics is direct** (no dynamics needed): (a) gauge invariance
`A_μ → A_μ + ∂_μχ` is the unobservable choice of carrier phase reference; (b) the **homogeneous**
Maxwell equations are the Bianchi identity `dF = 0`, automatic because `F = dA`. What is **not**
derived is the **dynamics** — the **inhomogeneous** equations `d⋆F = J` — which hold only if
integrating out the substrate around the carrier band yields the Maxwell action
`S_eff[A] = −¼ ∫ F∧⋆F`, requiring: (i) the right **form** (`∝ F²`, not another functional of `F`);
(ii) **masslessness** (no `m²A²`; the `U(1)` must be an exact flat direction); (iii) the emergent
**Hodge `⋆`** = isotropic light-cone (ties A4a / `[111]`); (iv) the **coupling** normalization
(undetermined at linear order, D2 — lives in the nonlinear fibre metric). Also open: what
dynamically forces a **nonzero, propagating** `F` (the phase must be non-integrable — the carrier
state must rotate in `ℂ^N`, `N≥2`, not merely carry an overall phase; a pure gradient is pure
gauge, `F=0`).
**Why it matters.** Closes the EM sector from "gauge kinematics" (done) to "Maxwell dynamics"
(the field equations + the photon). Until derived, the paper claims only the kinematic/gauge
structure, not full Maxwell dynamics (§2 non-claims). Recurs across §5.6, D2, and the complex-`U(1)`
discussion. The complex structure / the `i` is the carrier's rotation along the timelike worldtube
axis (`[[project_complex_u1_from_time]]`), so `A_μ` is intrinsically a 4-covector on `(x,t)`.
**Candidate approach.** Effective-field-theory reduction of the substrate action about the carrier
band → `S_eff[A]`; verify the `F²` form + masslessness; read the coupling from the fibre metric
(D2). Owner: physics-derivation + berry-validator. See `paper/05` §5.2.

### D4. EM = `U(1)` vortex / color = `SU(3)` texture; do they bind, and does color stabilize the electron? — `open` (active experiment, 2026-06-06)
**Statement.** The `U(3)=U(1)×SU(3)` field carries two topologically distinct particle
sectors (BACKBONE #25): the EM/electron-like object is a `π₁(U(1))=ℤ` carrier-phase
**vortex** (codim-2 worldsheet, donut cross-section, spin-½ from the `2π` axis tumble,
charge from the temporal carrier); the color/baryon is a `π₃(SU(3))=ℤ` **texture**
(`π₁(SU(3))=0` — no vortex). They are one field, so **protection ≠ participation**:
`SU(3)` is dynamically coexcited around a `U(1)` vortex and may be load-bearing even at
net-zero color (a color-singlet electron `SU(3)`-stabilized).
**Open questions.** (i) Does the substrate **bind** the `U(1)` vortex worldtube — which,
with `SU(3)` unfrozen, is a **semilocal vortex** (simply-connected full vacuum, *not*
topologically protected), so binding is a **dynamical** condition (the `β<1` /
hardening-bound-mode regime; Vachaspati–Achucarro), on a clean periodic-clamped `r_t>0`
eigenstate (supersedes the retracted C2-class attempts)? (ii) Does the `SU(3)` sector
coexcite and **stabilize** the color-neutral electron (its non-Abelian/semilocal core
moduli)? (iii) **Scale hierarchy:** does the size ratio of the `U(1)`
vortex vs the `SU(3)` texture track `f(α)` — the candidate origin of the QCD-vs-EM
scale split? *Caveat:* the canonical `U(1)/SU(3)` split is dynamical (running
couplings); a static `α` most naturally gives a coupling *ratio*, so a true *scale*
hierarchy is not guaranteed (no geometric-`α` precedent — see the prior-art note).
(iv) the worldsheet's periodic-closure / `(ω,T)` quantization condition. (v) EM =
the **trace `(1,1,1)`** lateral direction — a single-component seed is only `⅓` EM.
**Status.** This is the active single experiment — full spec in `EXPERIMENT.md`
(object, injection ansatz, one parameter set, measurement suite, scale ladder).
Owner: physics-derivation (ansatz/closure) + the build.

---

## E. Pointers to scoped non-claims (tracked elsewhere)

The following are bounded-scope limits already stated honestly in paper §2
"Non-claims"; listed here only so the central tracker is complete. They are
*scope statements*, not active derivation tasks, unless promoted above:
derived value of `G` (Newtonian limit of the contraction channel), full Maxwell
/ Yang–Mills dynamics in all regimes, a finished proton/neutron solution, the
full Standard-Model gauge group / generations / couplings, a closed
multi-particle scattering theory.
