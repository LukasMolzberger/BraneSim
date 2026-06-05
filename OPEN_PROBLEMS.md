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

### A4. Temporal-link form — model (a) vs model (b) — `resolved (model b, 2026-06-05)`
**Statement.** What is the precise form of the timelike link?
- **Model (a) — kinetic (baseline, = current code):** zero-rest-length
  quadratic `½ (m/Δt²) |ΔR|²`. EL = plain Newton = Verlet stencil.
- **Model (b) — central-force temporal spring:** `½ k_t (|ΔR| − r_t)²` with its
  own rest length `r_t` and prestress, Lorentzian sign. Fully symmetric
  4D-cubic; EL is non-Newtonian.

**Exact convergence.** At matched stiffness `k_t = m/Δt²` the two link energies
differ by a single closed-form term:
`E_(b) − E_(a) = −k_t r_t |ΔR| + ½ k_t r_t²`. The difference is `∝ r_t`, so the
models coincide **identically iff `r_t = 0`** (`α_t = 0`); no condition on the
held spacing `β Δt` is needed. See `discrete_4d_brane_action.md` §6.

**Why it matters / what `r_t` controls.** The distinguishing term carries the
Euclidean *norm* `|ΔR|` (not `|ΔR|²`), which generates a **geometric quartic**
in the time link with coefficient `∝ r_t` — the time-link analogue of the
Skyrme/contraction term (#17) and the natural home for backbone #22's unified
contraction (gravitational time dilation = spatial displacement stretching the
timelike link). Model (a) (`r_t = 0`) has no such quartic, so the convergence
point is exactly where model (b) loses its distinctive physics. `r_t` is the
dial between "matches Newton + existing code" and "a gravitating time link".
*(Note: the distinction is the norm-nonlinearity, not transverse stiffness —
`½ k_t|ΔR|²` already has isotropic Hessian.)*

**Resolution (2026-06-05, owner decision): model (b) chosen.** The timelike link is a
genuine central-force spring `½k_t(|ΔR|−r_t)²` with its own rest length `r_t ≠ 0`
(temporal prestress `α_t`). Rationale: (i) keeps the 4D lattice genuinely symmetric
rather than special-casing time; (ii) `r_t` is a real degree of freedom and should be
exposed, not hidden; (iii) the foundational block solver (root-find `∇S=0`) is required
regardless; (iv) backbone #22's unified contraction (gravitational time dilation) *needs*
the time-link geometric quartic, which only model (b) provides — so #22 was already
implicitly assuming it. **Consequences:** forward Verlet / IVP (model a) is now the
`r_t=0` *special case*, not the general dynamics — the block solver is mandatory; A2
well-posedness must be re-checked with the time-link quartic restored; and a new
sub-question opens (A4a). Code note: `ActionParams` keeps `temporal_model="a"` (r_t=0) as
the validated special case until the block solver implements model (b). Moved into
backbone #21/#22. Owner: contraction-channel (gravity derivation).

### A4a. Single 4D prestress `α_t = α` — `adopted (2026-06-05); verification open`
**Decision.** A **single** prestress `α` governs all four lattice directions: every
link's rest length is `α ×` its held spacing (spatial `ℓ₀ = αa`, temporal `r_t = α·βΔt`).
This is the symmetric reading of model (b) (A4) and ties soliton binding (spatial
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

---

## E. Pointers to scoped non-claims (tracked elsewhere)

The following are bounded-scope limits already stated honestly in paper §2
"Non-claims"; listed here only so the central tracker is complete. They are
*scope statements*, not active derivation tasks, unless promoted above:
derived value of `G` (Newtonian limit of the contraction channel), full Maxwell
/ Yang–Mills dynamics in all regimes, a finished proton/neutron solution, the
full Standard-Model gauge group / generations / couplings, a closed
multi-particle scattering theory.
