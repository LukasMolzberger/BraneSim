# Substrate Bridge — Status

## HAVE

- Discrete 4D brane action with codimension-zero embedding; Lorentzian sign as a property of the action (A5), not the ambient.
- Geometric nonlinearity derivation: the entire anharmonic sector is the norm term −k_s αa|ΔR|, exactly ∝ α; quadratic geometric stiffness ∝ (1−α); quartic coefficient +k_s α/(8a²) > 0 (hardening).
- Lattice-to-continuum reduction: acoustic tensor A_abcd = (k_s/a)[α Q_abcd + (1−α) δ_ac δ_bd]; wave speeds c_L² = k_s a²/m, c_T² = (1−α) k_s a²/m.
- Spring constitutive continuum limit: central-force spring ½k_s(|ΔR|−αa)² established as the constitutive law; StVK demoted to quadratic-order proxy (inverted quartic α-scaling).
- Cross-check of discrete and continuum wave speeds via conventions.py closed-form dispersion.
- Two-time (block) BVP well-posedness **resolved in the linear regime**: a two-past-slice chiral boundary condition gives an N-independent, bounded condition number (OPEN_PROBLEMS A2, `in-progress`, linear regime resolved 2026-05-30; stated in §substrate-model of the paper).
- See: discrete_4d_brane_action.md, geometric_nonlinearity_alpha_scaling.md, lattice_to_continuum.md, spring_constitutive_continuum_limit.md.

## MISSING

- Saddle-point solver for the full Lorentzian action (OPEN_PROBLEMS A1, `open`); the linear-regime well-posedness above does not yet extend to the nonlinear two-time problem with the time-link quartic restored.
- Proof that the 6-neighbor stencil with retuned shell weights produces an isotropic acoustic tensor; or explicit acceptance of the leading-order cubic anisotropy as a feature with quantified magnitude. (The paper takes the latter stance: the cubic anisotropy is the structural source of the gauge sector, not a defect to tune away.)
- Verification of α_t = α consistency for the canonical prestressed vacuum (OPEN_PROBLEMS A4a, `adopted; verification open`).
- Marsden–West reference still uncited: §substrate-model attributes the Störmer–Verlet = variational-integrator result to Marsden–West in prose, but references.bib has no matching entry.

## Open derivations

*Relocated from the former central `OPEN_PROBLEMS.md` (group A, foundational
solver — the 4D block-variational formulation). IDs retained so existing
cross-references (`OPEN_PROBLEMS.md A1`, etc.) resolve here.*

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
(The downstream Lorentz-invariance consequences are tracked in the Lorentz
bridge.)

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
**New test (2026-06-07, from D6).** The time-link binding derivation
(`paper/derivations/time_link_binding.md`) gives a *third* check: the same `r_t` that
sources gravitational time dilation (#22) also sources EM↔color binding (D6 Part 1), so
**binding strength and gravity strength must co-vary in `α`** (`∂ln κ_bind/∂ln α =
∂ln c_grav/∂ln α + 1`). If a simulation can tune binding independently of gravity,
`α_t=α` is falsified. Strengthens (ii). See D6 (color bridge) and the Gravity bridge.