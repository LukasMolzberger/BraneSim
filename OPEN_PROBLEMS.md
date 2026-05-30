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

### A2. Well-posedness of the two-time (block) BVP — `open`
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

### A4. Temporal-link form — model (a) vs model (b) — `parked`
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

**Status / resolution criterion.** Parked. Settled by the gravity/contraction
channel (layer 6): if backbone #22's time-dilation face requires the time-link
geometric quartic, model (b) with `r_t ≠ 0` is forced; otherwise keep (a).
Owner: contraction-channel.

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

## C. Pointers to scoped non-claims (tracked elsewhere)

The following are bounded-scope limits already stated honestly in paper §2
"Non-claims"; listed here only so the central tracker is complete. They are
*scope statements*, not active derivation tasks, unless promoted above:
derived value of `G` (Newtonian limit of the contraction channel), full Maxwell
/ Yang–Mills dynamics in all regimes, a finished proton/neutron solution, the
full Standard-Model gauge group / generations / couplings, a closed
multi-particle scattering theory.
