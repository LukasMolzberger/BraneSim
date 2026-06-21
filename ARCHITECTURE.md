# BraneSim — Experimental Code Architecture (block-solver-centric blueprint)

**Status:** blueprint for the `branesim/` reimplementation (now implemented; some sections
remain forward-looking). The former `components/` (forward-Verlet) code has been
**replaced by the `branesim/` package** described here.
This document is meant to be specific enough to reimplement from scratch.

**Scope note.** The *foundational solver algorithm* (root-finding the Lorentzian
action over a 4D block) is an **open research problem**
(`papers/core/derivations/status.md` §Open-derivations A1/A2). Where the algorithm is settled, this doc specifies it; where it is not,
it marks a **DESIGN DECISION** with candidate approaches and a recommended
default. Those decisions must be resolved (a short derivation/spike) before the
solver core is finalized — but the rest of the architecture does not depend on
which choice is made.

Canonical sources this blueprint is bound to: `BACKBONE.md` (#1, #15, #21,
#22), `PRINCIPLES.md` (§1.1, §1.2, §1.5, §6), `paper/derivations/discrete_4d_brane_action.md`
(the action and the EL=Verlet identity), `papers/core/derivations/status.md` §Open-derivations (group A).

---

## 0. Why a rewrite, and what changes

The legacy code is a **forward velocity-Verlet** initial-value marcher: prescribe
one slice of state `(R⁰, V⁰)`, step `l → l+1`. That is the *special case*. The
foundational object of the theory (backbone #1, #21) is the **full 4D
world-volume as a stationary point of the brane action `S[R]`**, found as a
**two-time boundary-value problem** (data on a past *and* a future spacelike
slice, interior root-found). This is the retrocausal/time-symmetric reading
(principles §1.5): causal direction is not built into the marcher but **selected**
by the solution (soliton chirality).

The rewrite makes the **block solver the core**. Forward Verlet survives only as
a degenerate IVP mode used for regression against the validated linear-dispersion
results.

**What is settled physics (the new code must reproduce it, bit-for-bit where
noted):**
- The discrete action §3.1 and the interior stationarity stencil §3.3.
- 6-neighbor axial-only spacelike stencil; central-force pair-spring energy.
- Conventions: `α := rest_length/spacing` (default 0.2); units `k_s=a=ρ=1`.
- Linear light-cone speeds `c_L²=k_s a²/m`, `c_T²=(1−α)k_s a²/m`; at α=0.2,
  `c_T/c_L=√0.8≈0.8944`. (Sprint 1 validated this to ~1e-5; it is the regression
  target.)

**What is genuinely new:** the block/BVP solver driver, the residual `∇S`
assembly, boundary-data specification on two slices, and block-specific
diagnostics (residual norm, BVP conditioning, chirality/causal measure).

---

## 1. The foundational object and the solve

### 1.1 State
A **4D world-volume**: a function on a 4D lattice
`(i, j, k, l) ↦ R ∈ ℝ⁴`, where `(i,j,k)` index the three **spacelike**
directions and `l ∈ {0,…,N}` indexes the **timelike** direction. `R_p^l ∈ ℝ⁴`
is the ambient position of node `p=(i,j,k)` on slice `l`. There is **no separate
amplitude field** — the 4th component of `R` is an ambient coordinate, and the
3+1 (gauge/gravity) split is an inside-observer interpretation, not a stored
quantity (backbone #22).

### 1.2 The action (settled — `discrete_4d_brane_action.md` §2)
```
S[R] = Σ_l Δt ( T^{l+½} − V^l )

V^l        = (k_s/4) Σ_p Σ_{δ∈𝒩_s} ( |R_{p+δ}^l − R_p^l| − α a )²      (spacelike, 6 links)
T^{l+½}    = Σ_p (m/2) ( (R_p^{l+1} − R_p^l)/Δt )²                       (temporal, 2 links)
```
- `𝒩_s = {±ê_x, ±ê_y, ±ê_z}` (6-neighbor axial-only, backbone #15).
- Spacelike links enter with `−` (potential), temporal with `+` (kinetic): this
  sign **is** the Lorentzian signature (backbone #21, A5). The ambient is
  Euclidean/symmetric; the signature lives in the action.
- The temporal link is a central-force spring `½ k_t(|ΔR| − r_t)²` like the
  spacelike links; its `r_t = 0` limit is the zero-rest-length kinetic increment
  (quadratic, no norm-then-subtract). See DESIGN DECISION D3.

### 1.3 The solve (the core problem)
Stationarity `δS/δR_p^l = 0` at every **interior** node gives the discrete
d'Alembertian, which is **term-for-term the Störmer–Verlet stencil**
(`discrete_4d_brane_action.md` §3):
```
m (R_p^{l+1} − 2 R_p^l + R_p^{l−1}) / Δt²  =  F_p^l  :=  −∂V^l/∂R_p^l
F_p^l = k_s Σ_{δ∈𝒩_s} ( |R_{p+δ}^l − R_p^l| − α a ) · (R_{p+δ}^l − R_p^l)/|·|
```
Define the **residual** at each interior node:
```
𝓡_p^l  :=  m (R_p^{l+1} − 2 R_p^l + R_p^{l−1}) / Δt²  −  F_p^l   =  −∇_{R_p^l} S
```
The world-volume is the configuration with `𝓡 = 0` at all interior nodes, subject
to boundary conditions on slices `l=0` and `l=N`.

**Critical:** `S` is **Lorentzian → a saddle, unbounded below** (`T` enters `+`,
`V` enters `−`). You **must not minimize `S`** (gradient descent diverges along
the kinetic direction and/or silently solves the Euclidean problem). The solver
**root-finds `𝓡 = 0`** (equivalently minimizes `‖𝓡‖²`). This is non-negotiable
(`papers/core/derivations/status.md` A1).

### 1.4 Two solve modes
- **IVP (forward, special case):** boundary data = two adjacent past slices
  `(R⁰, R¹)`; march `R^{l+1} = 2R^l − R^{l−1} + (Δt²/m) F^l`. Well-posed Cauchy
  problem. **Used only for regression** (must reproduce Sprint-1 dispersion).
- **BVP (block, foundational):** boundary data = a past slice `l=0` and a future
  slice `l=N` (plus chosen conditions, D2); root-find all interior slices
  `l=1..N−1` simultaneously. This is the retrocausal/time-symmetric solve.

---

## 2. Open design decisions (resolve before finalizing the solver core)

These are the deferred items. Each is a real research/engineering fork; the
architecture is built so the choice is localized to the solver module.

- **D1 — Root-finder for `𝓡=0` (A1).** Candidates: (a) **Newton–Krylov**
  (JFNK) on the residual `𝓡` — recommended default; the Jacobian-vector product
  is cheap (the residual is a sparse 4D stencil). (b) Nonlinear least-squares on
  `‖𝓡‖²` (Gauss–Newton/L-BFGS on the *squared* residual — never on `S`). (c)
  Picard/relaxation as a warm-start. Recommended default: **JFNK with a
  forward-Verlet IVP warm start.**
- **D2 — Boundary conditions for the two-time BVP (A2). RESOLVED in the linear
  regime (D2 spike, 2026-05-30; numerically confirmed).** Because `D(k)` is
  diagonal, the block BVP decouples per spatial mode into a scalar recurrence
  `a^{l+1} − 2cosθ(k) a^l + a^{l−1} = 0`, `θ(k)=arccos(1−Δt²ω²(k)/2)`.
  - **Naive Dirichlet two-time (fix `R⁰` and `Rᴺ`) is NOT buildable:** the modal
    operator determinant is `2i·sin(Nθ(k))`, singular at `Nθ=mπ` (normal-mode
    resonances), so the block is generically ill-conditioned
    (`κ ~ 1/min_k|sin(Nθ(k))|`; ~most time-extents fail at realistic mode counts).
  - **The fix is two-past-slice Cauchy data (verdict a, 2026-05-31 — implemented,
    SUPERSEDES the characteristic-future-condition sketch).** The characteristic
    "kill `a₋` per mode" idea is correct 2×2 algebra but **wrong for a real
    field**: reality couples `a₊(k)=conj(a₋(−k))`, so zeroing `a₋` per mode
    deletes both characteristics → non-real garbage. The correct, well-posed,
    reality-respecting chiral BC is simply **two adjacent past slices `(R⁰, R¹)`
    marched** — forward = matter/retarded, backward (reverse stencil) = antimatter.
    One slice can't encode direction; two can. No FFT, no future condition, no DC
    special case; reality is automatic (real stencil + real data). Condition is
    bounded `O(10)`, N-independent (incl. at Dirichlet resonances). This *is* the
    matter/antimatter worldtube orientation. Implemented in
    `solver/{boundary,bvp}.py`; `ChiralBC(R0,R1,chirality)`, `solve_block` routes
    it to `ivp.march` (no JFNK). Tested: recovers an exact forward eigenmode to
    8e-14 at a resonant N where Dirichlet is `κ≈1e14`.
  - **MAJOR IMPLICATION (reshapes step v / §14):** the well-posed chiral solve =
    Cauchy march ≈ IVP for free waves. So the **two-time Dirichlet BVP is NOT the
    vehicle for a bound particle.** A particle at rest is a **time-periodic
    (carrier+envelope, #18) + spatially-localized** worldtube → the particle-search
    block solver must be a **nonlinear eigen-BVP: periodic-in-time
    (`R^{l+P}=R^l`, period/frequency the eigenparameter) + localized-in-space**,
    not two-time Dirichlet and not IVP. Nonlinear global uniqueness stays open
    (`papers/core/derivations/status.md` A2).
- **D3 — Temporal-link form (A4, `discrete_4d_brane_action.md` §6).** The temporal
  link is a central-force spring with rest length `r_t`, the same law as the
  spacelike links — **one 4D-isotropic spring parameterized by `r_t`, not a fork**.
  `r_t = 0` is the linear/Verlet limit (zero-rest-length kinetic = plain Newton,
  matches the validated dispersion — the **default**); `r_t = α·β·dt` is the
  prestressed substrate, adding the time-link geometric quartic (backbone #22's
  time-dilation face). IMPLEMENTED 2026-06-05: `ActionParams.r_t`; the residual
  routes on `r_t` (0 → Verlet stencil, >0 → spring force), bit-identical at `r_t→0`.
  Newtonian-limit/isotropy of `α_t=α` settled via the gravity/contraction channel.
- **D4 — Solver state size.** A full block at `N` slices × `N_x³` nodes × 4
  components × float64 is large. DESIGN DECISION: matrix-free residual (never
  assemble the Jacobian); slab/streaming storage; optionally solve in a reduced
  (narrowband-envelope) subspace. Keep the residual operator the single source of
  truth so memory strategy is swappable.

---

## 3. Physics core (settled specification)

### 3.1 Lattice topology
- **Spacelike:** 6-neighbor axial-only on `(i,j,k)`. Offsets
  `{(±1,0,0),(0,±1,0),(0,0,±1)}`, all `|δ|=1`, equal weights (`axial_weight`,
  default 1). Per-axis periodic or open boundaries. Diagonal shells are
  **intentionally absent** — this is what makes the spatial dynamical matrix
  `D(k)` diagonal in the Cartesian basis at every `k`, `α` (backbone #15;
  Sprint-2 #9 certified).
- **Temporal:** 2-neighbor `{l±1}`, one link family, opposite sign.
- **Dimension-agnostic:** the spacelike stencil must generalize to 1D/2D/3D
  (`2·dim` axial neighbors) with **no hard-coded `3`** in core/forces/solver/
  diagnostics (principles §6). The "4" in `ℝ⁴` (ambient) and the timelike axis
  are separate from the spacelike `dim`.

### 3.2 Energy and residual
- Spacelike force / potential: central-force pair springs, exactly as §1.2/§1.3.
  The energy stays **bounded** as a link collapses (`→ ½k_s(αa)²`) — the StVK /
  central-force law is **non-coercive in compression** (`papers/matter_mass/derivations/matter/status.md` C1).
  Do **not** silently fix this; if a coercive (polyconvex) law is later adopted,
  it is a deliberate constitutive change with its own derivation.
- The residual operator `𝓡` (§1.3) is the **single shared primitive**: the solver
  drives it to zero; diagnostics report `‖𝓡‖`. Implement once.

### 3.3 Conventions & units (must match canon)
- `α := rest_length / spacing`; `α=1` ↔ no prestress, `α=0` ↔ max prestress;
  **default α=0.2**.
- Dimensionless units `k_s = a = ρ = 1` (`m = ρ a^{dim}`).
- Light-cone: `c_L² = k_s a²/m`, `c_T² = (1−α) k_s a²/m`. **Regression target:**
  `c_T/c_L = √(1−α)`; at α=0.2, `0.8944 ± <5%` at `|k|a=0.1` (Sprint-1 hit ~1e-5).

### 3.4 Non-negotiables (principles §; enforce in review)
- **No emergent back-reaction:** diagnostics/measurements never feed forces.
- **No clamps/cutoffs/hand thresholds:** any threshold must emerge from dynamics.
- **Geometric coupling only:** lateral↔amplitude coupling is via the induced
  metric (the 4-vector norm), never an added field.
- **English only.**

---

## 4. Component architecture

Four components, **file-mediated only** (separate processes, communicate via
files — no shared memory). This separation from the legacy design is kept; the
*simulation* component is replaced by a *solver* component.

```
[1] initialization  → boundary_problem.npz   (held lattice, params, slice data)
                          ↓
[2] solver          → worldvolume.zip         (solved 4D block; or IVP march)
                          ↓
        ┌─────────────────┴─────────────────┐
[3] diagnostics → diagnostics/*.{json,csv,npz}   [4] visualization → plots/*.mp4
```

- **[1] initialization** — generates **boundary data**: the held reference
  lattice, parameters (`α`, `k_s`, `m`, `Δt`, grid, `N` slices), and the
  prescribed slice configuration(s). For IVP mode: two past slices. For BVP mode:
  past slice `l=0` + future slice `l=N` (+ chirality condition, D2). Baryon/wave
  *seeds* live here (§10).
- **[2] solver** — the core. Reads `boundary_problem.npz`, runs IVP march or BVP
  root-find (D1), writes the full world-volume. Reports solver telemetry
  (residual history, iteration count, condition estimate).
- **[3] diagnostics** — reads the world-volume, computes metrics (§9). Read-only.
- **[4] visualization** — renders slices / world-volume to video.

---

## 5. State & data model

- **World-volume** = ordered stack of spacelike slices `R^l`, `l=0..N`, each
  `(N_x·N_y·N_z, 4)` float64. The legacy "trajectory" *is* a world-volume; for
  BVP it is the solved interior, for IVP it is the march.
- **Boundary data** = the subset of slices/components that are *fixed* during the
  solve, plus a mask identifying them. IVP: `{R⁰, R¹}` fixed. BVP: `{R⁰, R^N}`
  fixed (interior free).
- **Held reference** = the unstressed lattice `a·(i,j,k)` (+ held timelike
  spacing `βΔt`); displacements `u = R − R_ref` are derived, not stored
  separately.

---

## 6. Data contracts (file formats)

Carry over the proven, versioned formats (legacy `io.py`,
`FORMAT_VERSION` bumped to `branesim-block-v1`). All JSON blobs stored as
length-1 string arrays inside `.npz`; `allow_pickle=False` everywhere.

### 6.1 `boundary_problem.npz`
| key | type | meaning |
|---|---|---|
| `format_version` | str[1] | `branesim-block-v1` |
| `ref_positions` | f64 `(Nnodes,4)` | held reference lattice |
| `boundary_slices` | f64 `(Nb, Nnodes, 4)` | prescribed slice configs |
| `boundary_indices` | i64 `(Nb,)` | which `l` each boundary slice pins (e.g. `[0,1]` IVP, `[0,N]` BVP) |
| `boundary_mask_json` | str[1] | which components/nodes are fixed (for partial BCs / chirality) |
| `lattice_json` | str[1] | `{grid_shape, spacing, periodic_axes, axial_weight, dim}` |
| `action_json` | str[1] | `{k_s, m (or ρ), alpha, dt, n_slices N, r_t (0 = linear limit; α·β·dt = prestressed), k_t}` |
| `seed_json` | str[1] | seed/ansatz metadata (J, L, B_winding, …) |
| `metadata_json` | str[1] | provenance |

### 6.2 `worldvolume.zip`
- `manifest.json`: `{format_version, mode:"ivp"|"bvp", lattice, action, slices:[{index l, time, name}], solver_report}`
- `slices/slice_{l:06d}.npz`: `{positions (Nnodes,4)}` (velocities are derived as
  temporal differences; store only if model needs them).
- `solver_report`: `{residual_initial, residual_final, iterations, converged,
  condition_estimate, walltime_s}`.
- `aux/ref_positions.npy`.

### 6.3 diagnostics outputs
`diagnostics/metrics.csv` (per-slice rows), `diagnostics/summary.json`,
`diagnostics/*.npz` (e.g. berry series). Unchanged in spirit from legacy.

---

## 7. Proposed module layout

```
branesim/
  core/                      # dimension-agnostic, no I/O, the physics primitives
    lattice.py               #   4D lattice topology: spacelike 2·dim axial + temporal stencil
    action.py                #   V^l, T^{l+½}, S; the spacelike force F (central-force springs)
    residual.py              #   𝓡 = m·∂_τ²R − F  (the shared primitive; matrix-free)
    conventions.py           #   α, units, c_L/c_T helpers
  solver/                    # the foundational core (the new part)
    bvp.py                   #   block root-find of 𝓡=0  (D1: JFNK)  + boundary application (D2)
    ivp.py                   #   forward Verlet march (special case; regression)
    boundary.py              #   chirality / two-time BC handling (D2)
  initialization/            # boundary-data generators
    seeds/ (plane_wave, hedgehog, skyrme_twisted, trace_hedgehog, axis_triplet, spherical)
    build.py                 #   assemble boundary_problem.npz
  diagnostics/
    dispersion.py  berry.py  christoffel.py  confinement.py  block_solver.py(residual/condition/chirality)
  visualization/
    slice.py  volume.py
  io/
    contracts.py             #   §6 readers/writers, FORMAT_VERSION
  orchestration/
    (run via branesim/run_experiment.py — init+solve→worldvolume; diagnostics via
    #                              branesim/diagnostics/run_measurements.py)
    configs/*.json
  tests/                     # see §12 (this time: real coverage)
```

Key signatures (illustrative, dimension-agnostic):
```python
# core/residual.py
def residual(world: Array["L+1,N,4"], ref: Array, p: ActionParams) -> Array["L+1,N,4"]:
    """𝓡 at interior nodes; 0 on boundary slices. = −∇S. Matrix-free."""

# solver/bvp.py
def solve_block(problem: BoundaryProblem, opts: SolveOpts) -> WorldVolume:
    """Root-find residual()=0 over interior slices (JFNK). NEVER minimizes S."""

# solver/ivp.py
def march(problem: BoundaryProblem, opts: SolveOpts) -> WorldVolume:
    """Verlet special case: R^{l+1} = 2R^l − R^{l−1} + (Δt²/m) F(R^l)."""
```

---

## 8. Configuration schema (orchestration JSON)
```json
{
  "initialization": {"module": "...seeds.skyrme_twisted", "n_slices": 200,
                     "grid_size": 60, "spacing": 1.0, "alpha": 0.2,
                     "periodic_axes": [true,true,true], "seed_params": {...}},
  "solver": {"mode": "bvp", "r_t": 0.0,
             "method": "jfnk", "tol_residual": 1e-8, "max_iter": 200,
             "warm_start": "ivp", "dt": 0.05},
  "diagnostics": {"frame_stride": 5, "confinement_radius_factor": 0.5, ...},
  "visualizations": [{"mode": "slice", ...}]
}
```

---

## 9. Diagnostics catalog
**Carry-over (validated/used):**
- `dispersion` — plane-wave `ω(k)` fit; the IVP-mode regression vs `c_T/c_L=√(1−α)`.
- `christoffel` — closed-form `D(k)` reference (diagonal Cartesian eigenframe).
- `berry` — complex-envelope geometric phase/connection (per-wavepacket, backbone #18).
- `confinement` — `spread_ratio = radius_rms/box_fill_radius`, `confined_fraction`,
  `radius_growth`. **Do NOT use the legacy self-referential `leakage_fraction`
  for go/no-go** (it was misleading; see project history).

**New (block-specific):**
- `residual_norm` — `‖𝓡‖` over the world-volume (solve quality).
- `bvp_condition` — operator condition estimate vs `N` (resonance/well-posedness, D2).
- `chirality` — forward-vs-backward worldtube measure (causal-direction selection).

---

## 10. Initializers / boundary-data generators
Each declares its `(J, P; SU(3)-irrep; Q_U(1); B_winding)` label in metadata
(backbone #20) so diagnostics can project the right channel. Carry over:
`plane_wave` (dispersion), and the VSH menu `hedgehog (C1)`, `skyrme_twisted (C2)`,
`trace_hedgehog (C3)`, `axis_triplet (C4, negative control)`, plus generic
spherical-harmonic. **Known result to honor:** small-amplitude wide hedgehogs
**disperse** (no wide soliton under the central-force law); seeds are for the
solver/diagnostic plumbing and lattice-scale or alternative-stabilizer studies,
not a presumed wide baryon.

---

## 11. Validation strategy (this time: real coverage)

The legacy code had ~1 unit test; that is a primary risk. The rewrite must ship
with:
1. **Regression (the gold check):** IVP mode reproduces Sprint-1 dispersion —
   `c_L=1`, `c_T/c_L=√(1−α)` within 5% at `|k|a=0.1`, and `D(k)` diagonal in the
   Cartesian basis (Sprint-2 #9). If the rewrite can't reproduce these, it's wrong.
2. **Solver consistency:** on an IVP problem, `solve_block` (BVP with `l=0,1`
   pinned and a long extent) must agree with `march` to solver tolerance.
3. **Residual=0 ⇔ Verlet:** a Verlet-marched world-volume has `‖𝓡‖ ≈ 0`
   (machine precision) at interior nodes — ties the two cores together.
4. **Energy/symplecticity:** IVP energy drift bounded over long runs (Verlet is
   symplectic); the nonlinear regime (untested in legacy) gets explicit checks.
5. **Action sign sanity:** confirm a descent on `S` *diverges* (saddle), and the
   root-finder on `𝓡` converges — guards against accidentally solving the
   Euclidean problem.

---

## 12. What carries over vs what is genuinely new

| Element | Carry over (as spec) | New |
|---|---|---|
| Spacelike energy / force | ✅ exact | — |
| 6-nn axial stencil, conventions, units | ✅ exact | extend to temporal links + dim-agnostic |
| File-mediated 4-component pipeline | ✅ pattern | solver replaces simulation |
| I/O contracts | ✅ adapted | world-volume / boundary-problem schemas |
| Diagnostics (dispersion/berry/christoffel/confinement) | ✅ | residual/condition/chirality |
| Seeds (plane wave, VSH menu) | ✅ | recast as boundary-data generators |
| Forward Verlet | ✅ as IVP special case / regression | — |
| **Block BVP solver (`𝓡=0` root-find)** | — | ✅ the core (D1) |
| **Two-time BCs + chirality selection** | — | ✅ (D2, open) |
| Test suite | ❌ (was absent) | ✅ §11 |

---

## 13. Risks & open problems (tracked in the per-paper bridge `status.md` files)
- **A1 (root-find, not minimize):** non-negotiable; mis-implementing it silently
  solves the wrong (Euclidean) problem. Guard with test §11.5.
- **A2 (BVP well-posedness / chirality):** the core is **resolved in the linear
  regime** (D2: two-past-slice chiral BC; tested to 8e-14 eigenmode recovery,
  verdict-a 2026-05-31 — see `papers/core/derivations/status.md` A2 and `LESSONS_LEARNED.md`).
  Remaining open: the chiral march in `solver/boundary.py` currently uses the
  `r_t=0` Verlet stencil, so `r_t>0` chiral solves need an `r_t`-aware step; and
  nonlinear global uniqueness stays open.
- **A4 (temporal link):** one 4D-isotropic spring parameterized by `r_t` (0 = linear
  limit, the default; `α·β·dt` = prestressed) — implemented; `α_t=α` Newtonian-limit
  /isotropy settled via the gravity channel.
- **C1 (StVK non-coercivity):** the energy doesn't penalize collapse; relevant if
  lattice-scale/large-amplitude solitons are pursued — possible polyconvex
  completion is a deliberate, separate change.
- **Scale/memory (D4):** a full 4D block is large; keep the residual matrix-free
  and storage swappable.

---

## 14. Program strategy — "all ingredients at once" (owner, 2026-05-31)

A stable particle is **not** expected from any single ingredient: the
elasticity-only IVP hedgehog provably disperses (`papers/matter_mass/derivations/matter/status.md` C1,
`[[project-c2-skyrme-no-confinement]]`). The particle is a stationary point of
the **full 4D action** that requires, *together*: **BVP not IVP** (a whole
time-coherent worldtube, not a radiating forward march); **colour from
axis-alignment** (SU(3) on the lateral triplet, #16/#19); **Berry/Wilczek–Zee
diagnostics** to see the gauge/binding holonomy; and **time-symmetry /
chirality** (which is also what makes the BVP well-posed — D2). The target is
likely a **time-periodic (carrier+envelope, #18) worldtube**, not a static blob.

So the search is **holistic**: compose the ingredients into one block-solve and
look for the particle there — do not expect it rung-by-rung. **Engineering
caveat:** each ingredient must still be **unit-tested correct in isolation**
before composition (an earlier broken chiral BC was the cautionary tale — a wrong
ingredient yields a residual-zero *garbage* solution that looks converged and is
undebuggable inside a combined run). "Verify the parts, then run them together."

Minimal coherent ingredient set for a first holistic attempt: working chiral BVP
+ a topological+colour seed (hedgehog/Skyrme with axis-colour) + Berry/WZ +
confinement diagnostics, as one solve.

## 15. Deployment

Runs target **AWS, high-memory CPU** (numpy/scipy; no GPU rewrite — backend
decision 2026-05-31). The block solve is memory-bound. Entry point:
`branesim/run_experiment.py` (`branesim-run`), config-driven, writes
`worldvolume.zip` + `summary.json`; cost-safe EC2 scaffolding in
`orchestration/aws/`. **See `DEPLOYMENT.md`** for the memory-sizing formula,
instance table, and the S3 launch/fetch flow. Packaging via `pyproject.toml`
(`pip install -e .`).

---

*Recommended sequence: (i) D2 BVP well-posedness — RESOLVED (linear); (ii) `core/`
+ `solver/ivp.py` + regression — DONE; (iii) `solver/bvp.py` JFNK + Dirichlet —
DONE (chiral BC resolved in the linear regime; `r_t>0` march still needs an `r_t`-aware step); (iv) port diagnostics + colour/Berry + seeds;
(v) compose the holistic particle-search block-solve (§14); deploy on AWS (§15).*
