# The explicit discrete 4D brane action (kinetic-temporal-link model)

Author: `physics-derivation` agent.
Companion to `lattice_to_continuum.md` (spatial-link algebra) and
`status.md` *Open derivations* (group A — solver caveats this derivation points to but does not
resolve).

> **What this formalizes.** In the working Verlet pipeline the timelike
> direction is *implicit*: it appears only as the kinetic term `½ m v²`. This
> document makes the timelike direction an **explicit lattice link** in the
> 4th direction, so the brane action becomes a genuine 4D-cubic link sum with
> `6 + 2` links per node, and the Lorentzian signature (backbone #21) is a
> single sign choice on the temporal link family. The local Euler–Lagrange
> stencil is shown to be *identical* to Störmer–Verlet — Verlet **is** the
> discrete variational integrator of this action — so nothing changes in the
> force code; only the global solution philosophy (forward IVP vs 4D block
> BVP) differs.

Conventions match `PRINCIPLES.md` §1.1/§1.1a, `BACKBONE.md` #15/#21/#22, and
`lattice_to_continuum.md` (6-neighbor axial-only spatial stencil,
`α := rest_length/spacing`, default `α = 0.2`).

---

## 1. Assumptions

- **A1.** 4D-in-4D ontology: each node carries `R_p^l ∈ ℝ⁴`; no separate
  amplitude field. `p = (i,j,k)` is the spacelike triple, `l ∈ ℤ` the timelike
  index. (backbone #1)
- **A2.** Reference = held lattice
  `R_p^l = a(i,j,k,0) + β Δt (0,0,0,l) + u_p^l`, `u` small.
- **A3.** 6-neighbor axial-only spacelike stencil
  `𝒩_s = {±ê_x, ±ê_y, ±ê_z}`; no diagonal shells ⇒ spatial `D(k)` diagonal in
  the Cartesian basis at every `α`. (backbone #15)
- **A4.** 2-neighbor temporal stencil `𝒩_t = {l±1}` (one timelike link family).
- **A5 (signature).** Spacelike links enter `S` with `−` (potential `V`);
  the temporal link enters with `+` (kinetic `T`). This *is* the Lorentzian
  signature — a property of the action, not the ambient. (backbone #21)
- **A6.** Long wavelength / slow time: `a k_wave ≪ 1`, `Δt ω ≪ 1`.
- **A7.** Linear elasticity for the dispersion result (geometric quartic /
  Skyrme term, backbone #17, is next order).
- **A8.** Inversion-symmetric stencil (each `δ` paired with `−δ`, equal weights).
- **A9 (constitutive).** Spacelike: central-force pair springs
  `U_link = ½ k_s (|R_{p+δ}^l − R_p^l| − α a)²`. Temporal: the same central-force
  spring `½ k_t (|R_p^{l+1} − R_p^l| − r_t)²` (§6); its `r_t = 0` linear limit is
  the zero-rest-length kinetic increment `½ m ((R_p^{l+1} − R_p^l)/Δt)²`.

---

## 2. The explicit 4D action

Spacelike potential on time-slice `l` (the inner ½ removes bond
double-counting):

$$
V^l = \frac{k_s}{4}\sum_{p}\sum_{\delta\in\mathcal N_s}\Big(\big|R_{p+\delta}^l - R_p^l\big| - \alpha a\Big)^2 .
$$

Temporal kinetic term on the link `l → l+1` (`(·)² ` = full 4D dot product):

$$
T^{l+\frac12} = \sum_{p}\frac{m}{2}\left(\frac{R_p^{l+1}-R_p^l}{\Delta t}\right)^2 .
$$

Discrete action (A5 signature):

$$
\boxed{\;
S[R] = \sum_l \Delta t\Big(T^{l+\frac12} - V^l\Big)
= \sum_l \Delta t\!\left[\sum_p \frac{m}{2}\!\left(\frac{R_p^{l+1}-R_p^l}{\Delta t}\right)^{\!2}
- \frac{k_s}{4}\sum_p\sum_{\delta\in\mathcal N_s}\!\big(|R_{p+\delta}^l-R_p^l|-\alpha a\big)^2\right].
\;}
$$

**Why the temporal link is "different math" from the spatial links:**

| | Spacelike links (`𝒩_s`, 6) | Temporal link (`𝒩_t`, 2) |
|---|---|---|
| Form | `(|ΔR| − αa)²` | `|ΔR|²` (no norm-then-subtract) |
| Rest length | `αa > 0` (prestress) | `0` |
| Nonlinearity | geometric (norm of 4-vector) | none — quadratic in raw increment |
| Sign in `S` | `−` (potential) | `+` (kinetic) |
| Stiffness | `k_s` | `m/Δt²` |

Because the temporal rest length is zero, `|ΔR|²` carries no `(1−α)`
geometric-stiffness piece (contrast `lattice_to_continuum.md` §3.1): it is
purely Hookean in the increment. The Lorentzian signature is just its opposite
sign relative to the six spacelike springs.

---

## 3. Stationarity = discrete d'Alembertian = Verlet stencil

Vary `S` w.r.t. an interior node `R_p^l` (it appears in `T^{l−½}`, `T^{l+½}`,
and `V^l`). The temporal contribution is

$$
\frac{\partial}{\partial R_p^l}\,\Delta t\big(T^{l+\frac12}+T^{l-\frac12}\big)
= -\frac{m}{\Delta t}\big(R_p^{l+1}-2R_p^l+R_p^{l-1}\big),
$$

and the spacelike contribution is the substrate force
`F_p^l := −∂V^l/∂R_p^l` (the principles §1.2 force `F = −∂U/∂R`):

$$
F_p^l = k_s\sum_{\delta\in\mathcal N_s}\big(|R_{p+\delta}^l-R_p^l|-\alpha a\big)\,\widehat{(R_{p+\delta}^l-R_p^l)} .
$$

Setting `δS/δR_p^l = 0`:

$$
\boxed{\;
m\,\frac{R_p^{l+1}-2R_p^l+R_p^{l-1}}{\Delta t^2} = F_p^l = -\frac{\partial V^l}{\partial R_p^l}.
\;}
$$

This is the **discrete d'Alembert/wave operator** (discrete `∂_t²` equated to
the spacelike elastic operator) and is **term-for-term identical to
Störmer–Verlet**:

$$
R_p^{l+1} = 2R_p^l - R_p^{l-1} + \frac{\Delta t^2}{m}\,F_p^l .
$$

This is the discrete-mechanics / variational-integrator statement (Marsden–West):
the Euler–Lagrange equations of the *discrete* action `S = Σ_l Δt(T−V)` are
exactly the Störmer–Verlet update, and the flow is symplectic. **Verlet is the
discrete variational (symplectic) integrator of this action.**

**Consequence (load-bearing).** The local stencil is unchanged whether the
model is read as a forward IVP or a 4D block BVP — both demand
`δS/δR_p^l = 0` at every interior node. Only the global solution philosophy
differs:

- **IVP / forward** (current pipeline): prescribe `(R⁰, R¹)`, march `l → l+1`.
  Well-posed Cauchy problem.
- **BVP / block** (foundational): prescribe spacelike-slice data on a *past*
  slice `l=0` and a *future* slice `l=N`, root-find the interior. This is the
  retrocausal / time-symmetric reading (principles §1.5, `status.md` A1/A2).

---

## 4. Continuum limit and light-cone speeds

With `R_p^l = R(an, lΔt)` smooth, the temporal difference gives `∂_τ²R` at
`O(Δt²)`, and the spacelike side reproduces the StVK acoustic operator of
`lattice_to_continuum.md` §3.4 restricted to shell I. With `ρ = m/a³`:

$$
\boxed{\;\rho\,\partial_\tau^2 R = \partial_b\,A_{abcd}\,\partial_d R_c,\qquad
A_{abcd} = \frac{k_s}{a}\big[\alpha\,Q_{abcd} + (1-\alpha)\,\delta_{bd}\,\delta_{ac}\big].\;}
$$

For `k = k ê_x` the Christoffel matrix is diagonal (A3), giving the lattice
light-cone speeds:

$$
\boxed{\;
c_L^2 = \frac{k_s\,a^2}{m},
\qquad
c_T^2 = (1-\alpha)\,\frac{k_s\,a^2}{m}.
\;}
$$

In dimensionless units `k_s = a = ρ = 1` (`m=1`) these reduce to the canonical
6-neighbor values `c_L = 1`, `c_T = √(1−α)` (principles §1.1a). At `α = 0.2`:
`c_T/c_L = √0.8 ≈ 0.894`.

**Tuning conditions for emergent isotropic Lorentz invariance** (stated, not
claimed — `status.md` A3). The discrete `□ = ∂_τ² − c²∇²` is
isotropic-Lorentzian only when:

1. **Cone-matching (CFL-like):** the temporal-to-spatial stiffness ratio
   `m/Δt²` vs `k_s` is tuned so the discrete cone matches the elastic cone;
   at the matched timestep `Δt = a/c_L` the temporal `Δt²/12 ∂_τ⁴` and spatial
   `O((ak)⁴)` truncations cancel to leading order ("magic timestep").
2. **Long wavelength:** `a k_wave ≪ 1`, so `O((ak)⁴)` anisotropy is below
   tolerance.
3. **L–T degeneracy:** `c_L = c_T` requires `α → 0`. At `α = 0.2` there is a
   residual `~10.6%` L–T gap along `[100]`, `0%` along `[111]` — **leading-order
   cubic anisotropy, not `O((ak)²)`**. Hence isotropic Lorentz invariance is
   *conjectured at the inside-observer level* (backbone #8), not present at the
   lab level. This derivation supplies the explicit gap to be closed by the
   dual-observer renormalization; it does not itself establish it.

---

## 5. Unboundedness and BVP ill-posedness (pointers, not re-derived)

`S = Σ_l Δt(T−V)` is Lorentzian (`T` enters `+`, `V` enters `−`) and therefore
**unbounded below — a saddle, not a minimum**. Consequences (full statements in
`status.md` *Open derivations*, group A):

- **A1.** The foundational solver must **root-find `∇S = 0`** (Newton–Krylov on
  the discrete d'Alembertian, or minimize `‖∇S‖²`), never gradient-descend `S`.
  Forward Verlet is the well-posed IVP special case.
- **A2.** The two-time (past+future) BVP is **not unconditionally well-posed**:
  at time-extents `NΔt` hitting a spatial normal-mode resonance the homogeneous
  problem has nontrivial solutions ⇒ non-uniqueness. The conjectured selection
  principle (soliton chirality, principles §1.5) is left open.

---

## 6. The temporal central-force spring and its linear limit

**The temporal link is a central-force spring** along time (own rest length
`r_t = α_t β Δt`, prestress), Lorentzian sign — the same law as the spacelike
links, so the substrate is one 4D-isotropic spring lattice. The `r_t = 0` limit
recovers the zero-rest-length kinetic increment of §2 (the linear/Verlet limit):

$$
T^{l+\frac12}_{\mathrm{time}} = \sum_p \frac{k_t}{2}\Big(\big|R_p^{l+1}-R_p^l\big| - \alpha_t\,\beta\,\Delta t\Big)^2 .
$$

All 8 links are then the same kind of central-force spring, differing only by
the sign of the time family. Its Euler–Lagrange time term is **not** plain
Newton (with `ΔR^± ≡ R_p^{l±1} − R_p^l`):

$$
-k_t\!\left[\big(|\Delta R^+|-\alpha_t\beta\Delta t\big)\widehat{\Delta R^+} + \big(|\Delta R^-|-\alpha_t\beta\Delta t\big)\widehat{\Delta R^-}\right] + \Delta t\,F_p^l = 0 .
$$

The bracket is a **nonlinear time-link force**.

**Exact reduction to baseline.** Write the temporal rest length as
`r_t := α_t β Δt`. Comparing the two link *energies* at matched stiffness
`k_t = m/Δt²` gives a single closed-form difference:

$$
E(r_t) - E(0) = \tfrac12 k_t\big(|\Delta R| - r_t\big)^2 - \tfrac12 k_t|\Delta R|^2
= -\,k_t\,r_t\,|\Delta R| + \tfrac12 k_t r_t^2 .
$$

The entire difference is the term `−k_t r_t |ΔR|` (plus an irrelevant additive
constant). It is **exactly proportional to `r_t`**, so the temporal spring and
its `r_t = 0` (kinetic) limit coincide
**identically (not approximately) iff `r_t = 0`** — i.e. `α_t = 0` (maximal
temporal prestress), with `k_t = m/Δt²`. No condition on `β Δt` is needed: the
identity `|ΔR|\widehat{ΔR} = ΔR` is exact, so at `r_t = 0` the EL force collapses
to `−k_t(R^{l+1} − 2R^l + R^{l−1})`, plain Newton. (The reduction needs no
`β Δt → 0` condition.)

**What `r_t` controls.** The distinguishing term `−k_t r_t |ΔR|` involves the
Euclidean *norm* `|ΔR|`, not the squared norm — the irreducibly geometric piece.
Expanding it generates both the non-Newtonian force above and a **geometric
quartic** (and higher norm-nonlinearities) in the time link whose coefficient is
`∝ r_t` (at `r_t = 0` there is no square root left to expand). That geometric
quartic is the time-link analogue of the Skyrme/contraction term (#17) and is
the natural home for backbone #22's unified contraction — gravitational time
dilation sourced by spacelike/kinetic displacement stretching the timelike link.
**Note:** the distinction is *not* transverse stiffness — `½ k_t |ΔR|²` already
has isotropic Hessian, so the `r_t = 0` limit's time link does respond to
transverse changes (that response is literally the kinetic energy). The genuine
distinction is the norm-nonlinearity / geometric quartic, present iff `r_t ≠ 0`.

So `r_t` dials between the **linear/Verlet limit** ("matches plain Newton and the
existing Verlet code", `r_t = 0`) and the **canonical prestressed substrate** ("a
non-Newtonian time link that carries gravity's second face", `r_t = α·β·dt`); the
two limits coincide only at `r_t = 0`. This is **one model parameterized by `r_t`** —
resolved and implemented (`ActionParams.r_t`; `status.md` A4); the `α_t = α`
consistency is the open verification (A4a).

---

## 7. Regime of validity

- `a k_wave ≪ 1`, `Δt ω ≪ 1`; neglected: `O((ak)⁴)` spatial, `Δt²/12 ∂_τ⁴R`
  temporal.
- Linear elasticity (`|∇u| ≪ 1`); geometric quartic / Skyrme excluded.
- `c_T² ≥ 0` requires `α ≤ 1` (at `α = 1` transverse modes are zero-frequency).
- Isotropic-Lorentzian behavior only in the tuned window of §4 — **not claimed**.

---

## 8. Falsifiable numerical prediction

**Light-cone L–T speed ratio from the temporal-link stiffness.** With
`k_s = a = ρ = 1` (`m=1`), independent of IVP-vs-BVP:

$$
\boxed{\;\frac{c_T}{c_L} = \sqrt{1-\alpha}\;}\quad\Rightarrow\quad
\frac{c_T}{c_L}\Big|_{\alpha=0.2} = \sqrt{0.8} = 0.8944,\quad c_L = 1.
$$

**Procedure** (dispersion-analyst layer): periodic `32³` spacelike lattice,
seed `u_a(x,0) = ε p̂_a cos(k·x)`, `ε = 10⁻³ a`, `k = k ê_x`,
`|k|a ∈ {0.1, 0.2, 0.4}`; run `p̂ = ê_x` (L) and `p̂ = ê_y` (T) with the §3
Verlet stencil; fit `ω(k)`, extrapolate `c²(k) → c²(0)`. **Stronger test of the
explicit-temporal-link claim:** vary `m/Δt²` (change `Δt` at fixed `m, k_s`) and
confirm both `c_L, c_T` are **independent of `Δt`** in the continuum
extrapolation (they depend only on `k_s a²/m`), while the *discrete* cone
`cΔt/a` scales linearly with `Δt`. **Failure thresholds:** `|c_L − 1| > 5%` or
`|c_T/c_L − √(1−α)| > 5%` at `|k|a = 0.1` falsifies the derivation or the force
code; `Δt`-dependence after `k → 0` falsifies the zero-rest-length temporal-link
assumption (A9) or the integrator's variational property.

---

## 9. What is derived vs. posited

1. **Derived** (from A1–A9): the action §2, the EL = Verlet identity §3, the
   continuum operator and `c_L² = k_s a²/m`, `c_T² = (1−α)k_s a²/m` §4.
2. **Posited (A5):** the Lorentzian sign. The ambient is symmetric; the
   signature is an input (backbone #21).
3. **Posited (A2):** the held timelike spacing `β Δt`; whether `β` is physical
   or a discretization gauge is not fixed in the `r_t = 0` limit (the time link
   has zero rest length). It becomes physical when `r_t > 0`.
4. **Not derived:** emergent isotropic Lorentz invariance (§4 gives only the
   tuning condition); inside-observer blindness to the `√(1−α)` anisotropy is
   the load-bearing conjecture (backbone #8, `status.md` A3).
5. **Pointed to, not derived:** unboundedness / BVP well-posedness (§5,
   `status.md` A1/A2).
6. **Resolved + implemented:** the temporal link is a central-force spring
   parameterized by `r_t` (§6); `r_t = 0` is the linear/Verlet limit. The open
   piece is the `α_t = α` consistency (A4a), not the link form.
7. **Constitutive dependence:** the EL = Verlet identity (§3) and unboundedness
   (§5) are constitutive-law-independent (survive any isotropic hyperelastic
   `V`); the specific `c_L, c_T` values (§4) are StVK/central-force-specific.
