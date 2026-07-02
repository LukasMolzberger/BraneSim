# Field-Strength Tensors from the Spring Stencil (Paper VII) — Status

Scaffolding stage, second attempt. Statuses marked `closed` / `in-ansatz` / `open` per project
convention. The defining upgrade over the first attempt: the complex carrier is *earned* from phase
space, and `α` is handled consistently — one parameter activating the colour sector in one direction.

## HAVE (imported / established — close as restatement)

- **8-link 4D stencil.** 6 spacelike axial + 2 temporal axial links per node; central-force springs;
  single anharmonicity `−k_sαa|ΔR|` `∝α`. From Papers I, III. See `link_holonomy.md`.
- **BZ curvature vanishes ∀α.** `D(k)` real symmetric ⇒ `k`-space Berry/WZ curvature `≡0`. The field
  strength lives over physical `(x,t)`, not the BZ. From Paper III (BACKBONE #16). Restated from the
  plaquette side in `berry_reconciliation.md`.
- **Carrier phase = time rotation.** The complex `U(1)` "`i`" is carrier rotation along the timelike
  link (`[[project_complex_u1_from_time]]`); two time slices minimal. *Why* the temporal links are
  mandatory for `F_{0i}` (the `E`-field).
- **Sector split.** `u(3)=u(1)⊕su(3)`: trace = EM, traceless = colour (Paper III). Inherited as the
  split of the single plaquette into the two field tensors.

## TARGET (this paper's new content)

- **`carrier_construction.md`** — `ψ^a=δR^a+iδṘ^a/ω` from phase space; `i` grounded as the time-link
  rotation; real degeneracy `O(3)` ⇒ complex carrier `U(3)`; single-vector `U(1)` phase vs transported
  three-frame `U(3)` connection (no `ℂ³`-vector / `u(3)`-matrix conflation); dynamic phase = raw
  material, not yet curvature. **The load-bearing new content.**
- **`link_holonomy.md`** — `U_μ(n)` from the carrier eigenframe; gauge law; plaquette `□_μν`;
  `log □_μν=ia²F_μν+O(a⁴)`; antisymmetry from orientation reversal.
- **`field_tensors.md`** — both sectors from the one `U(3)` plaquette: Faraday (trace `U(1)`, `E`/`B`,
  homogeneous Maxwell = exact Bianchi) and Yang–Mills (`G^a_μν`, `f^{abc}` from non-commuting links,
  non-Abelian Bianchi); **the single-direction `α`-activation of the colour sector** (`U(1)³` at `α=0`
  → full `U(3)` as `α→1`; weight trace `∝(3−2α)`, traceless `∝α`) — one monotonic activation, no tension.
- **`berry_reconciliation.md`** — continuum limit of the `(x,t)` plaquette holonomy `=` Berry / WZ
  curvature of Paper III; the two constructions are one object; BZ-vanishing restated.

## OUT OF SCOPE (deferred — do not claim)

- **Dynamics (D3).** `−¼F²`/`−¼G²` action, inhomogeneous equations `d⋆F=J`, masslessness, propagation
  speed `c`. Open in Paper III; this paper builds the tensor, not its action.
- **Coupling / α (D2).** `α` undetermined at linear order; no `1/137`, no `α_s`. Even with the two
  roles separated, `α` fixes the *relation*, not `α_EM/α_s`.
  `[[project_alpha_undetermined_at_linear_order]]`.
- **Confinement.** Archived 2026-06-28; not referenced as a result.
- **Soliton-layer objects.** Vortex/texture cores, spin-½ holonomy — Paper IV / matter sector.

## Named results (to promote to proposition in the paper)

- **Prop 0 (carrier).** `ψ=δR+iδṘ/ω` is the phase-space carrier; the real degenerate triplet carries
  only `O(3)` frame freedom, and complexification promotes it to `U(3)`. The `i` is the time-link
  rotation.
- **Prop 1 (link variable & gauge law).** `U_μ(n)` is well-defined off band-crossings and transforms
  as a lattice connection under carrier rephasing `u→Vu`, `V∈U(d)`.
- **Prop 2 (curvature & antisymmetry).** `□_μν=exp(ia²F_μν+O(a⁴))`, `F_μν=−F_νμ`, antisymmetry
  inherited from plaquette orientation — the structural origin of the rank-2 antisymmetric tensor.
- **Prop 3 (Faraday + Bianchi).** Trace sector ⇒ Faraday tensor; homogeneous Maxwell = exact lattice
  Bianchi identity.
- **Prop 4 (Yang–Mills field strength).** Traceless sector ⇒ `G^a_μν` with `f^{abc}` self-coupling
  from ordered non-commuting links; non-Abelian Bianchi.
- **Prop 5 (`α` activates the colour sector).** The single parameter `α` turns the non-Abelian sector
  on monotonically: `α=0` (linear) ⇒ Abelian `U(1)³`, traceless weight `∝α → 0`; `α→1` ⇒ degenerate
  triplet, full `U(3)`, traceless `SU(3)` on, weight `(3−2α):α`. Degeneracy/frame-freedom and
  curvature-weight increase together — no tension, no interior sweet spot.
- **Theorem (equivalence).** Continuum limit of the `(x,t)` plaquette holonomy equals the Berry / WZ
  curvature — the plaquette construction *is* Paper III's identification, now grounded in the stencil.

## Open derivations

*Per-paper bridge entries (`[[project_open_problems_tracker]]`). This bridge = `field-strength`.*

1. **Link-variable definedness across band-crossings.** `U_μ` needs band isolation; characterize where
   the carrier band is gapped enough that `U_μ` is single-valued (ties to Paper III adiabaticity).
2. **`O(a⁴)` lattice-artifact terms.** Confirm the subleading plaquette expansion is a genuine
   higher-derivative correction (no spurious symmetric or parity-odd piece) — clover/improvement check.
3. **Non-Abelian ordering convergence.** Verify the `f^{abc}` term is the unique non-commuting remnant
   and that the symmetric/path-ordering choice does not shift the continuum `G^a_μν`.
4. **Numerical equivalence check.** Transport a band-isolated `ℂ³` carrier around a fixed `(x,t)`
   plaquette and confirm `log □ = ia²F` matches `branesim/diagnostics/berry_holonomy.py` (the same loop
   the D2 falsifier `R(0.5)/R(0.2)=0.3077` uses).

## Self-flagged risks (carry into critique when reviewed)

1. **Band-isolation is load-bearing.** `U_μ` is single-valued only off band-crossings; the colour `d=3`
   subspace is degenerate by construction and `g([111])=0` means the colour-coherent direction has no
   defined holonomy. The construction is valid *where the subspace is band-isolated* — for colour,
   essentially the soliton layer (Paper IV), not a free linear direction.
2. **"Field-strength tensor" ≠ "the gauge theory."** Constructing `G^a_μν` is not deriving QCD. Keep the
   action, running coupling, and confinement (archived) explicitly out.
3. **The complex structure is earned.** The `i` is the timelike-link rotation, now realized concretely as
   the phase-space pairing `ψ=δR+iδṘ/ω` (`carrier_construction.md`). Always trace it to the temporal
   spring; cite `[[project_complex_u1_from_time]]`.
4. **No `ℂ³`-as-`u(3)`-matrix slip.** The split is of the *connection* valued in `u(3)`, not of the
   carrier vector. Single vector ⇒ `U(1)` phase; transported three-frame ⇒ `U(3)` connection. Keep them
   distinct (`carrier_construction.md` §5).
5. **`O(a⁴)` honesty.** Lattice plaquettes carry higher-derivative artifacts; the claim is the *leading*
   term is `a²F_μν`, pending the clover/improvement check (Open derivation 2).
6. **Keep `α` consistent and monotonic.** Every effect of `α` points one way: the colour/non-Abelian
   sector turns on with `α` and vanishes at the linear limit `α=0`. Do not reintroduce a
   degeneracy-vs-weight "tension" or an interior sweet-spot `α` — those came from stating the degeneracy
   direction backwards against the anchored weighting `traceless ∝α`. State the two clean limits, not an
   experiment-fitted value.

## Honesty headline

This paper closes the **kinematic construction** of `F_μν` and `G^a_μν` from the phase-space carrier on
the 8-link spring stencil and proves it is the Berry/WZ object of Paper III. It does **not** close the
**dynamics** (no action, no field equations, no propagation) or the **coupling** (`α` undetermined — it
fixes the *relation/ratio* between sectors, not `α_EM/α_s`). Say so plainly in the scope section.