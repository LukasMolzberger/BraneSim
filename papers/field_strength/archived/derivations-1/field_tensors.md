# Faraday and Yang–Mills Tensors from the One `U(3)` Plaquette

Both field strengths come from the *single* carrier plaquette of `link_holonomy.md`. The `U(3)`
connection splits as `u(3) = u(1) ⊕ su(3)` (Paper III): the trace part is electromagnetism, the
traceless part is colour. This note builds both, then separates the two distinct jobs the rest-length
parameter `α` does.

## 1. The sector split of the connection

Project the `u(3)`-valued connection of `link_holonomy.md` §4 onto the trace and Gell-Mann directions:

```
    A_μ = A^0_μ (𝟙/√6) + A^a_μ T^a,    a=1..8,    T^a = ½λ^a,    tr(T^aT^b)=½δ^{ab}.
```

`A^0_μ` is the trace `U(1)` (the overall carrier phase of `carrier_construction.md` §5 — one vector,
one phase); the eight `A^a_μ` are the traceless `su(3)` colour potentials (the relative orientation of
the transported three-frame). It is the *same* plaquette feeding both; the carrier vector is never the
connection (critique #4).

## 2. Faraday tensor (trace `U(1)`)

The trace direction is the carrier-phase trace `Ψ_tr = (ψ_x+ψ_y+ψ_z)/√3`; its link variable is the
scalar phase `U_μ(n)=e^{−iaA_μ}`, `A_μ∈ℝ`. The complex phase is the timelike-link rotation, so `A_μ` is
intrinsically a 4-covector on `(x,t)` (`[[project_complex_u1_from_time]]`). The Abelian plaquette
(`[A_μ,A_ν]=0`) gives the **Faraday tensor**

```
    F_μν = ∂_μ A_ν − ∂_ν A_μ,        F_μν = −F_νμ,
    E_i = F_{0i}   (time–space plaquette, uses a temporal link),
    B_i = ½ ε_{ijk} F_{jk}   (space–space plaquette, spatial links only).
```

The split `E ↔ F_{0i}`, `B ↔ F_{ij}` is not imposed by hand: `E` is literally the holonomy of a
plaquette that steps once in time, `B` the holonomy of a purely spatial plaquette. The electric field is
the curvature that requires the temporal spring; the magnetic field is the curvature of the spatial
stencil — the concrete payoff of the 2 temporal links.

**Homogeneous Maxwell = exact lattice Bianchi.** The product of the six face plaquettes around an
elementary 3-cube spanned by `(ê_λ,ê_μ,ê_ν)` is the identity (each link traversed once each way), so

```
    ∏_{faces} □ = 1   ⇒   ∂_{[λ} F_{μν]} = 0   exactly on the lattice (no a→0 needed).
```

These are `∇·B=0` (spatial cube) and Faraday's law `∇×E+∂_tB=0` (cubes with one temporal edge):
automatic and exact, the statement that `F=dA` ⇒ `dF=0`. Pure kinematics; no dynamics invoked.

## 3. QCD field-strength tensor (traceless `SU(3)`)

For the near-degenerate triplet the link variable is the matrix `U_μ(n)∈U(3)` (`link_holonomy.md` §2).
The non-Abelian plaquette does **not** drop its commutator: adjacent link matrices fail to commute, and
the BCH expansion leaves

```
    G^a_μν = ∂_μ A^a_ν − ∂_ν A^a_μ + f^{abc} A^b_μ A^c_ν,
```

with `f^{abc}` the `su(3)` structure constants `[T^b,T^c]=if^{abc}T^a`. The `f^{abc}A^bA^c` term — the
gluon self-interaction — is **exactly the non-commutativity of the ordered link product** (`−i[A_μ,A_ν]`
of `link_holonomy.md` §4, projected to the adjoint). The Abelian sector has no such term because `U(1)`
link phases commute. So the structural difference between the Faraday and QCD tensors is precisely
whether the elementary plaquette's link matrices commute — a statement about the carrier subspace
dimension (`d=1` vs `d=3`), nothing more.

The same closed-cube argument, now with path-ordered products, gives the covariant Bianchi identity

```
    D_{[λ} G_{μν]} = 0,    D_μ • = ∂_μ • − i[A_μ, •] ,
```

automatic and exact on the lattice; the homogeneous Yang–Mills equations are kinematic.

## 4. The rest length `α` activates the colour sector (one direction)

The single parameter `α` controls the colour sector, and **every effect it has points the same way**:
the non-Abelian structure turns on monotonically with `α` and switches off in the linear limit `α=0`.
Two descriptions of the same activation:

**Which group (degeneracy / frame freedom).**

```
    α → 0 :  linear limit; the three lateral carriers decouple  ⇒  Abelian U(1)³  (three phases, no mixing)
    α → 1 :  carriers lock into the degenerate triplet          ⇒  full    U(3)   (frame freedom for SU(3))
```

The traceless `su(3)` generators are the off-diagonal rotations of the degenerate triplet; they exist
only to the extent the triplet is locked, i.e. only for `α>0` (`carrier_construction.md` §4).

**How strongly (curvature weight).** The two sectors of the *one* plaquette are weighted (Paper III §3):

```
    trace holonomy      ∝ (3 − 2α)
    traceless holonomy  ∝ α
    ratio               R(α) = (3 − 2α)/α          [D2 falsifier: R(0.5)/R(0.2) = 0.3077]
```

with the cross-sector quartic vertex `∝ α |Ψ_tr|² |Ψ_⊥|²` (also vanishing at `α=0`).

**The two descriptions agree.** Degeneracy/frame-freedom and traceless curvature-weight both increase
with `α` and both vanish at `α=0`: the colour sector is exactly the part of the connection switched on
by the anharmonicity `−k_sαa|ΔR|` `∝α` (`[[project_geometric_nonlinearity_alpha_scaling]]`). There is no
competition between them — they are one monotonic activation read two ways — and therefore **no preferred
interior value of `α`**. The clean theoretical statement is the two limits (`α=0`: Abelian `U(1)³`,
colour off; `α=1`: full `U(3)`, colour maximal), not a fitted intermediate point.

**Caveat (D2).** `α` fixes the *relation/ratio* between the sectors, **not** the physical `α_EM/α_s`
hierarchy — that is normalization-undetermined at linear order and lives in the nonlinear fibre metric
(`[[project_alpha_undetermined_at_linear_order]]`). No derived coupling hierarchy is claimed.

## 5. Scope discipline (sharp)

- **These are the field-strength tensors only — NOT the gauge theories.** No `−¼F²`/`−¼G²` action, no
  inhomogeneous equations `d⋆F=J`, no running coupling, asymptotic freedom, string tension, confinement
  (archived), mass, screening, or `1/r`-vs-exponential range. Those are dynamics (D3) or coupling (D2).
- **The colour identification is the linear/kinematic one of Paper III.** The known tension — no single
  linear direction is both coherent and colour-active (`g([111])=0`), so colour curvature is
  intrinsically a soliton-layer object (`[[project_spin_half_is_soliton_layer]]`) — is inherited and
  flagged, not resolved. The non-Abelian *construction* is valid wherever the `d=3` subspace is
  band-isolated; whether a propagating colour excitation exists is a matter-sector question.

## 6. Status

- `closed`: `F_μν` (`E`/`B` split, homogeneous Maxwell as exact lattice Bianchi) as trace-`U(1)`
  plaquette curvature; `G^a_μν` (with `f^{abc}` self-coupling, non-Abelian Bianchi) as traceless-`SU(3)`
  plaquette curvature. (Modulo band-isolation and the `[111]` coherence/colour tension.)
- `closed`: the single-direction activation of the colour sector by `α` (`U(1)³` at `α=0`, full `U(3)`
  as `α→1`; traceless weight `∝α`, trace `∝(3−2α)`) — one monotonic activation, no interior sweet spot.
- `open`: ordering-convergence check (Open derivation 3); all dynamics (D3); coupling (D2).