# Faraday and Yang–Mills Tensors from the One `U(3)` Plaquette

Both field strengths come from the *single* carrier plaquette of `link_holonomy.md`. The `U(3)`
connection splits as `u(3) = u(1) ⊕ su(3)` (Paper III): the trace part is electromagnetism, the traceless
part is colour. This note builds both, then separates the two jobs the rest-length parameter `α` does.

## 1. The sector split of the connection

Project the `u(3)`-valued connection of `link_holonomy.md` §4 onto the trace and Gell-Mann directions:

```
    A_μ = A^0_μ (𝟙/√6) + A^a_μ T^a,    a = 1..8,    T^a = ½λ^a,    tr(T^aT^b) = ½δ^{ab}.
```

`A^0_μ` is the trace `U(1)` (the overall carrier phase of `carrier_construction.md` §5 — one vector, one
phase); the eight `A^a_μ` are the traceless `su(3)` colour potentials (the relative orientation of the
transported three-frame). It is the *same* plaquette feeding both; the carrier vector is never the
connection.

## 2. Faraday tensor (trace `U(1)`)

The trace direction is the carrier-phase trace `Ψ_tr = (ψ_x + ψ_y + ψ_z)/√3`; its link variable is the
scalar phase `U_μ(n) = e^{−iaA_μ}`, `A_μ ∈ ℝ`. The complex phase is the timelike-link rotation, so `A_μ`
is intrinsically a 4-covector on `(x,t)` (`[[project_complex_u1_from_time]]`). The Abelian plaquette
(`[A_μ,A_ν] = 0`) gives the **Faraday tensor**

```
    F_μν = ∂_μ A_ν − ∂_ν A_μ,        F_μν = −F_νμ,
    E_i = F_{0i}   (time–space plaquette, uses a temporal link),
    B_i = ½ ε_{ijk} F_{jk}   (space–space plaquette, spatial links only).
```

The split `E ↔ F_{0i}`, `B ↔ F_{ij}` is not imposed by hand: `E` is literally the holonomy of a plaquette
that steps once in time, `B` the holonomy of a purely spatial plaquette. The electric field is the
curvature that requires the temporal spring; the magnetic field is the curvature of the spatial stencil —
the concrete payoff of the 2 temporal links.

**Homogeneous Maxwell = exact lattice Bianchi.** The product of the six face plaquettes around an
elementary 3-cube spanned by `(ê_λ,ê_μ,ê_ν)` is the identity (each link traversed once each way), so

```
    ∏_{faces} □ = 1   ⇒   ∂_{[λ} F_{μν]} = 0   exactly on the lattice (no a→0 needed).
```

These are `∇·B = 0` (spatial cube) and Faraday's law `∇×E + ∂_tB = 0` (cubes with one temporal edge):
automatic and exact, the statement that `F = dA ⇒ dF = 0`. Pure kinematics; no dynamics invoked.

## 3. QCD field-strength tensor (traceless `SU(3)`)

For the near-degenerate triplet the link variable is the matrix `U_μ(n) ∈ U(3)` (`link_holonomy.md` §2).
The non-Abelian plaquette does **not** drop its commutator: adjacent link matrices fail to commute, and
the BCH expansion leaves

```
    G^a_μν = ∂_μ A^a_ν − ∂_ν A^a_μ + f^{abc} A^b_μ A^c_ν,
```

with `f^{abc}` the `su(3)` structure constants `[T^b,T^c] = if^{abc}T^a`. The `f^{abc}A^bA^c` term — the
gluon self-interaction — is **exactly the non-commutativity of the ordered link product** (`−i[A_μ,A_ν]`
of `link_holonomy.md` §4, projected to the adjoint). The Abelian sector has no such term because `U(1)`
link phases commute. So the structural difference between the Faraday and QCD tensors is precisely whether
the elementary plaquette's link matrices commute — a statement about the carrier subspace dimension
(`d=1` vs `d=3`), nothing more.

The same closed-cube argument, now with path-ordered products, gives the covariant Bianchi identity

```
    D_{[λ} G_{μν]} = 0,    D_μ • = ∂_μ • − i[A_μ, •] ,
```

automatic and exact on the lattice; the homogeneous Yang–Mills equations are kinematic.

## 4. The rest length `α`: colour curvature grows with the splitting it induces

The single parameter `α` controls the colour sector through the *spectrum* of the lateral triplet. Build
the long-wavelength Bloch operator from the bond matrices of `carrier_construction.md` §1a. On the axial
stencil `K_i^{AB} = κ[(1−α)δ^{AB} + α(ê_i)^A(ê_i)^B]`, so `D^{AB}(k) = Σ_i k_i² K_i^{AB}` is diagonal with
the three polarization eigenvalues

```
    λ_A(k) = κ[ (1 − α)|k|² + α k_A² ] ,     A ∈ {x,y,z}.
```

Split into trace (mean) and traceless (splitting) parts:

```
    λ̄        = (1/3)Σ_A λ_A = κ|k|²(1 − 2α/3)   ∝ (3 − 2α)     (trace,     U(1)/EM)
    λ_A − λ̄  = α κ|k|²(k̂_A² − 1/3)              ∝ α · g(k̂)     (traceless, SU(3)/colour)
```

with `g(k̂) = √(Σ_a(k̂_a² − 1/3)²)`, `g([111]) = 0`, `g([100]) = √(2/3)`.

**Degeneracy and colour move in *opposite* directions.**

```
    α = 0 : λ_A = κ|k|² ∀A  → EXACT degeneracy.  U(3) frame freedom is available, but the traceless
                              splitting is zero, so the colour curvature is FLAT. EM (trace) on, colour OFF.
    α → 1 : λ_A = κ k_A²    → MAXIMAL splitting.  The triplet is non-degenerate; the traceless SU(3)
                              curvature is maximal (∝ α). Colour ON.
```

So the exact `U(3)` frame freedom lives at `α=0` (perfect degeneracy) and is *destroyed* as `α` grows; the
colour field strength is generated precisely by the splitting `∝ α` that breaks it. **This is a genuine
near-degeneracy tension, not a smooth co-activation:** a well-defined transportable colour band needs
approximate degeneracy (coherence), while a nonzero traceless curvature needs the splitting that spoils
it. The two requirements fight — this is the `[111]` coherence-vs-colour tension (`g([111]) = 0`), the
kinematic-confinement mechanism of the colour bridge (`[[project_alpha_undetermined_at_linear_order]]`).
The Wilczek–Zee construction is therefore the standard *near-degenerate* one: the triplet is treated as
approximately degenerate (so the `U(3)`/WZ frame is defined), with a small splitting `∝ α` supplying the
traceless curvature. The one thing that is monotonic and unambiguous is the **colour curvature weight**,
`∝ α` off `[111]`, vanishing at `α=0` — do not read this as the degeneracy also growing with `α`.

**The holonomy ratio, derived here from our own `K_μ`.** Both sectors' `(x,t)` connection magnitudes are
built from the frequency matrix `diag(ω_A(k))`, `ω_A = √λ_A`; its trace part scales as `∝ (3−2α)` and its
traceless part as `∝ α g(k̂)` (the `√` convention cancels in the ratio, since `ω̄/δω = 2λ̄/δλ`). Hence the
off-`[111]` fibre-holonomy ratio

```
    R(α) = (C/g(k̂)) · (3 − 2α)/(√3 α) ,     R(α₂)/R(α₁) = [(3−2α₂)/α₂]/[(3−2α₁)/α₁],
    R(0.5)/R(0.2) = (2/0.5)/(2.6/0.2) = 4/13 = 0.3077.
```

The ratio `R(α₂)/R(α₁)` is **normalization-independent** — the scale
`C = tr(P_{U(1)}𝒢)/tr(P_{SU(3)}𝒢)` cancels — so `0.3077` is a clean falsifiable number that falls straight
out of the substrate's own stiffness matrix (the same `(1−α)` transverse weight that sets `c_T²`). This is
a *derived result of this paper*, and it reproduces the colour-bridge prediction of Paper III.

**Caveat (spectral vs coupling — D2).** `R(α)` is a ratio of eigenvalue *splittings* (dispersion). The
step to a physical coupling ratio `α_EM/α_s` needs the connection normalization `C`, which is
undetermined at linear order (`𝒢 ∝ 𝟙`, flat ℓ² envelope): `α` fixes the *relation/ratio* between sectors,
not the physical coupling hierarchy (`[[project_alpha_undetermined_at_linear_order]]`). And over a closed
loop the uniform clock cancels (`carrier_construction.md` §6), so the traceless holonomy is sourced by the
`k`-gradient of the differing channel frequencies `ω_A(k)`, not by any rotation of the (real) eigenframe.
The number is a genuine prediction; its identification as a coupling ratio is deferred.

## 5. Scope discipline (sharp)

- **These are the field-strength tensors only — NOT the gauge theories.** No `−¼F²`/`−¼G²` action, no
  inhomogeneous equations `d⋆F = J`, no running coupling, asymptotic freedom, string tension, confinement
  (archived), mass, screening, or `1/r`-vs-exponential range. Those are dynamics (D3) or coupling (D2).
- **The colour identification is the linear/kinematic one of Paper III.** The known tension — no single
  linear direction is both coherent and colour-active (`g([111]) = 0`), so colour curvature is
  intrinsically a soliton-layer object (`[[project_spin_half_is_soliton_layer]]`) — is inherited and
  flagged, not resolved. The non-Abelian *construction* is valid wherever the `d=3` subspace is
  band-isolated; whether a propagating colour excitation exists is a matter-sector question.

## 6. Status

- `closed`: `F_μν` (`E`/`B` split, homogeneous Maxwell as exact lattice Bianchi) as trace-`U(1)`
  plaquette curvature; `G^a_μν` (with `f^{abc}` self-coupling, non-Abelian Bianchi) as traceless-`SU(3)`
  plaquette curvature. (Modulo band-isolation and the `[111]` coherence/colour tension.)
- `closed` (derived here): the triplet spectrum `λ_A = κ[(1−α)|k|² + α k_A²]`; trace `∝ (3−2α)`, traceless
  splitting `∝ α g(k̂)`; the colour curvature weight is `∝ α` off `[111]` and vanishes at `α=0`. **Exact
  degeneracy (and full `U(3)` frame freedom) is the `α=0` limit; splitting/colour grows with `α` — the two
  move in opposite directions (near-degenerate WZ / `[111]` coherence-vs-colour tension).**
- `closed` (derived here): the holonomy ratio `R(α) ∝ (3−2α)/α`, normalization-independent, giving
  `R(0.5)/R(0.2) = 0.3077` — a genuine prediction of this paper, reproducing the colour bridge.
- `open`: identification of `R(α)` as a physical coupling ratio (D2, needs `C`); ordering-convergence
  check (Open derivation 3); all dynamics (D3).
