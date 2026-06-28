> **Archived 2026-06-28:** superseded by `derivations/audit_2026-06-28.md`; retained for provenance only.

# Connections, parallel transport, and holonomy: the unifying spine

**Status:** `closed` as a conceptual/structural framing that organizes the established results;
introduces no new dynamics. The individual holonomy facts are derived in the sibling notes.
**Sources:** `lorentz_gravity` (induced metric / geometry), `gauge_color/sections/02_gauge.tex`
(`eq:berry-connection`, `eq:wz-connection`, `subsec:spinorial-holonomy`),
`gauge_color/sections/03_color.tex` (`subsec:color-confinement`, `subsec:u3-decomp`),
`matter_mass/sections/02_matter.tex` (`subsec:spherical-material-coords`, `eq:hedgehog`).

> **Thesis.** Every structure in the model is a *connection* on the single `ℂ³` carrier bundle over
> the brane worldvolume, and every confinement statement is a constraint on a *holonomy* (the
> path-dependence of parallel transport). The prestress `α` controls the **reducibility** of the
> internal connection (Abelian `U(1)³` ↔ non-Abelian `U(3)`); the lattice **regularizes** every
> holonomy integral at the core. This is the differential-geometric backbone of the paper.

---

## 1. Curvilinear coordinates and bases (the setup)

The intrinsic brane domain `Ω` carries material coordinates; the embedding `X^A:Ω×ℝ→ℝ⁴` induces the
metric `g_{ij}=∂_iX^A ∂_jX^A` (Paper I). A **coordinate basis** of the tangent space is
`e_i = ∂_iX` with dual `e^i`; `g_{ij}=e_i·e_j`. The lateral triplet `ξ` (and its complex promotion
`Ψ∈ℂ³`) is a vector field carrying an index in this space.

**Ambient symmetry; the asymmetry is a brane property.** In the Cartesian description all four
ambient components `X^A` (`A=1,…,4`) are geometrically equivalent — the ambient `ℝ⁴` has no preferred
direction. The asymmetry appears only because the brane action gives one lattice direction Lorentzian
sign structure; the timelike direction is selected by the **brane/worldvolume**, not by the Euclidean
ambient (Paper I). On a spacelike slice the fourth component appears as a transverse amplitude
direction `u=X⁴`.

**Spherical coordinates are a defect-adapted chart, not a second ontology.** For particle-like states
we write `X^A = X^A(r,θ,φ,t)`, where `(r,θ,φ)` organize the *three spatial material coordinates* and
`u=X⁴` is a **scalar / fiber-like component over that spherical base**. The angular structure rotates
the *spatial* basis vectors `ê_r,ê_θ,ê_φ` (and, in the hedgehog, the locked lateral triplet), but it
does **not** rotate the fourth component `u`. The underlying dynamics remains the coordinate-free
embedding `X:Ω×ℝ→ℝ⁴` (Paper IV `subsec:spherical-material-coords`); the chart is a convenience for
localized excitations, not a new layer of reality.

**Why a spatial defect lifts to a 4-D object.** Because the timelike direction `u` is *not* part of
the spatial angular frame, a defect that is point-, line-, or toroidal-like in a 3-D slice is swept
along the fourth/timelike direction into a 4-D object: a **worldline** for a point texture (`SU(3)`,
`su3_texture.md`), a **worldsheet** for a vortex centerline (`U(1)`, `u1_vortex_core.md`), and a
**worldtube** for a closed toroidal mode (the electron-as-torus; `spin_chirality.md`). The spatial
chart sets the cross-section; the timelike sweep makes it a worldvolume object.

Two bases are in play and the transition between them is the first geometric object of the paper:

- **Cartesian / lattice-aligned** `{ê_x,ê_y,ê_z}` — constant in space; the natural basis of the
  `O_h` lattice and of the `U(3)` triplet `(ψ_x,ψ_y,ψ_z)`.
- **Curvilinear** — spherical `(r,θ,φ)` for the point texture, cylindrical `(ρ,φ,z)` for the line
  vortex. Here the basis vectors **vary with position** (`∂_a ê_b ≠ 0`).

The change Cartesian↔curvilinear is a **position-dependent frame rotation** `R(x)∈SO(3)`.

**Category note (important).** Expressing a vector field in a basis that turns as you move makes the
naive derivative non-tensorial; the correction is a connection. This is **moving-frame
covariantization**, *not* a new constitutive nonlinearity on the same footing as the norm law. The
norm term `−k_sαa|ΔR|` is a genuine dynamical interaction (`nonlinearity.md`); the `SH↔Cartesian`
connection terms are a *change of description* (covariant transport in a rotating frame). Both are
"geometric," but only the first changes the energy; the second reorganizes how derivatives are
written. Keeping this distinction avoids a category error.

## 2. The covariant infinitesimal derivative and the Christoffel symbols

Because `∂_a ê_b ≠ 0`, the ordinary partial derivative of a vector field mixes the true field change
with the basis rotation. The **covariant derivative** restores tensoriality:

    ∇_a V^b = ∂_a V^b + Γ^b_{ac} V^c,
    Γ^b_{ac} = ½ g^{bd}( ∂_a g_{cd} + ∂_c g_{ad} − ∂_d g_{ac} ),

the **Levi-Civita connection** of `g` (metric-compatible `∇g=0`, torsion-free `Γ^b_{ac}=Γ^b_{ca}`).
The infinitesimal **parallel transport** of `V` along `dx^a` is `δV^b = −Γ^b_{ac}V^c dx^a`, and the
failure to return after a closed loop is the **holonomy**, measured by the Riemann curvature
`R^b_{cad} = ∂_a Γ^b_{dc} − ∂_d Γ^b_{ac} + Γ^b_{ae}Γ^e_{dc} − Γ^b_{de}Γ^e_{ac}`.

**Concrete payoff — the radial-ODE barriers ARE Christoffel terms.** The "geometric nonlinearity of
the coordinate transformation" is not abstract: the connection coefficients are exactly the
angular-momentum barrier terms in the dimensionally-reduced equations.

- *Line vortex, cylindrical.* For `Ψ_tr=f(ρ)e^{inφ}`, `|∇Ψ_tr|² = f'² + (n²/ρ²)f²`; the `n²/ρ²`
  is the connection/centrifugal term from the `φ`-coordinate metric `g_{φφ}=ρ²` (`u1_vortex_core.md`).
- *Point texture, spherical.* For the hedgehog `ξ^i=f(r)x̂^i`,
  `∂_iξ_j = f' x̂_ix̂_j + (f/r)(δ_{ij}−x̂_ix̂_j)`; the `f/r` pieces are the Christoffel terms
  (the `x̂` basis rotates), and the radial operator `f'' + 2f'/r − 2f/r²` is the `L=1` **covariant
  vector Laplacian** — the `−2f/r²` is purely the Levi-Civita barrier (`su3_texture.md`).

So the SH/curvilinear reduction that makes both defects tractable *generates* the connection terms;
they are the price (and the tool) of trading 2-D/3-D Cartesian problems for 1-D radial ones.

## 3. Three connections on one bundle

The full covariant derivative of the carrier triplet combines a **spatial-frame** connection and an
**internal-frame** connection:

    D_a Ψ^i = ∂_a Ψ^i  +  Γ^i_{ac} Ψ^c        (Levi-Civita, spatial/lateral frame)
                          +  ( i a_a 𝟙 + 𝒜_a )^i_{\ j} Ψ^j   (gauge, internal carrier frame).

The three connections and their parallel-transport meaning:

| Connection | Symbol | Transports | Sector | Status |
|---|---|---|---|---|
| Levi-Civita | `Γ^k_{ij}` of `g_{ij}` | displacement vectors over the brane | geometry / gravity (Paper II) | **dynamical** |
| Berry | `a_μ=i⟨u\|∂_μu⟩` | the carrier phase | `U(1)` / EM (Paper III) | diagnostic |
| Wilczek–Zee | `𝒜_μ∈u(3)` | the colour frame | `SU(3)` / colour (Paper III) | diagnostic |

**Important distinction (kept throughout).** The Levi-Civita connection is the *actual substrate
geometry* and back-reacts (it is the elastic/gravitational dynamics). The Berry/WZ connections are
*diagnostics* read out of the state and never fed back as forces (No-Back-Reaction, PRINCIPLES §2).
They unify *conceptually* as connections-with-holonomy; they do **not** have equal dynamical status.
This is the one place the paper must not overclaim a "geometry = gauge" identity.

## 4. α controls the reducibility of the internal connection

The single dial `α` sets whether the internal (WZ) connection is reducible:

- `α→1` (no prestress): branches maximally split, `D(k)` maximally non-degenerate ⇒ `𝒜_μ` reduces
  to block-diagonal `u(1)³` — Abelian, three decoupled axis phases. Holonomy group `→ U(1)³`.
- `α→0` (max prestress): `D(k)∝I`, branches degenerate ⇒ the transported subspace is genuinely 3-D
  and `𝒜_μ∈u(3)` is irreducible non-Abelian. Holonomy group `→ U(3)`.

The off-diagonal coupling that makes the connection irreducible is the cross-vertex
`g=Θk_sα/a²` (`nonlinearity.md` §3) — `∝α`, vanishing exactly at `α=0`'s linear... (note the dual
role: `α→0` *degenerates the spectrum* but also *kills the anharmonic coupling*; the physical
soliton sits at intermediate `α≈0.5–0.8` where both the splitting and the coupling are `O(1)`).
**The U(1)↔SU(3) transition is therefore a connection-reducibility crossover mediated by the
prestress.**

## 5. Confinement as holonomy constraints

Each confinement result restated as a holonomy fact (derivations in the sibling notes):

| Phenomenon | Holonomy statement | Quantized by |
|---|---|---|
| `U(1)` charge | `∮ a = 2πn` around the centerline | `π₁(U(1))=ℤ` |
| Colour confinement | WZ holonomy is **well-defined only off `[111]`** (non-degenerate) but the triplet is **coherent only on `[111]`** (degenerate) — mutually exclusive ⇒ no colored holonomy survives to infinity | `subsec:color-confinement` |
| `ℤ₃` triality | single-valuedness ⇒ trivial **total** `U(3)` holonomy ⇒ `q≡−t (mod 3)` | `U(3)=(U(1)×SU(3))/ℤ₃` |
| Spin-½ | `ℤ₂` holonomy under a `2π` rotation of the combined (spatial+internal) frame | `π₁(SO(3))=ℤ₂` |
| Gravity (context) | holonomy of `Γ` around a loop = Riemann curvature of the deformed brane | Paper II |

**The lattice regularizes every holonomy.** The minimal loop is the unit plaquette (perimeter `~8a`);
all curvatures (Riemann, Berry, WZ) are evaluated on finite plaquettes, so the holonomy integrals are
finite and there is no singular core (`nonlinearity.md` §5). The continuum singularity that would
make a connection ill-defined at a defect center is cut off at `ρ=a`.

## 6. Parallel transport unifies the field strengths (conceptually)

Transporting the lateral triplet around a small closed loop and reading the residual rotation gives,
in each bundle, the curvature:

- spatial frame → Riemann `R^i_{jab}` (gravity / induced-metric curvature);
- internal phase → Berry `f_{μν}` ≅ `F_{μν}` (EM field tensor);
- internal frame → WZ `𝒡_{μν}` ≅ the colour field tensor.

Same machinery (loop holonomy = curvature), three bundles. This is the conceptual unity the paper
draws on; it is *not* a claim that the gauge curvatures are metric curvatures (they live on different
bundles, and the gauge ones are diagnostic).

## 7. Ledger

| Item | Status |
|---|---|
| Curvilinear-basis / covariant-derivative / Christoffel framing | `closed` (standard diff-geo) |
| Radial barriers = Levi-Civita connection terms | `closed` |
| Three connections on one `ℂ³` bundle; combined `D_aΨ` | `closed` (organizational) |
| Levi-Civita dynamical vs Berry/WZ diagnostic | `closed` (PRINCIPLES) |
| `α` = reducibility dial `U(1)³↔U(3)` | `closed` |
| Confinement = holonomy constraints (the table) | `closed` (restatement) |
| "Geometry = gauge" full identification | **rejected** (different bundles, diagnostic vs dynamical) |
