> **Archived 2026-06-28:** no `SU(3)`-valued order parameter was constructed, so the winding assignment is unavailable.

# The SU(3) point texture: winding, not harmonics; Derrick + lattice anti-collapse

**Status:** the SH-vs-winding dichotomy, kinematic confinement, and the `ℤ₃` triality lock are
`closed`; the size law is a **Derrick/lattice scaling within the hedgehog reduction** (not an
existence claim); the existence of a converged `B≠0` eigenstate is `open` (numerical). *Per the
external critique: topology classifies the sector, but stabilization in 3-D is a substantive
dynamical problem — the scaling law is plausible, the eigenstate is not yet demonstrated.*
**Sources:** `gauge_color/sections/03_color.tex` (`subsec:topology`, `subsec:color-confinement`,
`eq:holonomy-ratio`, `eq:hedgehog`), `matter_mass/sections/02_matter.tex` (`eq:vsh-decomposition`,
`eq:skyrme-twist`, `eq:derrick-radius`), `matter_mass/derivations/matter/status.md` (C3 scaling leg).

---

## 1. Object: a point defect, not a line defect

The colour sector carries different topology from the `U(1)` sector
(`subsec:topology`):

- `π₁(SU(3)) = 0` — **no vortices** in the colour sector;
- `π₃(SU(3)) = ℤ` — **textures** (Skyrmion/instanton class), **point-like** in a 3-D slice, a
  codim-3 **worldline** in 4-D.

So the two confined objects are topologically distinct: the `U(1)` charge is a codim-2 worldsheet
(`u1_vortex_core.md`); the colour/baryon is a codim-3 worldline. They are the same field `Ψ∈ℂ³`,
different homotopy classes of it.

## 2. Why spherical harmonics are NOT the SU(3) language — `closed`

This answers the owner's intuition that "for SU(3) we might need another framework than SH."

- The `U(1)` vortex is a codim-2 object with an **axis**; its transverse profile separates into a
  radial function times an angular harmonic (`Y₁¹`), so SH is exactly the right basis — the angular
  index labels the multipole content of the phase around the centerline.
- The `SU(3)` texture has **no axis** and is classified by a **degree/winding map**. The relevant
  invariant is the integral

      B = (1/24π²) ∫ ε^{ijk} tr(L_i L_j L_k) d³x,   L_i = U⁻¹ ∂_i U,   U ∈ SU(3),

  i.e. the homotopy degree of the map `S³ → SU(3)` (boundary-compactified slice → group). A
  multipole/SH expansion of `U(x)` does **not** see this integer — winding is not a harmonic. SH is
  the `U(1)` language; the winding integral is the `SU(3)` language. Stating this dichotomy crisply
  is one deliverable of the paper.

## 3. Reducing the texture to a 1-D radial ODE: the combined-SO(3) hedgehog — `closed`

Although the texture is not described by a single SH, its *symmetry* still gives a 1-D reduction —
not via the spatial-`SO(3)` SH but via the **combined** (diagonal) `SO(3)` that rotates spatial and
internal/colour indices together. The hedgehog (`eq:hedgehog`, `eq:vsh-decomposition` with
`J=0,L=1`):

    ξ^i(x) = f(r) x̂^i      (internal colour index locked to spatial direction).

Invariance under the combined `SO(3)` collapses the field to a single radial profile `f(r)` (and,
with the `X⁴` Skyrme twist `eq:skyrme-twist`, a single Skyrme angle `F(r)`, `F(0)=π`, `F(∞)=0`,
giving `B=1`). The 3-D texture problem becomes a **1-D radial ODE** — the analogue, in the colour
sector, of the `U(1)` radial reduction.

The hedgehog is a *section* `x̂^i` that locks the internal colour frame to the spatial frame, so it
rotates with position; differentiating it consistently needs the **combined** connection (spatial
Levi-Civita + internal gauge rotation), and the "combined `SO(3)`" is precisely parallel transport
under that diagonal connection (`connections_holonomy.md` §2–3). Concretely
`∂_iξ_j = f' x̂_ix̂_j + (f/r)(δ_{ij}−x̂_ix̂_j)`: the `f/r` pieces are the Christoffel terms of the
rotating `x̂` basis, and the `W₂`-radial operator `f'' + 2f'/r − 2f/r²` is the covariant `L=1` vector
Laplacian — the `−2f/r²` is purely the Levi-Civita barrier (the colour-sector counterpart of the
vortex's `n²/ρ²`). The trace part `f(r)(x̂¹+x̂²+x̂³)` averages to zero on the
sphere ⇒ trace-neutral far field (neutron-like); an `L=0` admixture shifts it to proton-like.

## 4. Closed size and lattice anti-collapse — `closed`

The texture has the same Derrick balance as the soliton of Paper IV. On the central-force lattice
there is **no independent Faddeev–Skyrme gradient-quartic**: the norm term depends only on the
scalar `|ΔR|²`, so the only 4-gradient operator is the **symmetric** geometric quartic
`(∂ξ·∂ξ)²`, never the antisymmetric `(∂ξ∧∂ξ)²` (matter-bridge C3 result). Therefore:

- the `σ`-model 2-gradient piece IS `E₂^⊥` (`∝A_⊥²|∂n̂|²`);
- the only `λ⁻¹` stabilizer is the geometric quartic `E₄` — `SU(3)` itself sets **no length**.

Derrick balance `λ* = √(E₄/E₂)` (`eq:derrick-balance`) gives

    R_h/a = κ (A/a) √(α/(1−α))        (κ = O(1)),

linear in amplitude, increasing in `α` — the **same** `√(α/(1−α))` as the `U(1)` healing length
(`u1_vortex_core.md` §4). **Anti-collapse** is the lattice UV cutoff: Derrick's continuous rescaling
`λ∈(0,∞)` is bounded below by `λ≥a`; energy concentrated at one node disperses immediately
(`subsec:derrick`). No separate quartic repulsion is needed to block collapse — the discreteness
does it. The geometric quartic supplies only the **anti-expansion** (IR) balance that sets `R_h`.

`SU(3)`'s job is therefore topological **existence** + **kinematic confinement**, not length-setting.

## 5. Kinematic colour confinement — `closed` (restated)

A free linear triplet wavepacket must be simultaneously **coherent** (the three branch frequencies
degenerate, `g(k̂)=0`, which holds **only** along `[111]`) and **colour-active** (branches
non-degenerate, `g(k̂)≠0`, i.e. **off** `[111]`, so the traceless holonomy is defined). No linear
direction satisfies both: on `[111]` colour is undefined; off `[111]` the triplet dephases over long
range (`subsec:color-confinement`). Hence:

> **No free colored asymptotic state.** Colour charge has support only in a nonlinear, phase-locked
> bound state held together by the geometric coupling, not by spectral degeneracy.

This is a statement about asymptotic states, not a linear inter-quark potential. Falsifiable handle:
the off-`[111]` holonomy ratio `R(α₂)/R(α₁) = [(3−2α₂)/α₂]/[(3−2α₁)/α₁]` = `0.3077` for `(0.2,0.5)`
(`eq:holonomy-ratio`).

## 6. ℤ₃ triality charge-lock — `closed` (the firmest confinement result)

Because `U(3) = (U(1)×SU(3))/ℤ₃` (the centre `ℤ₃⊂SU(3)` identified with cube-roots of unity in
`U(1)`), single-valuedness forces the trace charge `q` and the `SU(3)` triality `t∈{0,1,2}`
(`𝟙→0, 𝟑→1, 𝟑̄→2, 𝟖→0`) to satisfy

    q ≡ −t  (mod 3).

A colour **triplet** (single quark, `t=1`) is forced to carry fractional (`⅓`-integer) charge; a
colour **singlet** (`t=0`) carries integer charge. This is exact, dynamics-free, and given only that
the field is genuinely `U(3)`. It is a *representation*-level confinement: an isolated fractional
charge cannot be a triality-0 asymptotic state, so quarks combine into singlets (`B=3` baryon,
`q∈ℤ`). It does **not** by itself fix the spatial co-location of charge and colour (that is the
binding question, `binding.md`).

## 7. Ledger

| Item | Status |
|---|---|
| Point texture, `π₃(SU(3))=ℤ`, codim-3 worldline | `closed` |
| SH-vs-winding dichotomy (SH is the `U(1)` language) | `closed` |
| Combined-`SO(3)` hedgehog → 1-D radial ODE | `closed` |
| `R_h/a=κ(A/a)√(α/(1−α))`; `SU(3)` sets no length | `closed-in-ansatz` (hedgehog reduction) |
| Lattice anti-collapse (UV cutoff `λ≥a`) | `closed` |
| Kinematic colour confinement (no free colored state) | `closed` (asymptotic) |
| `ℤ₃` triality `q≡−t (mod 3)` | `closed` (dynamics-free) |
| Converged stable `B≠0` eigenstate | `open` (numerical) |
| `κ` prefactor; existence proof independent of Derrick | `open` |
