# Reconciliation: the `(x,t)` Plaquette Curvature IS the Berry / Wilczek–Zee Curvature

The central theorem. Paper III *identifies* `F_μν` with the Berry curvature `f_μν` and `G^a_μν` with the
Wilczek–Zee curvature `𝓕_μν`. This note shows the plaquette construction of the phase-space carrier
(`carrier_construction.md`, `link_holonomy.md`) is that same object — so the identification is grounded in
the stencil, not posited.

## 1. The two constructions

- **Paper III (differential):** connection `a_μ = i⟨u|∂_μu⟩`, curvature `f_μν = ∂_μa_ν − ∂_νa_μ`
  (rank-1); WZ connection `𝒜_μ = i⟨u_i|∂_μu_j⟩`, curvature
  `𝓕_μν = ∂_μ𝒜_ν − ∂_ν𝒜_μ − i[𝒜_μ,𝒜_ν]` (rank-d).
- **This paper (integral/lattice):** carrier eigenframe `|u(n)⟩` from `ψ = δR + iδṘ/ω`, link
  `U_μ(n) = ⟨u(n)|u(n+ê_μ)⟩/|·|`, plaquette `□_μν`, `log □_μν = −ia²F_μν + O(a⁴)`.

## 2. Equivalence

Expand the lattice overlap link to first order in `a`:

```
    U_μ(n) = ⟨u(n)|u(n+ê_μ)⟩/|·| = exp( −i a · i⟨u|∂_μu⟩ + O(a²) ) = exp( −i a · a_μ + O(a²) ).
```

So the lattice link variable is the exponentiated Berry connection sampled on the link, and the plaquette
`□_μν` is the Wilson loop of `a_μ` around the elementary square. Its log is, by Stokes / the non-Abelian
Stokes theorem, the curvature flux through the plaquette:

```
    log □_μν = −i a² f_μν + O(a⁴)   (Abelian),     −i a² 𝓕_μν + O(a⁴)   (non-Abelian).
```

**Therefore `F_μν = f_μν` and `G_μν = 𝓕_μν`** — the plaquette field strength and the Berry/WZ curvature
are the same tensor, one computed by an integral (lattice loop), the other by a derivative. This is the
standard lattice-gauge ↔ continuum-connection correspondence; the paper's content is that the *physical*
connection is the phase-space carrier overlap on the spring stencil.
[Open derivation 4: numerically confirm `log □ = −ia²f` on a fixed `(x,t)` loop against
`branesim/diagnostics/berry_holonomy.py`.]

## 3. Why the Brillouin-zone curvature vanishes (restated from this side)

A natural worry: if `F` is a plaquette curvature, why doesn't the *displacement* lattice — or the
Brillouin zone — already carry it? Because the dynamical matrix `D(k)` is **real symmetric**
(`carrier_construction.md` §1a), its eigenvectors `|u(k)⟩` can be chosen real, so the `k`-space overlap
`⟨u(k)|u(k+δk)⟩` is real and every `k`-space plaquette `□^{BZ}` is real-positive: `log □^{BZ} = 0`, i.e.
the **BZ Berry/WZ curvature is identically zero `∀α`** (BACKBONE #16;
`[[project_spin_half_is_soliton_layer]]`).

This holds *even when the rest length makes `u(k)` genuinely `k`-dependent.* The transverse-stiffness knob
`(1−α)` in the bond matrices `K_μ` can rotate the polarization eigenframe across `k`
(`carrier_construction.md` §1a), but real-symmetric `D(k)` keeps that frame real, so a `k`-dependent-but-
real `u(k)` still gives `A_μ = i⟨u|∂u⟩ = 0`. **`k`-mixing is necessary but not sufficient; it can seed at
most a `ℤ₂` (Zak) holonomy, not a `U(1)` curvature.** (That Zak `ℤ₂` is the spin-1 envelope's — not
spin-½, which is a soliton-layer real-space holonomy; `carrier_construction.md` §1b.)

The nonzero plaquette lives over **physical `(x,t)`**, where the carrier state genuinely rotates in `ℂ³`
— and it rotates because the complex part of `ψ = δR + iδṘ/ω` is the carrier's advance along the timelike
worldtube axis (`[[project_complex_u1_from_time]]`). A purely spatial snapshot is a real configuration
with `F = 0` in its spatial plaquettes; the field strength is switched on by the temporal link. This is
the same fact Paper III states ("the gauge structure does not live in the Brillouin zone"), now seen as:
*real links ⇒ trivial plaquette; the complex link is the time-rotation, so two time slices are minimal.*

## 4. Consequence for the paper's logic

The reconciliation lets the paper present **one** construction (the phase-space carrier plaquette) and
harvest both descriptions: the lattice plaquette gives the exact Bianchi identities and the
antisymmetry-from-orientation argument cleanly; the differential Berry/WZ form connects to Paper III and to
the holonomy-ratio prediction `R(0.5)/R(0.2) = 0.3077` — derived here from the triplet spectrum
`λ_A = κ[(1−α)|k|² + α k_A²]` (`field_tensors.md` §4), on the same `(x,t)` loop. The plaquette/Wilson-loop
objects are also the diagnostics of `lgt_diagnostics.md` — computed on the *same* substrate graph, not a
second lattice.

## 5. Status

- `closed` (as derivation): plaquette curvature = Berry/WZ curvature via the link expansion + Stokes.
- `closed` (restatement): BZ-vanishing from real-symmetric `D(k)`, including the `k`-dependent-but-real
  case.
- `open`: numerical loop-equivalence check (Open derivation 4).
