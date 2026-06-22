# Bell Bridge — Pancharatnam–Berry holonomy and the form of E(a,b)

Working note. **Derives the functional FORM `E(a,b) = −cos 2(a−b)` (incl. the
factor of 2), not the Bell violation.** Honest debt; not asserted in the paper.
Companion to `vvertex_energy_transport.md` (§"Phase language").

## Scope

`−cos 2(a−b)` is the Pancharatnam–Berry (Poincaré-sphere) geometry of
polarization. This note derives it cleanly and states exactly what the substrate
worldtube transport (B5) must reproduce. It is the QM-side *target form*; the
part that exceeds the local CHSH bound `S = 2` is supplied by the projective
threshold (B4) + retro-coupling at the vertex (B5), **not** by the geometry
(see Caveat). Algebra below verified numerically (`E = −cos2(a−b)`, CHSH
`|S| = 2√2`).

## 1. Polarization space and the Poincaré double angle

Linear polarization at physical angle θ:

    |θ⟩ = cos θ |H⟩ + sin θ |V⟩,    θ ∈ [0, π).

Amplitude overlap, and Malus's law as its modulus:

    ⟨a|θ⟩ = cos(θ − a),    P(click at analyzer a | θ) = |⟨a|θ⟩|² = cos²(θ − a).

Bloch / Poincaré vector  s(θ) = ⟨θ| σ |θ⟩  (σ = (σ_x, σ_y, σ_z)):

    s_x = ⟨σ_x⟩ = 2 cosθ sinθ = sin 2θ
    s_y = 0
    s_z = ⟨σ_z⟩ = cos²θ − sin²θ = cos 2θ
    ⇒  s(θ) = (sin 2θ, 0, cos 2θ).

So a physical angle θ embeds on the Poincaré sphere at azimuth **2θ** — the
Pancharatnam–Berry double angle. Orthogonal polarizations (Δθ = π/2) are
antipodal (Δ(2θ) = π); the linear-polarization state space is `RP¹` (angles mod
π) sitting as a great circle via θ ↦ 2θ. Frame inner products:

    s(θ₁)·s(θ₂) = cos 2(θ₁ − θ₂),
    |⟨θ₁|θ₂⟩|² = (1 + s₁·s₂)/2 = cos²(θ₁ − θ₂).   ✓ (consistent with Malus)

## 2. Measurement observable

Analyzer at `a` defines the ±1 observable

    A_a = |a⟩⟨a| − |a+π/2⟩⟨a+π/2| = s(a)·σ,   eigenvalues ±1, axis s(a) = (sin 2a, 0, cos 2a).

(One checks (s(a)·σ)|a⟩ = +|a⟩.) The threshold model (B4) supplies the ±1
outcome physically; `A_a` is its projective idealization.

## 3. The singlet correlation

The source vertex `S` emits the rotation-invariant singlet

    |Ψ⁻⟩ = (1/√2)( |H⟩_A|V⟩_B − |V⟩_A|H⟩_B ),

with the defining SU(2) identity  ⟨Ψ⁻| (m̂·σ) ⊗ (n̂·σ) |Ψ⁻⟩ = − m̂·n̂.  Hence

    E(a,b) = ⟨Ψ⁻| A_a ⊗ A_b |Ψ⁻⟩ = − s(a)·s(b) = − cos 2(a − b).    ∎

Three pieces, each with a distinct origin:

- **factor 2** — the Poincaré embedding θ ↦ 2θ (Pancharatnam–Berry geometry, §1);
- **cosine** — inner product of the two analyzer axes on the sphere;
- **minus sign** — singlet anti-alignment, locked at `S`.

CHSH at `(a,a',b,b') = (0, π/4, π/8, 3π/8)` gives `|S| = 2√2` (Tsirelson),
verified numerically.

## 4. The holonomy reading

The angle between the two analyzer axes on the sphere is

    ∠( s(a), s(b) ) = 2(a − b).

Parallel-transporting a polarization frame along the geodesic from `s(a)` to
`s(b)` rotates it by exactly this angle (Pancharatnam connection on `CP¹`); the
overlap of the transported frame with the target is its cosine. The singlet ties
the A- and B-frames anti-parallel at `S` (the minus), so the transported
correlation is `−cos 2(a−b)`.

Geometric-phase hallmark: for a closed geodesic circuit through `s(a)`, the
`S`-reference, and `s(b)`, the acquired phase is `−Ω/2` (half the enclosed solid
angle) — the Berry curvature of `CP¹`. The correlation is the real Pancharatnam
overlap around this circuit.

## 5. Substrate translation — what B5 must reproduce

In the worldtube picture the polarization frame is the local eigenvector
`|u(x)⟩` (Paper III), and the connection is the Berry connection
`a_μ = i⟨u|∂_μ u⟩` over `(x,t)`. The discrete holonomy of §4 is realized by:

- transport along arm `A` from analyzer `D_A` back to the vertex `S`,
- the singlet frame relation locked at `S` (the conserved polarization invariant
  of the V-vertex junction, `vvertex_energy_transport.md` §"Conservation laws"),
- transport along arm `B` from `S` forward to `D_B`.

**B5 claim:** the substrate's time-symmetric worldtube transport reproduces this
Poincaré-sphere holonomy — the accumulated frame rotation `A→S→B` equals
`2(a−b)`, and the junction supplies the singlet minus sign. The energy-language
picture (offer/confirmation flux) is conjugate to this phase-language picture.

## Caveat — geometry gives the form, not the violation

§3 uses the quantum projective rule: `A_a` is a difference of projectors and
`E = ⟨Ψ⁻| A_a ⊗ A_b |Ψ⁻⟩`. A **local** model in which each photon carries a
definite Poincaré point λ and outputs `sign( s(λ)·s(a) )` reproduces the
marginals and the single-arm `cos²` law, but yields the **triangular/sawtooth**
joint correlation, bounded by `|S| ≤ 2` — not `−cos 2(a−b)` (whose CHSH value is
`2√2`). So the Poincaré/Pancharatnam geometry fixes the functional form and the
factor 2; the excess over `S = 2` requires the measurement to **project**
(setting-chosen threshold, B4) fed by the **retro-coupling** at `S` (B5), not to
reveal a pre-existing λ. **Geometry necessary, not sufficient.**

## References (wire into `references.bib` if this graduates)

- Pancharatnam 1956; Berry 1984, 1987 — geometric phase of polarized light,
  Poincaré sphere, the `2θ` double angle.
- Standard singlet identity ⟨Ψ⁻|(m̂·σ)(n̂·σ)|Ψ⁻⟩ = −m̂·n̂ (any QM text).
