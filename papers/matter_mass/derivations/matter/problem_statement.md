# Matter Bridge — Problem Statement

Show that the substrate supports stable localized solitonic modes. The eigenproblem formulation is in place.

## Specific tasks

**(a) Existence of at least one stable localized solution.** Prove existence of at least one stable localized solution to the nonlinear substrate equations (with periodic boundary conditions and fixed prestretch). The decisive first test is the hedgehog/Skyrme-twisted baryon-like triplet with angular structure ξ^i = f(r) x̂^i (J=0, L=1 VSH). Periodic boundary conditions are mandatory: open/free BC lets the prestressed lattice relax (contract to αa) and dominates everything, collapsing even the A=0 vacuum.

**(b) Derrick's theorem.** Address Derrick's theorem: why does the soliton not collapse, and why does it not expand?

The two instability directions are blocked by different mechanisms:

- **Anti-collapse (UV):** Lattice discreteness. Derrick's scaling argument assumes continuous rescaling x → λx for all λ > 0. On a discrete lattice with spacing a, λ is bounded below: the soliton cannot shrink past one cell. A configuration with energy fully concentrated at a single node is not a local energy minimum — the spring network immediately disperses it. The lattice spacing a is the physical UV cutoff; no quartic term is needed or relevant for the collapse direction.

- **Anti-expansion (IR) — this *is* the confinement problem.** The IR instability direction (the localized mode spreading out / dispersing to box-fill) is not a separate problem: it is the confinement problem treated at length in the color bridge (kinematic color confinement, BACKBONE #24; the π₃(SU(3)) texture). It should not be framed or solved in isolation here. The two sectors do **distinct jobs** (scaling leg of C3):

  - *SU(3) mechanics — blocks dispersal (sets no length).* Topological existence (the π₃(SU(3)) winding, B = 1) plus kinematic confinement (coherence and color-activity are mutually exclusive on the cubic lattice → no free colored asymptotic state to disperse into). This is *why it stays a localized bound object*, but it is group-theoretic/topological and sets *no* equilibrium width. See the color bridge (`papers/gauge_color/derivations/color/status.md`, BACKBONE #24/#25).

  - *Geometric quartic — sets the width.* A direct consequence of the Pythagorean norm-term nonlinearity: the transverse quartic W₄ ∝ k_s α/a is the **only** Derrick λ⁻¹ term on the central-force substrate (the norm term depends only on the scalar |ΔR|², so it produces the symmetric (∂ξ·∂ξ)² contraction, never an antisymmetric Faddeev–Skyrme term). Under rescaling in 3D, W₂ → λ E₂ (gradient, favors expansion) is balanced against W₄ → λ⁻¹ E₄, with equilibrium at λ* = √(E₄/E₂), giving R_h/a = κ (A/a) √(α/(1−α)), linear in amplitude A and increasing in α (sweet spot α ≈ 0.5–0.8, A/a ≈ 10 → R_h/a ≈ 8).

  So anti-expansion is *not* purely the geometric quartic (the over-attribution this corrects): the quartic only sizes a bound object that the SU(3) topology/confinement is responsible for keeping localized at all — and the quartic-only width law is still unverified (the C2 Skyrme sweep with hardening disperses to box-fill, `[[project_c2_skyrme_no_confinement]]`). Discriminating test and the one axiom that could flip the verdict: open derivation C3 (`status.md`).

**(c) Finite-energy boundary conditions.** Establish finite-energy boundary conditions consistent with prestretch. The Skyrme twist F(r): [0,∞] → [π,0] gives winding B=1; the boundary conditions must be compatible with the prestressed vacuum (u → 0 as r → ∞ in the held frame, not the stress-free frame).