# Matter Bridge — Problem Statement

Show that the substrate supports stable localized solitonic modes. The eigenproblem formulation is in place.

## Specific tasks

**(a) Existence of at least one stable localized solution.** Prove existence of at least one stable localized solution to the nonlinear substrate equations (with periodic boundary conditions and fixed prestretch). The decisive first test is the hedgehog/Skyrme-twisted baryon-like triplet with angular structure ξ^i = f(r) x̂^i (J=0, L=1 VSH). Periodic boundary conditions are mandatory: open/free BC lets the prestressed lattice relax (contract to αa) and dominates everything, collapsing even the A=0 vacuum.

**(b) Derrick's theorem.** Address Derrick's theorem: why does the soliton not collapse, and why does it not expand?

The two instability directions are blocked by different mechanisms:

- **Anti-collapse (UV):** Lattice discreteness. Derrick's scaling argument assumes continuous rescaling x → λx for all λ > 0. On a discrete lattice with spacing a, λ is bounded below: the soliton cannot shrink past one cell. A configuration with energy fully concentrated at a single node is not a local energy minimum — the spring network immediately disperses it. The lattice spacing a is the physical UV cutoff; no quartic term is needed or relevant for the collapse direction.

- **Anti-expansion (IR):** Geometric quartic stabilization ∝ k_s α/a. Under Derrick rescaling in 3D: W₂ → λ E₂ (gradient, favors expansion) and W₄ → λ⁻¹ E₄ (quartic, favors confinement). Balance at λ* = √(E₄/E₂) gives the equilibrium soliton size. Make this rigorous: the Derrick balance gives R_h/a = κ (A/a) √(α/(1−α)), linear in amplitude A and increasing in α. Sweet spot: α ≈ 0.5–0.8, A/a ≈ 10 → R_h/a ≈ 8.

**(c) Finite-energy boundary conditions.** Establish finite-energy boundary conditions consistent with prestretch. The Skyrme twist F(r): [0,∞] → [π,0] gives winding B=1; the boundary conditions must be compatible with the prestressed vacuum (u → 0 as r → ∞ in the held frame, not the stress-free frame).