# Matter Bridge — Problem Statement

Show that the substrate supports stable localized solitonic modes. The eigenproblem formulation is in place.

## Specific tasks

**(a) Existence of at least one stable localized solution.** Prove existence of at least one stable localized solution to the nonlinear substrate equations (with periodic boundary conditions and fixed prestretch). The decisive first test is the hedgehog/Skyrme-twisted baryon-like triplet with angular structure ξ^i = f(r) x̂^i (J=0, L=1 VSH). Periodic boundary conditions are mandatory: open/free BC lets the prestressed lattice relax (contract to αa) and dominates everything, collapsing even the A=0 vacuum.

**(b) Derrick's theorem.** Address Derrick's theorem: why does the soliton not collapse? The answer must be the geometric quartic stabilization coefficient ∝ k_s α/a, which is the only term that prevents collapse (W₂ → λ⁺¹ E₂, W₄ → λ⁻¹ E₄ under Derrick rescaling; balance at λ* = √(E₄/E₂)). Make this rigorous: the Derrick balance gives size R_h/a = κ (A/a) √(α/(1−α)), linear in amplitude A and increasing in α. Sweet spot: α ≈ 0.5–0.8, A/a ≈ 10 → R_h/a ≈ 8.

**(c) Finite-energy boundary conditions.** Establish finite-energy boundary conditions consistent with prestretch. The Skyrme twist F(r): [0,∞] → [π,0] gives winding B=1; the boundary conditions must be compatible with the prestressed vacuum (u → 0 as r → ∞ in the held frame, not the stress-free frame).