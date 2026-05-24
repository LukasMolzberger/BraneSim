# backbone.md

This file is the short, non-negotiable backbone of the project.
It is intended to prevent core assumptions from being lost during paper edits, code changes, or simulation planning.

## Canonical project backbone

1. **Single substrate, substrate-only evolution**
   - Fundamental ontology: a 3D brane lattice embedded in 4D.
   - Evolve only microscopic substrate degrees of freedom: node positions and velocities.
   - Forces must come from the substrate energy only.

2. **Particles are not point-like**
   - Electrons, protons, neutrons, and other particles are modeled as extended standing-wave / solitonic patterns.
   - Point-like appearances are interpreted as interaction and localization effects, not as literal ontology.

3. **Primary nonlinearity is geometric**
   - The essential nonlinearity comes from full Euclidean distances in the 4D embedding and, in the continuum, from the induced metric.
   - Even if a spring law is linear in scalar extension, the map from coordinates to extension is nonlinear.
   - Do not replace this insight with claims that a separate ad hoc nonlinear force law is required.

4. **Pure geometric amplitude--lateral coupling**
   - Transverse amplitude couples to lateral contraction through geometry.
   - No additional coupling field is introduced by default.
   - Gravity-like behavior is expected to emerge from this contraction channel.

5. **No artificial clamps, cutoffs, or hand-imposed thresholds**
   - No hard amplitude clamps.
   - No piecewise “if energy > … then …” confinement rules.
   - Any threshold-like behavior must emerge smoothly from substrate dynamics.

6. **Prestretch parameter alpha is physically important**
   - The prestretch parameter is not cosmetic.
   - It controls the mismatch between held spacing and stress-free length.
   - It also controls inter-axis coupling / branch mixing in the cubic substrate picture.

7. **Cubic lattice is an ontic microphysical hypothesis**
   - Fine-grained rotational symmetry is not assumed.
   - Axis alignment is expected to matter microscopically and is the natural source of the triplet/color idea.

8. **Macroscopic isotropy must nevertheless emerge operationally**
   - The coarse-grained universe appears rotationally symmetric.
   - This can only be acceptable if rods, clocks, and signal propagation built from the substrate are renormalized by the same effective metric.
   - Residual anisotropy / birefringence is therefore a decisive diagnostic, not a detail to ignore.

9. **Gauge structure is emergent holonomy, not a fundamental extra field**
   - Berry and Wilczek--Zee connections are interpreted as effective gauge structure of an ordered carrier sector.
   - Axis-wise Berry phases are only meaningful in the decoupled limit.
   - With mixing, the correct invariant object is the non-Abelian transported subspace.

10. **Coherence of theory has priority over shortness of paper**
    - Necessary conceptual bridges should stay if they are required to make the theory work.
    - Redundancy should be reduced, but not at the cost of breaking the explanatory chain.

11. **The paper must present only the final theory, not the path of discussion**
    - Discussions between user and assistant, exploratory back-and-forth, and abandoned intermediate ideas are not part of the manuscript.
    - The paper should present one internally consistent formulation of the theory, not a reconstruction of how that formulation was reached.

12. **Simulation priority: decisive tests over elegant prose**
    - The next decisive numerical target is a baryon-like triplet mode (proton/neutron sector), not another purely electron-centered attempt.
    - Reason: the cubic substrate naturally singles out a three-channel internal structure.

13. **Dimensionality-agnostic core implementation**
    - Core solver, forces, initialization interfaces, and diagnostics should work for 1D, 2D, and 3D where possible.
    - Dimension-specific code is acceptable in geometry setup, visualization, and experiment orchestration only.

14. **English only policy**
    - Never ever use any other language than English in this project!

15. **Cauchy relation ≠ cubic isotropy (the difference matters)**
    - Central-force pair springs in a stress-free Bravais lattice automatically
      give the Cauchy relation `C_{1122} = C_{1212}` (`C_{12} = C_{44}` in
      Voigt notation). This reduces the cubic elastic tensor from 3 to 2
      independent constants.
    - This does **not** imply cubic isotropy. Isotropy requires the additional
      condition `C_{1111} − C_{1122} = 2 C_{1212}`, which depends on the shell
      weights. The `1/|δ|²` choice currently in `components/simulation/grid.py`
      does NOT satisfy this condition: it produces a leading-order static
      cubic anisotropy of roughly `−21%` at `α = 1` and `−2%` at `α = 0.2`.
    - See `paper-v4/derivations/lattice_to_continuum.md` for the closed form.
      The anisotropy is leading-order, not `O((ka)⁴)`. Shell-weight retuning
      to recover isotropy is an open project decision.

16. **Prestress `α` runs the U(1)³ → U(3) crossover**
    - At `α = 1` the three lateral channels decouple at the linear level
      → gauge group of the narrowband triplet is `U(1)³`.
    - At `α < 1` the diagonal-shell springs activate linear off-axis
      coupling → gauge group is the full `U(3) = U(1) × SU(3)`.
    - This is the operational meaning of `α` for the gauge sector and is
      directly testable by sweeping `α` and measuring the SU(3)-content of
      the Wilczek–Zee curvature.

17. **Geometric quartic provides Skyrme-class soliton stabilization**
    - The induced-metric correction `∂_i u ∂_j u` in StVK gives a `(|∇u|²)²`
      term that scales as `λ^{+1}` under Derrick scaling, balancing the
      `λ^{−1}` quadratic-gradient term. This defeats Derrick collapse at the
      continuum level without any added field.
    - **Target soliton width ≫ a** (Skyrme-stabilized regime), not `~ a`
      (lattice-stabilized regime). The latter has Peierls–Nabarro pinning
      that breaks emergent Lorentz invariance.
    - Absolute stability against unwinding to vacuum may still require
      topology (winding number, Hopf charge); existence + perturbative
      stability is not the same as topological stability.

18. **Narrowband is local per wavepacket, not a global postulate**
    - The complex-envelope description is required to promote the real
      lateral triplet to `Ψ ∈ ℂ³`. It is applied per band-isolated
      excitation, never as a globally coherent universe-wide carrier.
    - Therefore "narrowband ordered sector" as previously stated in v3 is
      replaced by "per-wavepacket band isolation". This removes the
      dependence on an unmotivated global oscillation.
