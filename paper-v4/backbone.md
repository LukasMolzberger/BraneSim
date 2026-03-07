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

11. **Simulation priority: decisive tests over elegant prose**
    - The next decisive numerical target is a baryon-like triplet mode (proton/neutron sector), not another purely electron-centered attempt.
    - Reason: the cubic substrate naturally singles out a three-channel internal structure.

12. **Dimensionality-agnostic core implementation**
    - Core solver, forces, initialization interfaces, and diagnostics should work for 1D, 2D, and 3D where possible.
    - Dimension-specific code is acceptable in geometry setup, visualization, and experiment orchestration only.
