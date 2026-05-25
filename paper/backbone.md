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

8. **Lab-frame anisotropy is real and acknowledged; inside-observer isotropy is conjectured (dual-observer split)**
   - The lab (outside) observer knows the lattice and sees direction-dependent
     wave speeds. On the minimal 6-neighbor axial-only stencil this is roughly
     `~7.5%` at `α = 0.2` for the longitudinal speed between `[100]` and `[111]`.
     This anisotropy is **not** retuned away; it is the structural feature that
     provides the basis for the emergent SU(3) gauge sector (see #16, #19).
   - The inside (soliton) observer is built from substrate excitations and is
     conjectured to be blind to this anisotropy because their rods, clocks, and
     signal cones renormalize coherently with the local wave speed. Their
     measured "speed of light" is therefore direction-independent in their own
     units. This is the **load-bearing assumption** for emergent Lorentz
     kinematics in this project.
   - Decisive empirical test: Sprint 3 subtask 11 (observer universality) in
     `paper/validation_roadmap.md`. If two solitons of different polarizations
     or orientations report different effective `g^eff_{μν}` after their own
     emergent rods/clocks are accounted for, this assumption fails.

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

15. **Minimal lattice: 6-neighbor axial-only**
    - The canonical substrate is the 6-neighbor cubic lattice (each node
      connected to its six nearest axial neighbors `±êᵢ`). Diagonal shells
      (face-diagonal, body-diagonal) are intentionally absent.
    - On this stencil the long-wavelength dispersion is cubically anisotropic:
      e.g. `c_L([100]) / c_L([111]) ≈ 1.075` at `α = 0.2`. This anisotropy is
      not retuned away — see #8 and #19 for the dual-observer framework that
      makes it load-bearing for the gauge sector rather than a defect.
    - Convention (matches code and derivations): `α := rest_length / spacing`,
      so `α = 1` ↔ no prestress and `α = 0` ↔ maximum prestress (rest length
      zero). The default operating point is `α = 0.2`.

16. **Prestress `α` and the U(1)³ → U(3) crossover (Mechanism ii)**
    - On the 6-neighbor axial-only lattice the dynamical matrix `D(k)` is
      *diagonal* in the axis-aligned basis at every `α`: there is no linear
      off-diagonal coupling between the three lateral channels `(ξ¹,ξ²,ξ³)`.
      What runs with `α` is the *eigenvalue spread* of `D(k)`, driven by the
      `(1 − α) |Δu|²` geometric stiffness from prestress.
    - At `α = 1` (no prestress) the energy is purely Hookean longitudinal:
      `D(k)` has maximally split eigenvalues and the transverse modes are at
      zero frequency. Channels are decoupled and the triplet is far from
      degenerate — closest to a U(1)³ structure.
    - At `α = 0` (maximum prestress) the energy is purely geometric and
      isotropic: `D(k) ∝ I`. The three lateral channels are fully degenerate
      and the Wilczek–Zee gauge group on the 3-dimensional retained subspace
      is the full `U(3) = U(1) × SU(3)`. The current default `α = 0.2` sits
      close to this end.
    - **What the cubic anisotropy provides is not the *existence* of U(3),
      which is a generic representation-theoretic fact about any 3-dim
      complex vector space.** It provides the *physically meaningful
      decomposition*: a preferred 3-dim subspace (the lateral triplet,
      distinct from `X⁴`), a preferred axis-aligned basis within it (which
      makes the SU(3) traceless generators carry operational "color"
      meaning), and the natural identification of the `U(1)` trace as the
      EM channel.
    - Diagnostic caveat: the Berry/WZ curvature on the *real* eigenframe of
      `D(k)` is identically zero (real symmetric → trivial real bundle). The
      meaningful gauge object lives on the **per-wavepacket complex envelope**
      `Ψ ∈ ℂ³` (per #18), not on the BZ link-variable construction in
      paper §5.6 as previously formulated. See
      `test-runs/sprint2_subtask8_u3_decomposition/` for the detailed result.

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

19. **Sector mapping: lab anisotropy is the structural source of three force
    sectors**
    - The 3-dim lateral triplet `(ξ¹,ξ²,ξ³)` carries the U(3) gauge sector
      of the per-wavepacket complex envelope. Its decomposition splits the
      Standard-Model-like sectors:
      - **U(1) trace ↔ electromagnetism.** The coherent in-phase oscillation
        of all three axis components, naturally identified with the EM
        channel.
      - **SU(3) traceless ↔ strong / colour.** Phase and amplitude
        relationships *between* the three axis components. "Colour" of a
        soliton is its axis-alignment, made operationally meaningful by the
        cubic anisotropy of the lab frame (#15).
    - The 4th embedding direction `X⁴` (transverse "amplitude") is **not** a
      fourth colour; it is the **gravity channel**. Transverse-amplitude
      gradients source in-brane strain via the Pythagorean `∂_i u ∂_j u`
      coupling. Gravity is therefore structurally distinct from gauge — the
      gauge sector lives in the lateral triplet, the gravity sector lives in
      the embedding-amplitude direction.
    - Inside (soliton) observers experience all three sectors as standard
      effective field theories: Lorentz-invariant EM with U(1) gauge,
      colour-charged matter with SU(3) gauge, and a long-range scalar/tensor
      gravitational potential. Lab-frame anisotropy is invisible to them
      under #8.
    - This mapping is what the validation roadmap tests sector by sector:
      Sprint 2 (gauge sector / SU(3)), Sprint 3 (observer universality /
      effective Lorentz), Sprint 4 (solitons / colour-charged matter),
      Sprint 5 (gravity channel).
