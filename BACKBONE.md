# backbone.md

This file is the short, non-negotiable backbone of the project.
It is intended to prevent core assumptions from being lost during paper edits, code changes, or simulation planning.

## Canonical project backbone

1. **Single substrate, substrate-only evolution**
   - Fundamental ontology: a **4D brane lattice embedded in 4D Euclidean
     ambient** (codimension 0). The "amplitude" direction of older 3D-in-4D
     formulations is absorbed into the time direction; the lattice IS the
     substrate. The four directions of the ambient are geometrically
     equivalent; their distinguished roles (time vs space, gauge vs gravity)
     are properties of the brane action and of the inside observer's frame,
     not of the ambient (see #21, #22).
   - The full 4D world-volume is a stationary point of the brane action
     `S[R]`; the foundational equations of motion come from the brane's
     elastic energy only. No external fields, no back-reaction from emergent
     diagnostics.
   - The current Verlet pipeline produces **Cauchy slices** through this 4D
     world-volume; it remains valid as a forward-evolution diagnostic of
     specific initial-condition problems, but it is not the foundational
     solver. A time-symmetric solver that finds full 4D stationary
     configurations under past+future boundary data is a separate
     development track (see project memory).

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
    - **Notation:** `α` is the **single 4D prestress** — one parameter governing all
      four lattice directions, spatial and temporal alike (rest length = `α ×` held
      spacing on every link; A4a, adopted 2026-06-05). It carries **no** spatial
      subscript. `α_t` appears only *provisionally* inside the A4a verification (testing
      `α_t = α`); a permanent split into `α_s`/`α_t` would be introduced only if that
      verification fails.

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
    - The induced-metric correction `∂_i u ∂_j u` gives a `(|∇u|²)²`
      term that scales as `λ^{+1}` under Derrick scaling, balancing the
      `λ^{−1}` quadratic-gradient term. This defeats Derrick collapse at the
      continuum level without any added field.
    - **Coefficient α-scaling (lattice-exact).** The quartic prefactor is
      `∝ k_s α/a` (vanishes at α→0), from the norm term `−k_s αa|ΔR|` of the
      central-force spring; a naive StVK identification `(μ/4ℓ₀⁴)` locks it to
      `μ ∝ (1−α)` and **inverts** the α-scaling. The Derrick `λ`-scaling above is
      α-independent and unaffected — only the *coefficient* was wrong. Radius
      grows with α (`R_h/a = κ(A/a)√(α/(1−α))`). See
      `paper/derivations/geometric_nonlinearity_alpha_scaling.md`.
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
    - **Origin of the complex structure (the `i`):** the `i` that promotes the
      real lateral triplet to `Ψ ∈ ℂ³` is not fundamental and is not an
      analytic-signal trick on a slice — it originates in the **timelike**
      direction of the worldtube. The carrier phase is the rotation of the real
      field as one advances along the time axis; the two real DOF spanning `ℂ`
      are the field and its time-quadrature (`Ψ ≈ ξ + (i/ω₀)ξ̇`, with `i` the
      time-evolution generator `i∂_tΨ = H_eff Ψ`). A single spacelike slice is a
      real snapshot with no phase; the `U(1)` exists only across the time extent.
      This is **why** the emergent gauge (EM/`U(1)`) connection's base space is
      `(x,t)`, not the Brillouin zone (paper §5.6 / backbone #9), and why two
      adjacent time slices are the minimum data that fix it (the block-solver
      chiral BC, OPEN_PROBLEMS A2). Consistent with the model-(b) timelike spring
      (#21): time is a genuine direction carrying real structure.

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
      the embedding-amplitude direction. **Note (per #22):** in the 4D-in-4D
      ontology of #1, `X⁴` is the inside observer's chosen timelike
      direction within the 4D ambient, and the lateral triplet is the
      perpendicular 3-dim spacelike subspace. The lateral/amplitude split is
      observer-relative but globally consistent under cosmological boundary
      conditions; the descriptions in this item remain valid as the inside
      observer's view.
    - Inside (soliton) observers experience all three sectors as standard
      effective field theories: Lorentz-invariant EM with U(1) gauge,
      colour-charged matter with SU(3) gauge, and a long-range scalar/tensor
      gravitational potential. Lab-frame anisotropy is invisible to them
      under #8.
    - This mapping is what the validation roadmap tests sector by sector:
      Sprint 2 (gauge sector / SU(3)), Sprint 3 (observer universality /
      effective Lorentz), Sprint 4 (solitons / colour-charged matter),
      Sprint 5 (gravity channel).

20. **Soliton-layer description protocol: vector spherical harmonics**
    - The natural description language for the angular structure of a
      localized 3D mode on the lateral triplet `ξⁱ ∈ ℝ³` is the
      **vector spherical harmonic** basis `Yⁱ_{JLM}(θ,φ)` paired with a
      radial profile:

        `ξⁱ(r,θ,φ) = Σ_{J,L,M} f_{JLM}(r) · Yⁱ_{JLM}(θ,φ)`

    - **This is a description-layer commitment, not a structural one.**
      The substrate is cubic (`O_h`, per #15); it is **not** `SO(3)`-symmetric.
      What spherical harmonics describe is the *emergent* approximate
      `SO(3)` at long wavelength (confirmed by Sprint 1 dispersion
      measurements), which is the regime where solitons live by the
      width-≫-a requirement of #17. Particles emerge from the lower layers;
      spherical harmonics organize how we *describe and classify* them —
      they do not determine which states exist. A `(J, L, M)` channel
      that does not support a stable stationary configuration of the
      substrate dynamics is simply absent from the physical spectrum.
    - **Cubic-anisotropy refinement at soliton scale.** `O_h` corrections
      appear as small splittings of `SO(3)` multiplets. The `L=1` rep maps
      intact to the `T_{1u}` irrep of `O_h` (still 3-dim — the hedgehog
      survives without restriction). Higher-`L` modes split (e.g.
      `L=2 → E_g ⊕ T_{2g}`) and should be parameterized in cubic
      harmonics when lattice corrections matter.
    - **Combined particle labels at the soliton layer:**
      `(J, P; SU(3)-irrep; Q_U(1); B_winding)`.
      - `J` = total angular momentum (from emergent `SO(3)` at scales ≫ a).
      - `P` = parity (discrete part of the emergent rotation group).
      - `SU(3)-irrep` and `Q_U(1)` from the U(3) decomposition of the
        lateral triplet (per #16, #19).
      - `B_winding` = topological winding number (per #17).
    - **Canonical baryon ansatz: hedgehog `ξⁱ = f(r) x̂ⁱ`** (`J=0, L=1`).
      The internal index `i` is locked to the spatial angular direction,
      putting the entire nontrivial structure in the SU(3) traceless
      sector; the U(1) trace averages to zero angularly (trace-neutral
      far field). Topological stabilization (winding `B=1`) is achieved by
      extending the hedgehog into a Skyrme-twisted configuration in which
      the embedding-amplitude direction `X⁴` (per #19) also rotates
      radially. Proton-vs-neutron distinction enters as a `J=0` trace
      admixture on top of the hedgehog, not as an unrelated seed.
    - Concrete ansatz menu used in the simulation search is recorded in
      `paper/baryon_simulation_roadmap.md` Phase 2.

21. **Asymmetry lives in the brane, not the ambient**
    - The 4D ambient (per #1) is fully symmetric Euclidean — just a stage
      for measuring distances (Pythagoras). No preferred direction at the
      ambient level.
    - The asymmetry that picks out one of the four lattice directions as
      **time** is a property of the brane's elastic action: the action has
      Lorentzian sign structure (standard `T − V` Lagrangian, recast on the
      lattice as link-energy terms with opposite sign for the timelike
      direction relative to the three spacelike directions). This is what
      makes one of the four lattice directions "time".
    - All previously stated anisotropy claims (cubic anisotropy at lab level
      per #8 and #15, the `U(1)³ → U(3)` crossover per #16) are
      consequently claims about the brane's structure on its 3D spacelike
      slice. The full 4D lattice has 4D-cubic symmetry at the geometric
      level; the 3D-cubic structure those items rely on is what a Cauchy
      slice perpendicular to the timelike direction exhibits.
    - The Verlet implementation already encodes this asymmetry implicitly
      via the kinetic term `½ m v²`: this is the discrete form of a
      timelike-direction "link" with an opposite-sign contribution to the
      action. The 4D-in-4D reading makes that structure explicit instead of
      implicit. The explicit discrete 4D action (6 spacelike central-force
      springs entering `−V`, 2 temporal links entering `+T`) is derived in
      `paper/derivations/discrete_4d_brane_action.md`. Two results from it:
      (a) its Euler–Lagrange stencil is *term-for-term* Störmer–Verlet —
      Verlet is the discrete variational integrator of the action, so the
      forward IVP and the 4D block BVP share the identical local stencil and
      differ only in boundary conditions; (b) the long-wavelength cone speeds
      are `c_L² = k_s a²/m` and `c_T² = (1−α) k_s a²/m`, fixing the lattice
      light-cone via the temporal-to-spacelike stiffness ratio.
    - **Temporal-link form (decided 2026-06-05): model (b).** The timelike link is
      a genuine central-force spring `½ k_t (|ΔR| − r_t)²` with its own rest length
      `r_t ≠ 0` (temporal prestress `α_t`) — fully symmetric 4D-cubic. This exposes a
      real model parameter rather than hiding it, and supplies the time-link geometric
      quartic that #22's unified contraction (gravitational time dilation) requires.
      Model (a) (zero-rest-length kinetic, plain Newton, current code) is the `r_t = 0`
      *special case*; the exact dynamics is non-Newtonian, so the block solver is the
      foundational integrator and forward Verlet is its `r_t = 0` IVP limit. A
      **single** prestress `α` governs all four directions (each link's rest length
      = `α ×` its held spacing; adopted 2026-06-05), tying soliton binding and
      gravitational time dilation to one dial; consistency (light-cone isotropy +
      Newtonian limit) is the open verification, `OPEN_PROBLEMS.md` A4a.
      See the derivation §6 and `OPEN_PROBLEMS.md` A4.
      The Lorentzian action is a saddle (unbounded below), so a foundational
      block solver must root-find `∇S = 0`, not minimize `S` (OPEN_PROBLEMS §A).

22. **The gauge / gravity split is observer-relative**
    - The 4-component displacement field `δX^μ` at each node has no
      privileged structural decomposition at the lab level. The split into
      "3 lateral channels carrying U(3) gauge + 1 timelike channel
      carrying gravity" (per #19) is performed by the **inside observer**
      after they pick a timelike direction.
    - Globally, the inside observers' choice is constrained to a single
      direction by **cosmological boundary conditions**: the Big Bang
      vertex emanates along one of the four lattice directions, and that
      direction is what every inside observer ends up calling "time". The
      split in #19 is therefore globally consistent in our universe even
      though it is not structurally built into the substrate.
    - Different inside observers (related by emergent Lorentz
      transformations on their spacelike slice) may disagree locally on
      what counts as the "amplitude" direction within `δX^μ`, but they
      agree on the global timelike axis by virtue of riding the same
      cosmological history.
    - **Unified contraction.** The geometric-quartic mechanism (#17) now
      acts symmetrically across all four lattice directions: a displacement
      perpendicular to any link stretches that link, regardless of which
      direction the perpendicular is. Consequence: gravitational length
      contraction (timelike displacement stretching spacelike links) and
      gravitational time dilation (spacelike displacement / kinetic energy
      stretching the timelike link) are two faces of the same rule. The
      two faces of gravity unify under one mechanism. This rests on the
      **model-(b)** timelike spring (#21 / A4, decided 2026-06-05): only a time
      link with its own rest length `r_t ≠ 0` carries the geometric quartic that
      lets a spacelike displacement stretch the timelike link. With the single 4D
      prestress `α_t = α` (A4a, adopted), the strength of both faces is set by the
      *same* `α` that governs soliton binding — gravity and confinement share one dial.

23. **Bell's theorem and the retrocausal worldtube interpretation**
    - BraneSim is manifestly local and deterministic — a classical lattice
      with action-based dynamics. Bell rules out theories that are
      simultaneously local, deterministic, and measurement-independent. The
      cleanest loophole consistent with this project's substrate is
      **retrocausality**: future measurement context propagates back along
      the particle's worldtube.
    - **Causality is a property of solitons, not of the substrate.** The
      brane's action is time-symmetric (Lorentzian signature does not
      single out a direction of time, only a *kind* of direction). The
      arrow of time is selected by the **chirality of soliton solutions**:
      matter solitons are forward-propagating worldtubes; antimatter
      solitons are backward-propagating worldtubes (Feynman–Stueckelberg-
      consistent).
    - **Entanglement is V-branching of one worldtube.** Spacelike
      correlations between entangled pairs are continuity of one extended
      4D object whose worldtube branches at a V-vertex in the shared past
      — not "spooky action at a distance".
    - This places BraneSim in the same interpretive family as Aharonov's
      two-state-vector formalism, Cramer's transactional interpretation,
      Price & Wharton's retrocausal models, and Sutherland's
      time-symmetric Bohmian model.
    - **Open derivations** (flagged as future work; do **not** claim in the
      paper as established, and do **not** present in the theory-structure
      diagram as theory): Tsirelson's bound (`2√2`), the no-signalling
      theorem, and the baryon-to-photon ratio. The live list is tracked
      centrally in `OPEN_PROBLEMS.md` §B, not duplicated in the manuscript.

24. **Color is kinematically confined: no coherent colored free mode**
    - The linear spectrum carries **no colored asymptotic free state**. Two
      requirements are mutually exclusive at linear order:
      - *Coherence* (a lateral-triplet wavepacket whose three components stay
        phase-aligned as it propagates) requires the three lateral branches of
        `D(k)` to be degenerate. On the 6-neighbor stencil this holds **only**
        on the body-diagonal `k ∥ [111]`, where the directional anisotropy
        `g(k̂)=0` (per #15, #16).
      - *Color-activity* (operationally meaningful `SU(3)`-traceless content; a
        finite fibre-internal holonomy) requires the three branches to be
        **non-degenerate**, i.e. `g(k̂) ≠ 0`, which holds only **off** `[111]`.
        Exactly on `[111]` the color holonomy is undefined (degenerate triplet,
        no preferred basis; the D2 ratio `R → ∞`).
    - Hence **no linear propagation direction is simultaneously coherent and
      color-active**: off `[111]` the unequal branch speeds dephase the triplet,
      on `[111]` there is no resolvable color. Color charge therefore cannot be
      carried to infinity by a free linear excitation. It has support only in a
      **nonlinear, phase-locked bound state** (a soliton), where the three
      components are held together by the geometric coupling rather than by
      spectral degeneracy. This is **kinematic confinement of color**.
    - **Scope caveat (load-bearing).** This is a statement about *asymptotic
      states* (there is no colored free particle), **not** a dynamical
      derivation of a linear potential / string tension / area law — those
      remain open (paper §2 non-claims). It is consistent with the identically
      vanishing k-space Berry/WZ curvature ∀α (#16): color is not a curvature on
      the Brillouin zone but a fibre-internal `(x,t)` object that only closes
      around a localized worldtube.
    - **Quantitative, falsifiable handle.** The off-`[111]` fibre-holonomy ratio
      scales as `R(α₂)/R(α₁) = [(3−2α₂)/α₂] / [(3−2α₁)/α₁]`; for `(α₁,α₂)=(0.2,0.5)`
      this is `0.3077` (`OPEN_PROBLEMS.md` D2;
      `paper/derivations/alpha_holonomy_estimator.md`). A `>10%` deviation
      falsifies the spectral-susceptibility factorization behind this picture.
