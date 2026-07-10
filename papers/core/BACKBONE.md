# BACKBONE.md

This file is the short, non-negotiable backbone of the project.
It is intended to prevent core assumptions from being lost during paper edits, code changes, or simulation planning.

## Canonical project backbone

1. **Single substrate, substrate-only evolution**
   - Fundamental ontology: a **static 4D hypercubic (quartic) brane lattice
     embedded in 4D Euclidean ambient** (codimension 0). The "amplitude"
     direction of older 3D-in-4D formulations is absorbed into the time
     direction; the lattice IS the substrate. The four directions of the ambient
     are geometrically equivalent, and the four internal lattice directions are
     treated symmetrically at the stencil/Bloch level. **Time is one of the four
     internal lattice directions, not an external evolution parameter.** Their
     distinguished roles (time vs space, gauge vs gravity) are properties of the
     brane action and of the inside observer's frame, not of the ambient.
   - The full 4D world-volume is a stationary point of the brane action
     `S[R]`; the foundational equations of motion come from the brane's
     elastic energy only. No external fields, no back-reaction from emergent
     diagnostics.

2. **Particles are not point-like**
   - Electrons, protons, neutrons, and other particles are modeled as extended standing-wave / solitonic patterns.
   - Point-like appearances are interpreted as interaction and localization effects, not as literal ontology.

3. **Primary nonlinearity is geometric**
   - The essential nonlinearity comes from full Euclidean distances in the 4D embedding and, in the continuum, from the induced metric.
   - Even if a spring law is linear in scalar extension, the map from coordinates to extension is nonlinear.
   - Do not replace this insight with claims that a separate ad hoc nonlinear force law is required.
   - **Corollary (no false no-go):** results obtained by linearizing this
     nonlinearity away — e.g. the bare-vacuum Bloch/Hessian — describe only the
     quadratic tangent operator about a straight background. They must never be
     stated as no-go theorems for the full theory. The relevant object is the
     fluctuation operator around a finite-amplitude nonlinear background, where the
     Pythagorean square root makes the Hessian a background-dependent frame tensor.
     Any negative conclusion that depends on dropping the Euclidean distance is an
     artifact of the approximation, not a property of the substrate. See
     LESSONS_LEARNED #3.

4. **Pure geometric amplitude--lateral coupling**
   - Transverse amplitude couples to lateral contraction through geometry.
   - No additional coupling field is introduced by default.
   - Gravity-like behavior is expected to emerge from this contraction channel.

5. **No artificial clamps, cutoffs, or hand-imposed thresholds**
   - No hard amplitude clamps.
   - No piecewise “if energy > … then …” confinement rules.
   - Any threshold-like behavior must emerge smoothly from substrate dynamics.

6. **Prestretch parameters are physically important (α_s, α_t, γ_t)**
   - The prestretch is not cosmetic. It is the mismatch between held spacing and
     stress-free length, and it controls inter-axis coupling / branch mixing.
   - **Sign and magnitude are separate.** `η_μ = (−1,−1,−1,+1)` carries only the
     *sign* (three spacelike one way, time-like opposite); the *magnitudes* of the
     directional prestress need not be equal in space and time. A single `α` was
     overloaded (magnitude, branch mixing, signature sign, temporal coupling all at
     once) and is now split:
     `r_i = α_s a` (i=1,2,3), `r_4 = α_t a`, with independent link stiffnesses
     `κ_i = κ_s`, `κ_4 = κ_t`. The signed prestress is
     `ρ_μ ~ η_μ κ_μ a (1 − α_μ)`. The **minimal dimensionless control set** is
     `{α_s, α_t, γ_t}` with `γ_t = κ_t/κ_s`; held spacing stays `a_s = a_t = a` to
     preserve the quartic/hypercubic substrate.
   - **α = 1 is the zero-mismatch (no-prestress) limit, not maximum prestress** —
     the transverse tangent stiffness is `∝ (1−α)`, so it *vanishes* at α = 1.
   - **`α_t ≠ 1` is required, and `α_t ≉ 0`.** `α_t = 1` removes the temporal link's
     transverse stiffness (`1−α_t = 0`), over-decoupling the time/amplitude direction
     from spatial geometry (kills the gravity contraction channel and the
     nonlinear-core mixing). `α_t ≈ 0` reverts to the old kinetic-limit notation
     where the time-like sign is no longer clearly a prestress. The working regime is
     `0 < α_t < 1`, likely `α_t = 1 − ε_t` with `0 < ε_t ≪ 1`. This replaces the
     earlier `r_t ≈ 0` kinetic limit and resolves that tension (T12).
   - **The Lorentzian signature is the sign pattern of the directional prestress,
     not an independent postulate** (`α → prestress state → ρ_μ = η_μ κ_μ a(1−α_μ)
     → Lorentzian-sign brane action`). Signature is *derived* from prestress — not
     the ambient metric, and not a special treatment of time.
   - **The fourth direction plays a dual role:** time-like internal direction in the
     Lorentzian action, and — viewed from a 3D spatial sub-slice — the transverse
     "amplitude" direction whose `∂_i X⁴ ∂_j X⁴` term drives the gravity channel.
     This dual reading is load-bearing and is why `α_t` must keep the fourth
     direction coupled (see #4).

7. **Hypercubic (quartic) lattice is an ontic microphysical hypothesis**
   - Fine-grained rotational symmetry is not assumed.
   - Axis alignment is expected to matter microscopically and is the natural source of the triplet/color idea.

8. **Lab-frame anisotropy is real and acknowledged; inside-observer isotropy is conjectured (dual-observer split)**
   - The lab (outside) observer knows the lattice and sees direction-dependent
     wave speeds. 
     This anisotropy is **not** retuned away; it is the structural feature that
     provides the basis for the emergent SU(3) gauge sector.
   - The inside (soliton) observer is built from substrate excitations and is
     conjectured to be blind to this anisotropy because their rods, clocks, and
     signal cones renormalize coherently with the local wave speed. Their
     measured "speed of light" is therefore direction-independent in their own
     units. This is the **load-bearing assumption** for emergent Lorentz
     kinematics in this project.

9. **Gauge structure is emergent holonomy, not a fundamental extra field**
   - Berry and Wilczek--Zee connections are interpreted as effective gauge structure of an ordered carrier sector.
   - Axis-wise Berry phases are only meaningful in the decoupled limit.
   - With mixing, the correct invariant object is the non-Abelian transported subspace.
   - **The carrier is a rank-3 internal space over a 4D base, giving U(3), not
     SU(4).** The connection is `𝒜_μ^a_b` with external Bloch/brane index
     `μ = 1..4` and internal carrier index `a,b = 1..3`; four-dimensionality of the
     base does *not* enlarge the gauge group. The rank-3 carrier is the complexified
     transverse fluctuation bundle of a nonlinear background (the moving `T⊥` of the
     link direction `Q̂`), **not** the fixed spatial axes. The fourth (time/amplitude)
     polarization mode is spectrally separated out of the near-degenerate triplet by
     `α_t` (and `γ_t`) — it stays in the substrate but leaves the carrier.
   - **Interpretation A is selected:** the trace of the U(3) connection is the
     electromagnetic U(1) (`a_μ = ⅓ Tr 𝒜_μ → f_μν`), the traceless part is the
     color SU(3) (`B_μ, Tr B_μ = 0 → 𝒢_μν`). This is a candidate identification
     until the Maxwell/Yang–Mills dynamics (T3/T4) are derived.
   - **Color is a universal vacuum gauge sector**, not a matter-bound holonomy: the
     U(1)×SU(3) structure exists over empty space, carried by the universal nonlinear
     periodic vacuum carrier (the gauge background `R̄`, *not* the bare straight
     lattice). Soliton/matter structure enters only at a much higher layer (#2,
     Layer 5), carrying color charge into the already-universal field; matter does
     not create color. Local nonlinear strain sets the SU(3) *coupling strength*
     (strong/short-range near matter), not color's existence.

10. **Coherence of theory has priority over shortness of paper**
    - Necessary conceptual bridges should stay if they are required to make the theory work.
    - Redundancy should be reduced, but not at the cost of breaking the explanatory chain.

11. **The paper must present only the final theory, not the path of discussion**
    - Discussions between user and assistant, exploratory back-and-forth, and abandoned intermediate ideas are not part of the manuscript.
    - The paper should present one internally consistent formulation of the theory, not a reconstruction of how that formulation was reached.

12. **Fundamental dynamics is deterministic; QM probability is emergent and threshold-triggered**
    - The substrate is **fully deterministic and ontic**. The only law is
      stationarity of the brane action (`∇S = 0` / the Verlet stencil); node
      trajectories `R_n^A(t)` are single-valued with no fundamental randomness,
      no wavefunction collapse, and no intrinsic probability.
    - The observed **probabilistic behavior of quantum mechanics is emergent, not
      fundamental**. It arises when the continuous, deterministic substrate state
      is read out through **classical thresholds**: a detection/localization event
      fires when a deterministic substrate quantity (e.g. accumulated energy,
      amplitude, or strain in a soliton–detector interaction) crosses a threshold.
    - Apparent randomness is **epistemic**, from unresolved/inaccessible substrate
      microstate and the sensitive, near-critical nature of the threshold trigger —
      not from any indeterminacy in the underlying law. The Born rule is a
      statistical regularity over these deterministic threshold crossings, to be
      *derived*, not postulated.
    - Discrete spectra likewise come from confined solitonic standing-wave modes,
      not from a quantization axiom (see #2). No hidden extra stochastic force is
      introduced; the trigger threshold must itself emerge smoothly from substrate
      dynamics (consistent with #5, no artificial clamps/cutoffs).
    - **Bell-theorem consistency** (`derivations/FOUNDATIONS_bell.md`): a deterministic,
      locally-mediated substrate must relax *measurement independence* to be
      compatible with quantum correlations. Here that is not a past conspiracy but the
      ordinary consequence of the all-at-once, time-symmetric variational formulation
      (Layer 0 saddle + periodic time BC) — the locally-mediated / time-symmetric class
      of Wharton & Argaman (RMP 92, 021002, 2020). This is a *non-exclusion* argument;
      actually reproducing Bell correlations and the Born rule remains the owed program
      above.

13. **Quantization originates in matter confinement; radiation inherits it**
    - Only **self-confined spin-1/2 solitons are intrinsically quantized**. Their
      discrete spectrum is the set of **confinement-equilibrium** standing-wave
      modes: the self-trapping balance selects discrete stable amplitudes and
      frequencies (rest mass, charge, spin). Quantization is a property of the
      bound equilibrium, not of the substrate field.
    - The free carried field is continuous and **not** intrinsically quantized.
      **Photons are inherently non-quantized** substrate excitations —
      continuous-amplitude waves with no fundamental quantum.
    - Photon quantization is **emergent and inherited**: it comes from interactions
      with spin-1/2 particles. A photon appears quantized only because its
      emission/absorption is gated by the discrete confinement-equilibrium
      transitions of the spin-1/2 solitons it couples to (consistent with the
      threshold-triggered readout of #12).
    - Consequence: energy exchanged in `hν`-like packets reflects the discrete
      level structure of the matter soliton, not a fundamental granularity of
      radiation.

14. **Dimensionality-agnostic core implementation**
    - Core solver, forces, initialization interfaces, and diagnostics should work for 1D, 2D, and 3D where possible.
    - Dimension-specific code is acceptable in geometry setup, visualization, and experiment orchestration only.

15. **English only policy**
    - Never ever use any other language than English in this project!

