> **Archived 2026-06-28:** this programme depended on invalid field-space and topology identifications.

# Confinement (Paper VI) — Problem Statement

Confinement is the central unsolved problem of the brane-substrate program. The pieces are
established across Papers I, III, IV but no single document gives a closed-form account of *why*
the only finite-energy, non-radiating excitations carrying `U(1)` charge and/or `SU(3)` colour are
**bound topological defects of the single `ℂ³` carrier envelope `Ψ=(ψ_x,ψ_y,ψ_z)`**, regulated by
the lattice. This paper assembles and extends those pieces into one closed treatment.

## The one field, the one nonlinearity

- EM/charge = the trace direction `Ψ_tr=(ψ_x+ψ_y+ψ_z)/√3` (a `π₁(U(1))=ℤ` line vortex).
- Colour = the traceless `su(3)` part `Ψ_⊥` (a `π₃(SU(3))=ℤ` point texture).
- The **only** anharmonicity is the Euclidean-norm term `−k_s α a |ΔR|`, exactly `∝α`. Everything
  below must vanish as `α→0` (the exactly-linear limit).

## Organizing theme: transitions, and the connection/holonomy spine

The paper is organized around **transitions** and the differential-geometric machinery that
describes them (`connections_holonomy.md`, `scale_transition.md`):

- **`U(1)↔SU(3)` mediated by the prestress `α`** — a *connection-reducibility* crossover (Abelian
  `U(1)³` ↔ non-Abelian `U(3)`), with the cross-vertex `g∝α` as the off-diagonal block.
- **Discrete lattice → continuous QM field** — a *nonlinearity along the scale axis*: the
  carrier–envelope coarse-graining (where the complex `i`/`U(1)` is born), the scale-dependent
  geometric quartic, and the emergence of quantization at the transition.
- **`SH↔Cartesian` curvilinear change of basis** — **moving-frame covariantization** (not a new
  constitutive nonlinearity): the connection coefficients (Christoffel/Levi-Civita) are exactly the
  angular-momentum barriers in the radial ODEs.

The unifying claim: every structure is a **connection** on the single `ℂ³` bundle (Levi-Civita /
Berry / Wilczek–Zee), and **confinement = constraints on holonomies**, regulated by the lattice.

## Four layers, then the confinement question (the paper's order)

Build the apparatus in four layers before asking the confinement question:

1. **Microscopic constitutive layer** — one norm law, one prestress `α`, the single-anharmonicity
   theorem (`∝α`, transverse, P-even).
2. **Coarse-graining layer** — the real lattice field becomes a complex carrier–envelope field only
   after temporal coarse-graining (where the complex `U(1)` description first becomes available);
   the bundle structure for Berry/WZ transport appears here. Emergent, not eliminativist.
3. **Reducibility layer** — `α` simultaneously controls branch splitting, triplet degeneracy, and
   the quartic trace↔traceless coupling: the `U(1)↔SU(3)` transition. Note the tension: `α→0`
   degenerates the spectrum but kills the coupling, so the physical regime is the **intermediate
   window `α≈0.5–0.8`** (a sweet spot that recurs across the notes — a strong internal-consistency
   signal).
4. **Moving-frame layer** — the defect-adapted reduction to cylindrical/spherical coordinates as
   **covariant transport** in a rotating frame, Christoffel terms furnishing the angular barriers.

Then the confinement answer splits into **three propositions**:

- **P1 (negative, closed).** Spatial modulus-only anharmonicity does **not** bind trace to traceless
  — the strongest result.
- **P2 (kinematic, closed).** No free colored asymptotic state (coherence ⊥ colour-activity, now
  with a closed dephasing length and an `O(1)` conjugacy bound) + `ℤ₃` triality lock.
- **P3 (positive, conditional).** Time-link/worldtube closure can flip the cross-vertex to attractive
  when `χ>1` — the conjectural completion, tied to the two-time boundary-value problem.

## Questions this paper must answer in closed form

1. **Where do the nonlinearities live, and how do they couple the two gauges?** (Norm term →
   transverse-only anharmonicity → trace↔traceless cross-vertex → `α`-reducibility of the internal
   connection.)
2. **The `U(1)` defect.** Closed core structure of the centerline (codim-2 in 3-D, a 2-D worldsheet
   in 4-D): radial profile, core energy, and the analytic statement that the **lattice forbids a
   true singularity** (the would-be `1/ρ` gradient saturates at `~1/a`).
3. **The `SU(3)` defect.** The point texture (codim-3 worldline in 4-D): why the natural language is
   the winding integral, not spherical harmonics; the combined-`SO(3)` hedgehog reduction; the
   Derrick size and lattice anti-collapse.
4. **The interaction.** Does the substrate confine charge to colour spatially? (Honest answer:
   spatial-only is repulsive; the binder is the time link / worldtube closure. Only the `ℤ₃`
   triality lock is a closed spatial statement, and it locks representation, not position.)
5. **Spin-½ and chirality.** The `ℤ₂` holonomy at the soliton layer; chirality as the arrow of time
   along the defect.

## What "closed" means here (and what stays open)

We deliver closed-form *structure* and closed-form *conditions*, marked `closed` / `conditional` /
`open` in `status.md`. The binding sign is `conditional` (rests on the worldtube closure axiom
`ωT=2πn`); the existence of a stable unit-winding `U(1)` eigenstate and a stable `B≠0` texture are
`open` (numerical, deferred). No new simulation is run; existing numerical results enter only as
falsifiers.

## Scope boundary

This is not a derivation of QCD: no asymptotic freedom, running coupling, string tension, area law,
or hadron spectrum is claimed. The confinement here is **kinematic + topological + lattice-regulated**
— a statement about which asymptotic states exist, not a dynamical inter-quark potential.
