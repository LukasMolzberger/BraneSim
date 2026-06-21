# Bell Bridge — Status

## HAVE

- Retrocausal worldtube interpretation: 4D brane lattice in 4D Euclidean ambient (codimension 0); time-symmetric dynamics; causal direction set by soliton chirality; entanglement = V-branching worldtube topology.
- Causality = soliton chirality: the arrow of time emerges from the chirality of the soliton (chiral characteristic boundary condition), not from the action's time-reversal properties (which are exactly time-symmetric).
- Entanglement = continuity of branching 4D worldtube: two entangled particles share a continuous 4D worldtube branch up to the branching event; Bell correlations are correlations within a single retrocausal worldtube, not across space.
- Placement in retrocausal completion family: TSVF (two-state vector formalism), transactional interpretation, time-symmetric Bohmian mechanics; all share the retrocausal structure forced by the substrate's time-symmetric BVP.
- The model commits explicitly to D = true and MI = true (no superdeterminism), placing it in the retrocausal quadrant of the Bell-option space.
- Non-probabilistic premise stated explicitly (`02_bell_constraint.tex` §"The probabilistic premise is not fundamental"): Bell constrains probability distributions, but the substrate is non-probabilistic at its core (deterministic field equations, no Born measure, no collapse). The model breaks with probability-fundamental interpretations (Copenhagen, objective collapse); probability/quantization are emergent. This does NOT exempt the theory from Bell — emergent correlations still relax MI via retrocausality; it relocates the probabilistic description one level up. Companion claim: quantization = emergent from soliton confinement + VSH geometry (Paper IV `02_matter.tex` §soliton-labels; Paper I `04_wave_structure.tex` §continuous-waves), with the Fraunhofer/cavity-mode analogy. Born-rule recovery remains the one open probabilistic debt (B-group below).

## MISSING

- Explicit formal argument that D + MI + local substrate → retrocausality uniquely; the "theorem" form with explicit premises, proof of exclusion of Bohmian alternatives, and proof of conclusion.
- No-signalling demonstration: show that retrocausal worldtube evolution preserves operational no-signalling (no faster-than-light information transfer in the effective theory).
- Born rule recovery: derive the Born probability rule from the deterministic substrate dynamics (or give an honest statement that this remains open and what the best current approach is — e.g., typicality argument analogous to Boltzmann).

## Open derivations

Drawn from `PRINCIPLES.md` §1.5 and `BACKBONE.md` #23. These are
constraints the retrocausal worldtube interpretation places on the program —
honest debts, not results. **Do not claim as established in the paper.**

### B1. Tsirelson's bound `2√2` — `open`
Reproduce the exact CHSH maximum from V-branching soliton correlations on the
substrate. Retrocausal models notoriously under-/over-shoot; deriving the exact
bound (not merely nonzero correlation) is the most direct empirical handle on
the stance.

### B2. No-signalling theorem — `open`
Derive no-signalling as a *theorem* about branching-worldtube dynamics (no
one-branch observable depends on the other branch's measurement context),
rather than assuming it. This is what prevents the retrocausal channel from
sending a macroscopic bit backward in time.

### B3. Baryon-to-photon ratio `η ≈ 6×10⁻¹⁰` — `open`
Show whether the geometric-vertex picture (matter/antimatter expanding from the
Big Bang vertex in opposite time directions) can match `η` quantitatively,
rather than merely permitting opposite-direction expansion in principle.

### B4. Threshold detection model — `open`
Model the electron–photon interaction at the detectors as a threshold event: the
electron soliton fires when the analyzer-projected carrier flux `∝ cos²(θ−a)`
crosses its transition gap, plus a hidden detector phase. Show this yields
single-arm Malus `cos²(θ−a)` for a uniform detector phase. This is the
*probabilistic ↔ energy* bridge and the concrete form of the Born-rule debt
(reduces it to "why is the detector phase uniform?", a typicality question).
Caveat: this component is local and, alone, Bell-bounded — it must NOT stand
alone or it rebuilds a detection-loophole model (the loophole Giustina 2015 /
Shalm 2015 closed). Supplies the *statistics*, not the Bell-violating angular
law. Feeds B1/B2.

### B5. V-vertex bidirectional energy transport — `open`
Specify the conserved-flux junction condition at the source vertex `S` and show
its time-symmetric (retarded offer + advanced confirmation) solution transports
the analyzer-setting information backward to `S` and forward to the other arm,
composing the two projections into `−cos 2(a−b)`. Formulated as a two-time BVP
on the branching worldtube (past BC = pump + vertex coupling; future BCs = the
two B4 threshold events). Must conserve energy, the singlet polarization
invariant, and carry no marginal advanced energy. This is the mechanism behind
"the V-vertex supplies the correlation"; it is the constructive route to B1
(exact Tsirelson) and B2 (no-signalling). Full note: `vvertex_energy_transport.md`.