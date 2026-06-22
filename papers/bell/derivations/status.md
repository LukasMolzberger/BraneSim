# Bell Bridge — Status

## HAVE

- Retrocausal worldtube interpretation: 4D brane lattice in 4D Euclidean ambient (codimension 0); time-symmetric dynamics with no intrinsic arrow; entanglement = V-branching worldtube topology.
- Arrow of time vs chirality — two distinct questions, answered separately (§3.2 `03_worldtube.tex`). (a) The **experienced** arrow (irreversibility, records of the past not the future, the thermodynamic gradient) is the *thermodynamic* arrow: it rests on a low-entropy boundary condition at one temporal end (the past-hypothesis / Big-Bang-vertex datum, = the single residual tuning named in §2.5), inheriting the standard time-symmetric account. It is NOT produced by soliton chirality, and the substrate action is exactly time-symmetric. (b) Soliton chirality answers a *different* question: it fixes matter vs antimatter and the local time-orientation of a given worldtube (matter = forward-propagating, antimatter = backward-propagating), but does not by itself supply the experienced arrow.
- Emission/absorption asymmetry (`emission_absorption_asymmetry.md`): one mechanism with a spatial leg (emitter releases a localized packet from its center; absorber is *grazed*, intercepting only a small cross-section) and a polarization leg (emitter has exact polarization; absorber is *misaligned*, intercepting only `cos²(θ−a)`). Both = absorption gathering an energy *deficit* (`sin²` of the misalignment / the un-intercepted flux) from the local vacuum against a probabilistic threshold — *identical* to B4's vacuum-borrowing, with the deficit supplied by B5's advanced confirmation. Grounds "randomness at absorption, definiteness at emission" in the past hypothesis (emitter prepared → phase fixed; absorber vacuum phase typical), so the asymmetry is emergent from the time-symmetric substrate, not assumed. Caveat (verified): the deviation must be the soft amplitude projection (→ Born); a hard local axis-deviation is a local HV model (CHSH `|S|=1.96≤2`, non-Malus). Falsifier: small detector-dependent deviation-from-Malus (polarization analog of the off-axis fidelity).
- Entanglement = continuity of branching 4D worldtube: two entangled particles share a continuous 4D worldtube branch up to the branching event; Bell correlations are correlations within a single retrocausal worldtube, not across space.
- Placement in retrocausal completion family (§2.7 `02_bell_constraint.tex`): Aharonov's two-state-vector formalism, Cramer's transactional interpretation, Price–Wharton retrocausal models, Sutherland's time-symmetric Bohmian model; all share the retrocausal structure forced by the substrate's time-symmetric BVP. Structural distinction: here retrocausality is *built into the time-symmetric action*, not imposed as an interpretive rule on top of forward-causal dynamics. Pilot-wave lineage: the model completes de Broglie's *double solution* (the particle is a soliton-singularity of a single **real wave in 3-space**, guidance = phase-locking theorem, local) — NOT Bohm's configuration-space pilot wave (point particle + 3N-dim guiding wave, irreducibly nonlocal), which is excluded on the same locality ground as every nonlocal completion (§2.6). de Broglie's unresolved many-body/entanglement gap is filled here by the 4D retrocausal worldtube rather than by config-space nonlocality. Full argument (incl. the guidance relation as a centroid *theorem* of the carrier, `dX/dt = ∇S/m_eff`) in `de_broglie_double_solution.md`.
- Bell-option commitment (§2.2 `02_bell_constraint.tex`): (D) deterministic and (L) local are non-negotiable substrate facts; therefore **(MI) measurement-independence is the assumption relaxed** — via retrocausality, NOT via superdeterminism (§2.5). The MI relaxation is sourced by a uniform local action coupling freely chosen future boundary data to the emission state, not by fine-tuned past data. This places the model in the retrocausal quadrant of the Bell-option space.
- Non-probabilistic premise stated explicitly (`02_bell_constraint.tex` §"The probabilistic premise is not fundamental"): Bell constrains probability distributions, but the substrate is non-probabilistic at its core (deterministic field equations, no Born measure, no collapse). The model breaks with probability-fundamental interpretations (Copenhagen, objective collapse); probability/quantization are emergent. This does NOT exempt the theory from Bell — emergent correlations still relax MI via retrocausality; it relocates the probabilistic description one level up. Companion claim: quantization = emergent from soliton confinement + VSH geometry (Paper IV `02_matter.tex` §soliton-labels; Paper I `04_wave_structure.tex` §continuous-waves), with the Fraunhofer/cavity-mode analogy. Born-rule recovery remains the one open probabilistic debt (B-group below).

## MISSING

- Fully formalized *theorem* form of the uniqueness argument: D + L + (one of D/L/MI must be relaxed) ⟹ relax MI via retrocausality, uniquely. The prose argument now exists in §2.6 `02_bell_constraint.tex` (four-option enumeration; nonlocal/Bohmian completions excluded on the locality ground; superdeterminism excluded as fine-tuning vs uniform-local-law in §2.5). Residual: restate it with explicit premises and a proof of exclusion + conclusion, rather than as structured prose.
- No-signalling demonstration: show that retrocausal worldtube evolution preserves operational no-signalling (no faster-than-light information transfer in the effective theory).
- Born rule recovery: *reduced* in `born_weight_typicality.md` to Liouville/microcanonical typicality on the complex carrier amplitude — the Born exponent `2` is the symplectic dimension of `ψ=q+ip` (`|amp|²` = phase-space area, `cos²` = projected area ratio), and the measure is the equal-a-priori postulate (uniform-on-disk ⇒ uniform intensity ⇒ Born; verified). This makes quantum probability the *same* as thermal probability, not a new kind. The only residual is the **universal** ergodic-relaxation debt (substrate vacuum samples the microcanonical measure) — identical to justifying the microcanonical ensemble in classical stat mech, not a quantum mystery. NOT a derivation ex nihilo (no approach is).

## Open derivations

Drawn from `PRINCIPLES.md` §1.5 and `BACKBONE.md` #23. These are
constraints the retrocausal worldtube interpretation places on the program —
honest debts, not results. **Do not claim as established in the paper.**

### B1. Tsirelson's bound `2√2` — `open` (reduced)
Reproduce the exact CHSH maximum from the V-branching correlations on the
substrate. Retrocausal models notoriously under-/over-shoot; deriving the exact
bound (not merely nonzero correlation) is the most direct empirical handle on
the stance. *Reduced* in `two_time_bvp_tsirelson.md` (= D3): the two-time
structure alone does not cap the correlation (a PR-box relaxes MI too, `S = 4`),
but the substrate's **linear complex `ℂ²` carrier** makes the two-time history
amplitudes genuine Hilbert overlaps — excluding PR-boxes and capping at
Tsirelson; with Born weights this gives *exactly* `±cos 2(a−b)`, `|S| = 2√2` (a
non-Born weight fails, verified). The Tsirelson cap and the functional form thus
follow; the **only** residual is the Born weight `|amp|²` = the joint two-time
`(τ_A,τ_B)` typicality, i.e. the same Born-rule debt as B4/below — localized, not
a separate miracle. That debt is itself reduced in `born_weight_typicality.md` to
Liouville/microcanonical typicality on the complex amplitude (exponent `2` =
symplectic area), leaving only the universal ergodic-relaxation residual.

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
*probabilistic ↔ energy* bridge and the concrete form of the Born-rule debt.
The hidden phase is identified concretely: **τ = the absorbing electron's local
substrate/vacuum phase** (absorption draws energy from carrier + ambient field;
Wheeler–Feynman / stochastic-electrodynamics reading), reducing the debt to a
typicality claim about vacuum-fluctuation statistics rather than "why uniform?".
Caveat: this component is local and, alone, Bell-bounded — it must NOT stand
alone or it rebuilds a detection-loophole model (the loophole Giustina 2015 /
Shalm 2015 closed). Supplies the *statistics*, not the Bell-violating angular
law. Feeds B1/B2. Formalized in `threshold_detection.md`: threshold ⇒ Malus
*iff* the effective threshold is uniform (the precise Born content, by a
converse argument), with single-arm no-signalling robust to the threshold
distribution; numerically checked. The uniformity itself is *grounded* in
`born_weight_typicality.md`: `τ̃` = vacuum intensity `|ψ_vac|²`, uniform under
Liouville/microcanonical typicality on the 2D amplitude phase space `(q,p)` ⇒
Born; the exponent `2` is the symplectic area.

### B5. V-vertex bidirectional energy transport — `open`
Specify the conserved-flux junction condition at the source vertex `S` and show
its time-symmetric (retarded offer + advanced confirmation) solution transports
the analyzer-setting information backward to `S` and forward to the other arm,
composing the two projections into `−cos 2(a−b)`. Formulated as a two-time BVP
on the branching worldtube (past BC = pump + vertex coupling; future BCs = the
two B4 threshold events). Must conserve energy, the singlet polarization
invariant, and carry no marginal advanced energy. This is the mechanism behind
"the V-vertex supplies the correlation"; it is the constructive route to B1
(exact Tsirelson) and B2 (no-signalling). Has a conjugate phase-language form:
the coupling as **Pancharatnam–Berry holonomy** of the polarization frame
transported `A→S→B` on the Poincaré sphere, which makes the `cos 2(a−b)` *form*
(incl. the double angle) geometric — but the *violation* above `S=2` still
requires the setting-chosen threshold projection (geometry necessary, not
sufficient). The form `E = −cos 2(a−b)` is derived (and numerically checked,
CHSH `|S| = 2√2`) in `pancharatnam_holonomy.md`; the energy-transport mechanism
is in `vvertex_energy_transport.md`. The junction condition at `S` is formalized
in `vvertex_junction_condition.md` (J1 energy, J2 singlet-invariant lock, J3
time-symmetric two-time closure via chiral-characteristic BC): it *encodes*
flat-marginal no-signalling (from J2) and MI-relaxation (from J3), and states the
conditional/collapse structure it must realize. Its own open debts: **D1** derive
J2 (singlet) from the substrate vertex nonlinearity; **D2** well-posedness /
uniqueness of the two-time branching BVP; **D3** prove the BVP yields *exactly*
`−cos 2(a−b)` (= B1 Tsirelson, the decisive step); **D4** no-signalling as a
theorem for the full nonlinear solution (= B2). Verified backbone: a *local*
junction is capped at `|S| = 2.000`, so the two-time BC is mandatory.
D1 is now *computed* in `bellstate_lock_from_vertex.md`: the geometric link
energy gives the cubic three-wave vertex `(ê·δu)|δu⊥|²`, whose symmetric
`δ`-contraction yields `M ∝ δ` ⇒ the emitted Bell state is **`|Φ⁺⟩`** (symmetric,
type-I-like), **not** the singlet — the real norm gives `δ`, not the
SU(2)-invariant `ε`. `|Φ⁺⟩` is local-unitary (waveplate) equivalent to the
singlet, gives `|S| = 2√2` and flat marginals (linear analyzers:
`E = +cos 2(a−b)`), so J2 is now derived from the Lagrangian. Residual (minimal):
confirm the longitudinal-pump three-wave channel + phase matching, and the
off-axis fidelity coefficient `F = 1 − O((Δ/Γ)²)`. (Corrects the earlier
`M ∝ ε`/singlet assumption.)