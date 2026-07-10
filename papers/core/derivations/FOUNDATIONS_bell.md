# Foundations note — Bell's theorem and the time-symmetric brane lattice

**Status: positioning / consistency argument, not a derived result.** This note
records *why Bell's theorem does not exclude* the deterministic brane-lattice class.
It does **not** claim the model reproduces quantum correlations, entanglement, or the
Born rule — that derivation is owed (BACKBONE #12: the Born rule is "a statistical
regularity ... to be derived, not postulated"). Keep this in a Discussion/Foundations
section, not among the T1–T12 derived results.

---

## 1. The argument

Bell's theorem excludes theories that *simultaneously* hold local causality **and**
measurement (statistical) independence. (Determinism is not itself required by Bell —
the theorem also constrains stochastic theories — but it is a feature of this
substrate.) The brane lattice keeps:

- **determinism** — the world-volume is fixed by the stationary local elastic action
  `∇S=0` (BACKBONE #1, #12); node trajectories are single-valued, no fundamental
  randomness;
- **locally-mediated causation** — the only microscopic couplings are nearest-neighbor
  links on the 4D lattice, so every mediating influence propagates along worldlines
  *within* the world-volume (Wharton–Argaman's "locally mediated" sense — slightly
  weaker than strict Bell local causality, but with no spacelike causal channels).

The assumption that must therefore fail is **measurement independence**, i.e.
`P(λ|a,b) ≠ P(λ)` for the effective hidden state `λ` of an experimental run given
settings `a,b`.

Crucially, in this model that failure is **not** a superdeterministic conspiracy
fixed in the remote past. Layer 0 is an **all-at-once, two-boundary variational
problem**: `S = T − V` is indefinite, so the fundamental formulation is a saddle-point
boundary-value problem ("never minimize"; ARCHITECTURE Layer 0), and nonlinear runs
use **periodic time boundary conditions** (LESSONS_LEARNED #1). Future measurement
settings therefore enter as boundary data of the *same* 4D solution, so `λ` is
generally conditioned on the full boundary context while every influence stays local
in the world-volume. This is the ordinary consequence of a time-symmetric variational
problem, not a past common cause.

This places the model in the class of **locally mediated, time-symmetric
reformulations** analyzed by Wharton & Argaman (2020), which are compatible with
Bell's theorem precisely because they relax measurement independence (an
*arrow-of-time* assumption) rather than introducing spacelike causal channels.

**Conclusion:** Bell's theorem *constrains* the brane-lattice theory — it requires the
emergent ensemble description to be measurement-dependent — but does **not** exclude
the time-symmetric brane-lattice class itself.

## 2. What this does and does not establish

- **Establishes (consistency):** the deterministic, locally-mediated, time-symmetric
  substrate is *not* ruled out by Bell; the required measurement-dependence is a
  natural feature of its all-at-once variational structure, not an ad hoc conspiracy.
- **Does NOT establish (owed):** that the model actually *violates* a Bell inequality
  / reproduces the quantum correlations, nor the Born-rule weights. That is the
  BACKBONE #12 program (deterministic microstate → classical threshold readout →
  emergent, epistemic probability), still to be derived, and it lives at Layer 5
  (solitons/detectors), which the current T1–T12 chain has not yet built.

## 3. Commitments and caveats (state these if used in the manuscript)

- The argument **requires the all-at-once / two-boundary reading to be fundamental**,
  not merely a computational alternative to the Störmer–Verlet initial-value march. A
  pure initial-value formulation would fix `λ` from the past alone and force either
  Bell-exclusion or superdeterminism. The indefinite saddle (`S=T−V`) and periodic
  time BC support the boundary-value reading, but adopting it is a physical commitment.
- Relaxing measurement independence is a **substantive, debated** interpretational
  stance; the "time-symmetric, not conspiratorial" distinction (Wharton–Argaman) is
  its main defense and should be stated, not assumed.
- "Locally mediated" ≠ strict Bell local causality; use the former term.
- Time-symmetry of the specific action `S = ½Σ η_μκ_μ(L_{nμ}−r_μ)²` should be stated
  explicitly (it contains no odd-in-time terms; the block/periodic-time formulation is
  T-symmetric), since the whole argument rests on it.

## 4. Reference

- **K. B. Wharton & N. Argaman**, *Colloquium: Bell's theorem and locally mediated
  reformulations of quantum mechanics*, Rev. Mod. Phys. **92**, 021002 (2020);
  arXiv:1906.04313. *(Verified to exist.)* Also cross-listed in `RELATED_WORK.md` §6.

Related internal: BACKBONE #12 (deterministic substrate, emergent threshold-triggered
probability, Born rule to be derived); ARCHITECTURE Layer 0 (all-at-once saddle-point
formulation); LESSONS_LEARNED #1 (periodic BC).
