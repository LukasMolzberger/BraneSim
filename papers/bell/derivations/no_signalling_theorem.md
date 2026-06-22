# Bell Bridge — No-signalling as a theorem (B2)

Working note. Derives no-signalling — B's local marginal is independent of A's
distant setting `a` — as a theorem about the branching-worldtube measurement,
modulo two already-named inputs (completeness / unit-efficiency, and the Born
weight of D3). The structural core (completeness erases the distant setting under
the sum over unobserved outcomes) is rigorous and numerically verified. Honest
debt; not asserted in the paper.

## 0. Statement and the retrocausal worry

No-signalling: `P(B | a, b)` is independent of `a`. The worry is real: the
two-time / retrocausal structure makes the hidden state at `S` depend on `a`
(MI relaxed, D3) — so why doesn't `a` leak into B's *local* marginal, which would
be a backward-in-time signal? The theorem shows the `a`-dependence is confined to
the *joint correlation* and erased from the marginal.

## 1. The branching-worldtube measurement

A's arm terminates (absorption) in one of A's two analyzer-port eigenstates
`{|a⟩, |a+π/2⟩}`; B's in `{|b⟩, |b+π/2⟩}`. The joint history amplitude is the
Hilbert overlap `⟨a_A, b_B | Φ⁺⟩` (J2 Bell state at `S`, two-time consistency,
D3), with Born weight `|⟨a_A, b_B|Φ⁺⟩|²`. B is ignorant of A's outcome, so B's
local marginal sums over it:

    P(B | a, b) = Σ_{A ∈ {+,−}} |⟨a_A, b_B | Φ⁺⟩|² .

## 2. The theorem: completeness erases the distant setting

A's two analyzer-port eigenstates are an orthonormal basis, so for **every** angle
`a`,

    Σ_A |a_A⟩⟨a_A| = 𝟙        (completeness; verified `‖·−𝟙‖ ~ 1e-16` for all a).

Therefore, using `|⟨a_A,b_B|Φ⁺⟩|² = ⟨Φ⁺| (|a_A⟩⟨a_A| ⊗ |b_B⟩⟨b_B|) |Φ⁺⟩`,

    P(B | a,b) = ⟨Φ⁺| ( Σ_A |a_A⟩⟨a_A| ) ⊗ |b_B⟩⟨b_B| |Φ⁺⟩
               = ⟨Φ⁺| 𝟙 ⊗ |b_B⟩⟨b_B| |Φ⁺⟩
               = ⟨b_B| ( Tr_A |Φ⁺⟩⟨Φ⁺| ) |b_B⟩ = ⟨b_B| ρ_B |b_B⟩.

This contains **no reference to `a`** — the sum over A's complete outcome set
collapses `Σ_A|a_A⟩⟨a_A|` to `𝟙`, which is `a`-independent for any basis. With J2
(`ρ_B = ½𝟙`, rotationally invariant, `bellstate_lock_from_vertex.md`),
`P(B|a,b) = ½` — independent of both `a` and `b` (verified flat). ∎

The single step that does the work is the projector identity under the sum over A.
The Born weight `|amp|²` is what lets the sum factor through `Σ_A|a_A⟩⟨a_A|` (it
is a sum of `⟨·|a_A⟩⟨a_A|·⟩`); we do *not* claim Born is the only weight that can
be no-signalling (for the symmetric `|Φ⁺⟩` even a `p=1` weight stays flat — checked
— so no-signalling is not a uniqueness handle on Born), only that Born + completeness
gives the clean, general collapse to the reduced state.

## 3. Where the retrocausal `a`-dependence goes (and why it cannot signal)

The advanced leg does carry `a` back to `S` and forward to B — but it appears only
in the *individual* `(A,B)` histories, i.e. the joint correlation. Verified: the
`A=+`-only partial marginal `P(B=+, A=+)` runs `0.50 → 0.42 → 0.10` as `a` varies.
But B alone cannot see a partial sum — B is ignorant of A's outcome and sees only
the full sum `= ½`. To expose the `a`-dependence one must **compare** A's and B's
records (classical communication). So the retrocausal channel modulates
*correlations* (visible only on comparison) and never the *local marginal* — no
backward signal, no superluminal signal. This is exactly how QM no-signalling
coexists with Bell nonlocality; the substrate inherits it through completeness +
the Hilbert (overlap) structure of D3.

## 4. The role of completeness — and the loophole

Completeness (`Σ_A |a_A⟩⟨a_A| = 𝟙`) is **unit-efficiency / fair sampling**: every
photon at A is absorbed in *some* port and counted. It is essential — verified:
counting only the `+` port, B's partial marginal becomes `a`-dependent (§3). So:

- **Theorem holds for complete (unit-efficiency) detection** — the physical
  idealization, and the regime the loophole-free experiments (Giustina 2015,
  Shalm 2015) approached.
- With inefficient detection (`η < 1`) completeness fails and *coincidence*
  statistics can fake `a`-dependence — the **detection loophole** (the same one B4
  flags). This is not an actual signal (B still cannot read `a` without
  comparison); it is the fair-sampling caveat.

## 5. Substrate grounding

- **Completeness:** A's analyzer is a lossless polarization splitter — a unitary
  on the 2D polarization space — and both output ports terminate in absorption
  events. Two orthonormal port-eigenstates ⇒ projectors sum to `𝟙`. (Substrate
  fact: unitary polarization optics + complete absorption.)
- **`ρ_B` rotationally invariant:** J2 (`bellstate_lock_from_vertex.md`) — the
  emitted pair's one-party reduced state is maximally mixed.
- **Born weight / Hilbert overlaps:** D3 (`two_time_bvp_tsirelson.md`,
  `born_weight_typicality.md`).

## 6. Established vs open

**Established (rigorous / verified):**
- summing over A's complete outcome set erases `a` from B's marginal (the
  completeness collapse to `⟨b|ρ_B|b⟩`);
- with J2, B's marginal `= ½`, independent of `a` and `b`;
- the retrocausal `a`-dependence is confined to the joint correlation (partial
  sums), unreadable locally — no signal.

**Open (already-named, not new):**
- completeness `=` unit-efficiency / fair sampling (detection-loophole caveat, B4);
- the Born weights / Hilbert overlaps (D3, Born = typicality residual).

## Result

No-signalling is a theorem of the branching-worldtube measurement: completeness of
the local absorption (`Σ_A|a_A⟩⟨a_A| = 𝟙`) erases the distant setting under the
sum over unobserved outcomes, so `P(B|a,b) = ⟨b|ρ_B|b⟩ = ½` (J2), independent of
`a`. The retrocausal MI-relaxation lives entirely in the joint correlation, never
the local marginal — **Bell nonlocality without signalling**. The residual inputs
are the universal completeness (fair-sampling) and Born (typicality) conditions,
both already named — no new debt.

## References (wire into `references.bib` if this graduates)

- Standard QM no-signalling (reduced-state / local-operations argument) — the
  structure transplanted here.
- B4 `threshold_detection.md` (completeness / efficiency, detection loophole);
  D3 `two_time_bvp_tsirelson.md` (Born weights, two-time);
  `bellstate_lock_from_vertex.md` (J2, `ρ_B = ½𝟙`);
  `born_weight_typicality.md` (Born / typicality).
