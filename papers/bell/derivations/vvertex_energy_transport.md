# Bell Bridge — V-vertex bidirectional energy transport

Working note. **Honest debt, not a result.** Do not assert in the paper.

This is the constructive route to the gap named in `critique.md` (item 3): *no
actual probability model reproduces the quantum statistics while preserving
no-signalling.* It splits that gap into two separable components and specifies
the second.

## The measurement has two separable components

1. **Threshold detection events** (at detectors `D_A`, `D_B`) — the
   *probabilistic ↔ energy* bridge. Local; sits at the detectors; converts the
   continuous delivered energy into a discrete "click" via the electron
   soliton's transition gap. Yields single-arm Malus `cos²(θ−a)` from a
   threshold on the analyzer-projected flux plus a uniform detector phase.
   Tracked separately (see `threshold_detection.md` when written; not the
   subject here).
2. **Bidirectional energy transport through the V-vertex `S`** — what actually
   couples the two arms. *This note.* It is distinct from the slogan "the
   V-vertex supplies the correlation"; that slogan is a label, this is the
   mechanism.

Component (1) alone, if local, is Bell-bounded (`S ≤ 2`) at unit efficiency, or
reproduces a detection-loophole model (Pearle 1970; Gisin & Gisin 1999) — the
very loophole the committed experiments (Giustina 2015, Shalm 2015) closed. The
Bell-violating angular law must therefore come from component (2).

## The object to model

On the branching worldtube the carrier field splits into two time-symmetric
legs (substrate realization of Wheeler–Feynman absorber theory / Cramer's
transactional handshake, both already cited in the paper):

- **retarded (offer) leg** — carries the actual photon energy *forward*,
  `S → D_A` and `S → D_B`;
- **advanced (confirmation) leg** — carries the detector threshold response
  *backward*, `D_A → S` and `D_B → S`.

"Energy bidirectionally transmitted through `S`" = this full time-symmetric
flux. Net **real** energy delivery is forward (pump → two clicks); the
**backward** leg carries no marginal energy but transports the setting
information into `S` and thence onto the other arm. Both legs and the ledger
tying them must be modeled.

## Conservation laws at the vertex `S`

`S` is a nonlinear matter vertex (the down-conversion crystal). The junction
condition is a Kirchhoff-like flux balance that must conserve three things at
once:

- **Energy.** `ℏω_pump = E_A + E_B` on the forward leg; the advanced fluxes from
  the two detectors balance against the source's own past (laser / past
  absorber) on the backward leg.
- **Polarization invariant.** The singlet frame relation `θ_B = θ_A + π/2`,
  invariant under a common rotation of both frames. This is the conserved
  "charge" routed through `S` — the analogue of polarization / angular-momentum
  conservation (phase-matching) in real SPDC.
- **No net advanced energy.** The backward leg carries *correlation* without a
  *marginal* (the no-signalling constraint below).

## How the backward leg carries the correlation (energy-explicit)

1. Threshold at `D_A` (analyzer `a`) fixes the advanced amplitude on arm `A`
   `∝ cos(θ_A − a)` — the analyzer-projected component.
2. That advanced flux reaches `S`; the singlet invariant ties it to arm `B`'s
   frame.
3. It re-emerges forward toward `D_B`, modulating the energy available to fire
   `B`'s threshold `∝ cos²(θ_B − b)` with `θ_B` locked to the `A`-confirmation.

Composing the two analyzer projections through the invariant at `S` is what
should yield the `cos 2(a−b)` angular law — *if* the junction is built
correctly.

## The two make-or-break constraints

- **B2 — no-signalling (sharp constraint on the backward leg).** The advanced
  flux must couple into `B`'s *correlation* with `A` but average to zero in
  `B`'s *marginal*: energy delivered to `B`, integrated over the hidden phase,
  must be independent of `a`. This is the precise condition that stops the
  bidirectional energy from being a backward telegraph. Naive
  bidirectional-energy models fail here first.
- **B1 — exact Tsirelson.** The junction composition must land on
  `E(a,b) = −cos 2(a−b)` / `S = 2√2` — not the local sawtooth, not an
  over/undershoot. Reproducing the exact bound (not merely nonzero correlation)
  is the decisive empirical handle.

## Formulation: a two-time BVP on the branching worldtube

This is the committed solver paradigm (stationary point of the time-symmetric
action, root-find `‖R‖ = 0` on the 4D block), specialized to a branching domain:

- **past BC** at `S`: pump/source + the nonlinear vertex coupling (enforces the
  three conservation laws above);
- **future BCs** at `D_A`, `D_B`: the two detector threshold events from
  component (1);
- **solution:** the time-symmetric field on the V; its **advanced component is
  exactly the energy bidirectionally transmitted through `S`.**

So the missing piece is concretely: *specify the vertex junction condition (the
conserved-flux coupling at `S`) and show its advanced component satisfies the
no-signalling marginal (B2) while composing to `−cos 2(a−b)` (B1).* This is
well-posed and separable from the threshold model.

## Open sub-questions (before this is a derivation, not a sketch)

- Exact functional form of the junction / phase-matching coupling at `S` (which
  conserved current, which nonlinearity).
- Energy ledger for the advanced leg — how "borrowed" advanced flux is returned
  (global vs per-instant conservation).
- Proof that the marginal is `a`-independent (B2) *given* that the joint depends
  on `a`.

## References (wire into `references.bib` if this graduates to the paper)

- Wheeler & Feynman 1945/1949 — absorber theory (retarded + advanced).
- Cramer 1986 — transactional interpretation (offer / confirmation). *Already
  cited.*
- Lamb & Scully 1969 — photoelectric effect without photons (semiclassical
  detection, detector-quantized) — supports component (1).
- Pearle 1970; Gisin & Gisin 1999 — local threshold / detection-loophole models
  reproduce QM only below unit efficiency — the trap component (1) must avoid.