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

## The hidden variable: τ = the absorber's local vacuum phase

The threshold model (1) needs a hidden variable τ whose ignorance produces the
apparent randomness. Identify it concretely: **τ is the absorbing electron's
local environmental/substrate phase at the moment of interaction.** Absorption
is the incoming carrier driving the electron soliton over a barrier ("up the
hill") into the next stable well; in relaxing it draws energy from its immediate
surroundings, and since the substrate field is continuous and not pre-quantized
into lumps, that energy is sourced partly from the carrier and partly from the
ambient substrate (vacuum) fluctuation it sits in. Which side of threshold the
electron lands on is fixed by the local fluctuation phase τ.

This is the substrate reading of two established ideas — Wheeler–Feynman
absorber accounting (the absorber's advanced response completes the transfer,
books balancing globally to `ℏω`) and stochastic electrodynamics (the
zero-point field supplies the local energy budget). The payoff is that it
converts the Born-rule debt from "why is τ uniform?" into a concrete typicality
claim about substrate vacuum-fluctuation statistics (equipartition of the local
phase), which is checkable rather than postulated.

**Conservation caution.** "Energy from the vacuum" must not drift into free
energy: net transfer is still `ℏω`; the vacuum exchange is a *borrowing* that
the advanced leg settles. Global conservation, local borrowing.

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

## Phase language: the V-vertex as Pancharatnam–Berry holonomy (conjugate to the energy picture)

The photon is the `U(1)` phase / Berry connection of Paper III, so there is a
phase-language description of the same junction, conjugate to the energy one.
The geometric phase of polarization is the **Pancharatnam–Berry phase**, living
on the **Poincaré sphere**. Two facts come for free:

- A physical analyzer angle `a` maps to `2a` on the sphere (orthogonal linear
  polarizations are 90° apart physically but antipodal on the sphere). This is
  the origin of the factor of 2 in `cos 2(a−b)` — Poincaré-sphere geometry, not
  an accident. (Connects to the spin-1 polarization frame and the double-angle
  cousin of Paper III's spinorial 2π/4π holonomy, `subsec:spinorial-holonomy`.)
- The singlet correlation is the sphere inner product:
  `E(a,b) = −cos(2a − 2b) = −cos 2(a−b)`.

So if the substrate's polarization transport through `S` reproduces
Poincaré-sphere parallel transport — which a Berry connection should — the
`cos 2(a−b)` *functional form* is geometric. The two descriptions are one
object:

| | phase language | energy language (this note) |
|---|---|---|
| photon | Pancharatnam–Berry connection on Poincaré sphere | carrier flux |
| V-vertex coupling | holonomy of the frame transported `A→S→B` | bidirectional offer/confirmation flux |
| `cos 2(a−b)` | sphere inner product (double angle automatic) | composed analyzer projections |
| randomness | τ = absorber's local vacuum phase | threshold draw from environment |
| no-signalling | off-diagonal phase only; diagonal fixed | marginal flux `a`-independent |

**Crucial caveat — geometry gives the *form*, retro-projection gives the
*violation*.** The Poincaré geometry applied *locally* to a shared λ still only
yields the local sawtooth (`S ≤ 2`). What turns "reveal a pre-existing point on
the sphere" into "project onto a setting-chosen axis" — the genuinely
Bell-violating step — is the threshold projection (1) being *setting-chosen*,
fed by the retro-coupling (2). Berry/Pancharatnam supplies the curve; the
V-vertex retro-coupling supplies the violation above `S = 2`. Do not let the
geometry oversell: it is necessary, not sufficient.

## Emission/absorption asymmetry is the thermodynamic arrow, not a new one

The intuitive asymmetry — an emitter is a localized high-energy source with the
future photon energy in a *defined* location, whereas an absorber sits at the
convergence of a wave that "could come from anywhere" — is the retarded/advanced
asymmetry: emission = diverging retarded wave from a point source, absorption =
converging advanced wave into a sink. The substrate action is time-symmetric and
treats the two identically; the *observed* direction (emit-then-absorb, sources
pointlike, sinks delocalized) is the **thermodynamic arrow** grounded in the
low-entropy past hypothesis — the same arrow committed in `03_worldtube.tex`
`subsec:soliton-chirality`, not a new fundamental asymmetry. Useful consequence
for this model: **the randomness enters at absorption, not emission** — the
photon is emitted into a definite mode; the τ-phase that decides the click lives
at the detector. This is consistent with "quantization/randomness supplied by
the detector" already held in the paper.

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
- Pancharatnam 1956; Berry 1984, 1987 — geometric phase of polarization on the
  Poincaré sphere; source of the `2a` double angle — supports the phase-language
  route.
- de la Peña & Cetto 1996 (stochastic electrodynamics) — zero-point field as the
  local energy budget for absorption — supports τ = local vacuum phase.