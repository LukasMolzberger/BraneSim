# Bell Bridge — The emission/absorption asymmetry (spatial grazing + polarization deviation)

Working note. Unifies two faces of one asymmetry: an **emitter** releases a
definite, exact, localized packet (in space *and* polarization), while an
**absorber** gathers a *deficit* — spatial (it is grazed, intercepting a small
cross-section) and angular (its axis is misaligned, intercepting only the
projected amplitude) — from the local vacuum, against a probabilistic threshold.
This grounds "randomness lives at absorption," gives a physical cause for the
vacuum-borrowing of B4/B5, and is grounded in the time-symmetric substrate + past
hypothesis (so the asymmetry is *emergent*, not assumed). Honest debt; not
asserted in the paper. Key caveat verified numerically.

## 1. The asymmetry, stated

- **Emission** — the emitting electron is in a definite excited state; it releases
  a self-contained retarded packet *from its center*, with a definite polarization
  fixed by the transition. Localized, exact, deterministic.
- **Absorption** — the absorbing electron is a tiny target met by a broad,
  already-spread wavefront. It directly couples to only a fraction of what it
  needs, and must **gather the rest of the quantum `ℏω` from the local vacuum**.
  Delocalized, approximate, probabilistic.

The same contrast appears in two independent "directions":

| | **spatial leg** | **polarization leg** |
|---|---|---|
| emitter | localized point source (energy at the center) | definite polarization `θ` |
| absorber's deficit | grazed: intercepts a small spatial cross-section of the broad wavefront | misaligned: intercepts only the projected amplitude `cos(θ−a)` |
| energy directly coupled | the intercepted flux | `cos²(θ−a)` |
| deficit to reach the gap | the rest of `ℏω` | `sin²(θ−a)` |
| supplied by | local vacuum / advanced confirmation | local vacuum / advanced confirmation |
| completion probability | Born | Born `= cos²(θ−a)` |

## 2. Both legs are the *same* vacuum-deficit mechanism (B4)

Normalize the transition gap to 1. The electron fires iff *directly-coupled
energy + borrowed vacuum energy ≥ gap*:

    fire ⟺ cos²(θ−a) + ε_vac ≥ 1 ⟺ ε_vac ≥ sin²(θ−a).

So the **polarization misalignment creates an energy deficit `sin²(θ−a)`** that the
vacuum must supply — *identical* to the spatial grazing, where the un-intercepted
flux is the deficit. With `ε_vac` uniform on `[0,1]` (the Born/typicality input,
`born_weight_typicality.md`),

    P(fire) = P(ε_vac ≥ sin²(θ−a)) = cos²(θ−a) = Born.

This is exactly B4 (`ε_vac = 1 − τ̃`). The two legs are one statement: *the
misalignment/grazing deficit is supplied by the vacuum, and the completion
probability is Born.*

## 3. Connection to the existing notes (one coherent picture)

- **Spatial (cause):** the grazed electron intercepts only a sliver → an energy
  deficit. *(New here; the missing spatial leg.)*
- **Polarization (cause):** the misaligned axis intercepts only `cos²` → the
  deficit `sin²`. The vacuum 2D amplitude `(q,p)` carries both: its **magnitude**
  is the energy threshold (`born_weight_typicality.md`), its **direction** is the
  polarization-frame deviation `δ`.
- **Energetic (B4):** the deficit is the borrowed `ε_vac`; Born iff uniform.
- **Temporal (B5):** the borrowing *is* the advanced confirmation reaching back to
  complete the transaction and deliver the full `ℏω`. Grazing/misalignment is
  *why* the advanced channel is needed at absorption and not at emission (emission
  releases a self-contained packet — no deficit).

## 4. The caveat that decides it (verified)

The deviation/grazing must be the **soft (amplitude) projection** — graded
response `cos(θ−a)`, squared via the 2D vacuum phase space to Born — tied to the
two-time/projective structure. It must **not** be a **hard, local** axis
deviation (absorber has a definite vacuum-deviated axis `a+δ` and reveals whether
the photon lands in its acceptance cone). That hard-local version is a local
hidden-variable model: verified to give CHSH `|S| = 1.96 ≤ 2` (sawtooth, not
Tsirelson) *and* a non-Malus single-arm shape. So a hard local "deviation" fails
on both counts. Same line as throughout: the grazing/deviation supplies the
*form/deficit*; the projective two-time structure (B5) is what reaches Tsirelson.

## 5. Why the asymmetry is emergent, not assumed

The substrate action is exactly time-symmetric (`subsec:time-symmetry`), so
emission and absorption are dynamically interchangeable. The asymmetry is
**thermodynamic / epistemic**, from the low-entropy past hypothesis (the same
datum as the arrow of time, `subsec:soliton-chirality`):

- the **emitter** is *prepared* by the low-entropy past — its internal phase is
  fixed (the excited state is a definite, pre-arranged condition), so emission is
  definite;
- the **absorber's** local vacuum phase `τ` is *typical / unknown* (a generic
  microstate drawn from the equilibrium measure), so absorption is probabilistic
  with the Born/typicality weight.

So "randomness at absorption, definiteness at emission" is the past hypothesis
acting on the detector microstate — not a fundamental emission/absorption
asymmetry, and consistent with the time-symmetric substrate.

## 6. Falsifiable consequence

In the ideal limit the deviation *is* exactly Malus. A finite absorber tolerance
beyond ideal would appear as a small, detector-dependent **deviation from Malus**
— the polarization analog of the off-axis Bell-state fidelity
`F = 1 − O((Δ/Γ)²)` (`bellstate_lock_from_vertex.md`). Separating a genuine
deviation-from-Malus from ordinary detector inefficiency is a concrete handle.

## 7. Result and residual

The emission/absorption asymmetry is one mechanism with a spatial and a
polarization leg: absorption gathers a deficit (cross-section + misalignment) from
the vacuum against a probabilistic threshold; emission is the definite release.
This gives the physical cause of the B4 vacuum-borrowing and the B5 advanced
channel, and grounds the definite/probabilistic split in the past hypothesis.

**Residual:** a quantitative cross-section model for the spatial grazing (how much
flux a point absorber intercepts from a mode of given transverse extent, and the
matching advanced-confirmation budget), and the explicit past-hypothesis argument
that the emitter phase is fixed while the absorber phase is typical (the same
ergodic-relaxation debt as `born_weight_typicality.md` §6).

## References (wire into `references.bib` if this graduates)

- Wheeler & Feynman 1945/1949 (absorber); Cramer 1986 (transactional) — the
  advanced confirmation gathering the deficit.
- de la Peña & Cetto 1996 (SED) — the vacuum as the local energy budget.
- B4 `threshold_detection.md`; B5 `vvertex_energy_transport.md`;
  `born_weight_typicality.md` (the 2D vacuum amplitude).
