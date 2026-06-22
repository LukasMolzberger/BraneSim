# Bell Bridge — Threshold detection model (B4)

Working note. **Derives single-arm Malus `cos²(θ−a)` from a threshold, and
reduces the Born rule to one typicality (uniformity) statement.** LOCAL — it does
NOT produce the Bell violation (that needs B5). Honest debt; not asserted in the
paper. Companion to `pancharatnam_holonomy.md` (the joint form) and
`vvertex_energy_transport.md` (the vertex coupling). Algebra below numerically
verified.

## Scope

This is the *probabilistic ↔ energy* bridge: it turns the continuous substrate
carrier + the confined detector electron into a discrete, probabilistic "click."
It establishes (1) single-arm Malus, (2) that the Born content is exactly the
*uniformity of the detector's effective threshold* (a typicality claim, not a
postulate), (3) that single-arm no-signalling is robust, and (4) that the model
is local and Bell-bounded on its own.

## 1. Delivered energy = squared projection (Malus, geometric)

The analyzer at `a` is a polarization projector onto axis `â`. A photon carrier
with unit polarization `ê(θ)` delivers transmitted amplitude `ê·â = cos(θ−a)`,
hence energy flux (Poynting, normalize `I₀ = 1`)

    I(θ, a) = |ê·â|² = cos²(θ − a).

This is Malus's law for the delivered energy, and it is the *same* projection
geometry as the Poincaré-sphere overlap `|⟨a|θ⟩|² = cos²(θ−a)` of
`pancharatnam_holonomy.md` §1 — not an extra assumption.

## 2. The detector as a threshold device (hill-climbing)

The detector electron is a confined soliton with a metastable internal
coordinate `q` in a potential `U(q)`: a ground well separated from the next state
by a barrier of height `B`. Detection = the electron surmounts `B` and relaxes
(registers a click, absorbing one transition gap `ΔE`).

Two energy contributions drive the climb:
- the carrier flux `I(θ,a)` it can deliver;
- a fluctuating energy `ε(τ)` the electron borrows from its local substrate
  (vacuum) environment, set by a hidden phase `τ` (the electron's internal
  limit-cycle phase + the local field fluctuation at the interaction event).

Firing rule:

    fire  ⟺  I(θ,a) + ε(τ) ≥ B   ⟺   I(θ,a) ≥ B − ε(τ) ≡ τ̃ ,

i.e. the carrier lowers the effective barrier `τ̃ ∈ [0,1]` the vacuum must clear.
`τ̃` is the **effective threshold** (gap minus borrowed vacuum energy), the
single hidden variable of the model.

*Conservation note (from `vvertex_energy_transport.md`):* net absorbed energy is
one gap `ΔE = ℏω`; the borrowed `ε(τ)` is settled by the advanced leg. Global
conservation, local borrowing — not free energy.

## 3. Threshold ⇒ Malus, and the Born debt as uniformity

With `F` the CDF of the effective threshold `τ̃`,

    P(click | θ, a) = P(τ̃ ≤ I(θ,a)) = F(cos²(θ − a)).

**Theorem.** If `τ̃ ~ Uniform[0,1]` then `F(x) = x` and

    P(click | θ, a) = cos²(θ − a)   — single-photon Malus, all randomness in τ̃.   ∎

**Converse (the honest content).** `P(click|θ,a) = cos²(θ−a)` for all `(θ,a)`
**iff** `F(x) = x`, i.e. iff `τ̃` is uniform. So the threshold model does not
*assume* Born/Malus — it *reduces* it to a single statement: the detector's
effective threshold is uniformly distributed over its energy window. That is a
**typicality / equipartition** claim about substrate vacuum fluctuations
(Boltzmann-style; cf. quantum-equilibrium typicality), checkable rather than
postulated. This is the precise residual form of the Born-rule debt.

(Numerically: `τ̃ ~ U[0,1]` reproduces `cos²(θ−a)` to sampling error.)

## 4. Discreteness: continuous carrier → one click

The threshold fires once and the electron absorbs exactly one gap `ΔE`,
regardless of how much carrier energy was available; the surplus stays in the
field/vacuum. A continuous, unquantized carrier input yields a discrete output.
This is where the "single photon, one lump `ℏω`" is manufactured — at the
detector electron, not in the field — consistent with "quantization supplied by
the detector" already held in the paper (`02_bell_constraint.tex`
§"Where the randomness of a photon measurement actually lives").

## 5. Single-arm no-signalling is robust

For the entangled pair each photon's reduced state is unpolarized, so the local
polarization `θ` is effectively uniform on `[0,π)`. The marginal click rate

    P(click | a) = (1/π) ∫₀^π F(cos²(θ − a)) dθ

is **independent of `a`** by translation invariance of the period integral — for
*any* threshold CDF `F`, not only the uniform one. (Verified numerically: a
deliberately non-uniform `F(x)=√x` still gives a flat marginal `≈ 0.637`; the
uniform case gives `1/2`.) So one-arm no-signalling (part of B2) holds
structurally and does not lean on the Born assumption of §3.

## 6. What B4 does NOT do (the load-bearing caveat)

This is a **local, single-arm** model. For the joint correlation: if each photon
carries a definite shared `λ = θ` from the source and each side applies this
threshold rule with **independent** `τ̃_A, τ̃_B`, the joint correlation is the
local one — the triangular/sawtooth, bounded by `|S| ≤ 2` — **not**
`−cos 2(a−b)`. Recovering the Bell-violating curve requires the two detectors'
hidden variables (or the polarization frames) to be **setting-correlated through
the vertex `S`**, i.e. the retro-coupling of B5; equivalently, the threshold must
be *setting-chosen* (project), not *reveal* a pre-existing `λ`.

Standing alone with post-selection on non-detections, the model is exactly a
detection-loophole construction (Pearle 1970; Gisin & Gisin 1999) — the loophole
the committed experiments (Giustina 2015, Shalm 2015) closed. **B4 must therefore
be wired to B5; it supplies the statistics and the discreteness, never the
violation.**

## References (wire into `references.bib` if this graduates)

- Lamb & Scully 1969 — photoelectric effect without photons (field classical,
  detector quantized) — the semiclassical basis for §1–4.
- de la Peña & Cetto 1996 (stochastic electrodynamics) — zero-point field as the
  local energy budget `ε(τ)`.
- Boltzmann typicality; Dürr, Goldstein & Zanghì 1992 (quantum equilibrium) —
  frame for the uniformity-of-`τ̃` (Born) reduction in §3.
- Pearle 1970; Gisin & Gisin 1999 — local threshold / detection-loophole models;
  the trap §6 names.
