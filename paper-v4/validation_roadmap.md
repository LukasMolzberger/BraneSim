# Validation roadmap

This file lists the discrete subtasks needed to rigorously test the BraneSim
core claims. Each subtask has a clear input, output, falsifiable success
criterion, and an owner agent (see `.claude/agents/`). Subtasks are organized
into six dependency-respecting sprints — earlier sprints establish the
foundations the later ones depend on.

The rule: **sub-agents must not skip levels.** A claim at level N is only
meaningful if the levels below have passed.

---

## Layer map (recap)

```
L0 substrate (lattice, ℝ⁴ embedding)
    │  homogenization
    ▼
L1 hyperelastic continuum (induced metric, StVK)
    │  linearization                         \
    ▼                                          \  slow in-brane mode
L2 linear branches (T₁, T₂, L; ω(k) per branch) \
    │  band isolation                            ▼
    ▼                                          L8 contraction channel
L3 eigen-bundle / mode frame                       │  Newtonian limit
    │  geometric phase                            ▼
    ▼                                          gravity-like Φ
L4 Berry / Wilczek–Zee connection
    │  effective-action projection
    ▼
L5 slow-sector EFT (Ψ, 𝒜, ℱ²)
    │
    ▼
L6 effective metric / Lorentz
    │
    ▼
L7 localized solitons
```

---

## Sprint 1 — IR foundations (L0 → L2)

The Cauchy-isotropy argument predicts isotropic long-wavelength dispersion *by
construction*. This sprint converts that prediction into measurable numbers
and confirms or falsifies the IR foundation.

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 1 | Lattice-to-continuum elastic-constants derivation: closed-form `μ(k_δ, a)`, `λ(k_δ, a)` from the central-force pair springs with `1/|δ|²` shell weights. **Done** (`paper-v4/derivations/lattice_to_continuum.md`). Finding: standard Cauchy relation `C_{1122} = C_{1212}` holds at `α=1`, but cubic isotropy `C_{1111}−C_{1122} = 2 C_{1212}` is **violated** with the `1/|δ|²` weights, with leading-order anisotropy index `η_cub = −7α/[2(39−22α)]`. Predicts `c_T ≈ 1.912`, `c_L ≈ 1.989`, `c_L/c_T ≈ 1.040` at `α=0.2` (the local_extensive default). | `physics-derivation` | ✅ closed-form expressions + numerical predictions |
| 1b | Solve for shell weights that recover cubic isotropy: find the manifold of `(w_I, w_{II}, w_{III})` satisfying `2 w_I = w_{II} + (8/3) w_{III}` (the isotropy condition; derived in subtask 1's §3). Choose a representative point, document the trade-off (e.g. effective `c_L/c_T`, condition number), and propose retuning. | `physics-derivation` | closed-form weight family + recommended choice + impact estimate |
| 2 | Long-wavelength dispersion + isotropy: measure `ω([100])` vs `ω([110])` vs `ω([111])` at `|k|·a ∈ {0.05, 0.1, 0.2, 0.3}`. | `dispersion-analyst` | < 1% deviation at `|k|·a ≤ 0.1`; report scaling of deviation with `(ka)`; cross-check against subtask 1's `(ka)²` vanishing prediction |
| 3 | Branch-speed match: measure `c_T`, `c_L` directly and compare to subtask 1's prediction. | `dispersion-analyst` | < 5% deviation, else subtask 1 wrong or numerical issue |
| 4 | Band isolation: identify `(α, ω₀)` window with `Δω_gap / ω₀ > 0.1` and verify mode leakage decays exponentially. | `dispersion-analyst` + `berry-validator` | yes/no with measured leakage time-constant |

Pass criterion for sprint: subtask 1 reproduces the empirical `c_T`, `c_L` to
within 5%, and direction-dependent dispersion is below 1% at `|k|·a ≤ 0.1`.

---

## Sprint 2 — Gauge structure (L3 → L5)

Once a clean isolated band exists, test that Berry/WZ holonomy is well-defined,
gauge-invariant, and converges. Then check the U(1)³ → U(3) crossover via `α`.

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 5 | Berry holonomy gauge invariance: random local re-phasing leaves invariants unchanged | `berry-validator` | ≤ 1e-3 deviation, hard pass/fail |
| 6 | Berry holonomy convergence: halve `Δt` and lattice spacing → ≤ 1% drift | `berry-validator` | hard pass/fail |
| 7 | Adiabatic-breakdown demonstration: drive across the gap, observe loss of geometric robustness | `berry-validator` | qualitative pass |
| 8 | U(3) → U(1) ⊕ SU(3) decomposition as function of `α`: traceless-part magnitude scales with `(1 − α)` and vanishes at `α = 1` | `berry-validator` | scaling law in `(1 − α)` |
| 9 | Effective-action stiffness `κ`: derive from substrate parameters and compare to measured curvature variance | `physics-derivation` + `berry-validator` | predicted vs measured `κ` |

Pass criterion: Berry invariants are stable, the SU(3) sector vanishes at
`α = 1`, and a derived `κ` matches the measured one within an order of
magnitude.

---

## Sprint 3 — Lorentz / observer universality (L2 → L6)

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 10 | Wave-cone closure: argue or measure that the dominant branch's superluminal channels (if any) are not excitable by physical solitons | `physics-derivation` | analytic argument or numerical leakage bound |
| 11 | Observer universality: two solitons with different polarizations report identical effective `g^eff_{μν}` | `physics-derivation` + `soliton-hunter` | numerical equality of measured time-dilation |

Pass criterion: leakage to non-dominant branches < 1% under realistic soliton
amplitudes; observer time-dilation factors agree to ≤ 1%.

---

## Sprint 4 — Solitons (L7)

Critical decision point: does the substrate support stable, localized 3D modes
at width ≫ a (Skyrme-stabilized) or only at width ~ a (lattice-stabilized)?

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 12 | Soliton existence: stable localized standing-wave mode survives ≳ 10³ carrier periods without radiation; sweep over seed width `R_seed ∈ {2a, 4a, 8a, 16a, 32a}` to identify regime A vs regime B | `soliton-hunter` | confinement score above threshold; explicit width-to-spacing ratio in report |
| 13 | Soliton stability: linear stability eigenproblem about `X_0` has no growing mode | `soliton-hunter` | spectrum of `ℒ_0` |
| 14 | Triplet baryon search: parameter sweep `(α, k_1, k_2, k_3, R, A, θ_mix)` ranking by confinement quality | `soliton-hunter` | top-10 candidate configs |
| 15 | Spin-½ holonomy: `4π`-vs-`2π` behavior of a 2-level polarization subspace under spatial rotation | `berry-validator` + `soliton-hunter` | sign-flip at `2π`, return at `4π` |
| 16 | Fermion-statistics gap: explicit acknowledgement in the paper that exchange statistics is **not** addressed by the current theory | `paper-writer` | label change in paper |

Pass criterion for sprint: at least one configuration in regime B (width ≫ a)
that is stable for 10³ carrier periods, with measurable triplet U(3) holonomy.

---

## Sprint 5 — Gravity channel (L1 → L8)

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 17 | Contraction-field extraction: measure ξ(x,t) sourced by a localized u-bump | `contraction-channel` | qualitative agreement with `∂ᵢu ∂ⱼu` source |
| 18 | Newtonian-limit reduction: derive `∇²Φ = 4πG_eff ρ_eff`; measure `G_eff` | `physics-derivation` + `contraction-channel` | a single number with error bar |
| 19 | Equivalence-principle proxy: two test wavepackets of different polarization fall identically through the same Φ | `contraction-channel` | %-level identity |
| 20 | α-scaling of gravity coupling: `G_eff → 0` as `α → 1` | `contraction-channel` | scaling law |

Pass criterion: a single `G_eff` value reproducible across runs, and EP proxy
holding to ≤ 1%. EP failure here is a fatal result for the gravity story.

---

## Sprint 6 — Closure (L5 → L7 / global)

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 21 | Maxwell-term derivation: show `tr(ℱ²)` is the leading IR term in the slow-sector EFT under stated power counting | `physics-derivation` | small parameter + neglected-term bound |
| 22 | Charge identification + Coulomb law: define `A_μ^eff`, identify candidate charge density of a localized soliton, derive `1/r` far-field | `physics-derivation` + `berry-validator` | analytic or numerical `1/r` tail |
| 23 | Calibration / parameter count: list substrate parameters vs target observables `(c, ℏ, e, m_e, G, Λ_QCD)`; report identifiability | `physics-derivation` + `paper-writer` | parameter table |
| 24 | Quantum-statistics toy: pick ONE phenomenon (double-slit, Stern-Gerlach analog, photoelectric threshold) and try to reproduce the statistics | `physics-derivation` + `soliton-hunter` | matches Born-rule frequencies in toy regime, or fails clearly |

---

## Anti-goals (what we will NOT do to make a sprint pass)

- Add a confinement force, a damping term, an amplitude clamp, or a hand-tuned
  saturation rule.
- Tune shell weights to obtain isotropy after we already showed it is automatic
  by the Cauchy relation. If isotropy fails empirically with central-force
  springs, the failure is itself the result.
- Postulate a global narrowband carrier sector to force Berry diagnostics to
  work. Per principles §4.3, narrowband is per-wavepacket.
- Hide preferred-frame effects behind "operational isotropy" hand-waving when
  they appear above the empirical bounds.

---

## How to claim a subtask

Sub-agents claim a subtask by setting their owner on the corresponding task in
the task list (see `TaskList`). Deliverables go under `test-runs/<run-id>/`
with a one-page Markdown report named `<subtask>_report.md`.