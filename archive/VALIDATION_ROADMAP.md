# Validation roadmap

> **ARCHIVED / historical** — see `archive/README.md`. Predates the
> 6-neighbor-only commitment and the U(3) reframing; uses the old `L0–L8`
> layering and 26-neighbor shell-weight language. Snapshot, not the current plan.

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

The lattice-to-continuum calculation now predicts finite leading-order cubic
anisotropy for the current `1/|δ|²` shell weights. This sprint converts the
closed-form acoustic tensor into measurable numbers, verifies the force code,
and decides whether to retune shell weights or carry a quantified anisotropy
forward.

| # | Subtask | Owner | Falsifiable output |
|---|---|---|---|
| 1 | Closed-form `D(k)` for the canonical 6-neighbor axial-only lattice; verify that `D(k)` is diagonal in Cartesian, eigenframe is k-independent, and branch speeds along `[100]` are `c_L² = 1`, `c_T² = (1−α)`. **Done** (`paper/derivations/lattice_to_continuum.md` + `test-runs/sprint2_subtask9_d_of_k_diagonal/`). | `physics-derivation` | ✅ closed-form + 5 numerical checks pass at `< 1e-14` tolerance |
| 1b | *(retired)* Shell-weight retuning to recover cubic isotropy on the 26-neighbor stencil. The project has committed to the 6-neighbor axial-only lattice (backbone #15) and accepts lab-frame anisotropy as the structural source of the gauge sector under the dual-observer framework. | — | — |
| 2 | Long-wavelength dispersion: measure `ω([100])`, `ω([111])` at `|k|·a ∈ {0.05, 0.1, 0.2, 0.3}` on the canonical 6-neighbor lattice and compare to closed form. | `dispersion-analyst` | branch speeds match `c_L = 1`, `c_T = √(1−α)` along `[100]` and the degenerate triplet speed along `[111]` after `k→0` extrapolation |
| 3 | Lab-frame anisotropy quantification: measure the structural ratio `c_L([100]) / c([111])` at the operating point `α = 0.2` and confirm it matches `1/√((3−2α)/3) ≈ 1.075`. | `dispersion-analyst` | < 1% deviation from the analytic prediction |
| 4 | Band isolation on the canonical lattice: structurally, the L–T gap along `[100]` is exact and equal to `1 − √(1−α)` (10.6 % at `α = 0.2`); the lateral triplet along `[111]` is exactly degenerate. No numerical sweep needed unless characterizing non-linear leakage at finite amplitude. | `dispersion-analyst` (only if non-linear leakage probe is needed) | analytic / yes-no |

Pass criterion for sprint: dispersion measurements reproduce the closed-form
`c_L`, `c_T` formulas to within 1 % at the working point, after finite-`ka`
extrapolation.

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
| 14 | Triplet baryon search: parameter sweep over the partial-wave ansatz menu (hedgehog, Skyrme-twisted hedgehog, hedgehog + trace admixture, axis-triplet control) with width `R`, amplitude `A`, prestress `α`, and profile shape; see `baryon_simulation_roadmap.md` Phase 2 for the (J,L) content of each ansatz family and backbone #20 for the description protocol | `soliton-hunter` | top-10 candidate configs across the four ansatz families, with confinement score per (J,L) channel |
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
- Claim isotropy from the Cauchy relation alone. If shell weights are retuned,
  the retuning must be explicit and recorded as a model choice; if they are not,
  the finite anisotropy must be carried forward quantitatively.
- Postulate a global narrowband carrier sector to force Berry diagnostics to
  work. Per principles §4.3, narrowband is per-wavepacket.
- Hide preferred-frame effects behind "operational isotropy" hand-waving when
  they appear above the empirical bounds.

---

## How to claim a subtask

Sub-agents claim a subtask by setting their owner on the corresponding task in
the task list (see `TaskList`). Deliverables go under `test-runs/<run-id>/`
with a one-page Markdown report named `<subtask>_report.md`.
