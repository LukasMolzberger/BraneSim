# α-Separability Experiment — linear-layer spec

**Status:** Track A + P1 **DONE & confirmed to machine precision** (2026-06-03,
`results/`, `branesim/diagnostics/alpha_separability.py`, 33 passing tests).
Track B B.1 derivation DONE (`derivation_H_eff.md`); P2/P3 baseline checks and the
L5 spin-½ ℤ₂ target remain. **Layer:** linear substrate (L1/L3), prerequisite for
the L5 soliton-window cross-check.

**P1 caveat:** confirmed via the standing-wave ω_T/ω_L ratio (+ analytic group
velocity), which shares the closed-form D(k) used in the derivation — so it is a
tight self-consistency check, not yet an independent translating-envelope centroid
measurement (which would probe the M_jl dispersion term). Strong pass; flagged for
rigor.
**Owner routing:** Track A → `dispersion-analyst`; Track B derivation →
`physics-derivation`; Track B diagnostic → `berry-validator`.

Canonical bindings: `BACKBONE.md` #16 (U(1)³→U(3) crossover), #18 (per-wavepacket
complex envelope), #19 (sector mapping); `OPEN_PROBLEMS.md` (gauge thread).
Prior art this builds on: `test-runs/sprint2_subtask9_d_of_k_diagonal/`
(D(k) diagonality + identity holonomy on the real eigenframe, certified for all α).

---

## 0. Motivation — separate the observables, re-couple through α

The emergent properties (U(1)×SU(3) split, spin ½, self-confinement) cannot be
*tuned* independently — they all hang off the single prestress knob **α** and the
same lattice. But they can be *measured* independently because each lives at a
different order of the dynamics and on a different carrier. The discipline:

1. **Isolate** each observable at the minimal substrate that supports its carrier.
2. **Re-couple** through α as the consistency check: the single α that produces the
   observed color split must also be the α that produces the correct Berry
   holonomy *and* the Q>0 confinement window (L5). Non-overlapping α-windows ⇒
   the theory fails — a test invisible to any single-soliton run.

## 1. The premise correction — this is TWO objects, not one diagonalization

The naive plan ("diagonalize D(k) once, read off both color split and Berry
phase") is **wrong** and contradicted by certified repo results:

- `D(k)` is diagonal in the fixed Cartesian basis at every (k, α); eigenvectors
  are k-independent (`branesim/core/conventions.py:d_of_k_eigenvalues`,
  lines 175–211). Hence the Berry/WZ connection on its **real eigenframe is
  identically zero** for all α (sprint2_subtask9, structural — not numerical).

| Observable | Carrier | Status |
|---|---|---|
| U(1)×SU(3) split | eigenvalue *spread* of the real `D(k)` | closed-form / solved |
| spin ½ / Berry | complex envelope `Ψ∈ℂ³` eff. Hamiltonian `H_eff(k₀,α)` | **not derived, not coded** |

---

## 2. Track A — U(1)×SU(3) split (closed-form; confirmatory)

**Observable.** Continuous-in-α trace/traceless decomposition of the lateral 3×3
block, with explicit projection operators (sprint2 reported only static gap
ratios). With `λ_a(k)=(2k_s/ρ)[α h_a+(1−α)H]`, `h_a=1−cos(k_a a)`, `H=Σ_a h_a`:

```
trace (U(1), dilatational):    λ̄        = (2k_s/ρ) H (1 − 2α/3)
traceless (SU(3), shear):      λ_a − λ̄  = (2k_s/ρ) α (h_a − H/3)
```

The **traceless / SU(3) content is exactly ∝ α**; the trace / U(1) part ∝ (1−2α/3).
SU(3)-to-U(1) ratio:

```
ρ_SU3(k̂, α) = g(k̂) · √3 α / (3 − 2α),   g(k̂) = √(Σ_a (h_a − H/3)²) / H
```

- → 0 at α=0 (full U(3) degeneracy); → maximal at α=1 (decoupled U(1)³).
- `g([111]) = 0` (certified triplet degeneracy); `g([100])` maximal.
- At default α=0.2 ⇒ ρ_SU3 ≈ 0.13·g (≈13% along [100]).

**Method.** (i) Evaluate the closed form on an α-grid and a k̂-grid (map the
directional anisotropy g(k̂)). (ii) Numerical confirmation against the actual
lattice operator, reusing the sprint2_subtask9 verifier (the trace/traceless
projectors P_U1 = (1/3)1 1ᵀ, P_SU3 = I − P_U1 applied to the measured 3×3 block).
**Cost:** ~half a day. **Deliverable:** `α-curve(g)`, projection operators for
downstream color diagnostics.
**Falsifier:** any deviation of the measured traceless content from the linear-in-α
law, or nonzero g along [111], breaks #16.

---

## 3. Track B — Berry / spin (DERIVATION DONE 2026-06-03 — mostly analytic)

**Derivation gate resolved.** `derivation_H_eff.md` (this dir). `H_eff` is derived,
and its coefficient matrices (ΔΩ, V_j, M_jl) are **all diagonal in the same
k-independent Cartesian frame as D(k)**. This forces three conclusions that
largely *settle* Track B analytically rather than by a sweep:

### B.1 k-space holonomy is trivial for all α (consistency check, berry-validator)
The rotating-wave `i` does **not** create base-space curvature: `A(k₀)≡0`,
`F≡0` ∀α (inherited from sprint2_subtask9). **Prediction P2:** Fukui–Hatsugai
plaquette holonomy of the H_eff eigenbundle = identity to machine precision,
every α. Fail if `‖U_Γ−𝟙‖ > 1e-10`. (No α-window — this is a flat consistency
check, not a measurement.)

### B.2 The surviving linear gauge object — fibre-internal U(1)/SU(3) (α-dependent)
Non-triviality lives in the **ℂ³ fibre**, transporting in **physical (x,t)**, not
in k: the per-axis U(1) phase (magnetic-curl channel) + relative-phase SU(3).
This is the real α-dependent linear gauge content; its α-dependence enters through
the group velocities `v^{(a)}_j ∝ [α δ_{aj}+(1−α)]`. **Prediction P1 (cheap,
dispersion-analyst):** transverse/longitudinal envelope-drift ratio = `√(1−α)`
(= 0.894 at α=0.2) for a [100] carrier at k₀a=π/4. Fail if >5% off at W/a≥8.

### B.3 Spin ½ is NOT a linear-layer object — it migrates to L5
`Ψ∈ℂ³` carries the **vector (J=1, spin-1)** rep of SO(3): a 2π rotation → `+𝟙`
(phase 0), α-independent. No half-integer rep exists in a real displacement
triplet. **Prediction P3 (decisive, berry-validator):** rotating a band-isolated
wavepacket carrier through 2π gives geometric phase **0** at all α; any measured
π falsifies the derivation. Spin ½ requires a **soliton-layer π₁(SO(3))=ℤ₂
holonomy** (rigid hedgehog rotation, orientation–isospin locking — standard
Finkelstein–Rubinstein). This is the expected Skyrme structure (phonon triplet =
spin-1 meson; baryon spin ½ = topological), a coherence *feature*, and becomes a
**new L5 target**, not a linear-layer measurement.

---

## 4. The joint α-window test (the payoff)

**Revised after the derivation.** The "Berry curvature window" is gone — at the
linear layer the gauge curvature is flat in k (B.1) and spin is α-independent
(B.3). So the linear layer contributes ONE real α-window (the color split, A),
plus the soliton-layer spin ½ check (B.3 → L5). The joint test is therefore
between the color window, the confinement window (L5, `Q = 8k_sα/a² > 0`), and
the soliton ℤ₂ spin holonomy:

```
α:  0 ───────────────── 0.2 ───────────────── 1
  SU(3) content (A, linear):   small ↑──────────────↑ large   (∝ α)
  confinement Q>0 (L5):        [ grows with α ]
  spin-½ ℤ₂ holonomy (L5):     [ requires odd-winding soliton at this α ]
                                 ↑ must all be simultaneously satisfiable at one α (default 0.2)
```

The discipline is intact — but the derivation has shown that **two of the three
windows are soliton-layer (L5), not linear.** The linear layer's job is now: (i)
nail the color α-curve (A), (ii) confirm the flat-curvature / spin-1 baseline (P2,
P3) so any soliton-layer π holonomy is unambiguously topological, not a
linear-layer artifact.

---

## 5. Artifacts
- `diagnostics/alpha_separability.py` — returns both α-curves (Track A closed-form
  now; Track B once `H_eff` exists). Read-only, no back-reaction (principles §).
- `derivation_H_eff.md` (this dir) — output of the physics-derivation gate.
- `results/` — α-curves (gitignored heavy arrays per test-runs policy; text/JSON
  summary versioned).