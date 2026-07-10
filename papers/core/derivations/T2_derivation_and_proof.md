# T2 — Derivation and proof of the emergent Lorentz metric

**Status: the provable core of T2 is proved; the full dual-observer universality is
sharpened to a precise conjecture.** From the substrate alone:

- the emergent metric has **Lorentzian signature (3,1), derived from the prestress
  sign pattern `η`** — not imposed;
- the branch speeds are given in **closed form** (the calibration among
  `α_s, α_t, γ_t, κ, a` that T2 requested), verified against the exact lattice
  dispersion;
- the long-wavelength dispersion is a **massless, linear light cone**;
- each single-signal cone is **Lorentz-invariant**, and slower modes are **causal**
  (timelike, "massive/subluminal");

and — honestly — the lab frame is **genuinely anisotropic for `α_s>0`** (BACKBONE #8
says exactly this), so full cross-sector isotropy cannot come from parameter tuning
and remains the **load-bearing dual-observer conjecture**, now with a concrete
obstruction (anisotropy `∝ α_s`) and a clean `α_s→0` isotropic limit. Reproduced by
`derivations/t2_emergent_metric.py`.

---

## 1. What T2 requires

OPEN_TASKS: *"Show the inside observer measures a Minkowski cone despite lab
anisotropy and `α_t≠α_s`; derive the calibration among `α_s,α_t,γ_t,κ_s,κ_t,a`.
(Load-bearing dual-observer conjecture, BACKBONE #8.)"*

So T2 splits into a **derivable** part (the effective metric, its signature, the
calibration, the cone and its Lorentz symmetry) and a **conjectural** part (that
the *inside* observer measures a single isotropic `c` despite the lab anisotropy).
This document proves the former and precisely characterises the latter.

---

## 2. The substrate dispersion

Linearising the action about the tensioned vacuum gives the per-axis tangent
stiffness `M_μ = (1−α_μ)I + α_μ ê_μê_μᵀ` and `C_μ = η_μκ_μ M_μ` (Layer 1). With
time = direction 4, small fluctuations obey (4D block, `η_4=+1`):

```
κ_t · 2(1−cos ω) · M_4 ε = Σ_i κ_s · 2(1−cos q_i) · M_i ε,
```

a generalised eigenproblem `A(q) ε = λ M_4 ε`, `ω² = λ/κ_t`. Because all `M_μ` are
diagonal in the axis basis, the polarisations are the coordinate axes and the four
branch speeds follow in closed form (§4).

---

## 3. Lorentzian signature — derived from `η` (proved)

For any fixed polarisation `ε`, the small-`k` dispersion is the vanishing of the
**spacetime quadratic form** `k^μ H_{μν} k^ν` with

```
H_{μν} = diag( η_μ κ_μ m_μ ),   m_μ = εᵀ M_μ ε > 0.
```

Since `κ_μ, m_μ > 0`, the signature of `H` is exactly the **sign pattern of `η`**:

```
η = (−1,−1,−1,+1)  →  signs (−,−,−,+)  →  signature (3,1)  LORENTZIAN  →  a light cone exists.
η = (−1,−1,−1,−1)  →  signs (−,−,−,−)  →  signature (4,0)  EUCLIDEAN   →  definite form, no waves.
```

(Both verified numerically.) This is the central structural result: **the (3,1)
Minkowski signature is the prestress sign pattern**, not an independent postulate
and not an ambient metric — confirming ARCHITECTURE Layer 0 / BACKBONE #6. The
opposite sign of the temporal link is what makes the stationary equation hyperbolic
(wave-like) rather than elliptic.

---

## 4. Effective metric and calibration (derived, verified)

Solving the dispersion at long wavelength gives the four branch speeds (units
`a=1`, `γ_t=κ_t/κ_s`):

```
c_L² = 1 / [γ_t (1−α_t)]                (longitudinal / compression)
c_T² = (1−α_s) / [γ_t (1−α_t)]          (spatial-transverse, ×2 — the gauge/"light" modes)
c_4² = (1−α_s) / γ_t                    (e_4 amplitude/time-polarised mode)
```

Numerically (`α_s=0.6, α_t=0.9, γ_t=2.322`): `c_4=0.4151, c_T=1.3125 (×2),
c_L=2.0753`, matching the exact lattice generalised-eigenproblem to `1e-3`. This
is the **calibration among `α_s,α_t,γ_t,κ,a`** T2 asked for. For the light branch
the effective (mostly-minus) line element is

```
ds² = −dt² + c_T⁻² (dx²+dy²+dz²),   c_T² = (1−α_s)/[γ_t(1−α_t)].
```

**A notable exact fact:** the `e_4` amplitude mode is **perfectly isotropic**
(`c_4` identical in every spatial direction), because `e_4 ⊥` all spatial links, so
it couples only to the isotropic `(1−α_s)` part of the stiffness. This mode is the
natural proper-time / gravity-eigenstrain channel (ARCHITECTURE Layer −).

---

## 5. Massless linear cone (proved)

Each branch has `ω = c_b|q| + O(|q|³)`: linear (gapless) dispersion at long
wavelength — a **relativistic massless cone**. Computed for the light branch:
`ω/|q| = 1.3104, 1.3120, 1.3124` at `|q| = 0.2, 0.1, 0.05 → c_T=1.31252`. The
`O(|q|²)` piece is the lattice-curvature correction, suppressed by `(qa)²`.

---

## 6. Lorentz invariance of a cone + causality (proved)

A single isotropic quadratic cone `ω² = c²|q|²` has the full Lorentz group of its
metric as isometry. Verified: a boost `β=0.5` applied to 1000 null 4-vectors
preserves `ω'² − |q'|² = 0` to `1.3e-15`. Moreover a slower branch (`c_4/c_T =
0.316 < 1`) has `ω² − c_T²|q|² < 0` — it is **timelike / inside the light cone**,
hence **causal**: relative to the `c_T` signal cone the amplitude mode behaves as a
massive/subluminal excitation. So the model has a consistent causal structure with
`c_T` as the signal speed and slower modes as matter-like.

---

## 7. The genuine anisotropy, and the dual-observer conjecture

Here is the honest heart of T2. For `α_s>0` the lab dispersion is **genuinely
anisotropic** — as BACKBONE #8 explicitly acknowledges ("the lab observer … sees
direction-dependent wave speeds … this anisotropy is not retuned away"). Computed
(`α_s=0.6`):

```
n=(1,0,0):     0.415, 1.313, 1.313, 2.075
n=(1,1,0)/√2:  0.415, 1.313, 1.736, 1.736     ← birefringence (two transverse speeds)
n=(1,1,1)/√3:  0.415, 1.608, 1.608, 1.608
```

and the fractional directional spread of the fastest spatial branch:

```
α_s = 0.60 → 0.254 ;  0.40 → 0.155 ;  0.20 → 0.072 ;  0.05 → 0.017 ;  0.00 → 0.000
```

**The anisotropy vanishes only as `α_s→0`** — but `α_s→0` disables the transverse
stiffness that activates the gauge sector (T1 needs `α_s>0`). This is a real,
newly-sharp tension:

> Emergent spatial isotropy cannot be obtained by tuning, because the same `α_s`
> that would make the lab isotropic switches off the SU(3)/EM sector. Isotropy for
> the inside observer must therefore come from **rod/clock renormalisation**, not
> from the dispersion — which is exactly the load-bearing dual-observer conjecture.

**What is provable about the inside observer.** If the inside observer uses a single
signal mode for synchronisation and rulers (radar/Einstein synchronisation), the
measured speed of *that* mode is isotropic and constant by construction, and the
inferred metric is Minkowski (§6). This is the operational content of "rods and
clocks renormalise with the wave speed." The nontrivial, unproved requirement is
**cross-sector universality**: that all matter and all signals share one cone. In
the bare vacuum they do **not** (three distinct speeds `c_L, c_T, c_4`, plus
transverse birefringence). So:

- **Proved:** per-sector operational Lorentz invariance (one mode → isotropic
  metric), Lorentzian signature, calibration, massless cone, causality.
- **Conjectural (BACKBONE #8):** a *single universal* cone for all sectors for
  `α_s>0`. Candidate resolutions to pursue (not established here):
  1. ~~the gauge (Berry) modes of the nonlinear carrier are more isotropic than the
     bare acoustic branches~~ — **TESTED AND FALSIFIED**
     (`t2_gauge_mode_isotropy.py`, `GENUINE_TESTS_results.md`): the orientation-averaged
     gauge kinetic tensor has anisotropy ≈ 0.81, *larger* than the bare acoustic 0.254.
     The gauge modes are **not** more isotropic than sound; this route does not rescue
     emergent Lorentz invariance.
  2. self-confined Layer-5 solitons average over the lattice and couple to an
     effectively isotropic cone (isotropy emergent at the matter scale, not the
     lattice scale) — still open;
  3. the inside observer's full apparatus (built from one sector) renormalises all
     of that sector's anisotropy coherently, leaving no *intra*-sector detectable
     anisotropy, with cross-sector effects suppressed — still open.

  Net: with route 1 falsified, the emergent-Lorentz obstruction is confirmed to
  persist in the gauge sector; only the (harder) matter-scale (2) and observer-
  renormalisation (3) routes remain, and full universality is firmly still a conjecture.

---

## 8. Scope — proved vs owed

**Proved / derived (T2 core):**
- Lorentzian signature (3,1) from `η` (§3) — the metric structure is derived, not
  imposed.
- Effective metric and closed-form calibration `c_L,c_T,c_4(α_s,α_t,γ_t)` (§4),
  verified against the exact lattice dispersion; `e_4` mode exactly isotropic.
- Massless linear cone (§5); Lorentz-invariance of each cone and causal ordering of
  slower modes (§6).

**Conjectural / owed (the dual-observer part, BACKBONE #8):**
- a single universal isotropic `c` for all sectors at `α_s>0` — obstructed by the
  computed acoustic anisotropy/birefringence, removable only at `α_s→0`;
- the explicit inside-observer construction (rulers/clocks from one sector) that
  would realise operational isotropy, and the size of residual cross-sector
  Lorentz violation (a quantitative target, tied to the gauge-mode dispersion and
  to Layer-5 solitons).

T2 therefore delivers the effective **metric, its Lorentzian signature, and the
calibration** requested, proves per-sector Lorentz invariance and causality, and
turns the dual-observer split from a slogan into a sharp conjecture with a measured
obstruction and concrete resolution routes.

---

## 9. Relation to Weinberg–Witten — why approximate Lorentz is a feature

The approximate character of the emergent Lorentz invariance found here (§7) is not
only an honest limitation — it is a **necessary** ingredient for the rest of the
programme. The **Weinberg–Witten theorem** (Phys. Lett. B 96, 59, 1980) forbids
emergent massless spin-1 gauge bosons (and spin-2 gravitons) in a theory with an
*exactly* Lorentz-covariant conserved current. The emergent photon (T3) and gluons
(T4) are exactly such objects, so an *exactly* Lorentz-invariant substrate could not
produce them.

Because the emergent Lorentz symmetry here is only a **long-wavelength attractor**
— exact isotropy requires `α_s→0`, and the lab is genuinely anisotropic/birefringent
for the `α_s>0` needed by the gauge sector — the Weinberg–Witten premise fails, and
the emergent gauge bosons are allowed. This is the same escape used by Volovik, Wen,
and the analogue-gravity literature (see `RELATED_WORK.md`). So the T2 result should
be read two ways: (i) full cross-sector Lorentz universality is the load-bearing
dual-observer conjecture; but (ii) the very fact that Lorentz invariance is emergent
and approximate is what makes the emergent `U(1)×SU(3)` sector consistent with known
no-go theorems. The cross-link is developed on the T4 side (`T4_derivation_and_proof.md` §8).

---

## 10. Reproducibility

`derivations/t2_emergent_metric.py` — self-contained (numpy/scipy). Prints (A) the
signature-from-`η` comparison, (B) the calibrated branch speeds vs closed form,
(C) the massless-cone limit, (D) the anisotropy-vs-`α_s` table with birefringence,
(E) the boost-invariance and causality checks. All numbers above are its output.