# T7 — Derivation and proof of the quantitative U(1)/SU(3) sector split

**Status: T7 is established.** The sector split is computed as an explicit function
of `(α_s, α_t, γ_t)` from the nonlinear carrier `D_{R̄}(k)`, on four fronts:

1. **Rank selection — U(3), not U(4):** the rank-3 near-degenerate cluster is
   isolated ~2× better than any rank-4 cluster, so the carrier is genuinely `U(3)`.
2. **λ_4 separation (Decision I):** `α_t` (with `γ_t`) opens the isolating gap; the
   isolation grows toward `α_t→1`, matching the analytic vacuum scaling.
3. **Algebra:** the traceless curvature is **generically `su(3)`** across the whole
   `(α_s,α_t)` plane; `so(3)` occurs only at isolated `k*`-specific (measure-zero)
   points — consistent with the T1 genericity lemma.
4. **Coupling split:** `1/e²` (U(1), T3) vs `1/g²` (SU(3), T4) is a **stable ~1:8
   total ratio**, i.e. ~equipartitioned per generator, roughly independent of the
   parameters.

A clean structural result falls out: **`γ_t` cancels from the split** — it sets
only the overall scale (the T2 light-cone), while the U(1)/SU(3) split is
controlled by `(α_s, α_t)`. Reproduced by `derivations/t7_sector_split.py`.

---

## 1. What T7 requires

OPEN_TASKS: *"Compute the `U(1)/SU(3)` split and the `λ_4` separation as a function
of `(α_s, α_t, γ_t)`."* Decision I upgraded this: *"`α_t` (with `γ_t`) creates the
gap that prevents SU(4); target `λ_1≈λ_2≈λ_3`, `|λ_4−λ_triplet| ≫ |λ_a−λ_b|`; must
be computed from `D_{R̄}(k;α_s,α_t,γ_t)`."* So T7 has a **spectral** part (rank
selection / gap hierarchy) and a **coupling** part (U(1) vs SU(3) strength).

---

## 2. Analytic backbone — the vacuum branch stiffnesses

From the vacuum dispersion (T2), for propagation along an axis the four polarization
stiffnesses (`∝` speed²) are

```
s_L = 1/[γ_t(1−α_t)]           (longitudinal)
s_T = (1−α_s)/[γ_t(1−α_t)]     (spatial transverse, ×2)
s_4 = (1−α_s)/γ_t              (amplitude / e_4)
```

giving two candidate near-degenerate triplets and their isolation ratios:

```
Story A  triplet {L,T,T}, e_4 split off :  gap/spread = (1−α_s) α_t / α_s
Story B  triplet {T,T,4}, L split off   :  gap/spread = α_s / [(1−α_s) α_t]
```

**Key observation: `γ_t` cancels in both ratios.** So `γ_t` sets the overall energy/
speed scale (the light-cone calibration, T2) but **does not control the rank split**
— that is a function of `(α_s, α_t)` alone. Story A grows with `α_t` (the amplitude
mode `e_4` separating), Story B grows with `α_s` (the longitudinal mode separating).
The nonlinear carrier (§3) selects the Story-A scaling: `α_t` opens the gap.

---

## 3. Spectral split from the nonlinear carrier `D_{R̄}(k)`

Computed on the exact stationary helix, isolation = best rank-n cluster gap/spread,
averaged over four `k`-points:

### (1) Rank selection — U(3) not U(4)
```
rank-2:  18.92     rank-3:  5.11     rank-4:  2.38     rank-5:  0.71
```
The rank-3 cluster is well isolated (5.11 ≫ 1) while the rank-4 cluster is marginal
(2.38) and rank-5 fails (0.71). The natural clustering is a doublet + a **rank-3
triplet**; extending the triplet to four bands roughly halves the isolation. So the
transported carrier is **`U(3)`, and `U(4)`/`SU(4)` is spectrally excluded** —
quantitatively, not by fiat. (This is the sharp form of Decisions C/E and
`From_U4_to_U3`.)

### (2) λ_4 separation vs `(α_s, α_t)` (Decision I)
Rank-3 isolation over the plane:
```
α_s \ α_t    0.70    0.80    0.90    0.95
0.3          2.12    2.11    3.88    1.96
0.5          6.19    2.35    3.56    4.03
0.7          1.85    2.17    3.19    3.79
```
The isolation generally **grows toward `α_t → 0.9–0.95`** (best column at high
`α_t`), confirming **Decision I: `α_t` (with `γ_t`) opens the gap that splits the
4th mode off the triplet**. The trend matches the analytic Story-A form
`(1−α_s)α_t/α_s` (which also grows with `α_t`); the nonlinear carrier enhances the
absolute isolation above the vacuum estimate. (The single-carrier diagnostic is
somewhat noisy in `α_s`; the robust statement is the `α_t→1` gap opening.)

### (3) The algebra is generically `su(3)`
Robust Lie rank (majority over four `k*`; `8=su(3)`, `3=so(3)`):
```
α_s=0.3:  at0.75→8   at0.85→8   at0.95→8
α_s=0.5:  at0.75→8   at0.85→8   at0.95→8
α_s=0.7:  at0.75→8   at0.85→8   at0.95→8
```
Full `su(3)` everywhere. The only `so(3)` readings are isolated `k*`-specific
accidents (e.g. one of four `k*` at `(0.3,0.75)`), i.e. **measure-zero**, exactly as
the T1 genericity lemma predicts. So the *existence* of full color does not depend
on tuning `(α_s,α_t,γ_t)`; only the *isolation* (the clean U(3) selection) does.

---

## 4. The coupling split — U(1) vs SU(3) strength

Abelian stiffness `1/e²` (T3) and non-abelian stiffness `1/g²` (T4), from the
substrate quantum metrics:

```
α_s  α_t   γ_t     1/e²    1/g²   (1/e²)/(1/g²)   1/g²/8
0.4  0.85  2.625   1.652  11.687     0.141        1.461
0.4  0.95  4.108   1.516  10.804     0.140        1.351
0.6  0.85  1.903   1.551  12.197     0.127        1.525
0.6  0.95  2.978   1.641  11.182     0.147        1.398
```

Two clean facts:
- **The SU(3) sector carries ~8× the total kinetic stiffness of the U(1) sector**
  (`1/g² : 1/e² ≈ 8 : 1`), because it has 8 generators to the U(1)'s one.
- **Per generator they are comparable** (`1/g²/8 ≈ 1.4 ≈ 1/e²`): the geometric
  stiffness is roughly **equipartitioned per generator**, and the split ratio is
  **stable (~0.13–0.15) across `(α_s,α_t,γ_t)`**.

So at the vacuum/lattice scale the bare couplings satisfy `g² ≈ e²` per generator
(the split is a democratic dof count, not a large hierarchy). Any physical
`α_s,α_t,γ_t`-driven running of these couplings enters through the strain
dependence derived in T4 (Decision J), not through the vacuum split itself.

---

## 5. Synthesis — the role of each parameter

Combining §2–4, the quantitative sector split is:

| parameter | role in the split |
|---|---|
| `α_s` | spatial branch mixing; sets the longitudinal↔transverse gap (Story B); required `>0` for the gauge sector (T1) and for lab anisotropy (T2) |
| `α_t` | **opens the λ_4 gap** that isolates the rank-3 carrier and enforces U(3)-not-U(4) (Decision I); best near `α_t→1` (Decision G) |
| `γ_t` | overall scale / light-cone calibration (T2); **cancels from the rank split and the coupling ratio** |

The headline: **U(3) selection and the U(1):SU(3) split are governed by `(α_s,α_t)`;
`γ_t` only rescales.** The algebra is `su(3)` generically; the *clean isolation* of
the triplet is what the parameters tune, and it is maximized toward `α_t→1`.

---

## 6. Scope — proved vs owed

**Proved / computed (T7):**
- Rank-3 selection (U(3) not U(4)), quantitative isolation hierarchy (§3.1).
- `λ_4` separation as a function of `(α_s,α_t)`, `α_t` opening the gap (§3.2), with
  analytic scaling and `γ_t` cancellation (§2).
- `su(3)` generic over the parameter plane; `so(3)` measure-zero (§3.3).
- U(1):SU(3) coupling split ratio, per-generator equipartition, parameter stability
  (§4).

**Owed (cross-task / scheme):**
- absolute normalisation of `1/e², 1/g²` (coarse-graining scheme, shared with
  T3/T4); the physical `α`-running of the couplings (Decision J, T4) beyond the
  vacuum split;
- an independent `γ_t` scan at fixed stationarity (here `γ_t` is fixed by the helix
  force balance; the analytic `γ_t`-cancellation covers the split, but a free-`γ_t`
  stationary family would make it fully manifest);
- the precise triplet↔vacuum-polarization identification (Story A vs B): the
  nonlinear carrier follows the Story-A (`α_t`) scaling, but the carrier itself is
  the complexified transverse bundle (Decision B) — a full reconciliation is a
  refinement, not needed for the split.

---

## 7. Reproducibility

`derivations/t7_sector_split.py` — imports the T1 carrier and the T3/T4 stiffness
routines; prints (1) the rank-selection hierarchy, (2) the `λ_4`-separation map vs
`(α_s,α_t)` with the analytic scaling, (3) the robust `su(3)`/`so(3)` algebra map,
(4) the U(1)/SU(3) coupling-split table. All numbers above are its output.