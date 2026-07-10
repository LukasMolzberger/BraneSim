# OPEN_TASKS.md — where the theory stands and what is still owed

Rebuilt from the `chatgpt/` discussion (chronological through
`The_anisotropic_stiffness_tensor.md`). The previous OPEN_TASKS.md was deleted as
outdated. `BACKBONE.md` holds the non-negotiable principles; `ARCHITECTURE.md` holds
the derivational layer stack; this file holds the **live task/decision state**.

Central status:

> **T1 core is PROVED (D5 settled).** The decisive calculation of §5 has been done
> (`derivations/T1_derivation_and_proof.md`, `derivations/t1_su3_witness.py`): an
> exact finite-amplitude stationary periodic vacuum carrier (`‖∇S‖≈2e-16`) has an
> isolated near-degenerate rank-3 Bloch carrier whose gauge-invariant traceless
> Wilczek–Zee curvature spans **full `su(3)` (Lie rank 8)**, robustly across
> parameters. The mechanism is identified: real noncommuting stiffness frames give
> the 3 `so(3)` generators; complex Bloch phases supply the 5 symmetric-traceless +
> Cartan generators that complete `so(3)→su(3)`. A fair real-frame control returns
> exactly `so(3)` (rank 3), so the test discriminates. What remains for the full
> QCD claim is T4 (Yang–Mills action), T7 (quantitative split), and confirming the
> EM `U(1)` identification (T3) — see §4.

---

## 1. Decisions now SETTLED (reflected in ARCHITECTURE/BACKBONE)

- **A. The bare linear vacuum is not the place for SU(3).** The gauge sector lives in
  `D_{R̄}(k)`, the Bloch fluctuation operator about a *finite-amplitude nonlinear*
  background — never `D_vac`. (Not a no-go; LESSONS_LEARNED #3.)
- **B. The carrier is the complexified transverse fluctuation bundle**, not the fixed
  spatial axes. `ℂ⁴ = ℂ_∥ ⊕ ℂ³_⊥`; the moving `T⊥` of the link direction `Q̂` is the
  color-carrier candidate.
- **C. 4D base, rank-3 internal carrier.** `𝒜_μ^a_b`, `μ=1..4` (base), `a,b=1..3`
  (fiber); `𝒜 ∈ Ω¹(M₄, u(3))`. 4D base does not imply SU(4).
- **D. Interpretation A.** EM `U(1)` = trace `a_μ = ⅓ Tr 𝒜_μ`; color `SU(3)` =
  traceless `B_μ`. (Candidate identification until T3/T4 close.)
- **E. The fourth mode is excluded from the *carrier*, not the substrate.** It is the
  longitudinal/time-amplitude direction, gapped/frozen relative to the triplet.
- **F. Split the single `α`** → `α_s` (spatial), `α_t` (temporal), keeping
  `η = (−1,−1,−1,+1)` as sign only. `r_i = α_s a`, `r_4 = α_t a`.
- **G. `α_t ≠ 1`** (and `α_t ≉ 0`): use `0 < α_t < 1`, likely `α_t = 1 − ε_t`.
  `α_t = 1` kills temporal transverse stiffness and over-decouples time/amplitude.
- **H. Rest distance matters and is not prestress.** `ρ_μ ~ η_μ κ_μ a (1 − α_μ)`;
  prestress is the *mismatch* `a − r_μ`. Held spacing stays `a_s = a_t = a`.
- **Anisotropic stiffness tensor is the mechanical bridge.**
  `C_{nμ}^{AB} = η_μ κ_μ[(1−r_μ/L̄)δ^{AB} + (r_μ/L̄)Q̂^AQ̂^B]` becomes the
  matrix-valued Bloch hopping coefficient; its link-to-link noncommutativity is the
  seed of the non-Abelian sector.
- **L. Near-degenerate triplet, not exact degeneracy** (was Open Decision 2).
  Working condition `λ_1 ≈ λ_2 ≈ λ_3` with a large outside gap (`ε ≪ Δ_outside`); the
  composite rank-3 subspace, not a strictly degenerate eigenvalue, carries the U(3).
- **M. Color is UNIVERSAL** (was Open Decision 3). The `U(1)×SU(3)` sector is a gauge
  sector of the substrate **vacuum** — it exists over empty space, not only inside
  matter. Soliton/matter structure enters at a **much higher layer** (Layer 5),
  carrying color charge and coupling to the already-universal field; matter does not
  *create* color. Consequences that follow from this choice:
  - **⚠ QUALIFIED by the D7 test (`GENUINE_TESTS_results.md`).** The color-carrying
    helix is a stable but **higher-energy texture** (+42.5% vs straight vacuum), and
    the straight ground state has *no* color (T1). So "universal over empty space /
    ground-state Yang–Mills" is **not established**; the supported claim is "su(3) on a
    stable substrate texture." Revisit this decision's wording, or find a
    color-carrying background at/below the straight-vacuum energy.
  - The gauge background `R̄` is a **universal nonlinear periodic vacuum carrier**
    (the vacuum microtexture), **not** a soliton core and **not** the bare straight
    lattice. (Settles Decision K toward the "universal periodic carrier" option; the
    bare hypercubic straight lattice is *not* the gauge vacuum.)
  - The connection lives over the Brillouin `k`-space of that periodic carrier.
    (Largely settles the old `k`-space vs `(k,X,λ)` question — Decision 4 — toward
    `k`-space; the adiabatic `(k,X,λ)` extension was the soliton-local route and is
    deferred to the matter layer.)
  - T4 target becomes a **universal Yang–Mills field over empty space** (the strong
    claim), not a soliton-bound color-like sector.

## 2. Provisional (working hypotheses, still need proof)

- **I. CONFIRMED (T7): `α_t` (with `γ_t`) creates the gap that prevents SU(4).**
  Computed from `D_{R̄}(k;α_s,α_t,γ_t)`: rank-3 cluster isolation grows toward
  `α_t→1` and is ~2× that of any rank-4 cluster → `λ_1≈λ_2≈λ_3` with `λ_4`
  separated, i.e. **U(3) not U(4)**. Analytic vacuum scaling `(1−α_s)α_t/α_s`; `γ_t`
  cancels from the split (sets only the T2 scale). See `derivations/T7_*`.
- **J. Strain sets SU(3) *coupling strength*, not color's existence** (revised under
  Decision M). Color exists universally; what scales with local nonlinear strain is
  the *magnitude* of the traceless curvature. In weak, long-wavelength regions `∂_k
  P_3` is small and the SU(3) sector is soft/perturbative; where the link frames vary
  strongly (high energy density, short range, `|Δu|/a ~ 1`) it becomes large and
  nonperturbative. So strong-coupling/confinement-like dynamics is short-range (near
  matter) while color itself is universal — analogous to QCD asymptotic freedom. Must
  derive the coupling-vs-strain scale from `α_s, α_t, γ_t, a`.
- **K. SETTLED → see Decision M.** The gauge background is the universal nonlinear
  periodic vacuum carrier. Periodic Bloch theory is not merely a "toy first" step — it
  *is* the physical carrier. (What remains is the vacuum-carrier justification below.)

## 3. Decisions still OPEN

> Decisions 2 (exact vs near-degenerate), 3 (universal vs matter-bound), and 4
> (`k`-space vs `(k,X,λ)`) are now **settled** — see §1 L and M. Labels below are kept
> stable for reference.

- **D1 — Is `γ_t = κ_t/κ_s` needed, or does `γ_t = 1` suffice?** Keep `γ_t` available
  as a light-cone calibration parameter for now; test `γ_t = 1` first.
- **D5 — SETTLED: the `su(3)` curvature spans all 8 generators.** Proved by
  explicit witness + genericity (`derivations/T1_derivation_and_proof.md`). Real
  transverse frame alone gives `so(3)` (3, confirmed by control); the 5
  symmetric-traceless + Cartan directions come from the complex Bloch phases acting
  on the noncommuting real anisotropic stiffness frames — exactly as conjectured.
  Full `su(3)` (rank 8) is the generic outcome; `so(3)` is a measure-zero
  (real-frame) special case.
- **D6 — SETTLED (T12).** `η_4 = +1` reads as genuine prestress: the temporal link
  is prestressed (`ρ_4 = +κ_t a(1−α_t)`) exactly like the spatial links, and its
  continuum/kinetic limit is the kinetic term `½κ_t a²[(∂_t u^4)²+(1−α_t)Σ(∂_t u^i)²]`
  (emergent inertia `κ_t a²`). The `(1−α_t)` anisotropy distinguishes it from the
  `r_4≈0` kinetic-sign artifact; `η_4=+1` makes `S=T−V` (Verlet dynamics, waves at
  `c_T`). See `derivations/T12_*`.
- **D7 — What *is* the universal nonlinear periodic vacuum carrier?** (New, forced by
  Decision M.) If color is universal then the gauge vacuum is a finite-amplitude
  nonlinear periodic state of the lattice, not the straight tensioned lattice
  (`A=0`) of LESSONS_LEARNED #1. Owe: identify this carrier, its period/amplitude,
  and *why it is universally present* — is it the true ground state, or a persistent
  microtexture, and what selects it? Must not conflict with the periodic-BC stable
  vacuum result.
  - **PARTIAL RESULT (`derivations/d7_vacuum_selection.py`, `GENUINE_TESTS_results.md`):**
    the specific helix carrier is **linearly (elastically) stable** (all link brackets
    PD, `V(q)` PSD, `ω²≥0`) but stores **+42.5% more elastic energy** than the straight
    vacuum → it is a **stable finite-energy TEXTURE, not the ground state**. Since the
    straight ground state has no color (T1), this **qualifies Decision M**: the honest
    claim is "su(3) on a stable substrate texture," not "the ground-state vacuum is
    Yang–Mills." Genuine (Floquet) selection of *the* carrier remains open; the
    "persistent microtexture" option is the current honest reading.

## 4. Task ledger (paper-section burdens)

- **T1 — SU(3) Wilczek–Zee sector.** CORE PROVED (see central status, D5). Route
  `C_{nμ}[R̄] → D_{R̄}(k) →` isolated rank-3 `P_3(k) → 𝒜_μ → 𝒢_μν → su(3)` executed
  end-to-end on an exact stationary universal periodic carrier (Decision M); the
  gauge-invariant traceless curvature spans full `su(3)` (Lie rank 8). Remaining T1
  polish is downstream: T4 (action), T7 (split), T3 (EM identification).
- **T2 — Emergent Lorentz / effective metric.** CORE PROVED; universality sharpened
  (`derivations/T2_derivation_and_proof.md`, `derivations/t2_emergent_metric.py`).
  Derived: the effective metric has **Lorentzian signature (3,1) from the `η` sign
  pattern** (flip `η_4`→−1 → Euclidean, no waves); closed-form calibration
  `c_L²=1/[γ_t(1−α_t)]`, `c_T²=(1−α_s)/[γ_t(1−α_t)]`, `c_4²=(1−α_s)/γ_t` (verified
  vs exact lattice); massless linear cone; each cone Lorentz-invariant (boost test
  1e-15) with slower modes causal/subluminal; the `e_4` amplitude mode is exactly
  isotropic. **Honest scope:** the lab is genuinely anisotropic for `α_s>0`
  (birefringence; spread→0 only as `α_s→0`, which disables the gauge sector), so a
  single universal cone for all sectors is not achievable by tuning — it remains the
  load-bearing dual-observer conjecture (BACKBONE #8), now with a concrete
  obstruction and resolution routes (nonlinear-carrier gauge-mode dispersion;
  Layer-5 soliton-scale isotropy; single-sector rod/clock renormalization).
- **T3 — Faraday U(1).** PROVED (`derivations/T3_derivation_and_proof.md`,
  `derivations/t3_maxwell.py`). Abelian holonomy of the substrate carrier: `f_μν`
  and Bianchi automatic (verified); U(1) gauge invariance + locality + power
  counting force Maxwell `−¼e⁻²f²` (mass term forbidden → **massless photon**,
  pure-gauge `f=0`); `1/e²` = substrate abelian quantum metric, **positive &
  finite**; inhomogeneous `∂^μf_μν=e²J_ν` with `∂·J=0`. Key result: the pristine
  vacuum is **EM-flat** (all abelian Berry curvatures ≈0, a PT-like symmetry — empty
  space has no background EM field), and EM curvature **switches on under
  deformation/matter**; consistent with `su(3)≠0` (off-diagonal, non-abelian).
  **Charge = U(1) vortex winding** (`∮=2πn` verified, `π₁(U(1))=ℤ`), conserved;
  photon continuous/non-quantized (BACKBONE #13). Still owed: exact `1/e²` norm,
  the T2 metric contraction, and the Layer-5 soliton current `J^μ`.
- **T4 — Yang–Mills.** FORM PROVED, coefficient derived
  (`derivations/T4_derivation_and_proof.md`, `derivations/t4_yang_mills.py`). The
  emergent gauge invariance of T1 + locality + power counting force the leading
  coarse-grained action to be uniquely `−¼g⁻² G^a_μν G^{aμν}` (mass term forbidden
  → gluon massless, verified: pure-gauge `G=0`). The substrate supplies `1/g²` as
  the integrated non-Abelian quantum metric of the carrier — **positive, finite**,
  and **strain/scale-dependent** (derives Decision J). Self-interaction is the
  unique `su(3)` vertex set (structure constants extracted; Casimir
  `f^{acd}f^{bcd}=2N δ`, `N=3`). Universal field over empty space (Decision M).
  Still owed (nonperturbative/cross-task): exact `1/g²` normalization + β-function,
  confinement/mass gap, the T2 metric contraction, and Layer-5 colored-soliton
  sources.
- **T7 — Quantitative sector split.** ESTABLISHED (`derivations/T7_derivation_and_proof.md`,
  `derivations/t7_sector_split.py`). Computed from `D_{R̄}(k;α_s,α_t,γ_t)`: (1) **rank
  selection U(3) not U(4)** — rank-3 cluster isolation (~5) ≫ rank-4 (~2.4),
  quantitatively; (2) **`λ_4` separation opened by `α_t`** (Decision I confirmed;
  best toward `α_t→1`), matching the analytic vacuum scaling `(1−α_s)α_t/α_s`;
  (3) algebra **generically `su(3)`** across the `(α_s,α_t)` plane, `so(3)` only at
  measure-zero `k*`-accidents (consistent with T1 genericity); (4) **U(1):SU(3)
  coupling split ≈1:8 total** (`1/e²≈1.6`, `1/g²≈11`), ~equipartitioned per
  generator, stable in the parameters. Clean structural result: **`γ_t` cancels from
  the split** (sets only the T2 scale); the split is governed by `(α_s,α_t)`. Owed:
  absolute `1/e²,1/g²` normalization (scheme, shared with T3/T4) and a free-`γ_t`
  stationary scan.
- **T12 — Temporal prestress.** PROVED (`derivations/T12_derivation_and_proof.md`,
  `derivations/t12_temporal_prestress.py`). The `r_4 = α_t a` temporal link is a
  genuine prestressed spring whose continuum limit **is** the kinetic term
  `T = ½κ_t a²[(∂_t u^4)² + (1−α_t)Σ_i(∂_t u^i)²]` with **emergent inertia** `κ_t a²`
  (no postulated mass). The **`(1−α_t)` anisotropy** is the fingerprint separating
  genuine prestress from the old `r_4≈0` **kinetic-sign artifact** (isotropic).
  `ρ_4 = +κ_t a(1−α_t)` is a finite prestress symmetric with `ρ_i`; `η_4 = +1` makes
  `S = T − V`, so 4D-block stationarity is Störmer–Verlet evolution carrying waves at
  `c_T` (measured 0.910 vs T2 0.913); the all-minus (Euclidean) sign blows up. Owed:
  full nonlinear temporal dynamics; proper-time/inside-observer clock (ties to T2);
  matching `κ_t a²` to a mass scale (Layer 5).

## 5. The concrete calculation (DONE — decided T1)

**This calculation has been executed** — see `derivations/T1_derivation_and_proof.md`
and `derivations/t1_su3_witness.py`. Result: steps 1–8 give full `su(3)` (rank 8).
The background used is exactly the circularly-polarized helix below, specialized to
propagation in the `(1,4)` plane and polarization in the `(2,3)` plane so that all
link lengths are `n`-independent and the force balance closes in one scalar
equation (solved for `γ_t`), yielding an exact critical point (`‖∇S‖≈2e-16`). The
steps, kept for reference:

Do this on a minimal finite-amplitude **periodic** background — the candidate
universal vacuum carrier (Decision M), e.g. a 4D helical/twisted link-frame texture
`Q̄_{nμ} = a e_μ + A[cos(K·n) v_μ + sin(K·n) w_μ]`, `v_μ, w_μ ⊥ e_μ`, real, periodic,
finite-amplitude:

1. choose `R̄` (stationary, `∂S/∂R[R̄] = 0`); compute all `Q̂_{nμ}`;
2. build all anisotropic stiffness tensors `C_{nμ}^{AB}`;
3. assemble `D_{R̄}(k; α_s, α_t, γ_t)`;
4. test for an isolated rank-3 band cluster (`λ_1≈λ_2≈λ_3`, `λ_4` separated);
5. compute `P_3(k)` and check `∂_{k_μ} P_3 ≠ 0`;
6. compute `𝓕_{μν} = i P_3[∂_μ P_3, ∂_ν P_3] P_3`;
7. split `𝓕_{μν} = f_{μν} I_3 + 𝒢_{μν}`;
8. **test `rank Lie⟨𝒢_μν, [𝒢_μν, 𝒢_ρσ], …⟩`**: 0 = trivial, 3 = only `so(3)`,
   **8 = genuine `su(3)` candidate**;
9. only then attempt the effective Yang–Mills action (T4).

## 6. Consistency notes / do-not-regress

- Do **not** reintroduce the analytic-signal / worldtube complex structure
  `Ψ ≈ ξ + (i/ω₀)ξ̇` (appears in the early `Spin_One_Half.md`). The complex phase is
  the character of the 4D lattice translation group (ARCHITECTURE Layer 2). The ℤ₂
  spin-1/2 holonomy conclusion is kept; its stated *mechanism* is not.
- Keep periodic BC for every nonlinear run (LESSONS_LEARNED #1).
- Never state a consequence of linearizing the Pythagorean length as a no-go
  (LESSONS_LEARNED #3).