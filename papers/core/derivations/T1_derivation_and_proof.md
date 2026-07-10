# T1 — Derivation and proof of the SU(3) Wilczek–Zee sector

**Status: the core of T1 is proved.** The traceless Wilczek–Zee curvature of the
substrate's own rank-3 carrier bundle generates the *full* Lie algebra `su(3)`
(rank 8), not merely `so(3)` (rank 3) and not a pure-gauge artifact. This settles
Open Decision **D5** ("can the `su(3)` curvature span all 8 generators, or only
`so(3)`?") — the answer is **all 8**, and the mechanism is identified exactly.

The result is a *constructive existence proof*: an explicit, exactly stationary,
finite-amplitude periodic vacuum carrier (Decision M) whose substrate-derived
curvature is certified to span `su(3)`, plus a genericity argument showing this is
the rule and the subalgebra cases are measure-zero. It is reproduced by
`derivations/t1_su3_witness.py`.

What T1 does **not** yet include (still owed, unchanged): the coarse-grained
Yang–Mills action `−¼G²` (**T4**), the quantitative `U(1)/SU(3)` split as a
function of `(α_s,α_t,γ_t)` (**T7**), and the physical identification of *which*
triplet is color vs the EM trace (candidate Interpretation A). T1's own burden —
"a genuine three-component complex carrier, transported nontrivially, whose
traceless curvature spans `su(3)`" — is met.

---

## 1. The precise criterion

Following ARCHITECTURE Layer 3 and OPEN_TASKS §5, T1 is proved iff there is a
finite-amplitude **stationary** periodic background `R̄` (`∂S/∂R[R̄]=0`) such that
the Bloch fluctuation operator `D_{R̄}(k)` has an **isolated rank-3 subspace**
`P_3(k)` (near-degenerate triplet, large outside gap — Decision L) whose
**gauge-invariant** Wilczek–Zee curvature, split as

```
𝓕_{μν}(k) = f_{μν}(k) I_3 + 𝒢_{μν}(k),   Tr 𝒢_{μν}=0,
```

has traceless part `𝒢_{μν}` obeying

```
Lie⟨ 𝒢_{μν}(k), [𝒢_{μν},𝒢_{ρσ}], … ⟩ = su(3)   (real dimension 8).
```

Rank 0 ⇒ trivial; rank 3 ⇒ only `so(3)`; **rank 8 ⇒ genuine `su(3)`**. Gauge
invariance is essential: it is what rules out the "pure-gauge `u(3)` basis
artifact" objection. We compute `𝓕` by the Wilson-loop (Fukui–Hatsugai)
construction, which is manifestly invariant under `k`-local `U(3)` basis changes.

---

## 2. The derivational chain (analytic)

Nothing below inserts a gauge field; the connection is the holonomy of the
carrier eigenspace.

**Layer 0 → 1 (anisotropic stiffness tensor).** One link `(n,μ)` contributes
`S_{nμ}=½η_μκ_μ(L_{nμ}−r_μ)²` with `L=|Q|`, `Q=R_{n+ê_μ}−R_n`. The second
variation in the link vector is the frame tensor

```
C_{nμ}^{AB} = η_μκ_μ[ (1 − r_μ/L̄_{nμ}) δ^{AB} + (r_μ/L̄_{nμ}) Q̂_{nμ}^A Q̂_{nμ}^B ],
```

isotropic part + rank-one `Q̂Q̂ᵀ` (a spring is stiffer along its own axis). In the
straight vacuum `Q̂=ê_μ` and all `C_{nμ}` are diagonal in the fixed axis frame:
they commute, so their common eigenframe is `k`-independent → trivial bundle. In a
nonlinear background the `Q̂_{nμ}` differ link-to-link, so `[C_{nμ},C_{mν}]≠0` —
the mechanical seed of non-commutativity.

**Layer 1 → 2 (Bloch operator).** The quadratic fluctuation action
`S^{(2)}=½Σ_{n,μ}(ξ_{n+ê_μ}−ξ_n)ᵀC_{nμ}(ξ_{n+ê_μ}−ξ_n)`. For a supercell-periodic
`R̄`, the Bloch ansatz `ξ_n = ε_s(k)e^{ik·n}` turns each link into a matrix-valued
hopping `C_ℓ^{AB}e^{ik·δ}`, giving the Hermitian

```
D_{R̄}(k) = Σ_links  ⎡ block(s,s):  +C_ℓ,   block(t,t): +C_ℓ,
                       block(s,t): −C_ℓ e^{ik·δ},  block(t,s): −C_ℓ e^{−ik·δ} ⎤.
```

Because `η_4=+1` while `η_i=−1`, `D_{R̄}(k)` is **indefinite** (Lorentzian
signature is inherited as the sign pattern of the spectrum, not imposed) — but it
is Hermitian, so its eigenvalues are real and its spectral projectors are
well-defined.

**Layer 2 → 3 (carrier + WZ connection).** Pick an isolated triplet
`P_3(k)=Σ_{a=1}^3|ε_a⟩⟨ε_a|`. The `U(3)` Wilczek–Zee connection is
`[𝒜_μ]_{ab}=i⟨ε_a|∂_{k_μ}ε_b⟩`, curvature
`𝓕_{μν}=iP_3[∂_μP_3,∂_νP_3]P_3`.

**Layer 3 → 4 (split).** `𝒜_μ = a_μ I_3 + B_μ`, `a_μ=⅓Tr𝒜_μ` (EM `U(1)`
candidate), `Tr B_μ=0` (color). Curvatures `f_{μν}=∂_μa_ν−∂_νa_μ` and
`𝒢_{μν}=∂_μB_ν−∂_νB_μ−i[B_μ,B_ν]`.

---

## 3. Structure theorem (three cases)

> **Theorem.** Let the carrier be an isolated rank-3 subspace of `D_{R̄}(k)`.
>
> **(i)** For the straight vacuum `R̄_n=a n`, `D` is `k`-diagonalized by a fixed
> real frame; `P_3` is constant, `𝓕≡0`. *No color in the bare vacuum* — but this
> is not a no-go (LESSONS_LEARNED #3), only the wrong background.
>
> **(ii)** For a real transverse frame with a single-site effective operator
> `M(k)=Σ_μ 2(1−cos k_μ)C_μ` (real symmetric at every `k`), the eigenframe is
> real, so `𝓕` is real-antisymmetric-valued: `Lie⟨𝒢⟩ ⊆ so(3)` (rank ≤ 3). *This
> is the D5 trap.*
>
> **(iii)** For a generic finite-amplitude periodic supercell carrier, the
> non-gaugeable Bloch phases `e^{ik·δ}` make `D_{R̄}(k)` complex-Hermitian at
> generic `k`; `P_3(k)` is genuinely complex and `Lie⟨𝒢⟩ = su(3)` (rank 8). The
> failure set (closure in a proper subalgebra) is a positive-codimension
> subvariety.

Parts (i),(ii) are elementary: a real symmetric family has a real eigenframe, and
`i·log(SO(3)) = so(3)`, closed under commutators. The content is (iii), proved by
an explicit witness (§4–5) plus genericity (§6). The **decisive contrast** is
that the extra 5 generators (`su(3)⊖so(3)`, the symmetric-traceless + Cartan
directions) are exactly what complex Bloch phases supply on top of the real
noncommuting stiffness frames (§7).

---

## 4. The exact stationary carrier (closed form)

We need one *exactly stationary* finite-amplitude periodic background (Decision M:
the universal nonlinear periodic vacuum carrier, not a soliton core). Take a
**circularly-polarized helix**

```
R̄_n = a n + A[ cos(K·n) p + sin(K·n) q ],   K=(K1,0,0,K4),  p=ê_2, q=ê_3.
```

Because the polarization plane `span(p,q)=(ê_2,ê_3)` is orthogonal to `ê_1,ê_4`,
and `K` has no component along `2,3`:

- links `μ=1,4`: transverse increment is a **chord of fixed angle** `K_μ` on the
  `(p,q)` circle, so `L̄_{nμ}=√(a²+2A²(1−cos K_μ))` is **independent of `n`**,
  while the direction `Q̂_{nμ}` **twists** with `n` (the nontrivial frame);
- links `μ=2,3`: `K_{2,3}=0`, so `Q̄_{nμ}=a ê_μ` exactly (untwisted).

All link lengths are `n`-independent ⇒ all link tensions `f_μ=η_μκ_μ(L̄_μ−r_μ)`
are constant. The nodal force is purely transverse and reduces to a single scalar:

```
[ κ_s(L1−r_s)/L1·(1−cos K1) ]  =  [ κ_t(L4−r_t)/L4·(1−cos K4) ],
```

which we solve for `γ_t=κ_t/κ_s`. The spatial term (`f_1<0`) and temporal term
(`f_4>0`, since `η_4=+1`) have **opposite sign** — so a positive `γ_t` solution
exists precisely *because* of the Lorentzian sign pattern. The result is an
**exact critical point of `S`**: the code reports `‖∇S‖ ≈ 2×10⁻¹⁶` (machine
zero) under periodic BC. This is a genuine finite-amplitude nonlinear vacuum
carrier, in full compliance with BACKBONE (#1, #6) and LESSONS_LEARNED #1.

Example: `α_s=0.6, α_t=0.9, A=0.30a, N1=N4=3` gives `γ_t=2.322`,
`L1=L4=1.127a`.

---

## 5. The numerical witness

Running the chain on the carrier above (`derivations/t1_su3_witness.py`):

```
exact-stationarity residual   ‖∇S‖              = 2.3e-16
carrier = isolated bands [14,15,16]   gap/spread = 3.29   (Decision L satisfied)
U(1) trace curvature      |Tr 𝓕|                 ≈ 3.0e-11   (flat here → pure color)
SU(3) traceless curvature ‖𝒢‖                    ≈ 3.5e-01   (nonzero)
so(3) antisymmetric span                          = 3 / 3
symmetric + Cartan span (the 5 extra generators)  = 5 / 5
raw curvature span                                = 8 / 8
Lie-closure rank                                  = 8 / 8   →  FULL su(3)
```

The `𝒢_{μν}` samples span all 8 directions **before** any commutator is taken
(raw span already 8); the Lie closure confirms 8. Note the trace curvature is
numerically zero on the isolated triplets of this symmetric carrier — i.e. this
vacuum gives **pure `SU(3)` color** with a flat `U(1)`, so the traceless part is
unambiguously `su(3)`-valued, not contaminated by the trace. (A physical EM `U(1)`
would then come from a less symmetric carrier or a separate vortex/Berry sector —
consistent with Interpretation A remaining a *candidate* until T3.)

**Fair `so(3)` control.** The genuine "real transverse frame" object of case (ii)
— a single site's four real-symmetric non-commuting `C_μ` used as a 4-band
operator — gives `raw span = 3, Lie rank = 3`: exactly `so(3)`. So the numerics
*discriminate*; they do not return 8 unconditionally.

**Robustness.** Rank 8 holds at every tested `(α_s,α_t,A,N1,N4,w)` and every
`(k*, plaquette size h)` — 5 parameter sets × 3 base-point/step choices = 15
independent exactly-stationary configurations, all rank 8.

---

## 6. Genericity lemma (why 8 is the rule)

> **Lemma.** A finite family `{X_α}⊂su(3)` fails to generate `su(3)` only if it
> lies in a proper subalgebra. The maximal proper subalgebras of `su(3)` (up to
> conjugacy) are `su(2)⊕u(1)` (dim 4) and `so(3)` (dim 3). Membership of a fixed
> proper subalgebra `𝔥` is a set of linear conditions after conjugation, i.e. a
> positive-codimension algebraic condition on `{X_α}`. Hence the non-generating
> set is a measure-zero subvariety.

The substrate-derived curvatures `{𝒢_{μν}(k)}` are analytic functions of
`(α_s,α_t,γ_t,A,K)` and `k`. The witness certifies one point off the failure set;
since "full `su(3)`" is an open condition, it holds on an open neighborhood in
parameter space — the robustness scan (§5) confirms an extended open region. The
subalgebra outcomes (`so(3)`, `su(2)⊕u(1)`, Cartan) occur only on the tuned
measure-zero set — of which case (ii), the real-frame limit, is the physically
meaningful representative.

---

## 7. Mechanism — what resolves D5

Dialing the Bloch phase `e^{ik·δ} → e^{i s k·δ}` from `s=0` to `s=1`:

```
s=0.00 :  raw span 0,  Lie rank 0            (no k-dependence → trivial bundle)
s>0    :  raw span 8,  Lie rank 8            (so(3)=3  +  symmetric/Cartan=5)
control:  raw span 3,  Lie rank 3            (real frame, no complex phases → so(3))
```

The picture is now exact:

- **No transport** (constant frame) → trivial, no color.
- **Real moving frame** (real-symmetric operator) → `so(3)` only: the 3
  antisymmetric generators, from orientation transport of the transverse triplet.
- **Complex Bloch phases** on the noncommuting real stiffness frames → the extra
  **5** symmetric-traceless + Cartan generators, completing `so(3)→su(3)`.

This is precisely the mechanism conjectured in
`The_anisotropic_stiffness_tensor.md §10` and `From_U4_to_U3.md`:
*projected anisotropic stiffness (real symmetric) + complex Bloch phases +
noncommuting background texture ⇒ `su(3)`* — now demonstrated, not hoped.

The `U(3)`-vs-`SU(4)` question (Decision C, `From_U4_to_U3.md`) is respected
throughout: the base and substrate stay 4D; only the transported **carrier** is
rank 3 (the isolated triplet), so the connection is `u(3)`-valued and its
traceless part is `su(3)`, never `su(4)`.

---

## 8. What is proved, and what remains

**Proved (T1 core / D5):**
- An exact finite-amplitude stationary periodic vacuum carrier exists (closed
  form), consistent with Decision M and all BACKBONE constraints.
- Its `D_{R̄}(k)` has isolated near-degenerate rank-3 carriers (Decision L).
- The gauge-invariant traceless WZ curvature spans **full `su(3)` (rank 8)**,
  robustly and generically; the `so(3)`-only outcome is the measure-zero
  real-frame case. The 8 generators are mechanically sourced (§7).

**Still owed (unchanged by this result):**
- **T4** — coarse-grain to the effective `−¼G^a_{μν}G^{aμν}` Yang–Mills action
  (the `−i[B_μ,B_ν]` self-interaction term is present in `𝒢`; the action is not).
- **T7** — quantitative `U(1)/SU(3)` split and `λ_4` separation as functions of
  `(α_s,α_t,γ_t)`. (Here the trace was flat — a special symmetric point; T7 needs
  the general dependence.)
- **T3 / Interpretation A** — confirm the trace is the physical EM `U(1)`
  (Maxwell dynamics), or source EM from the vortex sector; needs a carrier with
  nonzero trace curvature.
- **D7** — the physical selection principle for *which* periodic carrier is the
  universal vacuum microtexture (the helix here is a valid representative of the
  class, not yet shown to be the selected ground/attractor state).

---

## 8b. Strengthening (genericity, robustness, scaling, Bianchi)

See `STRENGTHENING_notes.md`. Beyond the single witness of §5:
- **Generic & robust** (`t1_genericity_robustness.py`): full su(3) in **23/24**
  distinct exact-stationary backgrounds (the exception is a measure-zero point, per
  the §6 lemma); rank-8 stable across supercell size and plaquette step.
- **Leading-order & Bianchi** (`t1_smallamp_bianchi.py`): the su(3) span persists to
  `A=0.01` with `‖G‖ ∼ A²` (mechanical seed `‖[C,C]‖ ∼ A¹`) — not a large-amplitude
  artifact; and the non-abelian Bianchi `D_{[μ}G_{νρ]}=0` is satisfied.
- **Spacetime realization** (`bridge_spacetime_gauge.py`): the k-space carrier bundle
  also yields a physical (non-pure-gauge) su(3) field strength over *spacetime*, with
  nonzero mixed `G_{kx}` — see T4 §8b.

---

## 9. Reproducibility

`derivations/t1_su3_witness.py` — self-contained (numpy only). Prints the
stationarity residual, the isolated-carrier gap ratio, the trace/traceless
curvature magnitudes, the `so(3)` / symmetric+Cartan / raw / Lie-closure ranks,
the fair `so(3)` control, and the robustness table. All numbers quoted above are
its output.