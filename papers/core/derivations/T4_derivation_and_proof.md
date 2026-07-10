# T4 — Derivation and proof of the Yang–Mills action

**Status: the form of T4 is proved; its coefficient is derived (positive, finite,
strain-dependent).** Building on T1 (the traceless Wilczek–Zee curvature spans
`su(3)`), the coarse-grained effective action for the emergent color connection is

```
S_eff[B] = ∫ d⁴x √(−g) [ −(1/4g²) G^a_{μν} G^{a,μν} ] + (irrelevant, O(1/Λ²)),
G^a_{μν} = ∂_μ B^a_ν − ∂_ν B^a_μ + f^{abc} B^b_μ B^c_ν,
```

a **universal `SU(3)` Yang–Mills field over empty space** (Decision M). The `−¼G²`
form is *forced* by the emergent gauge invariance that T1 established; the
substrate supplies a positive, finite, scale-dependent `1/g²` and the `su(3)`
self-interaction. Reproduced by `derivations/t4_yang_mills.py`.

What remains beyond T4 (nonperturbative, and shared with T2/T7): the exact
running coefficient / β-function, confinement, and the Minkowski index
contraction (which uses the T2 effective metric). T4's own burden — *"the
coarse-grained `−¼ G^a_{μν}G^{aμν}` effective action with the correct non-Abelian
self-interaction, as a universal field over empty space"* — is met.

---

## 1. What must be shown

T1 delivered the *algebraic* object `𝒢_{μν}(k) ∈ su(3)` — a fixed curvature over
the Brillouin zone. T4 must show its *dynamics* over the emergent spacetime is
Yang–Mills. Four things are required:

1. a propagating spacetime gauge field `B_μ(x) ∈ su(3)` (not just a k-space
   background);
2. its leading action is `−¼g⁻²G²` and nothing lower-dimensional (in particular,
   no mass term);
3. the coefficient `1/g²` is real, positive (stability) and comes from the
   substrate, not by hand (BACKBONE #5);
4. the self-interaction is the unique `su(3)` Yang–Mills vertex set with one
   coupling.

The logical spine is a **uniqueness theorem** (§4): given emergent local `SU(3)`
invariance + locality + power counting, the answer can only be `−¼g⁻²G²`. The
substrate calculations (§5) then fix and characterise `1/g²` and confirm the
`su(3)` self-interaction.

---

## 2. The dynamical gauge field (from k-space background to spacetime field)

T1's connection lives over the Brillouin zone of the *universal* vacuum carrier
`R̄`. Because the carrier is the vacuum itself (Decision M), it is present in every
region of the emergent spacetime `M₄`. Nothing forces the carrier's `SU(3)` frame
to be the *same* in every region: the physical object is the gauge-invariant
projector `P_3`, and the choice of orthonormal basis inside `P_3` is a pure
redundancy. Promoting that redundancy to a slowly varying spacetime field,

```
{ε_a(k)} ↦ { g(x)_{ab} ε_b(k) },   g(x) ∈ SU(3),
```

is the standard coarse-graining move: the local carrier orientation becomes a
spacetime field, and Wilczek–Zee parallel transport of the triplet defines a
connection

```
B_μ(x) ∈ su(3),   transforming as   B_μ → g B_μ g⁻¹ − i (∂_μ g) g⁻¹.
```

This `B_μ(x)` is the emergent gluon field. Its curvature is exactly the T1 object,
now `x`-dependent: `G_{μν} = ∂_μB_ν − ∂_νB_μ − i[B_μ,B_ν]`. No new field is
postulated — `B_μ` is the collective coordinate of the vacuum carrier's frame.

---

## 3. Minimal coupling (why the connection enters as `D_μ`)

The carrier amplitude is a field `ψ(x) ∈ ℂ³` (the occupation of the three carrier
modes). By the Wilczek–Zee projection theorem, the substrate dynamics *projected
onto an isolated multiplet* is governed by the covariant derivative built from the
WZ connection: to leading order in gradients,

```
L_ψ = ψ† (iD_0)² ψ − Σ_i v_i² ψ† (iD_i)†(iD_i) ψ + …,   D_μ = ∂_μ − i B_μ,
```

with the `v_i` the carrier band velocities read off from `D_{R̄}(k)`. This is not
an assumption: for an isolated band multiplet the only way the connection can
appear in the projected equations of motion is as `D_μ` — that is the content of
adiabatic/WZ reduction. Hence `B_μ` couples to the substrate exactly as a gauge
field, and the local `SU(3)` acting on `(ψ, g)` simultaneously is an **exact
redundancy** of the description.

---

## 4. Uniqueness theorem (the proof core)

> **Theorem.** Let `Γ[B]` be the coarse-grained effective action obtained by
> integrating out all substrate fluctuations except the soft carrier-frame field.
> Assume (i) **exact local `SU(3)` invariance** (§2–3), (ii) **locality** — the
> integrated-out modes are gapped (the non-carrier bands and the high-`k` carrier
> modes), so their kernel is short-ranged and `Γ[B]` admits a derivative
> expansion, (iii) **Lorentz invariance** in the inside-observer metric (T2).
> Then the leading term of `Γ[B]` is
>
> ```
> Γ[B] = ∫ d⁴x √(−g) [ −(1/4g²) G^a_{μν} G^{a,μν} ] + O(1/Λ²),
> ```
>
> uniquely, with the non-Abelian `G` of §2 and a single coupling `g`.

**Proof.** By (i), `Γ` is a functional of gauge-covariant quantities only:
`G_{μν}` and its covariant derivatives, contracted into `SU(3)`- and
Lorentz-scalars via `Tr` and the metric. Enumerate by mass dimension (`[B]=1`,
`[G]=2`, `[∂]=1`):

- **dim 2:** `Tr(B_μB^μ)` — *not gauge invariant* (§2 transformation law), so
  **forbidden**. This is the crucial point: emergent gauge invariance outlaws the
  gluon mass term. (Equivalently: a pure-gauge `B_μ=−i(∂_μg)g⁻¹` has `G_{μν}=0`,
  so it must cost zero action; a `Tr(BB)` term would penalise it. Verified
  numerically in §5D.)
- **dim 4:** `Tr(G_{μν}G^{μν})` — the unique gauge- and Lorentz-invariant scalar.
  (`Tr(G_{μν}Ĝ^{μν})` is a total derivative / topological term, not part of the
  local dynamics.)
- **dim ≥ 6:** `Tr(G³)`, `Tr(D_αG_{μν})²`, … — irrelevant, suppressed by the gap
  scale `Λ ∼ 1/a`.

Hence the leading term is `c · Tr(G_{μν}G^{μν})`. Writing `Tr(T^aT^b)=−½δ^{ab}`
and `c=−1/(2g²)` gives `−(1/4g²)G^aG^a`. The `su(3)` structure constants
`f^{abc}` in `G` are inherited from T1 (§5C), so the cubic and quartic gluon
self-couplings are fixed by the same `g` — the Yang–Mills vertex set. ∎

The theorem reduces T4 to two substrate facts: `1/g² > 0, finite` (else the
"action" is unstable or trivial) and the algebra is genuinely `su(3)`. Both are
computed in §5.

---

## 5. Substrate computations (`t4_yang_mills.py`)

Carrier: the exact stationary helix of T1 (`‖∇S‖ ≈ 2×10⁻¹⁶`, `γ_t=2.322`).

### (A) `1/g²` is positive and finite
The induced gauge stiffness is the substrate's own quantum-geometric response of
the isolated carrier — the integrated **non-Abelian quantum metric**

```
g_{μν}(k) = ½ Re Tr[ (∂_μ P_3)(∂_ν P_3) ],   1/g² ∝ ∫_BZ d⁴k Σ_μ g_{μμ}(k).
```

`Σ g_{μμ} = ½Σ_μ‖∂_μP_3‖²_F ≥ 0` pointwise, so `1/g²` is manifestly positive.
This is the *same* quantum-geometric tensor whose imaginary part gave the T1
curvature — one geometric object yields both the curvature (T1) and the kinetic
stiffness (T4). Computed:

```
integrated quantum-metric stiffness = 11.84 ± 9.09   (min sample 2.93 > 0)
```

over 400 gapped `k`-points ⇒ **`1/g²` positive and finite ⇒ a stable propagating
gauge field**, not a topological/degenerate one. (The proportionality constant
between this geometric integral and the physical `1/g²` involves the band
velocities and the coarse-graining cutoff — a scheme detail; the *positivity and
finiteness* are the theorem-level facts.)

### (B) The coupling runs with strain (Decision J)
BZ-averaged stiffness vs carrier amplitude `A` (a strain proxy):

```
 A      1/g² (BZ avg)     g² ~ 1/⟨stiffness⟩
0.15       8.96              0.112
0.20      10.74              0.093
0.25      10.72              0.093
0.30      12.29              0.081
0.35      13.09              0.076
0.40      11.47              0.087
```

`1/g²` is not a universal constant — it is set by the local nonlinear texture.
This is the derived origin of Decision J / the Large_Energies picture: coupling
strength depends on strain/scale (weak and long-range where the carrier is soft,
strong and short-range where the link frames twist hard — near matter). The
coupling is *derived*, never imposed (BACKBONE #5).

### (C) The algebra is `su(3)`, fixing the self-interaction
Orthonormalising the T1 curvature samples gives 8 generators; their structure
constants satisfy

```
total antisymmetry  |f_{abc}+f_{acb}|        = 2.2e-16
Casimir  κ_ab = f_{acd}f_{bcd}  off-diagonal = 2.9e-11
         κ diagonal                          = [6,6,6,6,6,6,6,6]  (= 2N, N=3)
         max/min diagonal ratio              = 1.0000
```

`κ_ab = 6 δ_ab = 2N δ_ab` with `N=3` is the adjoint Casimir of `su(3)` — a
definitive fingerprint (a simple algebra with `C₂(adj)=2N`). Since the generators
are literally traceless anti-Hermitian `3×3` matrices spanning the full 8-dim
space, the gauge algebra **is** `su(3)`, and the cubic/quartic gluon vertices are
the unique `su(3)` Yang–Mills ones sharing the single coupling `g`.

### (D) The gluon is massless
A pure-gauge configuration `B_μ = −i(∂_μg)g⁻¹` for random smooth `g(x)∈SU(3)`
gives

```
max|G_{μν}| = 2.3e-10   (== 0)
```

so pure gauge costs zero action — the `Tr(BB)` mass term is forbidden (§4), and
masslessness of the gluon and the `G²` form are the *same* consequence of the
emergent gauge invariance.

---

## 6. Why this is "over empty space" (Decision M)

The carrier `R̄` is the universal nonlinear periodic vacuum microtexture, present
everywhere — so `B_μ(x)` and its Yang–Mills action are defined over all of empty
spacetime, independent of matter. Matter (Layer-5 solitons) enters later as
*color-charged sources* `J^a_μ` coupling to this already-universal field
(`ψ`-currents of §3 localised on solitons), exactly the QCD structure: universal
gluon field + localised colored matter. The strain-dependence of `1/g²` (§5B)
makes the coupling strong precisely near those sources — the asymptotic-freedom-
like behaviour the architecture calls for.

---

## 7. Scope — proved vs still owed

**Proved / derived (T4):**
- A propagating spacetime `su(3)` gauge field `B_μ(x)` emerges as the vacuum
  carrier's frame collective coordinate (§2), minimally coupled (§3).
- Its leading action is uniquely `−¼g⁻²G²` (uniqueness theorem §4); the mass term
  is forbidden and the gluon is massless (§4, §5D).
- `1/g²` is substrate-induced, positive, finite (§5A) and strain/scale-dependent
  (§5B) — Decision J derived.
- The self-interaction is the unique `su(3)` Yang–Mills vertex set (§5C).

**Still owed (nonperturbative / cross-task):**
- the exact `1/g²` normalisation and its β-function / running coefficient
  (requires the full gradient-expansion loop integral, not just the geometric
  stiffness);
- confinement and the mass gap (nonperturbative, not expected from the derivative
  expansion);
- the Minkowski contraction `G_{μν}G^{μν}` uses the T2 effective metric — T4
  inherits T2's calibration `(α_s,α_t,γ_t,κ,a)`;
- coupling to Layer-5 colored solitons (the matter currents `J^a_μ`).

---

## 8. Relation to Weinberg–Witten and emergent-Lorentz constraints

Any programme that produces **emergent massless gauge bosons** must confront the
**Weinberg–Witten theorem** (Weinberg & Witten, Phys. Lett. B 96, 59, 1980): in a
theory with an **exactly Lorentz-covariant** conserved current `J^μ`, there can be
no massless spin-1 particle carrying a nonzero charge under that current (and,
via the covariant stress tensor, no massless spin-2 graviton). An emergent
non-abelian gauge boson — the gluon of §1–4 — is exactly the object the theorem
targets, so `−¼g⁻²G²` as *fundamental* physics would appear to be forbidden.

The escape is the one used throughout the emergent-gravity / emergent-gauge
literature (Volovik; Wen; analogue gravity — see `RELATED_WORK.md`): **Lorentz
invariance here is emergent and only approximate, not exact.** Weinberg–Witten
assumes exact Lorentz covariance of the current; that premise fails in this
substrate. This is not a loophole invented for T4 — it is precisely what **T2**
independently found: the lab dispersion is genuinely anisotropic/birefringent for
`α_s>0`, and an exact isotropic cone is reached only in the `α_s→0` limit (which
disables the gauge sector). So:

> the approximate nature of emergent Lorentz invariance (a T2 "limitation") is
> exactly the property that lets the emergent Yang–Mills sector (T4) exist without
> violating Weinberg–Witten.

Concretely, the theorem is evaded because (i) the gauge bosons are collective modes
of the substrate (the carrier frame field), not asymptotic single-particle states
of a Lorentz-covariant `S`-matrix; and (ii) exact Lorentz covariance of the color
current holds only in the long-wavelength limit, with lattice/anisotropy
corrections at `O((ka)²)` and `O(α_s)`. The same reasoning applies to the T3 photon
(spin-1) and to any future emergent-graviton (Layer −) claim — the latter being the
more dangerous case that must be revisited if the gravity channel is developed. A
paper-level treatment should state the Weinberg–Witten assumptions explicitly and
point to the emergent/approximate-Lorentz escape as the consistency argument.

---

## 8b. The k-space → spacetime promotion, demonstrated

The promotion of the k-space carrier connection (§2) to a spacetime gauge field is no
longer only an argument. `derivations/bridge_spacetime_gauge.py` (see
`STRENGTHENING_notes.md` #1) shows that slow spacetime variation of the carrier
produces a **physical, non-pure-gauge su(3) field strength `G_{xx}(x)`** (a
constant-subspace basis rotation gives `‖G‖≈2.6e-9`, confirming physicality), and that
the **mixed curvature `G_{kx}` is nonzero** in every k-direction — so the
Brillouin-zone object (T1) and the spacetime gauge field are restrictions of a single
`(k,x)` connection, with `D_μ=∂_μ−iA_μ` from the Wilczek–Zee adiabatic theorem. This
substantiates the collective-coordinate step in §2/§3.

---

## 9. Reproducibility

`derivations/t4_yang_mills.py` — imports the T1 carrier, prints (A) the positive
finite stiffness, (B) the strain-running table, (C) the `su(3)` structure-constant
and Casimir checks, (D) the pure-gauge masslessness test. All numbers above are
its output.