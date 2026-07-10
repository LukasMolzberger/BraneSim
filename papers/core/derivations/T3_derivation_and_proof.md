# T3 — Derivation and proof of the Faraday U(1) (Maxwell) sector

**Status: T3 is proved.** The electromagnetic sector is the *abelian* holonomy of
the same substrate carrier that gave T1/T4. Its effective dynamics is Maxwell,

```
S_EM = ∫ d⁴x √(−g) [ −(1/4e²) f_{μν} f^{μν} + a_μ J^μ ],
f_{μν} = ∂_μ a_ν − ∂_ν a_μ,   a_μ = ⅓ Tr 𝒜_μ  (or any non-degenerate Berry line bundle),
```

with: the Bianchi identity (homogeneous Maxwell) automatic; the inhomogeneous
Maxwell equation `∂^μ f_{μν} = e² J_ν`; a **massless photon**; a substrate-induced
**positive, finite** `1/e²`; and **charge = U(1) vortex winding** in `π₁(U(1))=ℤ`.
A clean, physically decisive feature emerges: the pristine vacuum is **EM-flat**
(no background electromagnetic field), and EM curvature switches on only where
matter/deformation breaks the vacuum's reality symmetry. Reproduced by
`derivations/t3_maxwell.py`.

What remains beyond T3 (shared with T2/L5): the exact `1/e²` normalisation, the
Minkowski contraction (T2 metric), and the explicit soliton current `J^μ` (Layer
5). T3's own burden — *"Maxwell dynamics + sourced equation + photon propagation;
charge = U(1) vortex winding"* — is met.

---

## 1. What must be shown

OPEN_TASKS: *"Algebraic `f_{μν}` and Bianchi are automatic; still owe Maxwell
dynamics + sourced equation + photon propagation. Charge = U(1) vortex winding."*
So T3 must deliver:

1. the algebraic `f_{μν}` and its Bianchi identity (homogeneous Maxwell);
2. Maxwell dynamics `−¼e⁻²f²` with a substrate-derived, positive `1/e²`;
3. the sourced (inhomogeneous) equation and a propagating, **massless** photon;
4. charge quantisation as U(1) vortex winding, and current conservation.

The logic mirrors T4 but is **abelian**, hence more robust: U(1) needs no
degeneracy, so it does not even rely on the hard `su(3)` result — only on the
generic emergence of a Berry line bundle plus U(1) gauge invariance.

---

## 2. The abelian connection and its index structure

Two equivalent sources give the same U(1):

- **Trace of the color `U(3)`** (Interpretation A): `a_μ = ⅓ Tr 𝒜_μ`, the common
  phase of the carrier triplet, embedded internally as `f_{μν} I_3` — a 4D
  antisymmetric tensor with *trivial* internal action (it moves the common phase,
  never mixes the triplet). Index structure (settled in `Trace_vs_Traceless.md`):
  `𝒜_μ{}^a{}_b = a_μ δ^a{}_b + B_μ{}^a{}_b`, `μ=1..4` base, `a,b=1..3` fiber.
- **Any non-degenerate band** of `D_{R̄}(k)`: an ordinary Berry connection
  `a_μ(k) = i⟨u|∂_{k_μ}u⟩`.

Promoted to spacetime exactly as in T4 (§T4.2): the local U(1) phase becomes a
field `a_μ(x)`, transforming as `a_μ → a_μ + ∂_μ λ`. Its curvature is the Faraday
tensor `f_{μν}` — a 4D antisymmetric tensor, 6 independent components (`E`, `B`).

---

## 3. Bianchi identity — homogeneous Maxwell (automatic)

Because `f = da`, `df ≡ 0` identically:

```
∂_λ f_{μν} + ∂_μ f_{νλ} + ∂_ν f_{λμ} = 0.
```

Verified numerically on a random smooth potential: `max|∂_[λ f_{μν]}| = 1.1e-13`.
These are Faraday's law and the absence of magnetic charge — half of Maxwell,
free. (A magnetic source would require `a_μ` to be globally ill-defined, i.e. a
Dirac string / monopole — not present for a single-valued Berry phase.)

---

## 4. Maxwell dynamics — abelian uniqueness theorem

> **Theorem.** The coarse-grained effective action `Γ[a]` for the soft U(1) mode,
> assuming (i) exact U(1) gauge invariance, (ii) locality (gapped modes → short-
> range kernel), (iii) Lorentz invariance (T2 metric), has leading term
>
> ```
> Γ[a] = ∫ d⁴x √(−g) [ −(1/4e²) f_{μν} f^{μν} ] + O(1/Λ²),
> ```
>
> uniquely. **Proof.** By (i), `Γ` depends only on `f_{μν}` and derivatives. By
> dimension: `a_μ a^μ` (dim 2) is **not gauge invariant → forbidden** (⇒ massless
> photon); `f_{μν}f^{μν}` (dim 4) is the unique invariant (`f f̃` is a total
> derivative); higher terms are irrelevant. ∎

This is the abelian shadow of the T4 theorem. Adding the matter coupling `a_μ J^μ`
and varying gives the **inhomogeneous Maxwell equation**

```
∂^μ f_{μν} = e² J_ν,
```

while gauge invariance of `a_μ J^μ` forces **current conservation** `∂^ν J_ν = 0`
(Noether). Together with §3 (Bianchi) this is the full set of Maxwell equations.

---

## 5. Substrate coupling `1/e²` — positive and finite

`1/e²` is the substrate's **abelian quantum metric** of a U(1) band,

```
g_{μν}(k) = ½ Re Tr[(∂_μ P_1)(∂_ν P_1)],   1/e² ∝ ∫_BZ d⁴k Σ_μ g_{μμ}(k) ≥ 0.
```

Computed over isolated-band `k`-points:

```
integrated U(1) stiffness 1/e² = 2.62 ± 7.83   (min sample 1.34 > 0)
```

⇒ **positive, finite `1/e²` ⇒ a stable, propagating photon**. Crucially, `1/e²`
is the *symmetric* (metric) part of the quantum geometric tensor; it is nonzero
**even though the antisymmetric (curvature) part vanishes in vacuum** (§6). So the
photon has a finite kinetic stiffness and propagates, while the vacuum carries no
background field. (The proportionality constant to the physical `1/e²` is
coarse-graining-scheme dependent, as in T4; positivity/finiteness are the
theorem-level facts.)

---

## 6. The vacuum is EM-flat; EM is matter-sourced (the key physical result)

A surprising, clean feature — first flagged by the flat color-trace in T1 —
turned out to be a theorem-quality statement:

```
|abelian Berry curvature|  pristine vacuum   = 4.7e-10   (~0, EM-flat)
|abelian Berry curvature|  symmetry-broken   = 3.4e+01   (turns ON)
color trace |Tr F_color|                     = 3.0e-11   (~0)
color traceless ||G_color||                  = 3.5e-01   (nonzero)
```

**All abelian (U(1)) Berry curvatures vanish pointwise in the symmetric vacuum**
(≈1e-10 by the accurate Wilson-loop measure), enforced by a PT-like reality
symmetry of the helical carrier. Physically: *empty space carries no background
electromagnetic field* — exactly as it must. This is fully consistent with the
nonzero `su(3)` of T1, because color lives in the **off-diagonal, non-abelian**
sector, untouched by the vanishing of the abelian diagonal/trace.

When the reality symmetry is broken — a generic non-symmetric texture, i.e. a
deformation or matter region — the abelian curvature switches on (`~34`). So EM is
**entirely matter-sourced**: the vacuum is neutral and flat, and electromagnetic
fields appear where matter deforms the carrier. This is the Faraday counterpart of
Decision J and dovetails with §7 (charge = winding) and BACKBONE #13 (photon
quantisation inherited from matter, not intrinsic).

---

## 7. Charge = U(1) vortex winding

Electric charge is topological: `π₁(U(1)) = ℤ`. A carrier phase field with winding
`n` around a loop has

```
∮ ∂_μ θ dx^μ = 2π n     (verified: n = 0,±1,2,3 → flux = 2πn exactly),
```

so by Gauss's law `∮ f = Q`, the enclosed charge is quantised in integer units.
Because the winding is a topological invariant, the associated current is
conserved (`∂^μ J_μ = 0`) — charge cannot leak. This is the substrate origin of
`charge = U(1) vortex winding` (BACKBONE, Layer 5): a Layer-5 soliton with phase
winding `n` is a source of `n` quanta of electric charge for the field of §4.

**The photon itself is continuous and non-quantised** (BACKBONE #13): `f_{μν}` is a
smooth classical field with continuous amplitude; the wave equation `∂^μf_{μν}=0`
(source-free §4) gives propagating photons with no intrinsic quantum. All apparent
`hν` granularity enters only through emission/absorption gated by the discrete
Layer-5 soliton transitions — i.e. through the *charged matter*, never through
`f_{μν}`.

---

## 8. Scope — proved vs still owed

**Proved / derived (T3):**
- `f_{μν}` and Bianchi/homogeneous Maxwell — automatic (§3).
- Maxwell dynamics `−¼e⁻²f²` unique (§4); inhomogeneous `∂^μf_{μν}=e²J_ν` and
  `∂·J=0`.
- `1/e²` substrate-induced, positive, finite (§5); photon massless (§4) and
  propagating; continuous / non-quantised (§7, BACKBONE #13).
- Vacuum EM-flat, EM matter-sourced (§6).
- Charge = quantised U(1) vortex winding, conserved (§7).

**Still owed (cross-task / matter):**
- exact `1/e²` normalisation and its scale dependence (as with T4's `1/g²`);
- the Minkowski contraction `f_{μν}f^{μν}` uses the **T2** effective metric;
- the explicit soliton current `J^μ` and the electron/proton charge spectrum
  (**Layer 5**); the `U(1)` mixing with the color trace in matter regions (a
  quantitative **T7** question — here vacuum-separated by symmetry).

---

## 9. Reproducibility

`derivations/t3_maxwell.py` — imports the T1 carrier and prints: (A) the Bianchi
residual; (C) the EM-flat vacuum vs symmetry-broken abelian curvature and the
color comparison; (B) the positive finite `1/e²`; (D) the pure-gauge masslessness
check; (E) the winding-quantised charge. All numbers above are its output.