# Bell Bridge — Deriving the singlet lock J2 from the substrate vertex nonlinearity (D1)

Working note. Reduces J2 (the source emits the polarization singlet) from an
*imported* conservation law to a *consequence of the vertex's polarization
isotropy*, with a definite validity regime, a controlled breaking parameter, and
a falsifiable fidelity prediction. Honest debt; the residual tensor computation
is flagged. Closes part of D1 of `vvertex_junction_condition.md`. Group theory
below numerically verified.

## What this closes

J2 was stated as "polarization angular-momentum conservation ⇒ singlet" — an
import. This note shows the singlet is the **unique** output of *any*
SU(2)-polarization-invariant pair-creation vertex, then grounds that invariance
(and its breaking) in the substrate. It does not yet compute the vertex tensor
from first principles (residual, §6).

## 1. The substrate vertex coupling

`S` is a nonlinear region; the carrier's geometric nonlinearity (`∝ k_s α`,
Paper I `eq:W4` and its cubic cousin) provides a three-wave coupling
`pump → A + B`. The pair-creation operator on the transverse polarization
doublet (`α,β ∈ {H,V}`, `ℂ²`) is

    C† = Σ_{αβ} M_{αβ} a†_{A,α} a†_{B,β},

with `M` the vertex polarization tensor (set by the nonlinearity's tensor
structure and phase matching). The produced pair is `|Ψ_pair⟩ ∝ M_{αβ} |α⟩_A |β⟩_B`.

## 2. The group-theory lock (rigorous)

Split `M = M_sym + M_anti`. For two doublets, `2 ⊗ 2 = 3 ⊕ 1` (symmetric triplet
⊕ antisymmetric singlet). Under the diagonal polarization rotation
`U_A = U_B = U ∈ SU(2)`:

- `M_anti ∝ ε_{αβ}` is the **unique invariant** (`U ε Uᵀ = det(U) ε = ε`); it
  creates the singlet `|Ψ⁻⟩`.
- `M_sym` transforms as the triplet (a vector) — **not** invariant; it creates
  triplet Bell states, each with a preferred polarization axis.

Verified numerically: `|Ψ⁻⟩` is `U⊗U`-invariant (overlap 1 for random SU(2));
`|Φ⁺⟩`, `|Ψ⁺⟩` are not; the `ε`-contraction operator yields `|Ψ⁻⟩` at fidelity 1.

**So: if the vertex coupling is invariant under diagonal SU(2) polarization
rotations, then `M = M_anti ∝ ε` uniquely and the emitted pair is the singlet.**
This is J2 — derived modulo the symmetry (§3).

## 3. Where the substrate HAS the symmetry — and where it breaks

The diagonal SU(2) acts on the transverse polarization doublet; the substrate
realizes it exactly when the two transverse polarizations are **degenerate**.

- **On-axis, long wavelength.** For propagation along a cubic axis (`[100]`
  etc.) the two transverse directions are exchanged by the 4-fold point symmetry
  ⇒ the doublet is degenerate ⇒ its 2D subspace carries a full SU(2) at the
  linear level. The geometric nonlinearity is isotropic in the transverse
  displacement (built from rotational invariants of `ΔR`), so the vertex coupling
  is naturally SU(2)-invariant there ⇒ `M = M_anti ∝ ε` ⇒ singlet. ✓
- **Breaking.** Off-axis propagation and finite `k·a` (lattice anisotropy /
  birefringence) split the transverse doublet, reducing SU(2) → SO(2) (about the
  axis) or lower. `M` then acquires a symmetric (triplet) part `∝ (anisotropy)`.

**Why SO(2) is not enough (sharp point).** Helicity (about-axis) conservation
alone leaves the `J_z = 0` sector 2-dimensional: *both* `ε` (singlet) and `δ`
(the `|Φ⁺⟩` triplet member) are SO(2)-invariant (`R ε Rᵀ = ε`, `R δ Rᵀ = δ` for
real `R`). Only the **full** SU(2) (which includes the relative-phase /
birefringent generators) kills `δ`: `U δ Uᵀ ≠ δ` for complex `U`. So the singlet
lock is **exactly as good as the transverse degeneracy** — degeneracy is the
physical content of J2.

## 4. Falsifiable consequence: singlet fidelity vs anisotropy

With `Δ` the transverse-doublet splitting (set by off-axis angle and `k·a`) and
`Γ` the vertex bandwidth, the symmetric admixture gives singlet fidelity

    F = 1 − O((Δ/Γ)²).

Prediction: Bell-state purity degrades with off-axis angle and toward short
wavelength, and is exact only on-axis in the long-wavelength limit. This is the
substrate's version of the birefringence/walk-off entanglement degradation seen
in real SPDC — a concrete, falsifiable handle on J2.

## 5. Bose consistency

`ε` is antisymmetric in the `(A,B)` polarization indices. Photons are bosons, so
the total two-photon state must be symmetric; the daughters occupy distinct
directional modes (`A → Alice`, `B → Bob`), and the antisymmetric polarization
pairs with the antisymmetric two-mode assignment to keep the full state
symmetric. Consistent with a two-arm (type-II-like) emission geometry — and it is
why the antisymmetric **singlet**, not a symmetric triplet, is the natural
two-mode output.

## 6. Residual debt (still open under D1)

- Compute `M` explicitly from the substrate geometric-nonlinearity tensor for a
  concrete on-axis vertex geometry: confirm `M_sym = 0` on-axis and obtain the
  off-axis admixture coefficient (the `O((Δ/Γ)²)` of §4).
- Confirm the cubic three-wave term (`pump → A + B`) is actually present in the
  carrier's geometric nonlinearity with the assumed structure (the nonlinearity's
  leading anharmonic term is quartic `∝ k_s α/a`; identify the cubic vertex or the
  effective three-wave coupling it induces in a pumped background).
- Relate `Δ` quantitatively to off-axis angle and `k·a` from the Paper I
  transverse-branch dispersion.

## Result

J2 is the unique SU(2)-polarization-invariant vertex output; the invariance is
exact for on-axis long-wavelength propagation (cubic-symmetry degeneracy of the
transverse doublet) and broken controllably by anisotropy. This moves J2 from
"imported conservation law" to "consequence of substrate transverse-doublet
degeneracy," with a falsifiable singlet-fidelity prediction — modulo the explicit
vertex-tensor computation (§6).

## References (wire into `references.bib` if this graduates)

- `2 ⊗ 2 = 3 ⊕ 1`; `ε` the unique SU(2)-invariant rank-2 tensor (any group-theory
  / QM text).
- Paper I — geometric nonlinearity `∝ k_s α`, transverse-branch degeneracy and
  dispersion.
- Standard SPDC: birefringence / walk-off degrades polarization entanglement (the
  empirical analogue of §4).
