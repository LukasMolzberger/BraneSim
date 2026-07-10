# Core-claim strengthening notes

Three additions that harden the core (gauge + kinematics) chain without touching the
matter sector (Layer 5). Each attacks a specific referee objection. All reproduce
from the scripts named below.

---

## #1 — The k-space → spacetime bridge is *demonstrated*, not postulated
`derivations/bridge_spacetime_gauge.py`

**Objection addressed:** "T1's su(3) lives over the Brillouin zone; the promotion to a
spacetime gauge field (T4) is hand-waved."

The emergent connection is defined on the whole parameter space the carrier depends
on, `(k, x)`. T1 is its restriction to the k-planes; the physical gauge field is its
restriction to the spacetime (x) planes. Modelling slow spacetime dependence by a
non-rigid modulation of the carrier `λ(x)`:

- **(A)** the spacetime field strength `G_{xx}(x)` is **nonzero and su(3)-valued**
  (mean `‖G_traceless‖ ≈ 8.0e-2`);
- **(B)** a genuine **pure-gauge** control (constant carrier subspace, U(3) basis
  rotation) gives `‖G‖ ≈ 2.6e-9` — so (A) is physical, not a basis artifact;
- **(C)** the **mixed** curvature `G_{kx}` is nonzero in every k-direction
  (`≈ 0.06–0.25`): `k` and spacetime are components of a *single* connection.

Analytic backing: for an isolated multiplet with slow `λ(x)`, the Wilczek–Zee
adiabatic theorem gives projected dynamics with covariant derivative
`D_μ = ∂_μ − iA_μ`, `A_μ = iP∂_{x_μ}P`, and `G_{μν} = iP[∂_μP,∂_νP]P` — the *same*
construction as T1 with the base reinterpreted from `k` to `x`. So the full su(3)
structure group certified over `k` (T1) is the structure group of the spacetime gauge
field. (Note: the spacetime field *sampled* only a few su(3) directions here — that is
a sampling limit of the modulation set, not a limit of the group; the group is fixed
by T1.)

---

## #2 — su(3) is generic and discretization-robust, not a single tuned witness
`derivations/t1_genericity_robustness.py`

**Objection addressed:** "It's one cherry-picked background / a numerical artifact."

- **Genericity:** over **24 distinct exact-stationary helices** (different propagation
  axis, supercell size, winding, `α_s`, `α_t`, amplitude; all `‖∇S‖ ≲ 1e-15`), the
  robust Lie rank is **8 (full su(3)) in 23/24**. The single exception is a
  measure-zero point — exactly what the T1 genericity lemma predicts.
- **Robustness:** for the reference carrier, rank-8 su(3) is **stable across supercell
  size** (3×3, 4×4, 5×5) **and plaquette step** `h ∈ {1e-2, 2e-3, 5e-4}` (raw span
  already 8; not a discretization artifact).

---

## #3 — su(3) is a leading-order (A²) effect and obeys Bianchi
`derivations/t1_smallamp_bianchi.py`

**Objection addressed:** "It's a large-amplitude / strong-coupling artifact, and it's
all numerics."

- **Small-amplitude scaling:** the su(3) span (8 directions) **persists down to
  A = 0.01**; tracking a fixed band/plaquette gives **‖G‖ ∼ A^2.07 (= A², leading
  order)**, while the mechanical seed — stiffness-matrix noncommutativity — scales as
  **‖[C,C]‖ ∼ A^0.90 (= A¹)**. This matches the analytic mechanism (`∂_kP ∼ O(A)` ⇒
  `F ∼ (∂P)² ∼ O(A²)`), with the complex Bloch phases (O(1)) supplying the full 8
  directions already at that order. So su(3) is present at *arbitrarily small* finite
  amplitude and vanishes continuously in the flat-vacuum limit `A→0`.
- **Non-abelian Bianchi:** in a parallel-transport gauge, the relative residual
  `‖D_[μG_νρ]‖/‖G‖` falls as `1.1e-4 → 6.8e-6` with the step `h` (≈ O(h²)) — the
  identity `D_{[μ}G_{νρ]}=0` holds, confirming the T4 field strength is a consistent
  curvature.

---

## Net effect on the paper

These convert the three most exposed claims from "asserted / single-witness / numeric"
into "demonstrated / generic-and-robust / analytically-scaled + identity-satisfying."
Cross-referenced from `T1_derivation_and_proof.md` and `T4_derivation_and_proof.md`.
The genuine-test items (D7 vacuum selection; T2 gauge-mode isotropy) are *not* covered
here — they remain the honest limitations, to be run before submission if desired.
