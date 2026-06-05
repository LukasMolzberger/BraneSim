# Where the geometric nonlinearity enters, and its exact α-scaling

**Claim (exact, lattice ground truth).** In the central-force spring model the
*entire* anharmonic sector is carried by a single term, `−k_s α a |ΔR|`, and is
therefore **exactly proportional to α** (the rest-length ratio). The harmonic
sector is α-independent. Consequently the Skyrme-class geometric quartic scales as
**∝ α** (vanishing at α→0), opposite to the quadratic "geometric stiffness" which
scales as ∝ (1−α).

Date: 2026-06-04. Owner: physics-derivation / dispersion-analyst. Primary model:
the lattice spring energy (NOT the continuum StVK approximation).

---

## 1. Exact decomposition

One link, held spacing `a`, rest length `ℓ₀ = αa`, bond unit vector `δ̂`,
relative displacement `Δu = u_{n+δ} − u_n`, `Δu_∥ = δ̂·Δu`,
`|Δu_⊥|² = |Δu|² − Δu_∥²`:

    E_link = ½ k_s (|ΔR| − αa)²,   |ΔR| = √((a + Δu_∥)² + |Δu_⊥|²).

Expand the square algebraically (exact, no Taylor truncation):

    E_link = ½ k_s |ΔR|²  −  k_s αa |ΔR|  +  ½ k_s α²a².

Now use the *exact* identity `|ΔR|² = a² + 2a Δu_∥ + |Δu|²` (the squared Euclidean
norm is a degree-2 polynomial in the displacements):

    ½ k_s |ΔR|² = ½ k_s a²  +  k_s a Δu_∥  +  ½ k_s |Δu|²     ← EXACTLY quadratic.

Therefore:

- The squared-norm term is harmonic + linear + constant. **It contains no cubic or
  higher term and no α.**
- The **norm** term `−k_s αa |ΔR|` is the *only* source of the square-root
  nonlinearity, and it appears multiplied by `α`. **Every anharmonic (cubic+) term
  in the model is exactly ∝ α.**
- `½ k_s α²a²` is a constant.

**α is the coupling constant of the geometric nonlinearity.** At `α = 0` (zero rest
length / max prestress) `E_link = ½ k_s |ΔR|²` is exactly harmonic and isotropic —
the dynamics is *linear* and supports no solitons. Nonlinearity turns on linearly
with α.

## 2. The nonlinearity is intrinsically transverse

For pure longitudinal motion (`Δu_⊥ = 0`): `|ΔR| = |a + Δu_∥|` is **linear** in
`Δu_∥` — axial stretching of a central-force spring is exactly Hookean. The square
root only bites through `|Δu_⊥|`. Hence the geometric nonlinearity is a
transverse-self / longitudinal↔transverse coupling, never pure-longitudinal. In the
4D embedding "transverse" includes the `X⁴` amplitude direction, so the
cross-term in `|Δu_⊥|⁴ ⊃ 2|Δu_⊥^{lat}|² (Δu_{X⁴})²` **is** the lateral↔amplitude
geometric quartic of backbone #17.

## 3. Quadratic spectrum (recovers the locked result)

Hessian at `Δu = 0` (from `½k_s|ΔR|²` minus the quadratic part of `k_sαa|ΔR|`):

- longitudinal stiffness: `k_s`              (from `½k_s|ΔR|²` only; α-independent)
- transverse stiffness:  `k_s − k_s α = k_s(1−α)`
  ( `+k_s` from `½k_s|ΔR|²`, `−k_s α` from the quadratic part of `−k_sαa|ΔR|`,
    namely `−½k_s α |Δu_⊥|²` )

⇒ `c_L² = k_s a²/m`, `c_T² = (1−α) k_s a²/m`. Matches `conventions.py`,
`lattice_to_continuum.md` §3.1, backbone #16. The "geometric stiffness ∝ (1−α)"
and the "geometric nonlinearity ∝ α" are *both* faces of the norm term
`−k_s αa |ΔR|` — but at different orders, with opposite α-scaling.

## 4. Quartic coefficient (the correction)

Pure-transverse expansion (`Δu_∥ = 0`, `q ≡ |Δu_⊥|`):

    |ΔR| = √(a² + q²) = a + q²/(2a) − q⁴/(8a³) + …
    ⇒ −k_s αa|ΔR| ⊃ +(k_s α/8a²) q⁴.

Discrete transverse quartic coefficient = **+k_s α/(8a²) > 0** (hardening, ∝ α).
Continuum (`Δu_⊥ ≈ a ∂_δ u_⊥`, density per cell `a³`):

    W_4 ∝ (k_s α / a) |∂u_⊥|⁴      (∝ α; vanishes at α → 0).

**Consistency:** sign and α-scaling match the breather result `Q ∝ k_s α/a² > 0`
(`[[project_pythagorean_breather_go]]`).

**Correction:** `vsh_channel_decomposition.md` §2.5 takes the quartic prefactor as
`~μ/(4ℓ₀⁴)` with `ℓ₀ = αa`, i.e. **∝ 1/α⁴ (diverging as α→0)** — the α-dependence
is inverted. The hedgehog-radius scaling `R_h/a ~ (ℓ₀/u₀)^{2/3}` derived from it
must be re-derived with `W_4 ∝ k_s α/a`. (Backbone #17 claims only a *Derrick*
scaling `λ^{+1}`, not an α-scaling, so it is silent rather than wrong; the α-scaling
here is a refinement.)

## 5. Why the StVK identification can disagree

The lattice energy `½k_s(|ΔR|−αa)²` is hyperelastic (depends only on the
frame-indifferent length `|ΔR|`) but is **quadratic in the stretch**, not in the
Green–Lagrange strain `E`. True StVK (`W = μ tr E² + (λ/2)(tr E)²`) locks its
quartic to `μ` by the `E ⊃ ½∂uᵀ∂u` structure; the central-force spring does not
obey that locking at quartic order. So a naive "quartic ∝ μ ∝ (1−α)" from assuming
StVK is **not** the lattice truth. The exact spring quartic is ∝ α (§4). This is
relevant to `OPEN_PROBLEMS.md` C1 (StVK vs the actual constitutive law).

## 6. Consequences

- **Soliton sweet spot α ≈ 0.5–0.8 is derived, not heuristic.** Binding needs the
  nonlinearity ∝ α (large α); confinement needs positive transverse stiffness
  ∝ (1−α) (not α→1). The product `α(1−α)` peaks at α=0.5; the useful window is the
  broad plateau α≈0.5–0.8.
- **α-estimator (D2) gains explicit α-dependence in the nonlinear normalization.**
  The fibre-metric correction `C` that the linear theory could not fix
  (`alpha_holonomy_estimator.md`) is sourced by this quartic, hence `C = C₀ + O(α·
  amplitude²)` — so the nonlinear computation can in principle determine α.

---

**Files:** `branesim/core/conventions.py`; `paper/derivations/lattice_to_continuum.md`
§3.1; `BACKBONE.md` #16/#17; `paper/derivations/vsh_channel_decomposition.md` §2.5
(correction target); `paper/derivations/alpha_holonomy_estimator.md` §5/§7;
`OPEN_PROBLEMS.md` C1, D2.