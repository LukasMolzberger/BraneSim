# Continuum limit of the EXACT central-force spring — the constitutive law

**Purpose.** Establish the central-force spring `½k_s(|ΔR|−αa)²` as **the** constitutive
law of the substrate (it is exact, frame-indifferent, and is what the simulation kernel
`branesim/core/action.py` actually computes), and derive its continuum expansion
`W = W₂ + W₃ + W₄ + …`. **St. Venant–Kirchhoff (StVK) is demoted to a quadratic-order
proxy:** it agrees with the spring on linear elasticity (wave speeds, Lamé) but inverts
the α-scaling of the geometric quartic, so it must not be used for soliton sizing.

Consolidation only — no new physics. Grounded in `branesim/core/action.py`,
`branesim/core/conventions.py` (closed-form dispersion), `paper/derivations/lattice_to_continuum.md`
§3.1, `paper/derivations/geometric_nonlinearity_alpha_scaling.md`, `principles.md` §1.1a.
Stencil: 6-neighbor axial-only `N_6` (the locked minimal model; the 26-neighbor sums in
`lattice_to_continuum.md` §3.3+ are superseded here).

---

## 1. Constitutive law and coarse-graining

    V = (k_s/4) Σ_p Σ_{δ∈N_6} ( |R_{p+δ} − R_p| − α a )²        (action.py)

Frame-indifferent (depends only on link length). Map: R_p = X_p + u(X_p), X_p = ap; bond
along ê_i, Δu = u(X+aê_i) − u(X) ≈ a(ê_i·∇)u. Expand about the **held** spacing a (not the
stress-free αa), so the (1−α) geometric stiffness survives. `k_s/4` = directed-bond
double-count convention.

## 2. Exact link decomposition (no truncation)

With Δu_∥ = δ̂·Δu, |Δu_⊥|² = |Δu|² − Δu_∥²:

    E_link = ½k_s(|ΔR|−αa)² = ½k_s|ΔR|²  −  k_s αa|ΔR|  +  ½k_sα²a².

Since |ΔR|² = a² + 2a Δu_∥ + |Δu|² is a degree-2 polynomial:
- `½k_s|ΔR|²` is **exactly quadratic** (harmonic + linear + const) — no α, no cubic-or-higher.
- the **entire anharmonic sector** is the Euclidean-norm term `−k_s αa|ΔR|`, carrying an
  explicit factor **α**. ⇒ every cubic+ term in the model is ∝ α.

## 3. Quadratic continuum energy W₂ (6-neighbor)

Per-link quadratic strain (lattice_to_continuum §3.1, stencil-independent):

    (s_δ²)⁽²⁾ = α(δ̂·Δu)² + (1−α)|Δu|².

Continuum (insert Δu_a ≈ a δ̂_b ∂_b u_a, ×k_s/4 per directed bond, ÷a³ cell volume), with
the axial-only sums `T⁽²⁾_ab = 2δ_ab`, `T⁽⁴⁾_abcd = 2Q_abcd` (Q=1 iff a=b=c=d):

    A_abcd = (k_s/a)[ α Q_abcd + (1−α) δ_ac δ_bd ].

Christoffel `M_ac(k)=A_abcd k_b k_d` along k = k ê₁:

    ρ_m c_L² = A_1111 = k_s/a  →  c_L² = k_s a²/m,
    ρ_m c_T² = A_2121 = (1−α)k_s/a  →  c_T² = (1−α) k_s a²/m.

**Exact cross-check:** `conventions.py` closed form ω_a²(k)=(2k_s/ρ)[α h_a+(1−α)Σ_b h_b],
h_i=1−cos k_i a, small-k → exactly the above. ✓

Effective Lamé on N_6 (axial-only): **μ(α) = (1−α)k_s/a**, **λ = C_1122 = 0**,
C_1111 = k_s/a. The stencil is structurally cubic-anisotropic
(C_1111 − C_1122 − 2C_1212 = (2α−1)k_s/a ≠ 0) — a feature (backbone #8/#16/#19), not a bug.

## 4. Cubic / quartic continuum energy

All of W₃, W₄, … descend from `−k_s αa|ΔR|`. Pure-transverse (Δu_∥=0, q≡|Δu_⊥|):

    |ΔR| = √(a²+q²) = a + q²/2a − q⁴/8a³ + …
    ⇒  −k_s αa|ΔR| ⊃ + (k_s α/8a²) q⁴.

- **W₃ = 0** on longitudinal modes (pure longitudinal |ΔR|=|a+Δu_∥| is exactly Hookean →
  no pure-longitudinal cubic) and on pure-transverse modes (even in q). The geometric
  nonlinearity is **intrinsically transverse**.
- Continuum transverse quartic (Δu_⊥ ≈ a ∂_δ u_⊥, ÷a³):

      W₄ = (k_s α / 8a) |∂u_⊥|⁴      (∝ α, vanishes at α→0; hardening, coeff > 0).

  In 4D, |Δu_⊥|² ⊃ |Δu_⊥^lat|² + (Δu_{X⁴})² ⇒ W₄ ⊃ (k_sα/4a)|∂u^lat|²(∂u_{X⁴})² — the
  lateral↔amplitude geometric quartic (backbone #17). Sign/α-scaling match the breather
  hardening Q = 8k_sα/a² > 0.

## 5. Spring ↔ StVK (spring primary, StVK demoted)

| | spring (constitutive law) | StVK (proxy) |
|---|---|---|
| quadratic in | the **stretch** `|ΔR|−αa` | the **Green–Lagrange strain** `E` |
| W₂ / wave speeds / μ,λ | — | **identical** to the spring |
| geometric quartic W₄ | **∝ α** (vanishes at α→0) | locked to μ ∝ **(1−α)** — *inverted*; naive `μ/(4ℓ₀⁴)` sizing even gives ∝1/α⁴ |

> **Crisp framing for §3.2:** StVK is the *quadratic-order proxy* of the central-force
> spring. They agree on linear elasticity (wave speeds, μ, λ) but the spring's geometric
> quartic is ∝ α whereas StVK's is ∝ (1−α). Soliton stabilization is set by the quartic, so
> the **spring law — not StVK — must be used for soliton sizing.**

Robustness: any isotropic `f(|ΔR|)` with `f(αa)=f'(αa)=0` reproduces the same W₂ (hence the
same c_L, c_T, μ, λ); the quartic coefficient depends on `f'''(αa)` and the `√(a²+q²)`
geometry, so `W₄ ∝ α` is specific to the spring and does **not** transfer to StVK.

## 6. Falsifiable check (discriminates spring ∝α from StVK ∝(1−α))

Transverse standing-wave hardening `ω²(ε) = ω₀²[1 + β(α)(ε/a)² + …]`. Predicted (spring):

    β(α) ∝ W₄/W₂^⊥ = (k_sα/8a)/((1−α)k_s/2a) = α / [4(1−α)]  — increasing in α, →0 as α→0.

StVK alternative: β = const (or ∝1/α⁴, decreasing). **Failure threshold:** if β does not
vanish as α→0 (e.g. β(0.1)/β(0.6) > 0.3 vs predicted ≈0.074), the spring-quartic derivation
is wrong. Procedure: periodic 32³, transverse seed kα·a=0.2, amplitudes ε/a∈{0.02..0.15},
α∈{0.1,0.2,0.4,0.6,0.8}, fit β(α).

## 7. Regime / open

- Valid for a|k|≪1 and |∇u|≪1; leading neglected: O((ak)⁴) dispersion, W₆ sextic.
- μ>0 needs α<1; α=1 → transverse zero modes (U(1)³); α=0 → D(k)∝I (U(3)).
- **Open (ansatz):** W₄ computed on the pure-transverse slice; the full mixed
  longitudinal–transverse quartic tensor is not assembled (expected still ∝α — same
  `−k_sαa|ΔR|` origin). General-direction W₃ from L–T cross terms claimed sub-leading,
  not fully assembled. Both fold into the §6 falsifiability run.
