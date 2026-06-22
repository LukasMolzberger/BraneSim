# Bell Bridge — The vertex tensor M and the Bell-state lock J2 (D1)

Working note. Computes the pair-creation tensor `M` explicitly from the substrate
geometric nonlinearity. **Result: `M ∝ δ` (symmetric) — the vertex emits the
symmetric Bell state `|Φ⁺⟩`, not the antisymmetric singlet `|Ψ⁻⟩`.** The two are
local-unitary (waveplate) equivalent and give identical Tsirelson violation, but
the substrate fixes the convention. This *computes* J2 (the D1 residual) from the
Lagrangian rather than importing it. All algebra below numerically verified.

Supersedes the `M ∝ ε`/singlet assumption used earlier: the explicit geometric
vertex is built from a real symmetric norm and produces `δ`, not `ε` (§3).

## 1. The geometric link energy and its cubic vertex

A link along unit vector `ê` carries relative displacement `δu = u_j − u_i`.
With `x = (ê·δu)/a` (longitudinal strain) and `y² = |δu⊥|²/a²` (transverse),

    |ΔR|/a = √((1+x)² + y²) = 1 + x + y²/2 − x y²/2 − y⁴/8 + …

The anharmonic (geometric-nonlinear) energy is the norm term `−k_s α a |ΔR|`
(`∝ k_s α`; `α=0` is exactly linear). Its leading cubic piece is the **three-wave
vertex**

    V₃ = +(k_s α / 2a) · Σ_links (ê·δu) |δu⊥|² .

It couples one **longitudinal** leg `(ê·δu)` to **two transverse** legs `|δu⊥|²`
— the down-conversion vertex `pump → A + B`, with coefficient `∝ k_s α/a`
(consistent with the geometric scaling; the quartic `y⁴` term is the next,
four-wave, order). *Feature/caveat:* the three-wave vertex requires a
**longitudinal pump leg** (there is no `y³` all-transverse cubic term — the
expansion is even in `y` at `x=0`); a transverse pump couples only through its
longitudinal projection, or via the quartic (four-wave) process.

## 2. The vertex tensor M

Extract the resonant pair-creation term (`pump → A + B`, both daughters
transverse). The transverse factor is

    |δu⊥|² = δ_{ab} δu⊥^a δu⊥^b ,

so with the longitudinal pump amplitude `ψ_p` the pair-creation operator is

    C† = M_{ab} a†_{A,a} a†_{B,b} ,   M_{ab} = (k_s α / 2a) · ψ_p · δ_{ab} .

Thus **`M ∝ δ_{ab}` (symmetric, identity tensor); `M_anti = 0`.**

## 3. Which Bell state: |Φ⁺⟩, not the singlet

The `δ`-contraction creates

    M_{ab}|a⟩_A|b⟩_B ∝ |H⟩_A|H⟩_B + |V⟩_A|V⟩_B = |Φ⁺⟩   (symmetric) ,

not the antisymmetric singlet `|Ψ⁻⟩ ∝ |HV⟩−|VH⟩` (which is the `ε`-contraction).
Verified numerically.

**Why `δ` and not `ε`.** The geometric coupling is the *real symmetric norm*
`|δu⊥|²`. It is invariant under SO(2) real transverse rotations — under which
*both* `δ` and `ε` are invariant — but the norm specifically produces `δ`. The
*full* SU(2)-Jones invariance that would force `ε` (the unique SU(2)-invariant) is
**not** a symmetry of the real geometric vertex; only the real SO(2) is. So the
naive "polarization-isotropy ⇒ singlet" expectation fails at the explicit level:
the real norm gives the symmetric `δ` ⇒ `|Φ⁺⟩`.

## 4. This is fine — |Φ⁺⟩ gives Tsirelson and is LU-equivalent to the singlet

- **Maximally entangled / no-signalling:** `ρ_A = ρ_B = ½𝟙`, so the marginal is
  flat (the B2 marginal). ✓
- **Linear analyzers:** `E(a,b) = +cos 2(a−b)` (vs the singlet's `−cos 2(a−b)`);
  CHSH `|S| = 2√2`. ✓ (verified)
- **Local-unitary equivalence:** all four Bell states are related by a local
  unitary on one arm — a fixed **waveplate** at one detector maps `|Φ⁺⟩ ↔ |Ψ⁻⟩`.
  So the Bell-test violation is identical; "singlet vs `|Φ⁺⟩`" is a convention,
  and the substrate fixes it to `|Φ⁺⟩`.
- **Bose-natural:** `|Φ⁺⟩` is *symmetric* under `A↔B`, pairing with a symmetric
  two-mode emission — the natural bosonic output (type-I-like). (The antisymmetric
  singlet would instead need an antisymmetric mode assignment.)

So replacing the singlet by `|Φ⁺⟩` changes nothing observable in the standard
linear-polarizer Bell test (`|S| = 2√2`), and the downstream machinery (B4
threshold, the junction's no-signalling, the `cos 2(a−b)` form) carries over with
`+cos 2(a−b)` in place of `−cos 2(a−b)`.

## 5. Group-theory tool (which symmetry ↔ which Bell state)

For two doublets `2 ⊗ 2 = 3 ⊕ 1`:

- `ε` (antisymmetric) = the unique **SU(2)**-invariant ⇒ singlet `|Ψ⁻⟩`;
- `δ` (symmetric) is **SO(2)**-invariant (real rotations) ⇒ symmetric Bell state
  `|Φ⁺⟩` (a triplet member).

The vertex's actual symmetry decides which. The explicit geometric vertex has the
real SO(2) and the symmetric norm, so it picks `δ ⇒ |Φ⁺⟩` (§3).

## 6. Validity regime and falsifiable fidelity

On-axis (`[100]`) long-wavelength, the transverse doublet is degenerate (cubic
4-fold symmetry), so the `δ`-coupling is isotropic among the two transverse
directions ⇒ a clean `|Φ⁺⟩`. Off-axis propagation and finite `k·a` split the
doublet, making the coupling anisotropic ⇒ Bell-state fidelity

    F = 1 − O((Δ/Γ)²) ,   Δ = transverse splitting, Γ = vertex bandwidth.

Falsifiable: entanglement degrades off-axis / toward short wavelength (the SPDC
birefringence/walk-off analogue), exact only on-axis long-wavelength.

## 7. Residual debt (now minimal under D1)

- Confirm the **longitudinal-pump** three-wave channel and its phase matching
  (or the transverse-pump longitudinal projection), tying `ω_p = ω_A + ω_B`,
  `k_p = k_A + k_B` to the Paper I dispersion.
- Compute the off-axis anisotropy coefficient in `F` (§6) from the transverse
  branch-splitting `Δ(k·a, angle)`.

## Result

`M` computed from the geometric cubic vertex `(ê·δu)|δu⊥|²`: `M ∝ δ` ⇒ the vertex
emits the symmetric Bell state `|Φ⁺⟩`, local-unitary-equivalent to the singlet,
giving `|S| = 2√2` and flat marginals. **J2 is now derived from the substrate
Lagrangian** (not imported), with the corrected Bell-state identity `|Φ⁺⟩`.

## References (wire into `references.bib` if this graduates)

- `2 ⊗ 2 = 3 ⊕ 1`; `ε` unique SU(2)-invariant, `δ` symmetric/SO(2)-invariant
  (any group-theory / QM text).
- Bell states are local-unitary equivalent (Nielsen & Chuang).
- Paper I — geometric nonlinearity `∝ k_s α`, transverse-branch degeneracy and
  dispersion.
- Standard SPDC: type-I (symmetric, `|Φ⟩`-like) vs type-II (`|Ψ⟩`-like); the
  geometric vertex is the type-I-like symmetric case. Birefringence/walk-off
  degrades entanglement (the empirical analogue of §6).
