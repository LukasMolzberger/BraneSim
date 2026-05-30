# Sprint 2 subtask 9 — D(k) is diagonal in Cartesian on the 6-neighbor lattice

**Status:** closed analytically; verified numerically.
**Scope:** structural certificate for backbone #15/#16. Closes the diagnostic
caveat that was added to paper §5.6 (BZ-link-variable Berry on real D(k)
eigenframes is identity-by-construction on the canonical lattice).

## What is claimed

On the canonical 6-neighbor axial-only cubic lattice (bonds `±êₓ, ±ê_y, ±ê_z`,
each with `|δ|=1`), the dynamical matrix is

    D_{ac}(k) · ρ / k₀  =  2 Σᵢ hᵢ · [ α · δ_{a,i} δ_{c,i}  +  (1 − α) · δ_{ac} ]

with `hᵢ ≡ 1 − cos(kᵢ a)`. This evaluates to

    D(k) = (2 k₀ / ρ) · diag(  α hₓ + (1−α) (hₓ+h_y+h_z),
                                α h_y + (1−α) (hₓ+h_y+h_z),
                                α h_z + (1−α) (hₓ+h_y+h_z)  )

i.e. **diagonal in the Cartesian basis at every k and every α**, with
eigenvectors `ê_x, ê_y, ê_z` k-independent.

Derivation in 3 lines: the contribution of a bond `δ = ±ê_i` to the bond
tensor is `α δ̂_a δ̂_c + (1−α) δ_{ac} = α δ_{a,i} δ_{c,i} + (1−α) δ_{ac}`,
which is the Kronecker pattern that has `δ_{i}` on one slot. The off-diagonal
entry (a ≠ c) requires `δ_{a,i} δ_{c,i}` to be nonzero for some axial `i`,
which forces `a = c = i` and contradicts `a ≠ c`. The second piece is
diagonal by `δ_{ac}`. Hence `D_{ac}(k) = 0` for `a ≠ c`. The diagonal entries
follow immediately.

## Numerical verification

`verify.py` reproduces `D(k)` two ways — once by explicit summation over the
six bonds, once by the closed-form expression in
`components/diagnostics/christoffel_6nn.py` — and runs five checks (33 k
samples × 6 α values for T1–T2; analytic anchors for T3–T5). All five pass:

| Test | What it checks | Tolerance | Result |
|---|---|---|---|
| T1 | `D(k)` off-diagonal entries vanish | `< 1e-14` | **PASS** (max `|D_off|` ≈ 0) |
| T2 | Closed form matches the 6-bond sum | `< 1e-14` rel | **PASS** |
| T3 | At `k ∥ [111]` the three eigenvalues are equal | `< 1e-14` rel spread | **PASS** |
| T4 | At α = 0, `D(k) ∝ I` for every k | `< 1e-14` rel spread | **PASS** |
| T5 | Fukui–Hatsugai plaquette holonomy at α = 0.20 around `k ∥ [111]` is identity (per-band and rank-3) | `< 1e-12` | **PASS** |

Eigenvalue gap along [100] and [111] (ω² convention):

| α    | gap_{100} = (ω²_L − ω²_T) / ω²_L | gap_{111} | √(1−α) = ω_T / ω_L |
|------|----------------------------------|-----------|--------------------|
| 0.00 | 0.000 | 0 | 1.000 |
| 0.10 | 0.100 | 0 | 0.949 |
| 0.20 | 0.200 | 0 | 0.894 |
| 0.50 | 0.500 | 0 | 0.707 |
| 0.80 | 0.800 | 0 | 0.447 |
| 1.00 | 1.000 | 0 | 0.000 |

Read across: along [100] the gap in `ω²` is **exactly** `α`; equivalently the
frequency ratio is `ω_T/ω_L = √(1−α)`, so the relative *frequency* gap is
`1 − √(1−α)` (10.6 % at α = 0.20). Along [111] the triplet is exactly
degenerate at every α, by the analytic isotropy of the bracket above.

## Implications for the gauge sector

- **The eigenframe is k-independent.** Therefore the Wilczek–Zee connection
  built from `D(k)` eigenvectors has zero curvature and identity holonomy
  on any k-space plaquette, *not just on small ones* — this is a structural
  statement, not a numerical accident, and it holds for every α.
- **What runs with α is the eigenvalue spread**, not the eigenframe. At
  α = 0 the spread vanishes everywhere (full U(3)); at α = 1 it is maximum
  (U(1)³ with one zero-frequency direction at each k); at α = 0.2 the
  triplet is near-degenerate (gap ≈ 10.6 % at [100], 0 % at [111]).
- The meaningful gauge object remains the per-wavepacket complex envelope
  `Ψ ∈ ℂ³` (backbone #18). Its effective Hamiltonian `H_eff(k₀, α)` enters
  via the rotating-wave / multiscale reduction, where the imaginary unit
  introduces structure not present in `D(k)` itself. That derivation is the
  next step on the gauge thread; it is *not* done here.

## Files

- `components/diagnostics/christoffel_6nn.py` — closed-form `D(k)` and
  per-axis eigenvalue accessor for the canonical lattice.
- `test-runs/sprint2_subtask9_d_of_k_diagonal/verify.py` — bond-sum
  reference + 5 checks.
- `test-runs/sprint2_subtask9_d_of_k_diagonal/results.json` — full results.

## What this does and does not establish

- **Establishes:** the BZ link-variable construction in paper §5.6 returns
  identity for *every* α on the canonical lattice; revisiting paper §5.6
  on this stencil cannot, in principle, distinguish U(1)³ from U(3).
- **Establishes:** the U(3) → U(1)³ crossover backbone #16 describes is an
  eigenvalue-spread phenomenon, not an eigenframe-rotation phenomenon.
- **Does not establish:** that U(3) gauge structure emerges in the complex
  envelope. That requires deriving `H_eff(k₀, α)` and computing WZ on its
  eigenframe — a separate piece of work (physics-derivation).
- **Does not establish:** anything about non-linear sectors (solitons,
  contraction field). This is the linear-regime structural certificate only.