> **Archived 2026-06-28:** coarse-graining and classical mode discreteness do not derive quantum mechanics or confinement.

# The discrete→continuum transition as a nonlinearity along the scale axis

**Status:** `closed` as a framing of established results (carrier–envelope coarse-graining,
scale-dependent geometric quartic, emergent quantization); does **not** claim a full renormalization
-group derivation.
**Sources:** `core/sections/04_wave_structure.tex` (`subsec:continuous-waves`),
`core/sections/03_substrate_model.tex` (`eq:W4`, lattice cutoff),
`gauge_color/derivations/color/vsh_channel_decomposition.md` (§3 multiplet splitting `α(a/w)²`),
`[[project_complex_u1_from_time]]`, `[[project_non_probabilistic_quantization_emergent]]`.

> **Thesis.** The map from the discrete real lattice to the continuous complex `ℂ³` field of quantum
> mechanics is itself a transition along a **scale axis**: the coarse-graining is a modulus/phase
> (carrier–envelope) decomposition, the only coupling it generates is scale-dependent (`∝1/a`), and
> the **discrete quantum numbers emerge only at the bounded-eigenproblem end of this transition**,
> not at the field level. (Carefully: this is a statement about *where the descriptions become
> available and discrete*, not an eliminativist claim that EM or QM are "mere bookkeeping.")

---

## 1. The coarse-graining map (where QM sits on top of the lattice)

The substrate state is discrete, real node displacements `ξ_p∈ℝ⁴`. The quantum-mechanical field is
the smooth complex envelope `Ψ(x,t)∈ℂ³`, obtained by the carrier–envelope split

    ξ(x,t) ≈ Re[ Ψ(x,t) · ε · e^{−iω₀t} ]     (band-isolated carrier ω₀).

Two things happen *at* this map, not before it:

1. **The complex `i` is born.** The imaginary unit is the carrier quadrature
   (`Ψ ≈ ξ + (i/ω₀)ξ̇`), i.e. rotation along the timelike worldtube axis
   (`[[project_complex_u1_from_time]]`). The substrate slice is real; the *complex `U(1)` carrier
   description* (and with it the gauge phase) first becomes **available** only after temporal
   coarse-graining — it is not present on a single real slice. This is a statement about when the
   description becomes available, not a claim that electromagnetism is an illusion.
2. **The map is nonlinear.** `Ψ` is a modulus/phase reorganization of node data; for `O(1)`
   amplitudes the envelope is a nonlinear functional of the discrete field (the carrier and envelope
   are not independent). This is the "nonlinearity along the scale dimension."

## 2. Scale as an axis; the geometric quartic as a scale-dependent coupling

Treat the resolution length `ℓ` as a coordinate. The theory's anharmonicity flows along it:

- **UV end (`ℓ→a`).** Lattice spacing `a`, Brillouin cutoff `|k|≤π/a`. The norm
  `|ΔR|=√(a²+|Δu|²)` is analytic at every node, so the would-be defect singularity is cut off — the
  gradient saturates at `~1/a`, core energies are finite (`nonlinearity.md` §5). This is the
  *anti-collapse* end.
- **IR end (`ℓ≫a`).** Solitons self-select a width `R_h≫a` by Derrick balance
  (`su3_texture.md`); the continuum/StVK description is valid there. This is the *anti-expansion* end.
- **The coupling carries the scale.** `W₄ = k_s α/(8a)|∂u_⊥|⁴` has an explicit `1/a`; the
  dimensionless anharmonicity at scale `ℓ` is

      W₄/W₂ ~ α (a/ℓ)²   (e.g. the multiplet splitting Δω/ω ~ α(a/w)², VSH note §3).

  Anharmonic (and `O_h`-anisotropy) effects **grow toward the UV** and are negligible in the deep IR.
  The single dial `α` sets the strength of this scale-flow, and is the *same* dial that controls the
  `U(1)↔SU(3)` reducibility (`connections_holonomy.md` §4) — the two transitions share a parameter.

This is a coarse-graining flow with the geometric quartic as the leading generated operator. We do
**not** claim a closed RG `β`-function; we claim the scale-dependence is explicit and `∝α`, with the
lattice as the hard UV fixed length.

## 3. Quantization emerges at the transition

The continuous field carries a **continuum** of modes up to `|k|~π/a` — nothing is discrete at the
field level (`subsec:continuous-waves`). Discreteness appears only through **particle formation**:

1. confinement (kinematic + topological + lattice-regulated, this paper) imposes boundary conditions
   on the continuous medium;
2. the self-adjoint linearized eigenproblem on the bounded soliton then has a **discrete** spectrum;
3. the angular geometry admits only integer-labeled VSH multiplets, and the topology only integer
   winding `B∈π₃(SU(3))=ℤ`, triality-locked charge `q≡−t (mod 3)`.

So the familiar quantum numbers `(J,P; SU(3)\text{-irrep}; Q_{U(1)}; B)` are the **discreteness of the
bound modes that live at the discrete→continuum→bound-state ladder**, not a primitive of the
substrate (`[[project_non_probabilistic_quantization_emergent]]`). The discreteness appears at the
**bounded-eigenproblem end** of the scale transition: below `a` there is only the lattice; far above
`a` there is a continuous classical field with a *continuum* of modes; the discrete quantum spectrum
is the bound-state structure that confinement carves out in between. (We claim emergence of the
discrete labels, not that "quantization is created" ex nihilo.)

## 4. The transitions are coupled

The three transitions are not independent:

- The branch-splitting `g(k̂)` that gates colour-activity (the `U(1)↔SU(3)` transition) is itself
  `k`-dependent — a *scale* statement; the soliton lives at `k~1/R_h`, so "coherence vs
  colour-activity" is evaluated at a particular point on the scale axis.
- The curvilinear (`SH↔Cartesian`) reduction is what makes the bound-state eigenproblem of §3
  tractable; its Christoffel terms (`connections_holonomy.md` §2) set the angular barriers that
  quantize `J,L`.

All three meet in the bound soliton: it is the object where the prestress sets the gauge content, the
scale flow sets the size and the validity of the continuum, and the curvilinear connection sets the
angular quantum numbers.

## 5. Ledger

| Item | Status |
|---|---|
| Coarse-graining = carrier–envelope (modulus/phase) map; `i` born here | `closed` |
| Map is nonlinear at `O(1)` amplitude | `closed` (qualitative) |
| Geometric quartic is scale-dependent `∝α(a/ℓ)²` | `closed` |
| Lattice `a` = hard UV length (anti-collapse); `R_h` = IR (anti-expansion) | `closed` |
| Quantization emerges at the transition (boundary-condition mechanism) | `closed` |
| Three transitions share `α` and meet in the soliton | `closed` (organizational) |
| Closed RG `β`-function / fixed-point flow | **not claimed** (`open`) |
