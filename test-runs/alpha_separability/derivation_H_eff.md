# Derivation: complex-envelope H_eff(k₀,α), its connection, and where spin lives

**Source:** physics-derivation agent, 2026-06-03. Gate for Track B of the
α-separability experiment (`SPEC.md`). Bindings: `branesim/core/conventions.py`
(175–211, closed-form ω_a²), `branesim/core/action.py` (141–205, central-force
EOM), `test-runs/sprint2_subtask9_d_of_k_diagonal/report.md` (certified
zero-curvature on the real D(k) eigenframe), `BACKBONE.md` #16/#18/#19/#20,
`paper/05_geometric_phase_and_gauge_diagnostics.tex` (§5.6 caveat closed here).

---

## Headline result

The rotating-wave reduction yields H_eff whose coefficient matrices (carrier
detuning ΔΩ, group velocity Vⱼ, dispersion Mⱼₗ) are **all diagonal in the same
k-independent Cartesian frame** as D(k). Therefore:

1. The **k-space (Brillouin-zone) Berry/Wilczek–Zee connection is identically
   zero for all α** — the complex `i` of the rotating-wave step does NOT create
   base-space curvature. Inherited directly from the sprint2_subtask9 certificate.
2. Non-triviality is **fibre-internal**: per-axis U(1) phase (φ_a = arg Ψ_a) and
   relative-phase SU(3), transporting in **physical (x,t)**, not in k.
3. **Spin ½ is NOT a linear-layer object.** Ψ∈ℂ³ carries the **vector (J=1,
   spin-1)** rep of SO(3); a 2π rotation gives +𝟙 (phase 0), α-independent. The
   π phase of spin ½ requires a soliton-layer **π₁(SO(3))=ℤ₂** rotation holonomy
   (hedgehog/Skyrme, orientation–isospin locking) — the standard Finkelstein–
   Rubinstein structure.

---

## Part 1 — The envelope equation

Map (BACKBONE #18): ξⁱ(x,t) = Re[Ψⁱ(x,t) e^{i(k₀·x−ω₀t)}], Ψ∈ℂ³ slowly varying,
band-limited |q|≪|k₀|,π/a.

Assumptions: (1) linear regime (drop geometric quartic #17); (2) scale separation
ε=|∇Ψ|/(k₀|Ψ|)≪1, width W≫a; (3) D(k) diagonal in Cartesian basis — exact for the
6-neighbor axial stencil, not an approximation; (4) single dispersion sheet per
axis, no inter-channel linear coupling; (5) rotating-wave / +ω analytic-signal
branch (this introduces the `i`); (6) narrowband Taylor to O(q²).

Linearized EOM: ∂²_t ξ̂_a = −ω_a²(k) ξ̂_a, with
ω_a²(k) = (2k_s/ρ)[α h_a + (1−α)Σ_b h_b], h_b = 1−cos(k_b a).

Result:
```
i ∂_t Ψ = H_eff(k₀,α; −i∇) Ψ
H_eff = ΔΩ − i V_j ∂_{x_j} − ½ M_{jl} ∂_{x_j}∂_{x_l}
ΔΩ = diag(ω_a(k₀) − ω₀),  V_j = diag(v^{(a)}_j),  M_{jl} = diag(M^{(a)}_{jl})
v^{(a)}_j = k_s a [α δ_{aj} + (1−α)] sin(k_{0,j} a) / (ρ ω_a(k₀))
```
All three coefficient objects are diagonal 3×3 matrices in (ê_x,ê_y,ê_z).

Regime: ε≪1, W≫a, k₀a ≲ O(1). Leading neglected term O(ε³ω₀) + dropped quartic.

**Open caveat:** a single common ω₀ across three non-degenerate branches is an
ansatz away from [111]. Off [111] the ω_a(k₀) differ (L–T gap), ΔΩ≠0, channels
dephase at rate |ω_a−ω_b|. Coherent Ψ∈ℂ³ requires k₀∥[111] (exact degeneracy) or
α→0 (D∝I) or a beat time longer than the experiment.

### Falsifiable prediction P1 (cheap, dispersion-analyst)
Carrier along [100] at k₀a=π/4: group-velocity ratio of transverse to longitudinal
envelope drift = √(1−α) = **0.894 at α=0.2**. Launch a band-isolated wavepacket,
track per-Cartesian-channel centroid drift, fit speed. Fail if measured ratio
differs from √(1−α) by >5% at W/a≥8.

---

## Part 2 — The connection (k-space)

Eigenvectors of H_eff = simultaneous eigenvectors of ΔΩ,V_j,M_{jl}, all diagonal in
the fixed Cartesian basis {e_a} with ∂_{k₀}e_a=0. Hence
A_{ab}(k₀) = i⟨e_a|∇_{k₀}|e_b⟩ = 0, F = dA + A∧A = 0 identically, all α.

**The k-space connection is trivial** — exactly the sprint2_subtask9 result, now
inherited by H_eff. The `i` provides a per-channel U(1) phase that evolves in
physical (x,t) and a relative phase between axes (SU(3)); the complex structure
lives in the ℂ³ fibre, not the k₀ base curvature.

Validity: exact on the canonical 6-neighbor axial stencil for all k₀,α. Adding
diagonal-shell bonds would make D(k) off-diagonal and could give nonzero k-curvature.

### Falsifiable prediction P2 (consistency, berry-validator)
Fukui–Hatsugai plaquette holonomy of the H_eff eigenbundle on any k₀-plane loop =
identity to machine precision, every α∈[0,1]. Fail if ‖U_Γ − 𝟙‖ > 1e-10.

---

## Part 3 — Which loop carries holonomy: k-space vs SO(3)

Under R∈SO(3): Ψ'ⁱ = Rⁱ_j Ψʲ and k₀→Rk₀ (ξ is a genuine displacement 3-vector;
the 𝟑/L=1/T_{1u} irrep, BACKBONE #20).

- **(a) k-loop:** A≡0 ⇒ holonomy = 𝟙, all α. No spin from the BZ.
- **(b) SO(3) 2π loop:** D^{(1)}(R_n(2π)) = exp(−i 2π n·L^{(1)}) = +𝟙.
  Carrier-triplet holonomy = +𝟙: **phase 0, spin-1, α-independent.**

To get −𝟙 (π, spin ½) needs a spinor (J=½, the 𝟐 of SU(2)); no half-integer rep
exists in the linear lateral triplet (ℂ³ = J=1 vector ⊕ J=0 trace). Spin ½ must
come from the **soliton layer**: hedgehog/Skyrme winding (ξⁱ=f(r)x̂ⁱ) whose
orientation–isospin locking gives a 2π→−1 holonomy via π₁(SO(3))=ℤ₂. The linear
layer **sets up the bundle (vector rep)** but cannot produce π.

```
Hol_{k-loop} = 𝟙;   Hol_{SO(3),2π} = D^{(1)}(2π) = +𝟙  (spin-1, α-independent).
Spin ½ (phase π) is NOT realized at the linear layer — it is a soliton-layer
π₁(SO(3))=ℤ₂ rotation holonomy.
```

This is the standard Skyrme structure (pion = spin-0/1 meson; baryon spin ½ from
Finkelstein–Rubinstein / Hopf quantization). The split is a coherence *feature*.

### Falsifiable prediction P3 (decisive, berry-validator)
Adiabatically rotate a band-isolated wavepacket's carrier orientation through 2π,
measure the accumulated geometric phase of the envelope. **Predicted: 0 mod 2π
(trivial, spin-1) at all α.** Fail-the-derivation: any measured π. Confirming
spin ½ requires the *soliton* experiment (rigidly rotate the hedgehog) — predicted
π there iff the configuration has odd Skyrme/winding parity. **This is a new L5
target.**

---

## Implications for the experiment program
- Track A (U(1)×SU(3) split, ∝ α) stands; add P1 (group-velocity √(1−α) test).
- Track B at the linear layer is now mostly **analytically settled**: k-holonomy
  trivial (P2 consistency check), SO(3) holonomy spin-1 (P3). The α-dependent
  linear gauge object that survives is the **fibre-internal U(1)/SU(3) phase
  transport in physical (x,t)** — the magnetic-curl channel.
- **Spin ½ migrates to L5**: rigid hedgehog rotation → check for ℤ₂ (π) holonomy
  with odd winding parity. This replaces the "spin window" in the joint α-test.