# L5 target — spin ½ as a ℤ₂ soliton rotation holonomy

**Status:** scoped (2026-06-03), **gated** on a confined soliton to rotate.
**Layer:** 5 (soliton). **Owner:** soliton-hunter (produce the soliton) →
berry-validator (apply the rotation holonomy). Feeds `paper/baryon_simulation_roadmap.md` Phase 2.

## Why this target exists
The α-separability derivation (`derivation_H_eff.md`, Part 3) proved that spin ½
is **not** a linear-layer object: the linear envelope `Ψ∈ℂ³` carries the SO(3)
vector (J=1, spin-1) rep, so a 2π rotation gives `+𝟙` at every α (confirmed
numerically, P3). The π phase of spin ½ can only arise as a topological
`π₁(SO(3))=ℤ₂` holonomy of a **real-space soliton** under rigid rotation — the
standard Finkelstein–Rubinstein mechanism (a Skyrmion is quantized as a fermion
iff it carries odd winding). This is the expected Skyrme structure: phonon triplet
= spin-1 meson, baryon spin ½ = topological.

## The test
1. **Obtain a confined soliton** with odd winding parity (B=1 hedgehog/Skyrme,
   `ξⁱ = f(r) x̂ⁱ` + the X⁴ twist of BACKBONE #20). *This is the gate* — provably
   no wide elasticity-only hedgehog confines (`OPEN_PROBLEMS C1`,
   `[[project_c2_skyrme_no_confinement]]`); needs the breather/Skyrme-stabilized
   line to deliver one first.
2. **Rigidly rotate the field configuration by 2π** about a fixed axis. The
   hedgehog locks the spatial direction to the internal triplet index (`x̂ⁱ`), so
   a spatial rotation automatically co-rotates the internal/color frame — this
   orientation–isospin locking is what can make the holonomy ℤ₂ rather than the
   trivial spin-1 of the unlocked envelope.
3. **Compute the rotation holonomy** with the *already-built, already-validated*
   kernel `branesim/diagnostics/berry_holonomy.py` (the P3 diagnostic; its
   synthetic spin-½ control correctly returns −𝟙 at 2π, so it can detect a
   fermionic holonomy when one is present). Feed it the soliton's eigenframe /
   field and rotate.

## Predicted result & falsifier
- **Spin ½ (baryon):** 2π → `−𝟙` (phase π), 4π → `+𝟙`, iff the configuration has
  **odd** Skyrme/winding parity.
- **Falsifier:** a confined B=1 soliton that rotates 2π → `+𝟙` (spin-1) would
  contradict the baryon = spin-½ expectation and impugn the orientation–isospin
  locking story.

## CRITICAL implementation caveat (from the P3 build — do not skip)
The SO(3) rotation holonomy MUST be computed with the **open-path** Fukui–Hatsugai
link product (θ: 0 → 2π, **no closing link**). It directly yields `D^{(j)}(R(2π))`.
**A closed-loop FH product composes the holonomy with its own inverse and returns
the identity for BOTH spin-1 and spin-½** — i.e. a closed-loop construction gives a
silent false-negative and would hide a real fermionic soliton. (Contrast: the P2
*k-space* plaquette is a genuine closed loop because the k-eigenframes are truly
periodic; the orientation loop is not closed in the same sense.) This is encoded
in `berry_holonomy.py`; preserve it when applying to the soliton.

## Joint α-window connection
This supplies the "spin ½" window of the `SPEC.md` §4 joint test — which, post-
derivation, is a **soliton-layer** window (alongside confinement `Q>0`), not a
linear one. The linear layer's contribution is the color split (∝ α) plus the
flat-curvature / spin-1 baseline (P2/P3) that makes any soliton π holonomy
unambiguously topological.