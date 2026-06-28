> **Archived 2026-06-28:** the required soliton configuration space and quantization have not been constructed.

# Spin-½ and chirality of the confined defects

**Status:** `closed` for the topological framing (which homotopy group supplies spin-½);
`conditional`/`open` for its realization on a converged mode.
**Sources:** `gauge_color/sections/02_gauge.tex` (`subsec:spinorial-holonomy`),
`bell/sections/03_worldtube.tex` (`sec:worldtube`), `[[project_spin_half_is_soliton_layer]]`,
`[[project_retrocausal_worldtube_interpretation]]`.

---

## 1. Spin-½ is a soliton-layer ℤ₂ holonomy, not a linear envelope phase

At the **linear** level, the polarization eigenbundle of a 2-level subspace has holonomy in the
double cover of `SO(3)`, so a `2π` loop can flip the sign of the transported state and a `4π` loop
restores it (`subsec:spinorial-holonomy`). But this linear `4π` aspect is a property of the
**envelope**, and the established result is that the linear envelope is spin-1 (vector), with
identically zero Berry curvature in the BZ ∀α — spin-½ is **not** carried by the linear envelope
(`[[project_spin_half_is_soliton_layer]]`).

Spin-½ lives at the **soliton layer** as a `ℤ₂` holonomy:

    π₁(SO(3)) = ℤ₂   (Finkelstein–Rubinstein).

A confined defect carrying internal orientation picks up a sign under a `2π` spatial rotation —
fermionic statistics — because the configuration space of the localized mode has `π₁=ℤ₂`. For the
`U(1)` centerline this is the `2π` **tumble** of the donut/line cross-section about its axis as one
encircles it; for the `SU(3)` texture it is the spin-from-isospin of the iso-rotated hedgehog
(`matter_mass` D1). In both cases spin-½ is a property of the *bound* defect, not of any free linear
wave — consistent with quantization being emergent (`subsec:continuous-waves`).

This is `closed` as a topological statement and `open` as a computation: it must be verified that a
converged soliton actually carries the `ℤ₂` holonomy under a `2π` rotation (Paper IV open problem 5).

## 2. Chirality sets the arrow of time along the defect

The substrate dynamics is time-symmetric (the timelike link is a central-force spring like the
spatial ones; `[[project_temporal_link_4d_spring]]`, `[[feedback_time_symmetry_no_damping]]`). The
defects are 4-D objects, obtained by sweeping the spatial-slice cross-section along the timelike
direction (which is *not* part of the spatial angular frame, `connections_holonomy.md` §1): the
`U(1)` charge centerline lifts to a **worldsheet** (codim-2), the `SU(3)` point texture to a
**worldline** (codim-3), and a closed toroidal spatial mode (the electron-as-torus) to a
**worldtube**. The arrow of time is not put in by hand; it is set by the **chirality** of
the soliton's carrier rotation along its worldtube axis (the same temporal carrier rotation that is
the `U(1)` charge phase, `[[project_complex_u1_from_time]]`).

Matter and antimatter are **opposite-chirality worldtubes** of the same object (Bell paper
`sec:worldtube`, `[[project_retrocausal_worldtube_interpretation]]`). This ties the confinement
picture to the worldtube interpretation: a confined particle is a finite, closure-locked
(`ωT=2πn`) worldtube whose chirality fixes both its charge sign and its time orientation — and whose
closure condition is exactly the linchpin of the time-link binding (`binding.md` §2).

## 3. Ledger

| Item | Status |
|---|---|
| Spin-½ = `π₁(SO(3))=ℤ₂` soliton holonomy (not linear envelope) | `closed` (topological) |
| `ℤ₂` holonomy realized on a converged mode | `open` (numerical, Paper IV #5) |
| Chirality sets arrow of time along the worldtube | `conditional` (worldtube interpretation) |
| Matter/antimatter = opposite-chirality worldtubes | `conditional` (Bell paper) |
