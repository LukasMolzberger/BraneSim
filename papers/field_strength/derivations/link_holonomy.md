# Link Variables and Plaquette Curvature from the Carrier

Goal: turn the carrier of `carrier_construction.md` into a lattice connection whose plaquette curvature
is the field-strength tensor. This is the shared spine for both sectors in `field_tensors.md`.

## 1. The carrier eigenframe

From `carrier_construction.md`: node `n=(n_t,n_x,n_y,n_z)`, spacing `a`, axial links `±ê_μ`,
`μ∈{t,x,y,z}` — **8 signed links per node**. The carrier `ψ(n) = δR(n) + iδṘ(n)/ω ∈ ℂ³` is the
phase-space amplitude of the lateral triplet; its `i` is the timelike-link rotation
(`[[project_complex_u1_from_time]]`).

At each node the band-isolated carrier selects a **polarization eigenframe**:

- rank-1 `|u(n)⟩ = ψ(n)/|ψ(n)| ∈ ℂ³` for the trace/EM sector (one carrier vector);
- rank-`d` orthonormal frame `{u_1,…,u_d}(n)` spanning the (near-)degenerate carrier subspace for the
  non-Abelian sector (`d=3` for the colour triplet).

The physics is *transport of the eigenframe* from node to node; that transport carries the gauge
information — not the real displacement geometry, which is curvature-free in the BZ
(`berry_reconciliation.md`).

## 2. Link variable

Define the link variable on the directed link `n → n+ê_μ` as the (normalized) eigenframe overlap

```
    U_μ(n) = ⟨u(n) | u(n+ê_μ)⟩ / |⟨u(n) | u(n+ê_μ)⟩|              (rank-1, Abelian)
    U_μ(n) = polar-unitary part of  M_{ij} = ⟨u_i(n)|u_j(n+ê_μ)⟩  (rank-d, non-Abelian, ∈ U(d))
```

The non-Abelian `U_μ(n)∈U(d)` is the lattice Wilczek–Zee parallel-transport matrix (the unitary part of
the overlap `M`). The reverse link is `U_{-μ}(n+ê_μ)=U_μ(n)^†`. **Spacelike and timelike links enter on
the same footing** — `μ` ranges over all of `{t,x,y,z}`; the only difference is that the timelike
overlap steps in `n_t` and so probes the carrier's time-rotation directly (this is what makes `F_{0i}`
nonzero; see `field_tensors.md` §2).

**Gauge (rephasing) law.** A local change of carrier frame `|u(n)⟩→V(n)|u(n)⟩`, `V∈U(d)`, gives

```
    U_μ(n) → V(n) U_μ(n) V(n+ê_μ)^† .
```

This is exactly the lattice-gauge transformation of a connection; the continuum limit (§4) reproduces
`A_μ → V A_μ V^† − i(∂_μV)V^†`, i.e. `A_μ → A_μ + ∂_μχ` in the Abelian case. The freedom `V` is the
carrier-frame redundancy identified in `carrier_construction.md` §5 — the gauge symmetry *is* the
arbitrariness of the local `U(3)` eigenframe.
[Open derivation 1: `U_μ` single-valued only off band-crossings — band-isolation requirement.]

## 3. Plaquette

The elementary oriented loop in the `(μ,ν)` plane (corners `n, n+ê_μ, n+ê_μ+ê_ν, n+ê_ν`):

```
    □_μν(n) = U_μ(n) · U_ν(n+ê_μ) · U_μ(n+ê_ν)^† · U_ν(n)^† .
```

This uses **two** of the 8 links at each visited node. The time–space plaquettes `(0,i)` use the
temporal links and produce the `E`-field; the space–space plaquettes `(i,j)` produce `B`. Without the 2
temporal links there is no `F_{0i}` — the temporal spring is structurally mandatory.

## 4. Continuum limit → curvature

Writing `U_μ(n)=exp(−ia A_μ(n+½ê_μ))` (the link variable is the path-ordered exponential of a
connection along the link, sampled at the midpoint), the standard lattice expansion gives

```
    □_μν(n) = exp( −i a² F_μν(n) + O(a⁴) ),
    F_μν = ∂_μ A_ν − ∂_ν A_μ − i[A_μ, A_ν].
```

The Baker–Campbell–Hausdorff commutator `[A_μ,A_ν]` is `0` in the rank-1 Abelian case and supplies the
`f^{abc}` self-coupling in the non-Abelian case (`field_tensors.md` §3).
[Open derivation 2: verify the `O(a⁴)` remainder is a genuine higher-derivative term — clover/improved-
plaquette check, no spurious symmetric or parity-odd piece.]

## 5. Antisymmetry is structural, not assumed

Reversing the plaquette orientation, `(μ,ν)→(ν,μ)`, traverses the loop backwards:
`□_νμ(n)=□_μν(n)^†`. Taking logs, `F_νμ=−F_μν`. **The antisymmetry of the field-strength tensor is
inherited from the orientation of the elementary lattice loop** — a rank-2 antisymmetric tensor is the
only object a 2D plaquette can produce. This is the cleanest statement of *why* both the Faraday tensor
and the QCD field-strength tensor are antisymmetric: they are curvatures of oriented loops on the
stencil. (Bivector/2-form language: `F = ½F_μν dx^μ∧dx^ν`, the wedge enforcing antisymmetry.)

## 6. What is established vs open here

- **Established (restatement):** the lattice-gauge dictionary link→plaquette→curvature; the gauge law;
  the BCH expansion. Standard; the novelty is the *physical* `U_μ` from the phase-space carrier on the
  spring stencil (`carrier_construction.md`).
- **Open:** definedness across band-crossings (1); `O(a⁴)` artifact structure (2); feed `status.md`.