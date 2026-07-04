# Lattice-Gauge Diagnostics on the Same Substrate Graph

A candidate section of the paper (final placement — dedicated section vs. cross-reference — is left open).
The claim it defends: the field strength of `field_tensors.md` can be *measured* with the standard
lattice-gauge-theory (LGT) toolbox — Wilson loops, plaquette curvature, topological-charge extraction,
finite-size scaling — **without introducing a second lattice.** The LGT objects are derived diagnostics
living on the existing substrate graph, computed from the carrier eigenbundle.

## 1. The wrong picture and the right one

**Wrong:** take the substrate lattice, then lay a separate gauge-theory lattice on top. That is a second
discretization; it obscures the point of the model and imports LGT's ontology (a fundamental gauge field
sampled on links) that this model rejects.

**Right:** use the *same* substrate links, and define LGT-like effective observables on them *after*
extracting the Berry/Wilczek–Zee connection from the local wave-branch structure. The microscopic
variables stay the node embeddings and spring link lengths; the gauge object is a diagnostic on the
existing substrate graph, not a new microscopic layer:

```
    substrate node/link data
        → local wave-branch eigenframe  U(x) = (u_1(x),…,u_N(x))
        → Berry overlap on neighbouring substrate sites
        → LGT-style holonomy / plaquette / Wilson-loop diagnostics.
```

This is exactly the link/plaquette construction of `link_holonomy.md`, now read as an *observable* rather
than a definition.

## 2. Effective link holonomies on existing links

At each coarse cell or node `x`, extract an orthonormal local frame of the `N` relevant carrier branch
eigenvectors, `U(x) = (u_1(x),…,u_N(x))`. On an existing substrate link `x → x + aî`, form the overlap

```
    M_i(x) = U(x)† U(x + aî) .
```

Numerical truncation makes `M_i` not exactly unitary; take its polar/unitary part

```
    U_i^{eff}(x) = M_i(x) ( M_i†(x) M_i(x) )^{−1/2}  ∈ U(N) .
```

This `U_i^{eff}(x)` looks formally like a lattice gauge link, but it is **not fundamental** — it is the
Berry/WZ holonomy extracted from substrate waves (`link_holonomy.md` §2, same object).

- `N = 1` → an emergent `U(1)` phase link (the trace/EM sector).
- `N = 3` → an emergent `U(3)` link, split as `U(3) ≃ U(1) × SU(3)/ℤ₃`. The approximate `SU(3)` part is
  `Ũ_i(x) = U_i^{eff}(x) / det(U_i^{eff}(x))^{1/3}`, the trace part is the phase of the determinant.

## 3. The plaquette is a Berry-curvature diagnostic

With effective holonomies on the *same* links, form the plaquette product

```
    U_{ij}^{eff}(x) = U_i^{eff}(x) U_j^{eff}(x+aî) U_i^{eff}(x+aĵ)^† U_j^{eff}(x)^† .
```

In ordinary LGT this is the basic curvature object; here it measures the **Berry curvature of the local
wave-branch bundle**. The analogy is sharp:

```
    LGT plaquette curvature   ↔   emergent Berry curvature of the substrate eigenbundle.
```

Again — not a second lattice; an observable computed on the original substrate lattice. Closed-loop
Wilson loops `W(C) = Tr ∏_{ℓ∈C} U_ℓ^{eff}` then test whether the emergent connection behaves like a gauge
field: whether curvature localizes near solitons, and whether large loops show perimeter-like, area-like,
or topological behaviour.

## 4. Regularization: same artifacts, different ontology

Both theories carry the same lattice artifacts, for opposite reasons:

```
    Lattice QCD:  lattice = regulator (artifact, sent to zero as a → 0).
    This model:   lattice = ontology  (physical; a ≠ 0 is a real short-distance scale).
```

But the *technical* consequences of `a ≠ 0` are shared: a finite Brillouin zone `k_i ∈ [−π/a, π/a]`, a
short-wavelength cutoff, finite-difference dispersion, anisotropy artifacts, and finite-volume
corrections. So the LGT artifact-handling toolbox is relevant — not because the model *is* lattice gauge
theory, but because both must handle the same kinds of lattice artifacts. The crucial divergence: **`a→0`
is not mandatory here** — the continuum limit is a long-wavelength effective description, not a definition
of the theory.

## 5. What is reusable, and what is not

**Directly useful (diagnostics on the frozen substrate state):**

- **Plaquette / Wilson-loop diagnostics** — closed-loop holonomies of the emergent Berry connection;
  curvature localization near solitons; loop-scaling behaviour.
- **Smearing / cooling / gradient-flow analogues** — suppress UV noise while preserving long-distance
  structure. Safe on the *diagnostic eigenbundle*; **must not** be applied to the *physical* substrate
  state (it would change the actual soliton).
- **Blocking / RG thinking** — coarse-grain substrate cells; track how effective branch speeds, Berry
  curvature, soliton mass, and anisotropy run with scale. Especially relevant because the model has a
  *real* cutoff and the effective observer sees only long-wavelength physics.
- **Finite-volume / finite-spacing scaling** — run the same soliton/wave-packet experiment at several
  `N`, `a`, and boundary conditions; measure observable drift. Distinguishes genuine substrate
  predictions from simulation artifacts.
- **Topological-charge extraction** — near-integer diagnostics (Chern numbers, winding numbers,
  Hopf-like invariants, loop holonomies) for solitons and Berry bundles.
- **HPC infrastructure lessons** — domain decomposition, stencil kernels, sparse linear algebra, Krylov/
  multigrid solvers, GPU data-parallel lattice operations (architectural lessons from Chroma/QUDA, not
  their physics).

**Partially reusable (use with care):**

- **Hybrid Monte Carlo.** The core dynamics are deterministic, so HMC must **not** replace the substrate
  evolution. It could still sample *ensembles of initial/boundary conditions* or thermalized
  perturbations — a statistical layer on top, never the microscopic law.
- **Gauge fixing.** The gauge freedom here is not fundamental; it is the local eigenbasis rotation of a
  degenerate carrier subspace (`carrier_construction.md` §5). Gauge-fixing-like methods can pick smooth
  eigenframes, but the safer route is to compute **gauge-invariant** quantities: traces of loop
  holonomies, plaquette eigenvalues, Chern numbers.
- **Wilson-action fitting.** Do **not** replace the elastic spring action with a Wilson gauge action. But
  one may *ask* whether the ensemble of emergent Berry links is approximately described by an effective
  Wilson-like action at long wavelengths:
  `S_{substrate}[X] → S_{eff}[U^{eff}] ≈ β Σ_p (1 − (1/N)Re Tr U_p) + …` — a *derived* effective action,
  never a postulate. A positive result here would be a genuine payoff.

**Not reusable (misleading if imported literally):**

- a fundamental Haar-measure path integral over `SU(3)` links;
- the Euclidean lattice gauge action as the *microscopic* dynamics (the microscopic action is the
  deterministic elastic/geometric brane action; its nonlinearity is the Euclidean link-length term and
  the induced-metric coupling, not a primitive non-Abelian gauge field);
- quark fields inserted as independent fundamental fermion variables;
- `a → 0` as mandatory.

## 6. Research plan (diagnostic pipeline)

1. Keep the substrate simulation unchanged; primitives stay `R(n)`, `Ṙ(n)`.
2. Compute local wave-branch eigenvectors on the substrate background (band-isolated; `link_holonomy.md`
   §1).
3. Define effective Berry links `U_i^{eff}(x)` via overlap matrices on existing substrate links (§2).
4. Compute plaquettes, Wilson loops, curvature density, topological charge from those effective links.
5. Compare scaling under lattice refinement, coarse-graining, and different stencils.
6. Test whether the emergent `U(3)` sector separates into a long-range `U(1)`-like part and a
   shorter-range `SU(3)`-like part (`[[project_u3_topology_scale]]` — the α scale-split hypothesis).

The last point is the direct experimental face of `field_tensors.md` §4: whether the substrate parameters
split the effective `U(3)` behaviour into a largely independent `U(1)` and `SU(3)` at different length
scales.

## 7. Status

- `closed` (framing): the LGT toolbox applies as *diagnostics on the existing substrate graph*, not a
  second lattice; effective links are the polar-unitary part of the carrier eigenframe overlap
  (`link_holonomy.md` §2).
- `closed` (triage): which tools are directly / partially / not reusable, with the ontology-vs-regulator
  distinction.
- `open` (placement): whether this becomes a dedicated paper section or a cross-reference to a diagnostics
  companion is undecided.
- `open` (numerical): the §6 pipeline is a program, not yet run; ties to Open derivation 4
  (`berry_reconciliation.md`) and `branesim/diagnostics/berry_holonomy.py`.
