> **Archived 2026-06-28:** statuses in this ledger were invalidated by the audit.

# Confinement (Paper VI) — Status

## HAVE (`closed`)

- **Connection/holonomy spine.** Every structure is a connection on the single `ℂ³` bundle
  (Levi-Civita / Berry / Wilczek–Zee); confinement = holonomy constraints; the lattice regularizes
  every holonomy at the core; `α` is the `U(1)³↔U(3)` reducibility dial. Levi-Civita is dynamical,
  Berry/WZ are diagnostic (no "geometry=gauge" identity). See `connections_holonomy.md`.
- **Three transitions, one parameter.** `U(1)↔SU(3)` (`α`-reducibility), discrete→continuum
  (scale-nonlinearity, quantization-at-the-transition), `SH↔Cartesian` (Christoffel barriers). See
  `scale_transition.md`, `connections_holonomy.md` §2.
- **Single-nonlinearity theorem.** All anharmonicity = the norm term `−k_sαa|ΔR|`, exactly `∝α`,
  transverse-only, P-even (sources energy not winding). See `nonlinearity.md`.
- **Lattice = physical UV regulator.** `|ΔR|=√(a²+|Δu|²)` analytic; gradients cap at `~1/a`; defect
  cores are finite with no true singularity; collapse blocked by `λ≥a`. See `nonlinearity.md` §5.
- **`U(1)` centerline.** **No `|Ψ|⁴` potential, no VEV, no mass, no Meissner** (`closed`) — the trace
  `U(1)` is a massless, ungauged, *global* defect on a localized envelope. Core-energy *scaling*
  `E/L=π(1−α)(k_s/a)n²f₀²ln(R/a)+(π/16)(k_sα/a³)n⁴f₀⁴` with `a` the physical cutoff, healing length
  `ξ=(|n|f₀/2)√(α/(1−α))` — **derived within a localized-envelope ansatz** (`closed-in-ansatz`,
  *not* a stable-particle existence claim). See `u1_vortex_core.md`.
- **`SU(3)` point texture.** `π₃(SU(3))=ℤ` worldline; SH-vs-winding dichotomy; combined-`SO(3)`
  hedgehog → 1-D radial ODE (all `closed`); `R_h/a=κ(A/a)√(α/(1−α))` is a Derrick/lattice scaling
  *within the hedgehog reduction* (`closed-in-ansatz`); `SU(3)` sets no length; lattice anti-collapse.
  See `su3_texture.md`.
- **Kinematic colour confinement — now quantitative.** No free colored asymptotic state; the colour
  survives at most a closed dephasing length `L_deph=2(3−2α)/(√3αg(k̂)k₀)` (finite off `[111]`,
  diverges on it), and the figure of merit `ℳ=(colour holonomy)×(coherence length)=O(1)` with
  `α,g,k₀` cancelling — colour rotation and coherence are conjugate. See `kinematic_dephasing.md`.
- **`ℤ₃` triality lock.** `q≡−t (mod 3)` — the firmest, dynamics-free confinement statement (binds
  representation, not position).
- **Spatial-only binding is a closed NEGATIVE.** Three independent routes all fail; the P-even
  modulus structure forbids a spatial binder. See `binding.md` §1.
- **Spin-½ = `π₁(SO(3))=ℤ₂` soliton holonomy** (not the linear envelope). See `spin_chirality.md`.

## CONDITIONAL

- **Time-link binding** `χ=2α(ω/ω_*)²(A₀/a)²>1` (energetic) + soft winding-lock `Q=κ(A)B` — rests on
  the worldtube closure axiom `ωT=2πn`. Binding↔gravity co-vary in `α` (A4a test). See `binding.md` §2.
- **Semilocal `β`** of the `U(1)` vortex — a dimensional model (`β_eff=4/n²` at Derrick scale;
  `(α/(1−α))(f₀/a)²` lattice-pinned), not a derived fluctuation spectrum.
- **Chirality → arrow of time; matter/antimatter = opposite-chirality worldtubes** (Bell paper).

## MISSING / OPEN

- Converged stable unit-winding `U(1)` eigenstate (branch-selection for the semilocal vortex).
- Converged stable `B≠0` `SU(3)` texture (the hard gate; the topological binding leg is vacuous
  without it).
- Linearized scalar-vs-vector fluctuation spectrum about the `U(1)` core (the proper `β`).
- `ℤ₂` spin holonomy realized on a converged mode.
- `κ` (Derrick prefactor), `Θ` (carrier-phase correlation), and the `c₁`/`κ(A)` normalization.
- Genuine colour *dynamics* beyond kinematic confinement (running coupling, string tension, area
  law, hadron spectrum) — explicitly out of scope.

## Honesty headline

The paper closes the **structure** of confinement (one nonlinearity, two lattice-regulated defects,
kinematic + triality confinement, no spatial binder) and the binding **condition** (`χ>1` via the
time link). It does **not** close the binding **proof** (rests on the closure axiom) or the
**existence** of the two converged eigenstates (numerical, deferred). Say this plainly in §6.

## Named results (promote to theorem/conjecture in the paper)

- **Theorem 1 (single anharmonicity).** All anharmonicity descends from `−k_sαa|ΔR|`, exactly `∝α`,
  transverse, P-even.
- **Theorem 2 (lattice UV regularization).** `|ΔR|=√(a²+|Δu|²)` analytic ⇒ defect cores finite, no
  singularity, collapse blocked by `λ≥a`.
- **Theorem 3 (no spatial binder).** Modulus-only local anharmonicity gives no attractive
  trace↔traceless binder (three independent routes).
- **Theorem 4 (kinematic no-free-colour).** Coherence ⊥ colour-activity ⇒ finite `L_deph`, `ℳ=O(1)`
  conjugacy bound; plus the `ℤ₃` triality lock `q≡−t (mod 3)`.
- **Conjecture (time-link binder).** Worldtube closure `ωT=2πn` flips the cross-vertex to attractive
  when `χ=2α(ω/ω_*)²(A₀/a)²>1` — conditional on the two-time BVP producing the closed family.
- **Numerical targets.** Existence of a stable `n=1` trace defect and a stable `B≠0` texture.

## External anchors to cite (from the critique)

- **Berry / Wilczek–Zee / Wilson-loop** (connection language for transported subspaces, with the
  "transport connection ≠ dynamical gauge field" caveat) — Berry 1984; Wilczek–Zee 1984;
  arXiv:1910.13991 (Wilson loop & WZ phase). Supports the holonomy framing.
- **Semilocal / electroweak strings** — arXiv:hep-ph/9904229 (Achúcarro–Vachaspati). Constrains the
  `U(1)` defect story; aligns the `β<1`/ungauged caveat with the standard field.
- **Skyrmion / Finkelstein–Rubinstein / fermionic Hopf-soliton quantization** —
  arXiv:hep-th/0503067. Supports `π₃` winding (not harmonics) and soliton-layer spin-½.
- **Derrick's theorem** — Derrick 1964. Baseline obstruction; our distinctive twist is lattice-UV
  anti-collapse, kept distinct from the IR anti-expansion.

## Open derivations

*These are the per-paper bridge entries for the confinement bridge (`[[project_open_problems_tracker]]`).*

1. **Converged semilocal `U(1)` vortex** — branch-selection machinery; verify `E_core∝α` and the
   `β` verdict (small-amplitude `n=1` on-axis vs `A₀≳1` drift).
2. **Converged `B≠0` texture** — periodic-clamped eigenstate (not seed-and-watch); verify
   `R_h∝A√(α/(1−α))` and compute the `ℤ₂` spin holonomy + the soft winding-lock `Q=κ(A)B`.
3. **Time-link binding test** — `Δu_∥∝ω²` across the closure-index family; `χ>1` crossing; the
   binding↔gravity α-co-variation (A4a).
4. **Proper semilocal `β`** — linearized scalar-vs-vector mode spectrum about the relaxed `U(1)` core
   (the load-bearing fluctuation analysis the semilocal literature calls for).
5. **Dephasing-length / `ℳ` verification** — direction scan `L_deph(k̂)·g(k̂)=const` and the
   `ℳ=O(1)` conjugacy bound (tests the `κ_color∝k₀δω/ω` proportionality, the one ansatz in
   `kinematic_dephasing.md`).
6. **Two-time boundary-value solver** — the linchpin: make the worldtube-closure axiom `ωT=2πn`
   either *emerge* or *fail*. Without it, the time-link binder stays a conjecture (critique's
   single highest-value missing piece).
