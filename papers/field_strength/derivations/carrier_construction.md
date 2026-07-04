# The Complex Carrier from Bloch Modes and Phase Space

The load-bearing note. It earns the complex carrier `ψ` — and its `U(3)` frame freedom — from the real
substrate, so that `link_holonomy.md` can build a connection without assuming a complex structure by
fiat. The logical order is: Bloch mode → one `U(1)` phase → Berry connection from the *eigenvector's*
phase freedom → what the rest length does to that eigenvector → why that is necessary but not sufficient
→ the phase-space carrier that is sufficient → `O(3) → U(3)` → vector vs connection → dynamic phase is
not yet curvature.

## 1. Bloch modes carry one `U(1)` phase; the Berry connection is the eigenvector's phase freedom

Linearize the spring dynamics about a reference lattice. A perturbation is a Bloch mode

```
    δR_p^A(k) = u^A(k) · e^{ik_μ p^μ} = u^A(k) · e^{i(−ω p^0 + k_1 p^1 + k_2 p^2 + k_3 p^3)} ,
```

where `A = 0,1,2,3` labels the ambient displacement component and `μ = 0,1,2,3` labels the intrinsic
lattice direction, with `k_μ = (−ω, k_1, k_2, k_3)` (the `−ω` in the `s_0 = −1` timelike slot).

Two clarifications guard the entry:

- **One `U(1)`, not `U(1)⁴`.** The exponential `e^{ik_μ p^μ}` is a *single* complex phase `θ`. It has
  four directional gradients `∂_μ θ = k_μ`, but its value lives in *one* circle, `θ ∼ θ + 2π`. Bloch
  theory gives a `U(1)` phase over a 4D wave-vector space — not four independently gaugeable phases.
  Calling it `U(1)⁴` would require four separately conserved, independently gaugeable phases, which the
  single exponent does not provide. (Distinguish this also from the three *carrier* phases of the lateral
  triplet in the degenerate `α=0` limit — those are three clocks in `ℂ³`, not the four lattice-direction
  gradients of one plane wave; and the physical connection there is still the single trace `U(1)`, since
  the traceless curvature is flat at `α=0`, `field_tensors.md` §4.)

- **The Berry connection is not built from the raw plane-wave phase.** A pure scalar Bloch wave with a
  *constant* polarization `u(k)` has trivial curvature — `e^{ik_μ p^μ}` alone carries no geometry. The
  physically relevant `U(1)` is the phase ambiguity of the *eigenvector*:

  ```
      u(k) → e^{iχ(k)} u(k) .
  ```

  That local phase freedom is the gauge `U(1)`. The Berry connection and curvature are

  ```
      A_μ(k) = i⟨u(k)|∂_{k_μ}u(k)⟩ ,     F_{μν}(k) = ∂_{k_μ}A_ν − ∂_{k_ν}A_μ .
  ```

  So the connection becomes meaningful only when the eigenmode has nontrivial internal structure —
  `u^A(k)` a genuinely `k`-dependent polarization vector, or coupled branches whose mixture twists as `k`
  changes. Because `R_p ∈ ℝ⁴` is vector-valued, the perturbation *does* have several polarization
  components, so this bundle exists. The question is whether `u(k)` actually twists — which is what the
  rest length controls.

## 1a. What the rest length does to `u(k)` — the stiffness geometry

The rest length does **not** touch the `U(1)` phase directly. It changes the linearized *bond stiffness
matrices*, and those decide whether the Bloch eigenmodes factor into separate directional channels or
mix.

Start from a central-force spring between neighbours, `V_{p,μ} = (κ_μ/2)(|R_{p+ê_μ} − R_p| − ℓ_μ)²`, and
linearize about a background `R̄`. With background bond vector `B_μ^A = R̄_{p+ê_μ}^A − R̄_p^A`, length
`L_μ = |B_μ|`, and unit direction `n̂_μ^A = B_μ^A/L_μ`, the bond's stiffness matrix is a longitudinal
projector plus a transverse projector weighted by the **pre-tension factor** `(1 − ℓ_μ/L_μ)`:

```
    K_μ^{AB} = κ_μ [ n̂_μ^A n̂_μ^B  +  (1 − ℓ_μ/L_μ)(δ^{AB} − n̂_μ^A n̂_μ^B) ] .
```

This is the key object.

- At **rest**, `ℓ_μ = L_μ`, the transverse part vanishes: `K_μ^{AB} = κ_μ n̂_μ^A n̂_μ^B`. The spring
  resists only motion along its own bond direction. Displacements transverse to the bond cost nothing.
- Under **pre-tension/pre-compression**, `ℓ_μ ≠ L_μ`, a transverse stiffness
  `κ_μ(1 − ℓ_μ/L_μ)(δ^{AB} − n̂_μ^A n̂_μ^B)` appears: off-bond displacements now contribute to the
  restoring force. **This is where dimensional coupling enters.**

With `ℓ_μ/L_μ = α`, the transverse weight is exactly `(1 − α)` — the same `(1 − α)` that sets the
transverse cone speed `c_T² = (1 − α)k_s a²/m` in the core substrate model. The Bloch operator is

```
    D(k) = Σ_μ f_μ(k) K_μ ,     f_μ(k) = s_μ · 4 sin²(k_μ a/2)  →  s_μ k_μ²  (continuum) ,
```

with `s_μ = (−1, +1, +1, +1)` the intrinsic Lorentzian sign. The Bloch eigenproblem is
`Σ_B D^{AB}(k) u_n^B(k) = λ_n(k) u_n^A(k)`.

- If all `K_μ` are **simultaneously diagonalizable** (`[K_μ, K_ν] = 0`), the polarizations can be chosen
  almost constant — the directions stay separate phase channels and `u(k)` is essentially `k`-independent.
- If the rest lengths make the `K_μ` **non-commuting** (`[K_μ, K_ν] ≠ 0`), they cannot be diagonalized at
  once. As `k` moves, the weighted sum `Σ_μ f_μ(k)K_μ` rotates in polarization space, so the eigenvector
  `u_n(k) ≠ const`. **This `k`-dependent polarization mixing is exactly the ingredient a Berry connection
  needs.** So the rest length is the natural tuning parameter for moving between separated dimensional
  phases and mixed phase geometry.

**Stencil caveat.** On the *canonical 6-neighbour axial stencil* the `n̂_μ` are the coordinate axes, so
every `K_μ` is diagonal in that basis and they commute at *every* `α` — `u(k)` stays the fixed L/T
coordinate directions, exactly as the core substrate model records. Genuine `k`-dependent mixing needs
off-axis (diagonal-shell) bonds. Either way the conclusion of §1b is identical.

## 1b. Why that is necessary but not sufficient — the broken fourth link

The natural chain "rest length → stiffness geometry → `k`-dependent mixing → Berry connection" has three
sound links and a broken fourth one.

Every `K_μ` is real symmetric and every `f_μ(k)` is real, so `D(k)` is real symmetric and its
eigenvectors can be chosen **real** at every `k`. A real eigenvector — even a `k`-dependent one — has

```
    A_μ(k) = i⟨u(k)|∂_{k_μ}u(k)⟩ = (i/2)∂_{k_μ}⟨u|u⟩ = 0 ,
```

so the `U(1)` Berry curvature over the Brillouin zone is identically zero **for all `α`**, mixing or not
(BACKBONE #16; `berry_reconciliation.md` §3). `k`-dependence of `u(k)` is **necessary** for a nonzero
connection but not **sufficient**.

What a real `k`-bundle can still carry is a `ℤ₂` (Zak) holonomy — a real sign, not a `U(1)` curvature.
That `ℤ₂` is a property of the **linear envelope, which is a spin-1 vector**, and it must **not** be
mistaken for spin-½. Spin-½ is a *different* `ℤ₂`: a real-**space** rotational holonomy
(`π₁(SO(3)) = ℤ₂`) of a bound soliton's framing relative to the far-field lattice, living at the soliton
layer — not this `k`-space bundle (`[[project_spin_half_is_soliton_layer]]`; `matter_mass` D1; the gear
picture `χ = ½(θ_env − θ_lat)` is developed there, not here). **Three distinct objects that must not be
conflated:**

```
   U(1) core winding   π₁(U(1)) = ℤ    — charge, vortex core (soliton layer)
   k-space Zak ℤ₂      real sign        — spin-1 envelope holonomy (this linear layer)
   spin-½ ℤ₂           π₁(SO(3)) = ℤ₂  — soliton framing vs far-field (soliton layer)
```

None of the three is a `U(1)` field strength. The genuine `U(1)`/WZ curvature is supplied by the complex
phase-space carrier of §2 — the `i` of `ψ = δR + iδṘ/ω` is the time-link rotation, and it is that complex
link, not the real `k`-space mixing, that makes a plaquette holonomy nontrivial.

## 2. Position + velocity = a complex amplitude (phase space)

A single oscillator degree of freedom has, at each instant, a displacement and a velocity. For carrier
frequency `ω`, position and velocity are 90° out of phase over the cycle, and the standard,
non-arbitrary way to package them is the complex phase-space amplitude

```
    ψ = q + i q̇/ω ,     (real part = configuration, imaginary part = phase-shifted motion).
```

For the lattice carrier, subtract the static background and apply this per embedding direction `a`:

```
    ψ^a_n = δR^a_n + i δṘ^a_n / ω ,     ψ_n = (ψ^a_n)_a .
```

The phase of `ψ` is the oscillator's **phase in its cycle**; a free mode runs as `ψ(t) ~ e^{−iωt}`.

**This is where the `i` comes from, concretely.** The velocity `δṘ` is the time-derivative of the
configuration; pairing it with `δR` via `i` is exactly encoding the carrier's rotation along the timelike
worldtube axis (`[[project_complex_u1_from_time]]`). A purely spatial snapshot of `δR` alone is real; the
imaginary part is supplied by stepping in time. Hence two time slices are minimal, and the complex
carrier is not added by hand — it is the phase-space state of the spring oscillator.

## 3. The carrier space: which components, and why three

The relevant carrier is the **lateral triplet** — the three transverse carrier polarizations of Paper
III. Promoting each to its phase-space amplitude gives

```
    ψ_n = (ψ^x_n, ψ^y_n, ψ^z_n) ∈ ℂ³ .
```

(The trace `U(1)` and the photon use a rank-1 reduction of this same carrier; see `field_tensors.md`
§2. The linear envelope is a spin-1 vector — `[[project_spin_half_is_soliton_layer]]` — so `ℂ³` is the
natural carrier space, not an imposed `SU(3)` representation.)

## 4. Real degeneracy gives `O(3)`; complexification gives `U(3)`

This is the precise chain from real substrate to complex gauge group.

- **Real triplet, exact degeneracy** `ω_x = ω_y = ω_z`. The real quadratic energy
  `½(|δR|² + |δṘ|²/ω²)` is invariant under real orthogonal mixing of the three components. The frame
  freedom of the *real* degenerate triplet is therefore only `O(3)`.
- **Complex phase-space carrier.** Once each real oscillator is promoted to `ψ^a = δR^a + iδṘ^a/ω`, the
  triplet lives in `ℂ³`, and the transformations preserving the complex norm `|ψ|²` (the scaled energy)
  are `ψ → Uψ`, `U†U = 𝟙`, `U ∈ U(3)`.

The complexification is exactly what enlarges the real `O(3)` frame freedom to the complex `U(3)` (the
extra generators are the `i`-times-symmetric part). So the "promote position+velocity to a complex
amplitude" step of §2 is not cosmetic: it is the reason the carrier frame group is `U(3)` and not `O(3)`.
Without it there is no complex gauge structure at all.

## 5. Vector phase vs frame transport: `U(1)` and `U(3)` are different objects

A standing trap (critique #4): **do not identify the carrier vector `ψ ∈ ℂ³` with the connection matrix
in `u(3)`.** They are different objects.

- **Single carrier vector.** Normalizing one carrier state `|u(n)⟩ = ψ_n/|ψ_n| ∈ ℂ³` and transporting it
  node-to-node tracks one *overall phase*. The node-to-node overlap `⟨u(n)|u(n+ê_μ)⟩` is a scalar; its
  phase is a **`U(1)` Berry connection**. This is the trace/EM sector.
- **Local three-frame.** Transporting a full orthonormal frame `{u_1, u_2, u_3}` spanning the degenerate
  `ℂ³` subspace, the overlap is a `3×3` matrix `M_{ij}(n) = ⟨u_i(n)|u_j(n+ê_μ)⟩` whose unitary part lies
  in `U(3)`. The **connection is `u(3)`-valued**; the carrier remains a vector. This is the object that
  splits into the trace `U(1)` and the traceless `SU(3)` (`field_tensors.md`).

The EM `U(1)` is the overall phase of one carrier vector; the colour `SU(3)` is the relative orientation
of a transported three-frame. One object (the `U(3)` connection) carries both, but the carrier vector is
never the connection.

## 6. Dynamic phase is raw material, not yet curvature

A caveat that keeps the paper honest. The phase-space amplitude has a phase, but for a free mode that
phase is just the clock,

```
    ψ(t) ~ e^{−iωt}    ⇒    uniform phase advance    ⇒    plaquette holonomy = 1,    F = 0 .
```

The clock phase alone is **not** a gauge field. A nonzero field strength requires the *relative* carrier
phase/frame between neighbours to be **path-dependent** — a nonzero plaquette holonomy
`ψ → ψ e^{iθ_loop}`, `θ_loop ≠ 0`. The construction of §§2–5 supplies the *carrier* (a genuine `U(3)`
connection over `(x,t)`); whether its curvature is nonzero is the statement that the carrier frame
actually rotates across `(x,t)`, which it does precisely because the complex link is the time-rotation.

This is why the field strength lives over physical `(x,t)` and not the BZ: the `k`-space frame is real
(no rotation, `F ≡ 0`); the `(x,t)` frame rotates (the time-link phase), so its plaquette holonomy is
nontrivial.

## 7. Status

- `closed` (construction): `ψ^a = δR^a + iδṘ^a/ω` from phase space; `i` = time-link rotation; real
  degeneracy `O(3)` ⇒ complex carrier `U(3)`; single-vector `U(1)` phase vs transported three-frame
  `U(3)` connection kept distinct.
- `closed` (scaffold, §1): Bloch mode carries one `U(1)` phase (four gradients, one circle), not
  `U(1)⁴`; the Berry connection is the eigenvector's phase freedom `u → e^{iχ}u`, not the raw plane-wave
  phase.
- `closed` (guard, §§1a–1b): the rest length sets the transverse stiffness `(1−α)` of `K_μ` and hence
  whether `u(k)` mixes in `k`; but real-symmetric `D(k)` forces the `U(1)` BZ curvature to zero `∀α`
  regardless — `k`-mixing is necessary, not sufficient; the `U(1)`/WZ curvature is the complex `(x,t)`
  carrier, not the real `k`-space bundle. Three `ℤ₂`-adjacent objects kept apart.
- `closed` (spectrum, feeds `field_tensors.md` §4): the exact degenerate triplet (full `U(3)` frame
  freedom) is the `α=0` limit; the splitting `λ_A − λ̄ ∝ α g(k̂)` grows with `α` and *breaks* the
  degeneracy, and it is that splitting — not the frame freedom — that generates the traceless colour
  curvature. Degeneracy and colour move in opposite directions (near-degenerate WZ; `[111]`
  coherence-vs-colour tension).
- `open`: band-isolation of the three-frame off carrier band-crossings (feeds `link_holonomy.md` §2).
