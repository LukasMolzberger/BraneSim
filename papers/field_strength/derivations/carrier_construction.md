# The Complex Carrier from Phase Space

The load-bearing new content of the second attempt. It earns the complex carrier `ψ` — and the `U(3)`
frame freedom — from the real substrate, so that `link_holonomy.md` can build a connection without
assuming a complex structure by fiat.

## 1. The substrate state is real

Each node carries a real embedding coordinate and velocity,

```
    R_n(t) = (X_n, Y_n, Z_n, A_n) ∈ ℝ⁴,    Ṙ_n(t) ∈ ℝ⁴.
```

The linearized spring dynamics are `δR̈ = −D·δR`, with `D(k)` the Bloch dynamical matrix. Because the
springs are central-force and the energy is a real quadratic form, **`D(k)` is real symmetric**: its
normal-mode polarization vectors `u(k)` can be chosen real. This is the fact behind the vanishing
Brillouin-zone curvature (`berry_reconciliation.md` §3): nothing complex lives in the displacement
basis or over `k`-space.

So the complex structure the gauge fields need is **not** in `R` and it is **not** in the BZ. It has to
be constructed, and the construction must be physical.

## 1a. What the rest length does to `u(k)` — and why it is still not curvature

It is worth being explicit about what the rest length does to the polarization eigenvectors, because the
natural chain "rest length → stiffness geometry → `k`-dependent mixing → Berry connection" has three
sound links and a broken fourth one.

Linearize a bond `n→n+ê_μ` about the background link `B_μ = R̄_{n+ê_μ}−R̄_n`, length `L_μ=|B_μ|`, unit
direction `n̂_μ = B_μ/L_μ`. Its stiffness matrix is a longitudinal projector plus a transverse projector
weighted by the **pre-tension factor** `(1−ℓ_μ/L_μ)`:

```
    K_μ^{AB} = κ_μ [ n̂_μ^A n̂_μ^B + (1 − ℓ_μ/L_μ)(δ^{AB} − n̂_μ^A n̂_μ^B) ] .
```

With `ℓ_μ/L_μ = α` the transverse weight is exactly `(1−α)` — the same `(1−α)` that sets the transverse
cone speed `c_T² = (1−α)k_s a²/m` in the core substrate model. The Bloch operator is
`D(k) = Σ_μ f_μ(k) K_μ` with `f_μ(k) = s_μ · 4 sin²(k_μ a/2)`. So the rest length is indeed the knob
that decides whether the `K_μ` share an eigenbasis: commuting `K_μ` keep the directions as separate phase
channels with `k`-independent `u(k)`; non-commuting `K_μ` make `u(k)` rotate in polarization space as `k`
moves.

(On the *canonical 6-neighbour axial stencil* the `n̂_μ` are the coordinate axes, so every `K_μ` is
diagonal in that basis and they commute at *every* `α` — hence the Christoffel matrix is diagonal and the
polarizations are the fixed L/T coordinate directions, exactly as the core substrate model records.
Genuine `k`-dependent mixing needs off-axis diagonal-shell bonds. Either way the conclusion below is
identical.)

**The broken fourth link.** `k`-dependent polarization mixing does *not*, by itself, give a Berry
connection here. Every `K_μ` is real symmetric and every `f_μ(k)` is real, so `D(k)` is real symmetric
and its eigenvectors can be chosen **real** at every `k`. A real eigenvector — even a `k`-dependent one —
has

```
    A_μ(k) = i⟨u(k)|∂_{k_μ}u(k)⟩ = (i/2)∂_{k_μ}⟨u|u⟩ = 0 ,
```

so the `U(1)` Berry curvature over the Brillouin zone is identically zero **for all `α`**, mixing or not
(BACKBONE #16; `berry_reconciliation.md` §3). `k`-dependence of `u(k)` is *necessary* for a nonzero
connection but not *sufficient*: what a real bundle can still carry is a `ℤ₂`/real (Zak) holonomy — but
that is a property of the **linear envelope, which is spin-1 (a vector)**, and it must *not* be mistaken
for spin-½. Spin-½ is a *different* `ℤ₂` again: a real-**space** rotational holonomy (`π₁(SO(3))=ℤ₂`) of
the bound soliton's framing relative to the far-field lattice, living at the soliton layer — not this
`k`-space bundle (`[[project_spin_half_is_soliton_layer]]`; `matter_mass` D1). Neither `ℤ₂` is a `U(1)`
field strength. The genuine `U(1)`/WZ curvature is supplied by the complex phase-space carrier over
physical `(x,t)` constructed in §2: the `i` of `ψ = δR + iδṘ/ω` is the time-link rotation, and it is that
complex link — not the real `k`-space mixing — that makes a plaquette holonomy nontrivial.

**One notational guard.** The raw Bloch phase `e^{ik_μ p^μ}` is a *single* `U(1)`: its four directional
gradients are `∂_μθ = k_μ`, but its value lives in one circle, so it is not a `U(1)⁴`. Nor is it the
gauge `U(1)` — that is the rephasing freedom `u(k) → e^{iχ(k)}u(k)` of the *eigenvector*, not the
plane-wave exponential. (Both also differ from the `U(1)³` of the `α=0` limit in `field_tensors.md` §4,
which is the three *carrier* phases of the decoupled lateral triplet — three carriers, not four
lattice-direction phases.)

## 2. Position + velocity = a complex amplitude (phase space)

A single oscillator degree of freedom has, at each instant, a displacement and a velocity. For carrier
frequency `ω`, position and velocity are 90° out of phase over the cycle, and the standard, non-arbitrary
way to package them is the complex phase-space amplitude

```
    ψ = q + i q̇/ω,        (real part = configuration, imaginary part = phase-shifted motion).
```

For the lattice carrier we subtract the static background and apply this per embedding direction `a`:

```
    ψ^a_n = δR^a_n + i δṘ^a_n / ω,        ψ_n = (ψ^a_n)_a .
```

The phase of `ψ` is the oscillator's **phase in its cycle**; a free mode runs as `ψ(t) ~ e^{−iωt}`.

**This is where the `i` comes from, concretely.** The velocity `δṘ` is the time-derivative of the
configuration; pairing it with `δR` via `i` is exactly encoding the carrier's rotation along the
timelike worldtube axis (`[[project_complex_u1_from_time]]`). A purely spatial snapshot of `δR` alone is
real; the imaginary part is supplied by stepping in time. Hence two time slices are minimal, and the
complex carrier is not added by hand — it is the phase-space state of the spring oscillator.

## 3. The carrier space: which components, and why three

The relevant carrier is the **lateral triplet** — the three transverse carrier polarizations of Paper
III. Promoting each to its phase-space amplitude gives

```
    ψ_n = (ψ^x_n, ψ^y_n, ψ^z_n) ∈ ℂ³.
```

(The trace `U(1)` and the photon use a rank-1 reduction of this same carrier; see
`field_tensors.md` §2. The linear envelope is a spin-1 vector — `[[project_spin_half_is_soliton_layer]]`
— so `ℂ³` is the natural carrier space, not an imposed `SU(3)` representation.)

## 4. Real degeneracy gives `O(3)`; complexification gives `U(3)`

This is the precise chain, and the step that the first attempt skipped.

- **Real triplet, exact degeneracy** `ω_x = ω_y = ω_z`. The real quadratic energy
  `½(|δR|² + |δṘ|²/ω²)` is invariant under real orthogonal mixing of the three components. The frame
  freedom of the *real* degenerate triplet is therefore only

```
    O(3).
```

- **Complex phase-space carrier.** Once each real oscillator is promoted to `ψ^a = δR^a + iδṘ^a/ω`, the
  triplet lives in `ℂ³`, and the transformations preserving the complex norm `|ψ|²` (the scaled energy)
  are

```
    ψ → U ψ,    U†U = 𝟙,    U ∈ U(3).
```

  The complexification is exactly what enlarges the real `O(3)` frame freedom to the complex `U(3)`. The
  extra generators are the phases (`O(3) → U(3)` adds the `i`-times-symmetric part).

So the "promote position+velocity to a complex amplitude" step of §2 is not cosmetic: it is the reason
the carrier frame group is `U(3)` and not `O(3)`. Without it there is no complex gauge structure at all.

## 5. Vector phase vs frame transport: `U(1)` and `U(3)` are different objects

A standing trap (critique #4): **do not identify the carrier vector `ψ ∈ ℂ³` with the connection matrix
in `u(3)`.** They are different.

- **Single carrier vector.** Normalizing one carrier state `|u(n)⟩ = ψ_n/|ψ_n| ∈ ℂ³` and transporting
  it node-to-node tracks one *overall phase*. The node-to-node overlap `⟨u(n)|u(n+ê_μ)⟩` is a scalar;
  its phase is a **`U(1)` Berry connection**. This is the trace/EM sector.

- **Local three-frame.** Transporting a full orthonormal frame `{u_1, u_2, u_3}` spanning the degenerate
  `ℂ³` subspace, the overlap is a `3×3` matrix `M_{ij}(n) = ⟨u_i(n)|u_j(n+ê_μ)⟩` whose unitary part lies
  in `U(3)`. The **connection is `u(3)`-valued**; the carrier remains a vector. This is the object that
  splits into the trace `U(1)` and the traceless `SU(3)` (`field_tensors.md`).

The distinction matters: the EM `U(1)` is the overall phase of one carrier vector; the colour `SU(3)` is
the relative orientation of a transported three-frame. One object (the `U(3)` connection) carries both,
but the carrier vector is never the connection.

## 6. Dynamic phase is raw material, not yet curvature

A caveat that keeps the paper honest. The phase-space amplitude has a phase, but for a free mode that
phase is just the clock,

```
    ψ(t) ~ e^{−iωt}    ⇒    uniform phase advance    ⇒    plaquette holonomy = 1,    F = 0.
```

The clock phase alone is **not** a gauge field. A nonzero field strength requires the *relative* carrier
phase/frame between neighbours to be **path-dependent** — a nonzero plaquette holonomy
`ψ → ψ e^{iθ_loop}`, `θ_loop ≠ 0`. The construction of §§2–5 supplies the *carrier* (a genuine `U(3)`
connection over `(x,t)`); whether its curvature is nonzero is a statement about the carrier frame
actually rotating across `(x,t)`, which it does precisely because the complex link is the time-rotation.

This is why the field strength lives over physical `(x,t)` and not the BZ: the `k`-space frame is real
(no rotation, `F≡0`); the `(x,t)` frame rotates (the time-link phase), so its plaquette holonomy is
nontrivial.

## 7. Status

- `closed` (construction): `ψ^a = δR^a + iδṘ^a/ω` from phase space; `i` = time-link rotation; real
  degeneracy `O(3)` ⇒ complex carrier `U(3)`; vector-`U(1)` vs frame-`U(3)` kept distinct.
- `closed` (guard, §1a): the rest length sets the transverse stiffness `(1−α)` of `K_μ` and hence whether
  `u(k)` mixes in `k`, but real-symmetric `D(k)` forces the `U(1)` BZ curvature to zero ∀α regardless —
  `k`-mixing is necessary, not sufficient; the `U(1)`/WZ curvature is the complex `(x,t)` carrier, not the
  real `k`-space bundle. Raw Bloch phase = one `U(1)` (four gradients, one circle), not `U(1)⁴`.
- `in-ansatz`: the degenerate triplet (full `U(3)` frame freedom) is the `α→1` limit; the linear limit
  `α→0` decouples the triplet to `U(1)³` and switches the traceless sector off — one monotonic
  activation in `α`, no interior sweet spot (`field_tensors.md` §4).
- `open`: band-isolation of the three-frame off carrier band-crossings (feeds `link_holonomy.md` §2).