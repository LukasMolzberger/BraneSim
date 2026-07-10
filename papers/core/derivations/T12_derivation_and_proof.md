# T12 — Derivation and proof of the temporal-prestress continuum/kinetic limit

**Status: T12 is proved.** Adopting the Layer-0 temporal rest length `r_4 = α_t a`
(`0 < α_t < 1`), the temporal link term is shown to be a **genuine prestressed
spring** whose continuum limit **is** the kinetic energy, with:

- an **emergent inertia** `m = κ_t a²` (no separately-postulated mass);
- a **`(1−α_t)` anisotropy** that is the fingerprint distinguishing genuine
  temporal prestress from the old "kinetic-sign artifact" (`r_4 = 0`);
- a finite signed prestress `ρ_4 = +κ_t a(1−α_t)`, symmetric in structure with the
  spatial `ρ_i = −κ_s a(1−α_s)`;
- and the confirmation that **`η_4 = +1` makes `S = T − V`** (Hamilton's principle),
  so the static 4D-block stationarity is Störmer–Verlet time evolution carrying
  waves at the T2 speed `c_T`.

This resolves Decision D6 / the "T12 tension." Reproduced by
`derivations/t12_temporal_prestress.py`.

---

## 1. What T12 requires

OPEN_TASKS: *"Adopt `r_4 = α_t a` (`0<α_t<1`) as the new Layer-0 form; derive its
continuum/kinetic limit and confirm `η_4 = +1` reads as genuine prestress."*
Decision D6: *"How does `η_4 = +1` arise from temporal prestress, and what is the
continuum/kinetic limit of `r_4 = α_t a`?"*

The **tension** being resolved: in the old formulation `r_t ≈ 0`, the temporal term
`½κ_t L_4²` is a plain isotropic kinetic energy, so `η_4 = +1` looks like a mere
*kinetic sign*, not a prestress — indistinguishable from "time is just kinetic."
T12 must show the new `r_4 = α_t a` makes the temporal prestress **physical** and
that its continuum limit yields the kinetic term.

---

## 2. The temporal link expansion (derivation)

Temporal link vector about the vacuum `R̄_n = a n`, with fluctuation `u`:

```
Q_{n4} = a ê_4 + Δ_4 u,   Δ_4 u = u_{n+ê_4} − u_n,
L_{n4} = |Q_{n4}| = a + Δ_4 u^4 + (1/2a) Σ_i (Δ_4 u^i)² + O(u³).
```

Insert into `S_time = ½κ_t Σ_n (L_{n4} − α_t a)²` and keep the quadratic part
(`a − r_4 = a(1−α_t)`):

```
S_time^(2) = ½ κ_t Σ_n [ (Δ_4 u^4)² + (1−α_t) Σ_i (Δ_4 u^i)² ].
```

This is exactly `C_4 = η_4κ_t M_4` with `M_4 = (1−α_t)I + α_t ê_4ê_4ᵀ` (Layer 1,
used in T2). Verified to machine precision by the exact single-link second
derivative:

```
alpha_t   longitudinal u^4    transverse u^i
0.0       κ_t                 κ_t              (isotropic)
0.5       κ_t                 κ_t·0.5
0.9       κ_t                 κ_t·0.1
ratio(transverse/longitudinal) = (1−α_t) exactly.
```

---

## 3. The continuum limit *is* the kinetic term (proved)

With the temporal lattice spacing fixed (`a=1`, fundamental), the continuum/kinetic
limit is the **long-wavelength** limit `ω→0` (many steps per oscillation), *not*
`a→0`. The discrete temporal kinetic operator obeys

```
2 κ_t (1 − cos ω)  ⟶  κ_t ω²      (2(1−cos ω)/ω² → 1,  verified: 0.948, 0.987, 0.997, 0.999).
```

So `Δ_4 u → a ∂_t u` and

```
T = ½ κ_t a² [ (∂_t u^4)² + (1−α_t) Σ_i (∂_t u^i)² ],
```

a genuine kinetic energy with **emergent inertia**

```
m_long = κ_t a²   (u^4 mode),      m_tran = κ_t a²(1−α_t)   (u^i modes).
```

The mass is *derived* from the temporal link stiffness `κ_t` — it is not a separate
postulate. This is the "kinetic limit" T12 asked for.

---

## 4. `(1−α_t)` anisotropy = the fingerprint of genuine prestress

This is the crux of the D6 resolution. Two cases:

- **`r_4 = 0` (old):** `S_time^(2) = ½κ_t|Δ_4 u|²` — **isotropic** (`m_long = m_tran
  = κ_t`). Indistinguishable from a plain kinetic energy `½m|∂_t u|²`. The prestress
  is invisible: the "kinetic-sign artifact."
- **`r_4 = α_t a`, `0<α_t<1` (new):** the transverse temporal inertia
  `m_tran = κ_t a²(1−α_t) = ρ_4·a` is *reduced* by the prestress, so `m_tran ≠
  m_long` — **anisotropic**. A pure inertia cannot do this; only a *spring under
  tension* has distinct longitudinal and transverse responses. So the `(1−α_t)`
  anisotropy is the direct, observable fingerprint that `η_4 = +1` is a **genuine
  prestress**, not a kinetic sign.

The transverse temporal inertia is literally the prestress: `m_tran/a =
κ_t a(1−α_t) = ρ_4`. This same `(1−α_t)` factor is what set the amplitude-mode
speed `c_4` and the `α_t`-dependence used in T2/T7 — T12 is their Layer-0 ground.

---

## 5. Prestress reading and the allowed window (Decision G)

The signed prestress `ρ_μ = η_μ κ_μ a(1−α_μ)` gives, for time,
`ρ_4 = +κ_t a(1−α_t)` — finite and positive, structurally symmetric with the
spatial `ρ_i = −κ_s a(1−α_s)`; the **opposite sign is the Lorentzian signature**
(T2). Limits:

```
α_t → 0 :  ρ_4 maximal but effect ISOTROPIC — the kinetic-sign-artifact limit (avoid).
α_t = 1−ε: genuine anisotropic prestress — the Decision G working regime.
α_t → 1 :  ρ_4 → 0, transverse temporal stiffness → 0 — over-decoupled (FORBIDDEN).
```

The temporal prestress also has a kinematic reading: the temporal link "prefers"
length `α_t a`, i.e. a preferred advance rate `|ΔR|/step = α_t` through the ambient,
held at `a` by the periodic 4-block — genuine pretension, exactly as the spatial
links are held at `a` with rest length `α_s a`.

---

## 6. `η_4 = +1` ⇒ `S = T − V` ⇒ Verlet dynamics (proved)

Because `η_4 = +1` (opposite to `η_i = −1`), the action is
`S = S_time − |S_space| = T − V`: the temporal prestress term is the kinetic `T`,
the spatial prestress terms the potential `V`. Stationarity `∂S/∂R_n = 0` of the
static 4D block is then exactly the **Störmer–Verlet stencil** — Hamilton's
principle. Numerically marching it (linearised transverse field, 1D chain):

```
c_T predicted (T2)        = 0.9129
c_T measured (Verlet)     = 0.9097          ← wave propagates at c_T
η_4 = +1 : |u| bounded (0.98–1.00)          ← genuine wave dynamics
η_4 = −1 : |u| → 2×10⁶⁸ (elliptic blow-up)  ← no dynamics (Euclidean block)
```

So `η_4 = +1` is precisely what turns the static 4D block into causal time
evolution: the kinetic term from the temporal prestress carries waves; flipping to
the all-minus (Euclidean) sign destroys propagation. This is the dynamical, EOM-side
counterpart of the T2 signature result.

---

## 7. Scope — proved vs owed

**Proved / derived (T12):**
- The temporal link term's quadratic expansion `= ½κ_t[(Δ_4u^4)²+(1−α_t)Σ_i(Δ_4u^i)²]`
  (exact, §2); its continuum limit is the kinetic energy with emergent inertia
  `κ_t a²` (§3).
- The `(1−α_t)` anisotropy distinguishes genuine prestress from the `r_4=0`
  kinetic-sign artifact (§4) — resolving the D6/T12 tension.
- `ρ_4 = κ_t a(1−α_t)` genuine prestress, symmetric with space; the Decision-G
  window `0<α_t<1` (§5).
- `η_4 = +1` makes `S = T − V`; the 4D-block stationarity is Verlet evolution
  carrying waves at `c_T`; the all-minus sign gives no dynamics (§6).

**Owed / cross-task:**
- the *nonlinear* (finite-amplitude) temporal-prestress dynamics beyond the
  quadratic/continuum limit (the full Pythagorean temporal term enters the T1
  carrier's stationarity, already used there);
- the precise inside-observer proper-time identification (ties to the T2
  dual-observer conjecture — whose clock is built from which mode);
- matching `κ_t a²` to a physical particle mass scale (needs the Layer-5 soliton
  sector).

---

## 8. Reproducibility

`derivations/t12_temporal_prestress.py` — self-contained (numpy). Prints (A,B) the
exact link inertia coefficients and the `(1−α_t)` anisotropy, (B') the
`2(1−cos ω)→ω²` kinetic limit, (C) the prestress reading and its windows, (D) the
Verlet wave-propagation demonstration and the `η_4=−1` blow-up. All numbers above
are its output.
