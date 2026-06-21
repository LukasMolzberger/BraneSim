# Gauge Bridge — Status

## HAVE

- Berry connection derivation (rank-1 Abelian and non-Abelian Wilczek-Zee); gauge transformation law a_μ → a_μ + ∂_μχ.
- Matter coupling via envelope equation: covariant derivative (∂_μ − ia_μ)Ψ enters the narrowband envelope equation.
- Brillouin-zone link-variable method for lattice holonomy computation.
- Confirmation that BZ Berry curvature vanishes for real symmetric D(k) (gauge structure lives in (x,t) not BZ); k-space connection ≡ 0 ∀α proven.
- U(1)/SU(3) sector split: a_μ = trace part (EM) + traceless Gell-Mann part (color); gauge-invariant sector holonomies Φ_U(1), Φ_SU(3) defined.
- Holonomy ratio R(α) = [(3−2α)/α] · (C/g(k̂)) derived; falsifiable prediction R(0.5)/R(0.2) = 0.3077.
- See: alpha_holonomy_estimator.md.

## MISSING

- Effective action for a_μ fluctuations: integration over fast substrate modes to produce S_eff[a_μ].
- Maxwell dispersion from substrate spectrum: derivation that free gauge fluctuations propagate as ω = c|k|.
- Derivation that photon propagation emerges (kinematics → dynamics gap).

## Open derivations

*Relocated from the former central `OPEN_PROBLEMS.md` (group D, gauge sector —
the EM/`U(1)` items D2, D3, D5). IDs retained so existing cross-references
(`OPEN_PROBLEMS.md D2`, etc.) resolve here. The `U(1)`↔`SU(3)` binding items D4
and D6 live in the Color bridge.*

### D2. Prestress α from the EM/colour coupling ratio — `resolved (linear: undetermined)` · `open (nonlinear)`
**Statement.** Can the empirical EM-vs-strong coupling ratio fix the prestress α?
The linear *spectrum* cannot: the traceless/`SU(3)` content `∝ α` and the trace/
`U(1)` content `∝ (1−2α/3)` stay comparable (≈13% apart at α=0.2) and both vanish
together as α→0, so the spectrum produces **no** EM/colour hierarchy — the famous
`α_EM/α_s` must live in the *connection/curvature normalization*, not eigenvalue
magnitudes. Compounding this, the **k-space Berry/WZ curvature is identically zero
∀α** (real-symmetric `D(k)`; k-space connection ≡ 0 ∀α, BACKBONE #16), so any coupling ratio is a
*fibre-internal* `(x,t)` holonomy object.
**Why it matters.** A derived α would convert the project's one free linear
parameter into a prediction. A clean *negative* result (α undetermined at linear
order) is itself valuable: it localizes the EM/colour hierarchy in the nonlinear
normalization and rules out spectrum-ratio "derivations" of `1/137`.
**Lorentz-bound corollary (already usable).** Mutual photon/gluon Lorentz
invariance (no measured vacuum birefringence between EM and colour) does *not*
bound α small; it forces physical carriers onto the `[111]` isotropic locus where
all three lateral eigenvalues are exactly degenerate ∀α (`g([111])=0`). This
*frees* α to sit at O(0.1–0.5) and is the cleanest current "SM fact → lattice"
import.
**Resolution (linear order, 2026-06-04 — derivation; α NOT determined, by design).**
The fibre-internal holonomy ratio is `R(α) = (C/g(k̂))·(3−2α)/(√3 α)`, where
`g(k̂)=√(Σ_a(h_a−H/3)²)/H` is the directional anisotropy and `C = N_U(1)/N_SU(3)`
is the **relative gauge-kinetic normalization** = `tr(P_U(1)𝒢)/tr(P_SU(3)𝒢)` set by
the fibre metric `𝒢_ij=⟨∂_iΨ|∂_jΨ⟩`. At linear order `𝒢 ∝ 𝟙` (envelope norm is
flat ℓ² on ℂ³), so `C` reduces to pure dimension-counting (`1/√2`) with **no
α-dependence** — the EM/colour hierarchy is therefore *absent* from the linear
theory and the naive inversion gives `α≈1.2 > 1` (out of range, a self-flag that
the linear identification is invalid). **α is normalization-undetermined at linear
order.** The missing ingredient is the *sector-dependent* fibre metric supplied by
the geometric quartic (backbone #17), which couples the trace/dilatational
direction to the `X⁴` gravity channel differently from the traceless/shear
directions — a nonlinear computation. **Sharpest new gap surfaced:** the triplet is
coherent only for `k₀∥[111]`, but `g([111])=0` kills the colour content there
(`R→∞`); off-[111] the branches dephase. There is **no linear direction that is
both coherent and colour-active** — color curvature is intrinsically a
soliton-layer (nonlinear) object. **Falsifiable, normalization-independent
spin-off (decisive):** the α-scaling is `R(α₂)/R(α₁) = [(3−2α₂)/α₂]/[(3−2α₁)/α₁]`;
for `(α₁,α₂)=(0.2,0.5)` this is **0.3077** — measurable by extending
`branesim/diagnostics/berry_holonomy.py` to transport a band-isolated ℂ³ envelope
around a fixed `(x,t)` loop off-[111]; a >10% deviation falsifies the
spectral-susceptibility factorization. Full derivation:
`paper/derivations/alpha_holonomy_estimator.md`.

### D3. Dynamical Maxwell: effective action for the EM four-potential — `open`
**Statement.** Electromagnetism is the `U(1)`-trace Berry connection of the carrier — the
**electromagnetic four-potential** `A_μ = i⟨u|∂_μu⟩` — with **field tensor**
`F_μν = ∂_μA_ν − ∂_νA_μ` (`E_i = F_{0i}`, the time–space part; `B_i = ½ε_{ijk}F_{jk}`, the
space–space part). The **kinematics is direct** (no dynamics needed): (a) gauge invariance
`A_μ → A_μ + ∂_μχ` is the unobservable choice of carrier phase reference; (b) the **homogeneous**
Maxwell equations are the Bianchi identity `dF = 0`, automatic because `F = dA`. What is **not**
derived is the **dynamics** — the **inhomogeneous** equations `d⋆F = J` — which hold only if
integrating out the substrate around the carrier band yields the Maxwell action
`S_eff[A] = −¼ ∫ F∧⋆F`, requiring: (i) the right **form** (`∝ F²`, not another functional of `F`);
(ii) **masslessness** (no `m²A²`; the `U(1)` must be an exact flat direction); (iii) the emergent
**Hodge `⋆`** = isotropic light-cone (ties A4a / `[111]`); (iv) the **coupling** normalization
(undetermined at linear order, D2 — lives in the nonlinear fibre metric). Also open: what
dynamically forces a **nonzero, propagating** `F` (the phase must be non-integrable — the carrier
state must rotate in `ℂ^N`, `N≥2`, not merely carry an overall phase; a pure gradient is pure
gauge, `F=0`).
**Why it matters.** Closes the EM sector from "gauge kinematics" (done) to "Maxwell dynamics"
(the field equations + the photon). Until derived, the paper claims only the kinematic/gauge
structure, not full Maxwell dynamics (§2 non-claims). Recurs across §5.6, D2, and the complex-`U(1)`
discussion. The complex structure / the `i` is the carrier's rotation along the timelike worldtube
axis (`[[project_complex_u1_from_time]]`), so `A_μ` is intrinsically a 4-covector on `(x,t)`.
**Candidate approach.** Effective-field-theory reduction of the substrate action about the carrier
band → `S_eff[A]`; verify the `F²` form + masslessness; read the coupling from the fibre metric
(D2). Owner: physics-derivation + berry-validator. See `paper/05` §5.2.

### D5. Massive-vs-massless: the substrate must not Meissner-screen its own `U(1)` photon — `open` (2026-06-06)
**Ontology note (owner, 2026-06-07) — read first.** This entry tracks an *empirical* worry
about field **range** (does the trace field fall off as `1/r` or exponentially), **not** a
claim that mass originates from a Higgs mechanism. In this theory **mass is the energy of a
self-confined wave** — an excitation trapped in an infinite loop back on itself
(geometric-self-guidance / trapped-mode picture, paper §6.6) — *never* a property handed to a
gauge field by a symmetry-breaking condensate. The Anderson–Higgs/GL language below is the
**standard-physics expectation the substrate must NOT inherit**, recorded so the risk is
explicit; the resolution routes are exactly the ones that avoid it.
**Statement.** Standard Abelian-Higgs / Ginzburg–Landau would *naively* predict: a *condensed*
medium with a spontaneously fixed order-parameter amplitude `|Ψ|>0` makes its `U(1)` gauge
field **short-range / Meissner-screened** (penetration depth `λ<∞`). But the active `U(1)`
carrier-phase vortex (D4, `EXPERIMENT.md`) is meant to model the **electron/EM**, whose photon
is **free and long-range** (Coulomb `1/r`). A naïvely condensed prestressed substrate therefore
risks producing the *opposite* regime — a short-range, **self-trapped, non-radiative** trace
field where a free photon is required. The screening that GL hands you "for free" is exactly
wrong for the electron's long-range field, and the job here is to confirm the substrate does
**not** sit in that regime for the trace-`U(1)`.
**Why it matters.** This is load-bearing for the experiment being built now. If the
in-substrate vortex screens its own emergent gauge field, the object is not EM-like, and the
`F_μν → E,B` gauge-layer measurement would show exponential (not `1/r`) fall-off. It is the
condensed-matter shadow of the masslessness requirement already noted abstractly in **D3(ii)**
(`no m²A²`; the `U(1)` must be an exact flat direction) — D5 names the concrete mechanism that
would *violate* it and the regime that controls it.
**Candidate approach / resolution sketch.** (i) **Unbroken-combination route** (electroweak
analogy): the physical `U(1)_em` is the *surviving, unbroken* generator; the substrate may
break some `U(3)` directions while the trace-`U(1)_em` stays massless. (ii) **Not-a-condensate
route** (preferred, consistent with `[[project_complex_u1_from_time]]`): the carrier phase is a
*rotation along the timelike worldtube axis*, not a thermodynamic condensate amplitude that
Higgses the field — i.e. the substrate is **not in the broken phase** for the trace-`U(1)`, and
the vortex is a **semilocal** excitation on a simply-connected vacuum (the `β<1` regime of
D4(i), Vachaspati–Achucarro), where the gauge field is only *partially* Higgsed. Semilocal
vortices live precisely in the window where masslessness can survive. **Decisive diagnostic
(cheap, add to the active run):** extract the emergent gauge-boson screening length from the
`F_μν`/supercurrent profile (penetration depth `λ` vs core `ξ`) — long-range/unscreened
(`λ → ∞`, power-law `B`/current tail) confirms massless EM; exponential screening falsifies the
EM identification and flags a `W/Z`-class object instead. **Do not import `e* = 2e`** (Cooper
pairing) — the carrier is single-charge.
**Note (2026-06-07).** Kept out of the paper — this screening worry is an internal tracker
item only. (A `§6.5` "companion falsifiable check" paragraph was briefly added and then
**removed**: the owner does not want an Anderson–Higgs / screening discussion in the
manuscript. Mass in this theory is a self-confined wave, not a Higgs effect; see
[[feedback_mass_is_self_confinement]].)
**Status.** open. Owner: physics-derivation (which `U(1)` survives massless) + the active
build (screening-length diagnostic). Tightens D3(ii); shares the semilocal-binding regime with
D4(i) (Color bridge). Surfaced from R. Behiel's *Superconductivity and the Higgs Field* (Anderson–Higgs).