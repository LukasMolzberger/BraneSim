# Bell Bridge — The pilot-wave lineage: de Broglie's double solution, not Bohm

Working note. **This does NOT add new physics** — it places the substrate's
already-committed ontology (real field, extended soliton, local 4D action) inside
the pilot-wave family precisely, and shows that the *guidance relation* is a
**theorem** of the carrier structure, not an extra postulate. Net claim: the
substrate is the natural completion of de Broglie's *théorie de la double
solution* — particle = soliton-singularity of one real wave in physical space,
guidance derived — and is structurally **incompatible** with Bohm's
configuration-space pilot wave, which is excluded on the same locality ground as
every nonlocal completion (§2.5 `02_bell_constraint.tex`). The single-particle
guidance theorem is rigorous; the many-body extension is deferred to the
two-time worldtube (D3), and Born to typicality — both already open debts, not new
ones.

## 0. Honest framing

"de Broglie–Bohm" names **two different mathematical programs**, and the substrate
is compatible with one and hostile to the other. Conflating them is the standard
error. This note separates them, shows which relations transplant to the substrate
as theorems, and names exactly where each program's distinctive content lives.
Nothing here is asserted in the paper beyond the placement line already in
`status.md`.

## 1. The two programs are not the same theory

| | **Bohm pilot wave** (1952) | **de Broglie double solution** (1925–27, 1950s) |
|---|---|---|
| Ontology | point particle **+** separate guiding wave | particle is a **singularity/soliton of one real wave** |
| Wave lives on | configuration space `ℝ^{3N}` | physical space `ℝ³` (here: the 4D world-volume) |
| Guidance | **postulated** `v_k = ∇_k S / m_k` | a **theorem**: phase-locking of core to envelope |
| Locality (N≥2) | **nonlocal** (config-space `Q` couples all particles) | local by construction |
| Dualism | yes (two beables) | no (monistic — one field) |

Bohm's *distinctive* content — the part that makes it "Bohmian mechanics" rather
than "Madelung hydrodynamics" — is the `ℝ^{3N}` pilot wave and the resulting
nonlocality. That is exactly the content the substrate cannot host: there is no
`ℝ^{3N}` on a nearest-neighbour lattice, and nonlocal guidance is excluded by (L)
(§2.5). The substrate is therefore **not** a Bohmian theory.

de Broglie's double solution, by contrast, is almost a verbatim description of the
substrate: a real physical wave whose localized singularity *is* the particle,
with the statistical `ψ` a derived envelope. What de Broglie could not do — make
the double solution work for many-body entanglement in `ℝ³` — is what the 4D
retrocausal worldtube supplies (§5).

## 2. The substrate carrier admits a Madelung form

From Paper III, a narrowband excitation is a **real** substrate field with a fast
oscillation at `ω₀`; the slowly-varying complex envelope (here the *matter
excitation's* envelope — the de Broglie matter wave — not the EM/`U(1)` photon
carrier of the gauge bridge)

    ψ = ξ + (i/ω₀) ξ̇            (the `U(1)`-from-time amplitude, `ψ = q + i p`)

obeys an effective Schrödinger equation in the linear/narrowband regime

    i ℏ_eff ∂_t ψ = H_eff ψ,     H_eff = −(ℏ_eff²/2 m_eff) ∇² + V,

with `m_eff`, `ℏ_eff` fixed by the dispersion and carrier scale `ω₀` (Papers
I/III); the committed envelope form `i ∂_t Ψ = H_eff Ψ` is its `ℏ_eff → 1` case. Write the **Madelung
polar form**

    ψ = R e^{i S / ℏ_eff},       ρ ≡ R² = |ψ|² = q² + p²   (the §-born intensity).

Splitting `i ℏ_eff ∂_t ψ = H_eff ψ` into imaginary and real parts gives, with no
extra assumption:

    (continuity)   ∂_t ρ + ∇·(ρ v) = 0,            v ≡ ∇S / m_eff,
    (quantum HJ)   ∂_t S + (∇S)²/2m_eff + V + Q = 0,  Q ≡ −(ℏ_eff²/2 m_eff) ∇²R/R.

`v = ∇S/m_eff` is the local carrier momentum flux; `Q` is the envelope's
dispersive self-stress. **Both are ordinary real-space fields**, evaluated at a
point — for one excitation, configuration space = physical space, and the Bohmian
and double-solution descriptions coincide.

## 3. The guidance relation is a THEOREM, not a postulate (single particle)

This is the crux, and the point at which the substrate sides with de Broglie. A
soliton is a localized field configuration of finite width `ℓ` centred at

    X(t) = ∫ x ρ(x,t) d³x / ∫ ρ d³x.

Differentiate and use the continuity equation (integrate by parts, `ρ → 0` at
infinity):

    dX/dt = ∫ x ∂_t ρ / ∫ ρ = −∫ x ∇·(ρv) / ∫ ρ = ∫ ρ v / ∫ ρ = ⟨v⟩.

In the **narrow-soliton limit** `ℓ ≪ |∇S|/|∇²S|` (envelope width small against the
scale on which the phase gradient varies) — which the soliton's *nonlinear
self-confinement* (Papers III/IV) supplies; the centroid identity above uses only
linear continuity, the confinement is what holds `ℓ` small —

    dX/dt = ⟨∇S/m_eff⟩ → ∇S(X,t) / m_eff.        ★ (de Broglie guidance)

So the soliton centroid follows the local carrier phase gradient **as an
Ehrenfest/centroid theorem of the substrate field equations** — not as an added
beable law. This is precisely de Broglie's claim for the double solution
("the guidance formula is a consequence, the `u`-wave and `ψ`-wave share a phase"):
here the soliton core and its envelope are the *same* field, so the core is
necessarily carried along `∇S`. Phase-locking is automatic because there is only
one field.

Contrast Bohm, where `★` is *postulated* for an independent point particle riding
a separate wave. The substrate needs no such postulate and no second beable.

## 4. Why Bohm's program is excluded (locality), not merely declined

For `N ≥ 2` the Bohmian pilot wave is `Ψ(x_1,…,x_N,t)` on `ℝ^{3N}`, and

    v_k = ∇_k S(x_1,…,x_N,t) / m_k

makes particle `k`'s velocity depend **instantaneously on every other particle's
position**. This is irreducible nonlocality (it is *how* Bohm escapes Bell — as an
overtly nonlocal hidden-variable theory) and it requires a preferred foliation to
define the simultaneity in "instantaneously."

The substrate forbids both:

- **No `ℝ^{3N}`.** The only field is the one real world-volume field in 4D. There
  is nowhere for a configuration-space wave to live; it is not an object the
  ontology contains.
- **No nonlocal guidance, no preferred foliation.** (L) is a structural fact of
  the nearest-neighbour lattice (§2.2). Bohm's nonlocal guidance and preferred
  slicing are excluded *before* one even reaches the config-space problem — on the
  same single ground that excludes every nonlocal completion (§2.5).

So the exclusion of Bohm is not a stylistic preference; it is forced by the same
locality commitment that the whole Bell argument rests on.

## 5. The many-body gap: filled by the worldtube, not by config space

de Broglie abandoned the double solution largely because he could not represent
entangled `N`-body correlations with `N` solitons in a single `ℝ³` — the
correlations seemed to demand the `ℝ^{3N}` wave that Bohm then embraced. The
substrate takes the third road:

- Entangled particles are **one branching 4D worldtube** (HAVE; §3
  `03_worldtube.tex`), not two independent `ℝ³` solitons needing a config-space
  wave to correlate them.
- The correlation is carried by **worldtube continuity + a time-symmetric (future)
  boundary condition** — the two-time BVP — not by instantaneous config-space
  guidance.
- This relaxes (MI) via **retrocausality**, the substrate's chosen Bell escape,
  *instead of* relaxing (L) via config-space nonlocality, Bohm's escape.

Thus the many-body content that forced Bohm into configuration space is, in the
substrate, the content of D3 (`two_time_bvp_tsirelson.md`): the joint two-time
history amplitude `⟨a|⟨b|Φ⁺⟩` and its Born weight. **This note does not close that;
it identifies it as the precise locus where the double solution is completed by the
4D worldtube rather than by `ℝ^{3N}`.** The single-particle guidance theorem (§3)
is the `N=1` shadow of the same structure.

## 6. Born: two complementary typicality routes, one postulate

de Broglie–Bohm recover Born via *quantum equilibrium* — an equal-a-priori measure
in configuration space, equivariant under the guidance flow, yields `|ψ|²` on
outcomes (Dürr–Goldstein–Zanghì). The substrate carries **two** routes to Born,
at *different loci*, and they must agree:

- **Preparation / equilibrium side (this note).** The guidance theorem `★`
  supplies the equivariant flow — the `ρ = |ψ|²` continuity equation of §2 *is* the
  equivariance statement — so an equal-a-priori ensemble of soliton configurations
  stays `|ψ|²`-distributed (the DGZ mechanism).
- **Detection side (`born_weight_typicality.md`).** Equal-a-priori (Liouville)
  measure on the *detector's* vacuum amplitude `(q,p)` ⇒ uniform intensity ⇒ Born
  click probability, with the exponent `2` = the symplectic dimension.

These are **not** the same measure — one is on the guided particle's
configuration, the other on the absorber's vacuum phase — but they rest on the
**same equal-a-priori postulate** and must yield the same `|ψ|²`. Their
consistency (preparation statistics = detection statistics) is the substrate's
analog of the standard requirement that the quantum-equilibrium distribution match
the measurement rule; establishing it is *part of* the Born debt, not a separate
one. Either way the only genuinely remaining input is the universal
ergodic-relaxation problem, shared with classical stat-mech (B-group) — not a new
quantum postulate.

## 7. The quantum potential, read on the substrate

`Q = −(ℏ_eff²/2m_eff) ∇²R/R` is, in substrate terms, the **dispersive
self-reaction of the soliton envelope** — the term by which the spreading wave
shapes its own centroid trajectory. For a single excitation it is a local
real-space field (it is only in Bohm's `N≥2` config space that `Q` becomes the
nonlocal coupling). It is what makes `★` reproduce interference/tunnelling without
any nonlocal input, and it is the natural home for the soliton's internal-structure
back-reaction (Papers III/IV). No new postulate: `Q` is already present in the
real part of `i ℏ_eff ∂_t ψ = H_eff ψ`.

## What this note establishes vs leaves open

**Established (rigorous, modulo the narrowband regime of Paper III):**
- the carrier admits the Madelung continuity + quantum-HJ split (§2);
- the de Broglie guidance relation `★` is a **centroid theorem** of the substrate
  field equations, not a postulate (§3);
- Bohm's config-space pilot wave is **excluded** by (L)/no-`ℝ^{3N}` (§4), placing
  the substrate squarely in the double-solution lineage;
- guidance supplies the equivariant flow for the typicality/Born argument (§6).

**Open (already-named debts, not new):**
- the many-body / entangled guidance = the two-time worldtube amplitude and its
  Born weight (D3, `two_time_bvp_tsirelson.md`);
- the Born measure itself (B-group, `born_weight_typicality.md`);
- validity of `i ℏ_eff ∂_t ψ = H_eff ψ` beyond the linear narrowband regime
  (Papers I/III); `★` inherits that regime.

## Result

The substrate is the **completion of de Broglie's double solution**, not a Bohmian
theory: the particle is a soliton-singularity of one real wave in physical space,
the guidance relation `dX/dt = ∇S/m_eff` is a theorem of the carrier's continuity
equation, and Bohm's configuration-space pilot wave is excluded by the same
locality commitment that drives the whole Bell argument. The many-body content
that forced Bohm into `ℝ^{3N}` is supplied here by the 4D retrocausal worldtube
(D3), and Born by typicality with guidance furnishing equivariance. The placement
is thus not interpretive taste but a structural consequence of the substrate's
ontology.

## References (wire into `references.bib` if this graduates)

- de Broglie 1927; de Broglie 1956 (*Une tentative d'interprétation causale et non
  linéaire de la mécanique ondulatoire* — the double solution) — particle as
  wave-singularity, guidance as theorem.
- Bohm 1952 (I & II) — the configuration-space pilot wave (the program excluded
  here).
- Dürr, Goldstein & Zanghì 1992 — quantum equilibrium / typicality (§6).
- Couder & Fort 2006; Bush 2015 — hydrodynamic pilot-wave / wave-memory
  realization of the double solution (local, real-space; guidance dynamics
  template). Single-particle/local analogy only — it cannot itself violate Bell;
  the substrate's Bell content is the worldtube (§5), not the droplet-like guidance.
- Sutherland 2017 — time-symmetric (retrocausal) pilot wave; the bridge between
  the double solution and the worldtube (§5).
- D3 `two_time_bvp_tsirelson.md`; Born `born_weight_typicality.md`;
  §2.5–2.6 `02_bell_constraint.tex`.