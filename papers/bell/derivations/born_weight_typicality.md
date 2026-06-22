# Bell Bridge — The Born weight from τ-typicality (the last debt)

Working note. **This does NOT derive Born ex nihilo** — no approach does; every
one reduces Born to *some* postulate (quantum equilibrium, rationality,
envariance, non-contextuality). This note reduces it to the **same typicality
postulate as classical statistical mechanics** (Liouville / microcanonical / equal
a priori measure) and shows the Born *exponent* is geometric. Net claim: quantum
probability is stat-mech typicality on the carrier amplitude — not a new kind of
probability. Numerics verified. Honest debt; not asserted in the paper.

## 0. Honest framing

Born from nothing is unsolved. The honest target is to reduce it to a *standard,
non-quantum* assumption. This note reduces the Born weight to Liouville typicality
(the foundational postulate of all statistical mechanics), leaving as the only
residual the **universal** ergodic-relaxation problem — the same debt that
justifies the microcanonical ensemble classically, not a new quantum mystery.

## 1. The reduction target (from B4, D3)

- Single arm (B4): `P(click) = Born` ⟺ the effective threshold `τ̃` is uniform on
  `[0,1]`.
- Joint (D3): the joint weight `|⟨a|⟨b|Φ⁺⟩|²` ⟺ the joint two-time
  `(τ_A,τ_B)`-measure is Born.

So the task is: derive the uniformity of `τ̃` (single arm) / the `|amp|²` measure
(joint) from a typicality principle.

## 2. The complex carrier amplitude is a 2D symplectic phase space

The carrier amplitude is **complex**, `ψ = q + i p` (the `U(1)`-from-time
structure, Paper III). The two real quadratures `(q,p)` are a 2D **symplectic
phase space** with Liouville measure `dq dp`. The absorbing electron borrows
energy from the *local vacuum amplitude*; the effective threshold is the vacuum
**intensity**

    τ̃  ∝  I_vac = q² + p² = |ψ_vac|²  (the phase-space "energy").

## 3. Liouville typicality ⇒ uniform τ̃ ⇒ Born

The natural phase-space measure is Liouville `dq dp` (area). On the energetically
accessible region — intensity `≤ I_max`, a **disk** `q²+p² ≤ I_max` — the
microcanonical (equal a priori) measure is **uniform on the disk**. Then the
intensity is uniform:

    P(I_vac ≤ u) = area(q²+p² ≤ u) / area(disk) = u    ⇒   I_vac ~ Uniform[0, I_max].

(The linear-in-`u` law is the 2D area Jacobian `r dr → r²`.) Hence
`τ̃ ~ Uniform[0,1]`, and with the analyzer-projected signal intensity
`I_sig = cos²(θ−a)` (Malus, §projection),

    P(click) = P(I_sig ≥ I_vac) = cos²(θ−a) = Born.    ✓ (verified numerically)

**The Born exponent `2` is the symplectic dimension** of the complex amplitude:
`|amp|²` is phase-space area = energy, and the Malus `cos²` is the *projected area
ratio*. The quadratic is geometry, not an assumption.

## 4. Why this is non-circular

Uniform-on-phase-space (Liouville) is **not** Born-on-outcomes. The Born `|amp|²`
*emerges* as the area ratio between the equal-a-priori phase-space measure and the
projection. This is exactly the Bohmian quantum-equilibrium mechanism (equal a
priori measure in configuration/phase space → `|ψ|²` on outcomes) and the
classical stat-mech typicality mechanism. The *input* is the equal-a-priori
(microcanonical) measure; Born is the *output*.

## 5. Joint two-time case (the D3 residual, same postulate)

The joint weight `|⟨a|⟨b|Φ⁺⟩|²` is the squared modulus of the consistent-history
amplitude — again a phase-space area, under the **same** Liouville typicality
applied to the joint two-time amplitude (the two detectors' vacuum measures,
correlated by the two-time self-consistency at `S`, D3 §3). So the joint Born
weight needs no assumption beyond the single-arm one plus the already-established
two-time structure: same exponent (symplectic area), same measure (microcanonical).

## 6. The residual is universal, not quantum

The one remaining assumption is that the substrate vacuum **samples / relaxes to**
the microcanonical (Liouville) measure — ergodicity / mixing of the substrate
dynamics. This is *identical* to the debt that justifies the microcanonical
ensemble in classical statistical mechanics (the ergodic / typicality problem); it
is **not** a new quantum mystery. A substrate `H`-theorem (mixing of the vacuum
amplitude toward Liouville) would close it. That is the honest endpoint.

## Result

The Born weight reduces to **Liouville/microcanonical typicality on the complex
carrier amplitude**: the exponent `2` is the symplectic dimension (`|amp|²` =
phase-space area), the `cos²` is the projected area ratio, and the measure is the
equal-a-priori postulate. Quantum probability is thereby the **same** as thermal
probability — derived, not fundamental — with the only residual the **universal**
ergodic-relaxation problem, shared with classical statistical mechanics. This
closes the Bell-bridge probabilistic chain modulo that universal debt.

## References (wire into `references.bib` if this graduates)

- Dürr, Goldstein & Zanghì 1992 — quantum equilibrium / typicality (equal a priori
  configuration measure → `|ψ|²`).
- Liouville / microcanonical typicality and the ergodic problem (any
  statistical-mechanics text).
- B4 `threshold_detection.md`; D3 `two_time_bvp_tsirelson.md`.
