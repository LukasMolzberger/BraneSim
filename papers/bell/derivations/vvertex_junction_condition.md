# Bell Bridge — The junction condition at the vertex S (B5 core)

Working note. Formalizes the conserved-flux junction at the source vertex `S`,
states what it **encodes** and what structure it must **reproduce**
(`E = −cos 2(a−b)` + flat marginal), and delineates sharply what stays open.
**The junction does NOT derive Tsirelson from locality — Bell forbids that.**
Honest debt; not asserted in the paper. Core of B5; companion to
`vvertex_energy_transport.md`, `threshold_detection.md` (B4),
`pancharatnam_holonomy.md` (the joint form).

## 0. The Bell wall (stated up front)

Any junction that is **local-in-time and measurement-independent** gives
`|S| ≤ 2`. Verified: a local junction with a shared `λ` and no setting cross-talk
gives the triangular/sawtooth correlation, CHSH `|S| = 2.000` exactly (vs the QM
`2√2`). So the junction at `S` cannot be a forward emission rule; it must be a
**constraint on the full two-time worldtube boundary-value problem**, making the
shared state at `S` depend on *both* future settings (relaxing MI). This is not a
stylistic choice — it is what the entire bridge rests on, and it is supplied by
the substrate action's time-symmetry (`subsec:time-symmetry`).

## 1. Fields and currents on the worldtube

On each leg the carrier is a polarization 2-spinor field `ψ(x) ∈ ℂ²` (the
eigenbundle of Paper III) with Berry connection `a_μ = i⟨u|∂_μ u⟩`. Two
conserved Noether objects of the carrier:

- energy–momentum `T^μ_ν` (translation symmetry);
- the polarization current `J^μ_α` (`α = 1,2,3`), the Noether current of the
  `SU(2)` polarization-frame rotation symmetry — the spin angular momentum of the
  carrier.

Time-symmetric decomposition on each arm: `ψ = ψ_ret + ψ_adv` (retarded offer +
advanced confirmation, `vvertex_energy_transport.md`).

## 2. The junction conditions at S (3-valent node: pump, arm A, arm B)

Integrate the conservation laws over a small 4-ball around `S`:

- **J1 — energy–momentum conservation (phase matching).**
  `∮ T^μ_ν dΣ_μ = 0`  ⇒  `ℏω_p = ℏω_A + ℏω_B`, `k_p = k_A + k_B`.
- **J2 — polarization-invariant lock.**
  `∮ J^μ_α dΣ_μ = 0`  ⇒  the two-photon polarization state sits in the sector
  fixed by the pump's polarization angular momentum. For the rotationally
  unbiased configuration this is the `SU(2)` **singlet** `|Ψ⁻⟩` — the unique
  `U⊗U`-invariant in `ℂ²⊗ℂ²`. Rotationally covariant: no preferred polarization
  axis is introduced at `S`.
- **J3 — time-symmetric two-time closure.**
  The solution on the V is the stationary point of the time-symmetric action on
  the branching domain with **past BC at `S`** (J1+J2) and **future BCs at
  `D_A`, `D_B`** (analyzer axes `a`,`b` + the B4 thresholds). The advanced legs
  carry `a`,`b` back to `S`. Well-posedness requires a **chiral-characteristic
  BC**, not Dirichlet (two-time Dirichlet is ill-conditioned, `κ ~ 1/sin Nθ`;
  the chiral BC restores `κ = 1` — see the block-solver chirality result).

## 3. What the junction cleanly ENCODES (two consequences)

- **(a) Flat marginal / no-signalling, from J2 alone.** The singlet's one-party
  reduced state is `ρ_A = ρ_B = ½𝟙` (rotation-invariant). Hence `P(B | a)` is
  independent of `a` *regardless* of what the advanced channel carries — the
  backward leg can move joint correlations but cannot shift a one-party marginal.
  Verified: `P(B=+) = ½` for all settings. This is the marginal half of B2, and
  it follows structurally from J2.
- **(b) MI-relaxation, from J3.** The field configuration at `S` (the hidden
  variable `λ`) is fixed by the two-time BVP that contains *both* future BCs, so
  `λ = λ(a,b,τ_A,τ_B)` and `P(λ|a,b) ≠ P(λ)`. This is the precise formal content
  of "retrocausal coupling through `S`," and it is exactly the assumption Bell
  requires one to relax. It is *not* superdeterminism: `λ` is fixed by the
  uniform local action under freely posed future BCs, not by tuned past data
  (`subsec:no-superdeterminism`).

## 4. The structure the junction must realize → −cos 2(a−b)

Given J2 (singlet anti-alignment at `S`) and the B4 thresholds at both ends, the
target structure is the conditional/collapse one: A's setting+threshold fixes A's
frame; J2 ties B's forward frame to the anti-aligned partner; B's threshold fires
per Malus. Carrying this through (verified):

    P(A=+) = ½ ;  given A=+, B-frame = |a+π/2⟩ ;  P(B=+|A=+) = |⟨b|a+π/2⟩|² = sin²(a−b)
    ⇒  E(a,b) = sin²(a−b) − cos²(a−b) = −cos 2(a−b),   marginal P(B=+) = ½.

The three origins: **factor 2** = Poincaré geometry (`pancharatnam_holonomy.md`);
**minus** = J2 singlet; **Malus** = B4. CHSH `|S| = 2√2`.

**This section states the structure the junction must hit, not a proven output**
— see D3 below.

## 5. Open debts — what is NOT yet established (do not claim)

This note specifies the conditions and shows they are *consistent with* and
*encode* the QM structure. Four things remain genuinely open:

- **D1 — derive J2 from the substrate vertex nonlinearity.** *Partly closed* in
  `singlet_lock_from_vertex.md`: the singlet is the **unique** output of any
  SU(2)-polarization-invariant pair-creation vertex (`ε` is the unique invariant
  rank-2 tensor; numerically verified), and that invariance is exact on-axis in
  the long-wavelength limit (cubic-symmetry degeneracy of the transverse
  doublet), broken controllably off-axis (falsifiable singlet-fidelity
  `F = 1 − O((Δ/Γ)²)`). J2 is thus reduced from an import to "transverse-doublet
  degeneracy." **Residual:** compute the vertex tensor `M` explicitly from the
  geometric nonlinearity (confirm `M_sym = 0` on-axis) and identify the cubic
  three-wave term.
- **D2 — well-posedness + uniqueness of the two-time branching BVP.** With the
  chiral-characteristic BC (J3), prove a unique global solution exists for all
  `(a,b,τ_A,τ_B)` on the branching domain.
- **D3 — exactness / Tsirelson (B1).** Prove the BVP yields *exactly*
  `−cos 2(a−b)` (`S = 2√2`), not merely some nonlocal correlation in
  `(2, 2√2]` or an over/undershoot. §4 *assumes* the junction realizes the
  projective collapse rule; deriving that it does so exactly is the decisive step.
- **D4 — no-signalling as a theorem (B2).** Extend §3(a) from the marginal/linear
  argument to the full nonlinear two-time solution: no one-arm observable depends
  on the far setting.

## Caveat (one line)

The junction is a *constraint* encoding the singlet invariant (J2) + the two-time
BC (J3); it reproduces the QM structure but does not yet derive it from the
substrate Lagrangian. Geometry + statistics give the form; the two-time BC gives
the `>2` violation; **deriving J2 (D1) and proving exact Tsirelson (D3) are the
remaining debts.**

## References (wire into `references.bib` if this graduates)

- Wheeler & Feynman 1945/1949 (absorber); Cramer 1986 (transactional) — the
  retarded/advanced two-time structure of J3. *Cramer already cited.*
- Aharonov, Bergmann & Lebowitz 1964; Aharonov TSVF — pre/post-selected
  (two-time) states, the QM scaffold J3 mirrors. *Already cited.*
- Block-solver chirality result (internal) — chiral-characteristic BC for the
  well-posed two-time BVP (J3, D2).
