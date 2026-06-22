# Bell Bridge — D3: does the two-time BVP yield exactly ±cos 2(a−b)?

Working note. **This does NOT prove D3 outright** — that would mean deriving QM
(Born rule + Hilbert structure) from the substrate, the program's deepest open
problem. It does three rigorous things: (1) shows the two-time structure *alone*
does not pin the correlation (PR-boxes also relax MI — `S = 4`); (2) identifies
the substrate ingredient that excludes them and caps the correlation at Tsirelson;
(3) reduces D3 to a single sharp residual — the Born weight — which is the
program-wide typicality debt (B-group / B4), now localized to the joint case.
All numerics below verified. Honest debt; not asserted in the paper.

## 0. What "prove D3" can and cannot mean

Given `(a, b, τ_A, τ_B)` the two-time BVP is deterministic and outputs
`(A,B) ∈ {±1}²`; the correlation is `E(a,b) = ⟨A·B⟩_τ`. Because the hidden state
at `S` may depend on *both* settings (MI relaxed, J3), the model is *not*
Bell-bounded — so it *can* reach the QM value. But "can reach" is not "lands
exactly on." Proving it lands exactly on `±cos 2(a−b)` is the hard content, and
§2 shows why it does not follow from the two-time structure alone.

## 1. The two-time BVP, precisely

- **Past BC at `S`:** the Bell state `|Φ⁺⟩` (J2, computed in
  `bellstate_lock_from_vertex.md`).
- **Future BCs at `D_A`, `D_B`:** absorption terminates each arm in the analyzer
  eigenstate (`|a⟩` at `D_A`, `|b⟩` at `D_B`) — the B4 threshold event as a
  future boundary condition. This *is* the projective/post-selection structure
  (the TSVF insight: a future BC is the "collapse").
- **Solution:** the time-symmetric stationary field on the branching domain; the
  advanced legs carry `a, b` back to `S`.

A two-time *history* is a consistent assignment `(|Φ⁺⟩ at S, |a⟩ at D_A,
|b⟩ at D_B)`; its **amplitude** is `⟨a|⟨b|Φ⁺⟩`. The correlation is
`E = Σ_{A,B} (AB) · w(A,B|a,b)` over the four histories, with weights `w` to be
fixed.

## 2. Obstacle: MI-relaxation is necessary but NOT sufficient (rigorous)

A two-time BVP permits the shared state to depend on both settings. But so does a
**PR-box**: a no-signalling, deterministic-given-`λ(a,b)` correlation with
`S = 4` (verified) — well above Tsirelson `2√2`. The PR-box fits the "λ depends on
both settings" mold the two-time structure provides. **Therefore the two-time
structure alone permits everything up to `S = 4`; it pins neither `±cos 2(a−b)`
nor even the Tsirelson cap.** Any claim "two-time ⇒ QM" is false without more.

## 3. The substrate ingredient: a linear complex `ℂ²` carrier ⇒ Hilbert amplitudes

The polarization carrier is a **linear, complex `ℂ²` field** (the `U(1)`-from-time
structure, Paper III). Its mode overlaps *are* `ℂ²` inner products — the qubit
Hilbert space. Hence the two-time history amplitudes are **genuine Hilbert
overlaps** `⟨a|⟨b|Φ⁺⟩`, not arbitrary numbers. This is the substrate fact that a
PR-box lacks: a PR-box has no underlying Hilbert amplitudes. Linearity +
complex polarization ⇒ the amplitude structure is quantum.

## 4. Born weight ⇒ exactly ±cos 2(a−b), Tsirelson automatic

With Born weights `w = |⟨a|⟨b|Φ⁺⟩|²`:

    E(a,b) = Σ_{A,B} (AB) |⟨a_A|⟨b_B|Φ⁺⟩|² = +cos 2(a−b),   CHSH |S| = 2√2.  ✓ (verified)

The Tsirelson cap is then **automatic**: Tsirelson's theorem bounds *every* `ℂ²`
Hilbert correlation with `±1` observables by `2√2`, and `|Φ⁺⟩` saturates it at the
optimal angles. PR-boxes (`S = 4`) are excluded precisely because the substrate
amplitudes are Hilbert overlaps (§3), not arbitrary. So:

    [linear complex ℂ² carrier]  +  [two-time absorption BCs]  +  [Born weight]
        ⇒  exactly ±cos 2(a−b),  capped at Tsirelson.

**The Born weight is load-bearing and specific.** A non-Born weight `w ∝ |amp|¹`
(otherwise identical) gives `E ≠ cos 2(a−b)` and `S = 1.66 ≠ 2√2` (verified). So
it is the `|amp|²` rule itself — not merely "some weight" — that produces the
correlation.

## 5. The residual: the Born weight = joint τ-typicality (B4 extended)

The only ingredient not derived from substrate structure is the Born weight
`w = |amp|²`. This is exactly B4's claim — *uniform effective threshold `τ` ⇒
Born* — now in the **joint** two-time setting: the `(τ_A, τ_B)` measure,
correlated through the BVP self-consistency at `S`, must reproduce
`|⟨a|⟨b|Φ⁺⟩|²`. Single-arm, B4 establishes `uniform τ ⇒ Malus = |amp|²`; the
joint residual is that the two-time-correlated `(τ_A,τ_B)` measure reproduces the
*joint* `|amp|²`. This is the program-wide Born/typicality debt (B-group), now
sharply localized — D3 is reduced from "prove the whole correlation" to "prove
the joint two-time weight is Born."

## 6. What this note establishes vs leaves open

**Established (rigorous / verified):**
- two-time structure alone ⇏ `±cos 2(a−b)`, and does not even cap at Tsirelson
  (PR-box, `S = 4`);
- the substrate's linear complex `ℂ²` carrier makes the two-time history
  amplitudes genuine Hilbert overlaps — which *excludes* PR-boxes and *caps* at
  Tsirelson;
- Hilbert overlaps + Born weight ⇒ **exactly** `±cos 2(a−b)`, `|S| = 2√2`, and the
  Born weight is specifically required (a non-Born weight fails).

**Open (the residual D3 = B1/Born-rule debt):** that the joint two-time
`(τ_A, τ_B)` measure reproduces the Born weight `|amp|²`. This is the same
typicality debt as B4/B-group, localized to the joint case — *not* a separate
miracle.

## Result

D3 is **reduced, not closed**: given (i) the substrate's linear complex `ℂ²`
carrier (Hilbert amplitudes, which cap the correlation at Tsirelson and exclude
PR-boxes) and (ii) the Born weight (the joint extension of B4's τ-typicality), the
two-time BVP yields *exactly* `±cos 2(a−b)`. The Tsirelson cap and the functional
form follow; the one genuine remaining debt is the Born weight — the joint
two-time typicality, the deepest and last open problem of the Bell bridge.

## References (wire into `references.bib` if this graduates)

- Tsirelson 1980 — the `2√2` bound on `ℂ²` Hilbert correlations.
- Popescu & Rohrlich 1994 — the PR-box (no-signalling, `S = 4`); the obstacle of §2.
- Aharonov TSVF; Cramer 1986 — two-time / post-selection structure (J3, §1).
- B4 (`threshold_detection.md`); Dürr–Goldstein–Zanghì typicality — the Born-weight
  residual of §5.
