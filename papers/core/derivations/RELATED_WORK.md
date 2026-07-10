# Related Work

This work derives an emergent `U(1)×SU(3)` gauge sector (T1, T3, T4, T7), an
emergent Lorentzian metric (T2), and a temporal-prestress kinetic limit (T12) from a
single deterministic 4D brane-lattice substrate. It therefore sits at the
intersection of three established lines of research: (i) the geometry of Berry /
Wilczek–Zee phases, (ii) emergent gauge fields and emergent relativity from
condensed-matter / lattice substrates, and (iii) the no-go constraints on emergent
gauge bosons. The idea that gauge structure and Lorentz invariance can *emerge* from
holonomy and low-energy dynamics is well established; the novelty here is deriving
**both** an abelian Maxwell sector **and** a non-abelian `SU(3)` Yang–Mills sector
from the **same** carrier geometry of a **deterministic prestressed lattice**, and
verifying the `su(3)` closure and the sector split numerically.

> **Bibliographic note.** All entries below were **verified against authoritative
> sources** (publisher pages / NASA ADS / DOIs) — titles, authors, journal, volume,
> pages, and year confirmed. Two entries from an earlier auto-generated table were
> removed or flagged as unverifiable — see §7.

---

## 1. Geometric-phase foundations

- **M. V. Berry**, *Quantal phase factors accompanying adiabatic changes*, Proc. R.
  Soc. Lond. A **392**, 45–57 (1984). The abelian geometric phase. Foundation of **T3**:
  the emergent U(1) is the Berry connection of a non-degenerate band, equivalently the
  trace of the `U(3)` Wilczek–Zee connection.
- **B. Simon**, *Holonomy, the quantum adiabatic theorem, and Berry's phase*, PRL
  **51**, 2167–2170 (1983). The holonomy/fibre-bundle reading of Berry's phase — the
  language used throughout T1–T4.
- **F. Wilczek & A. Zee**, *Appearance of gauge structure in simple dynamical
  systems*, PRL **52**, 2111–2114 (1984). The non-abelian geometric phase: a degenerate
  rank-`N` eigenspace carries a `U(N)` connection. Theoretical basis of **T1/T4**
  (here `N=3`, giving a `U(3)` connection whose traceless part is the `SU(3)` sector).
- **A. Shapere & F. Wilczek (eds.)**, *Geometric Phases in Physics*, World Scientific
  (1989). Standard reprint volume establishing the general role of geometric phase in
  physics.

## 2. Berry curvature and quantum geometry in band theory (methods for T3/T4/T7)

- **J. P. Provost & G. Vallée**, *Riemannian structure on manifolds of quantum
  states*, Comm. Math. Phys. **76**, 289–301 (1980). The **quantum geometric tensor** /
  quantum metric (Fubini–Study). This is the object used as the induced coupling
  `1/e²` (T3) and `1/g²` (T4): the symmetric (metric) part of the same tensor whose
  antisymmetric part is the Berry curvature.
- **D. Xiao, M.-C. Chang, Q. Niu**, *Berry phase effects on electronic properties*,
  Rev. Mod. Phys. **82**, 1959–2007 (2010). Standard reference for Berry-curvature band
  theory, semiclassical transport, and the quantum metric — the backdrop for T3/T4/T7.
- **T. Fukui, Y. Hatsugai, H. Suzuki**, *Chern numbers in discretized Brillouin zone*,
  J. Phys. Soc. Jpn. **74**, 1674–1677 (2005). The gauge-invariant Wilson-loop lattice
  Berry-curvature method **used directly in the T1/T3/T4 code** (non-abelian
  generalisation for the rank-3 carrier). Should be cited in Methods.

## 3. Emergent gauge fields from lattice / many-body systems (peer group for T1/T4)

- **G. E. Volovik**, *The Universe in a Helium Droplet*, Oxford Univ. Press (2003).
  The closest relative: emergent Weyl fermions, emergent gauge fields, and emergent
  gravity from superfluid ³He. Establishes the whole paradigm of relativistic fields
  as low-energy modes of a non-relativistic substrate.
- **X.-G. Wen**, *Quantum Field Theory of Many-Body Systems*, Oxford Univ. Press
  (2004); **M. A. Levin & X.-G. Wen**, *String-net condensation: a physical mechanism
  for topological phases*, Phys. Rev. B **71**, 045110 (2005). Emergent gauge
  fields — **including non-abelian ones** — and even emergent fermions from purely
  bosonic lattice models. Direct precedent for "gluons/photons as lattice holonomy."
- *(genuine emergent-`SU(3)` example, condensed matter)* — e.g. emergent `SU(3)`
  flux on antiferromagnetic skyrmion lattices, Nat. Commun. (2024). A real neighbour
  to T1, in contrast to the fabricated entry removed in §7.

## 4. Emergent / analogue Lorentzian geometry (peer group for T2)

- **W. G. Unruh**, *Experimental black-hole evaporation?*, PRL **46**, 1351–1353 (1981).
  Phonons in a moving medium see an emergent Lorentzian ("acoustic") metric.
- **C. Barceló, S. Liberati, M. Visser**, *Analogue Gravity*, Living Rev. Relativity
  **8**, 12 (2005). Comprehensive review of emergent metrics from media. **T2's
  effective metric from the phonon dispersion is an acoustic metric** in this sense;
  the anisotropy/birefringence analysis of T2 is the analogue-gravity way of
  quantifying deviations from an exact cone.

## 5. Emergent gauge & Lorentz invariance as low-energy attractors (the philosophy)

- **D. Förster, H. B. Nielsen, M. Ninomiya**, *Dynamical stability of local gauge
  symmetry — Creation of light from chaos*, Phys. Lett. B **94**, 135–140 (1980). Gauge
  invariance as a low-energy attractor of a generic, non-gauge-invariant lattice theory
  ("random dynamics") — the philosophical ancestor of this programme.
- **J. D. Bjorken**, *A dynamical origin for the electromagnetic field*, Ann. Phys.
  **24**, 174–187 (1963). Early emergent-photon proposal (as a Goldstone mode).

## 6. Constraints and no-go theorems the programme must address

- **S. Weinberg & E. Witten**, *Limits on massless particles*, Phys. Lett. B **96**,
  59–62 (1980). Forbids emergent massless charged spin-1 (and spin-2) particles in an
  **exactly** Lorentz-covariant theory. Directly constrains **T3/T4** (and any future
  emergent graviton in Layer −). **Escape:** Lorentz invariance here is emergent and
  approximate (T2), so the exact-covariance premise fails — the same resolution used
  by Volovik/Wen/analogue-gravity. Developed in `T4_derivation_and_proof.md` §8 and
  `T2_derivation_and_proof.md` §9.
- **H. B. Nielsen & M. Ninomiya**, *Absence of neutrinos on a lattice: (I). Proof by
  homotopy theory*, Nucl. Phys. B **185**, 20–40 (1981). Fermion doubling on a lattice
  — a constraint the **Layer-5** matter/soliton sector must confront once fermionic
  matter is introduced.
- **J. S. Bell**, *On the Einstein Podolsky Rosen paradox*, Physics Physique Fizika
  **1**, 195–200 (1964). The foundational no-go: local causality + measurement
  independence ⇒ Bell
  inequalities. Constrains the **foundations** of the deterministic substrate (BACKBONE
  #12); the escape is measurement-dependence via the all-at-once time-symmetric
  structure (see next entry and `FOUNDATIONS_bell.md`).
- **K. B. Wharton & N. Argaman**, *Colloquium: Bell's theorem and locally mediated
  reformulations of quantum mechanics*, Rev. Mod. Phys. **92**, 021002 (2020);
  arXiv:1906.04313. *(Verified to exist.)* Locally mediated, time-symmetric
  reformulations that evade Bell by relaxing an arrow-of-time (measurement-
  independence) assumption rather than adding spacelike channels. The class the
  brane-lattice foundations belong to — see `FOUNDATIONS_bell.md`. **Consistency /
  positioning only**: it shows Bell does not *exclude* the class, not that the model
  reproduces quantum correlations (that is the owed BACKBONE #12 program).

## 7. Berry–Maxwell / parameter-space electromagnetism (context for T3)

- **G. Konstantinou**, *Time-Independent Parameters in Quantum Systems: Revisiting
  Berry Phase, Curvature and Gauge Connections*, arXiv:2507.23347 (2025). *(Verified
  to exist.)* Constructs "Berry–Maxwell equations" — Berry electric/magnetic fields
  from scalar/vector potentials of the quantum state. Supports the T3 reading of
  Berry curvature as an emergent electromagnetic structure in parameter space.
- *(more T2/T3-relevant)* arXiv:2308.00612, *Constructing Berry–Maxwell equations
  with Lorentz invariance and Gauss' law of Weyl monopoles in 4D energy-momentum
  space* (2023). Berry–Maxwell equations **with Lorentz invariance in 4D** — closer
  to the combined T2+T3 setting than the entry it supplements.

---

## Corrected comparison table

| Source | Main result | Relation to this work |
|---|---|---|
| Berry, Proc. R. Soc. A **392**, 45 (1984) | Abelian geometric phase | Basis of **T3** (U(1) = Berry connection / trace of `U(3)`) |
| Simon, PRL **51**, 2167 (1983) | Berry phase as holonomy | Bundle language of T1–T4 |
| Wilczek & Zee, PRL **52**, 2111 (1984) | Non-abelian `U(N)` geometric phase | Basis of **T1/T4** (`N=3` → `U(3)`, traceless = `SU(3)`) |
| Provost & Vallée, CMP **76**, 289 (1980) | Quantum geometric tensor / metric | The induced `1/e²`, `1/g²` (**T3/T4/T7**) |
| Xiao, Chang, Niu, RMP **82**, 1959 (2010) | Berry-curvature band theory | Backdrop for T3/T4/T7 |
| Fukui, Hatsugai, Suzuki, JPSJ **74**, 1674 (2005) | Lattice Wilson-loop Berry curvature | **Method used in the T1/T3/T4 code** |
| Volovik, *Helium Droplet* (2003) | Emergent Weyl fermions/gauge/gravity from ³He | Closest paradigm relative (**T1–T4, T2**) |
| Wen (2004); Levin & Wen, PRB **71**, 045110 (2005) | Emergent (non-abelian) gauge fields from lattices | Direct precedent for **T1/T4** |
| Unruh, PRL **46**, 1351 (1981); Barceló–Liberati–Visser, LRR **8**, 12 (2005) | Emergent/acoustic Lorentzian metric | **T2** is an acoustic-metric construction |
| Förster–Nielsen–Ninomiya, PLB **94**, 135 (1980) | Gauge invariance as low-energy attractor | Programme philosophy |
| Weinberg & Witten, PLB **96**, 59 (1980) | No emergent massless charged spin-1/spin-2 (exact Lorentz) | **Key constraint on T3/T4**; evaded via approximate Lorentz (T2) |
| Bell, Physics Physique Fizika **1**, 195 (1964) | Local causality + measurement independence ⇒ Bell inequalities | Constrains **foundations** (BACKBONE #12) |
| Wharton & Argaman, RMP **92**, 021002 (2020) | Locally-mediated, time-symmetric reformulations relax measurement independence | The brane-lattice foundations' class (`FOUNDATIONS_bell.md`) — consistency, not derivation |
| Konstantinou, arXiv:2507.23347 (2025) | Berry–Maxwell equations | Context for **T3** |

## Removed / flagged entries (from the earlier auto-generated table)

- **REMOVED — "C. Roberge, *Emergent SU(3) Gauge Geometry* (2026)."** Could not be
  located in any search (no arXiv ID, no venue). Its stated result matches T1 almost
  verbatim, and it shows the classic signature of an LLM-fabricated citation. **Do not
  cite** unless the actual paper is located and read.
- **FLAGGED — "Mukunda et al., Eur. J. Phys. 3, 132 (1982)."** Likely misattributed:
  geometric-phase work in 1982 (pre-Berry) in a physics-education journal is
  implausible; the well-known Mukunda–Simon geometric-phase papers are ~1993
  (Ann. Phys.). Verify or replace before citing.
- **FLAGGED — "O. Ulfbeck & A. Bohr (1992)."** Too vague to cite as-is; locate exact
  reference (the Bohr–Ulfbeck foundations papers I am aware of are ~2000) or drop.

---

## Statement of novelty

Emergent gauge structure from geometric phase/holonomy (Berry, Wilczek–Zee, Simon)
and emergent relativistic fields from condensed-matter substrates (Volovik, Wen,
analogue gravity) are both well established. The contribution here is to obtain, from
**one** deterministic prestressed 4D brane-lattice and its **single** nonlinear
periodic vacuum carrier:

1. a rank-3 Wilczek–Zee carrier whose traceless curvature is numerically certified to
   span **full `su(3)`** (T1), with an explicit `so(3)`-vs-`su(3)` mechanism (complex
   Bloch phases on noncommuting real frames);
2. the **`U(1)` Maxwell** (T3) and **`SU(3)` Yang–Mills** (T4) sectors from the *same*
   carrier geometry, with substrate-induced positive couplings;
3. an emergent **Lorentzian metric** with signature fixed by the prestress sign
   pattern (T2), and a quantitative **U(1)/SU(3) split** (T7);

and to situate all of this against the **Weinberg–Witten** constraint, with the
emergent/approximate nature of Lorentz invariance (T2) as the explicit escape.
