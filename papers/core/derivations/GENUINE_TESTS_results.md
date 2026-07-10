# Genuine tests (D7 and T2 gauge-mode) — results

These are the two "could-cut-either-way" tests. Both returned **cautionary** results
that bound the claims. Better surfaced by us than by a referee. Scripts:
`d7_vacuum_selection.py`, `t2_gauge_mode_isotropy.py`.

---

## D7 — Is the helix the selected vacuum? → No: a stable *texture*, not the ground state

- **Energetics:** the helix stores **+42.5% more elastic energy** per node than the
  straight tensioned vacuum (`E_helix = 0.359` vs `E_straight = 0.252`). So the helix
  is a **finite-energy texture**, not the elastic ground state.
- **Linear/elastic stability:** every link Hessian bracket is positive-definite
  (min `1−r/L = +0.20`, at the temporal link), so the spatial fluctuation Hessian
  `V(q)` is PSD (numerically min eigenvalue `+0.12` over sampled `q`), and with the
  temporal kinetic term PD all fluctuation `ω² ≥ 0`: the helix is **linearly stable**
  (no tachyonic/runaway elastic modes).

**Verdict:** the color-carrying helix is a **stable but non-ground-state microtexture.**
Full dynamical (Floquet) selection of *the* vacuum is not settled (the 4D action is an
indefinite saddle, so energy minimization doesn't decide it) — D7 stays open, but the
specific helix is now known to be a texture, not the ground state.

**Consequence for framing (important):** this *qualifies* Decision M ("color is a
universal vacuum gauge sector"). The straight ground state has **no** color (trivial
bundle, T1); the su(3) structure is carried by a higher-energy microtexture. So the
defensible claim is **"an su(3) gauge sector exists on a stable finite-amplitude
periodic texture of the substrate,"** not "the ground-state vacuum is a Yang–Mills
field." The "universal over empty space" language should be softened accordingly (or
D7 must exhibit a color-carrying background that is at or below the straight-vacuum
energy — not yet found).

---

## T2 gauge-mode isotropy — Are gauge modes more isotropic than sound? → No

The emergent gauge kinetic tensor is the spatial quantum metric
`g_ij(k) = ½ Re Tr[(∂_iP_3)(∂_jP_3)]`. Its eigenvalue spread sets the photon/gluon
cone anisotropy.

- single-helix carriers are strongly axis-anisotropic (spread `1.6–2.6`) — expected,
  since a helix picks its propagation axis;
- the **orientation-averaged** gauge tensor (fair measure) has anisotropy
  **≈ 0.81**, *larger* than the bare acoustic **0.254**.

**Verdict:** the gauge-mode cone is **not** more isotropic than the bare phonons — the
hope that "light" is closer to Lorentz-invariant than sound is **not** borne out; if
anything the gauge sector is more anisotropic. The **T2 dual-observer / emergent-Lorentz
obstruction therefore persists** and is not rescued by moving from phonons to gauge
modes. (Caveat: orientation averaging used only 3 axes with few isolated-carrier
k-points each; the qualitative conclusion is robust given the large single-orientation
values, but the exact 0.81 is a rough estimate.)

---

## Combined implication for the paper

Both tests point the same way and should reshape the framing of the two *weakest*
core claims (the gauge STRUCTURE results T1/T3/T4/T7 are unaffected and remain strong):

1. **Not a Lorentz-invariant ground-state vacuum.** Color lives on a stable,
   finite-energy, anisotropic **microtexture**, and its emergent metric is genuinely
   anisotropic in the gauge sector too. The honest headline is *"emergent U(1)×SU(3)
   gauge structure on a substrate texture,"* with emergent Lorentz invariance an
   explicit open conjecture (the dual-observer renormalization remains the only route,
   now with the obstruction quantified in both sectors).
2. **Decision M needs qualifying** (universal-vacuum-color) and **T2 universality stays
   conjectural** — recommend stating both plainly in the Limitations section.

These do not undermine the gauge-structure demonstrations; they correctly bound the
*vacuum* and *Lorentz-universality* interpretations around them.
