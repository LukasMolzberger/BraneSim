# Sprint 1 subtask 3 — Lab-frame anisotropy `c_L([100]) / c([111])` at `α = 0.2`

**Status:** **PASS** at `1.9 × 10⁻⁸` relative error against the analytic
prediction. The 1 % criterion in `paper/validation_roadmap.md` is met by
~6 orders of magnitude.

## What is claimed

On the canonical 6-neighbor axial-only lattice (backbone #15), the
long-wavelength branch speeds along `[100]` and `[111]` differ structurally:

  - `c_L([100])  = 1                    ` (longitudinal, α drops out)
  - `c   ([111]) = √((3 − 2α) / 3)      ` (Cartesian-pol eigenmodes, all three degenerate)

so the lab-frame anisotropy ratio is

  - `c_L([100]) / c([111])  =  √(3 / (3 − 2α))`

At `α = 0.2` this is `√(3/2.6) ≈ 1.07417`. This anisotropy is **not retuned
away** — backbone #8 and #19 carry it forward as the structural feature that
provides the operational basis for the SU(3) gauge sector under the
dual-observer framework.

## Measurement

The two branch speeds come from the same sweep that produced subtask 2 (see
`test-runs/sprint1_subtask2_dispersion/`). For each direction the speed
`c(|k|·a)` is measured at three `k`-modes (`n = 1, 2, 3` on `N = 32`) and
extrapolated to `k → 0` via a quadratic fit in `(|k|·a)²`:

| Direction | `|k|·a` samples | `c (k→0)` measured | `c (k→0)` analytic | rel err |
|---|---|---|---|---|
| `[100]` (L) | 0.196, 0.393, 0.589 | 0.99998810 | 1.00000000 | −1.19e-05 |
| `[111]` (T) | 0.340, 0.680, 1.020 | 0.93093828 | 0.93094934 | −1.19e-05 |

The same systematic −1.19e-5 bias appears on both branches; it is the
single-step integration error of `VelocityVerletSolver` integrated over ~6
oscillation periods. It is common-mode and largely cancels in the ratio:

| | measured | analytic |
|---|---|---|
| `c_L([100])`              | 0.99998810  | 1.00000000 |
| `c   ([111])`             | 0.93093828  | 0.93094934 |
| **ratio**                 | **1.07417229** | **1.07417231** |

**Ratio relative error: `−1.86 × 10⁻⁸`.**

## Files

- `run_dispersion.py`, `dispersion_results.json`, `dispersion_summary.json`,
  `dispersion_raw.npz` (in
  `test-runs/sprint1_subtask2_dispersion/`).
- The anisotropy ratio and its analytic prediction are stored under the
  `anisotropy` key of `dispersion_summary.json`.

## What this does and does not establish

- **Establishes:** at the canonical operating point, the lab-frame
  longitudinal anisotropy `c_L([100]) / c([111])` matches the closed-form
  prediction `√(3/(3−2α))` to numerical precision. The structural anisotropy
  is real and quantified.
- **Establishes:** common-mode integrator bias cancels in the ratio, so this
  measurement is more accurate than either branch speed in isolation.
- **Does not establish:** observer universality (`Sprint 3 #11`) — the
  inside-observer claim is that solitons of different polarization renormalize
  rods/clocks coherently and so report an isotropic effective metric. That is
  a separate sprint and depends on having a stable soliton (Sprint 4) before
  it can even be set up.
- **Does not establish:** that the lab-frame anisotropy is the *correct*
  structural source for the SU(3) sector (backbone #19) — that is the gauge
  thread (Sprint 2), which requires the complex-envelope reduction and is not
  performed here.