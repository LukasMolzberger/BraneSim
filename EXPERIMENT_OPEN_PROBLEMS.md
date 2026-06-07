# Experimental-Setup Open Problems — the U(1) vortex run

Tracker for problems with the **experimental apparatus** (seed → solve → render →
diagnose pipeline in `branesim/experiments/vortex_seed_render.py` + the AWS
orchestration), **not** with the physics derivations.

## Why this is separate from `OPEN_PROBLEMS.md`

`OPEN_PROBLEMS.md` tracks *theory* gaps — derivations that bridge physical layers
(solver formulation A, QM foundations B, soliton stability C, gauge sector D).
Those are answered with mathematics.

This file tracks *methodology / instrumentation* gaps — **does the apparatus
actually instantiate the object it claims to study, and does it measure that
object faithfully?** These are answered by fixing code, the seed, the solver
convergence criterion, or a diagnostic's pass/fail logic — no new physics. A run
can be physically sound in principle yet still produce a misleading verdict
because the instrument is mis-wired; that class of bug lives here.

Statuses: `open` · `in-progress` · `resolved` · `wip-expected`. Each entry names
the falsifier: the concrete check that would close it.

**`wip-expected`** marks items that are *not apparatus defects* — they are the
normal state of a work-in-progress run (the eigen-solve and the final seed are
not yet wired in). They resolve themselves once the build is complete; they are
listed only so a reader does not mistake the current run's output for a finished
result. The real punch list is the `open` structural items (E1, E2, E4, E5, E6, E7).

---

## Group I — Alignment: the apparatus does not yet instantiate the stated target

These are the gaps between EXPERIMENT.md's *target object* and what the current
seed/solve actually produces. While open, run verdicts about "the EM/electron
U(1) vortex" are about a **different object** than the one named.

### E1. The seed is NOT a pure EM/U(1) vortex — it is ⅓ EM + ⅔ color — `open`

The injector (`vortex_worldtube.py`) writes the carrier into a **single lateral
component** (`CARRIER_RE = 2`). By BACKBONE #25 and EXPERIMENT.md §"Injection
ansatz", a single-lateral-component vortex reads **⅓ trace (U(1)/EM) + ⅔
traceless (SU(3)/color)** under the U(3) projection. The D6 diagnostic confirms
this exactly: `u1_fraction_mean = 0.333`. So the headline object — "the EM /
electron-like sector" — is mislabeled; what is injected is dominantly color.
EXPERIMENT.md "Next (i): re-inject in the trace (EM) direction" is still pending.

**Falsifier / fix.** Seed the carrier in the symmetric `(1,1,1)/√3` trace
direction (all three lateral components sharing the winding phase). D6 should
then read `u1_fraction → 1` for the bare seed, with `su3_fraction` rising only if
the relaxation genuinely coexcites color. Until then, the run cannot answer
"does the EM vortex bind."

### E2. No tumble → spin-½ is not in the experiment at all — `open`

EXPERIMENT.md lists three deliverables; #3 is "the spin-½ ℤ₂ holonomy, measured
directly." The ingredient that produces it is the **2π SO(3) tumble of the axis
n̂(t) over the loop** (the twisted worldsheet → Finkelstein–Rubinstein ℤ₂). The
seed has a **fixed** axis (z), so the worldsheet is untwisted and there is no
spin holonomy to measure. There is also no diagnostic that would extract it.

**Falsifier / fix.** Add the axis tumble to `VortexParams`/the injector, and a
diagnostic that measures the ℤ₂ FR sign over the loop. Without it, drop spin-½
from the run's stated deliverables.

### E3. "Eigenstate, not seed-and-watch" is not realized — the solve does not converge — `wip-expected`

> **WIP note.** Expected: the converged eigen-solve is not yet wired in (the run
> still does a bounded fixed-iteration relaxation probe). Not a defect — resolves
> when the eigenstate solve lands. Listed so the current transient output is not
> read as a finished result.

EXPERIMENT.md §"Solve" commits to root-finding `‖𝓡‖=0` over the periodic loop.
In practice the run does a **fixed-iteration** JFNK relaxation
(`RELAX_ITERS = 15`, the eigensolve pretest used 20) that **does not converge**:
the last eigensolve reported `residual_initial = 78.2 → residual_final = 1.78`,
`tol = 1e-6`, **`converged = False`**. The module's own "5th-dimension" framing
concedes only the converged fixed point is physical and "iterates before it are
unphysical transients" — yet every diagnostic is run on exactly such a transient
(`iter_0015`/`iter_0020`). The measurements are therefore of a non-stationary
intermediate, not of an eigenstate.

**Falsifier / fix.** Either (a) drive the solve to `‖𝓡‖/DOF ≤ tol` and gate
diagnostics on `converged == True`, or (b) if the periodic operator cannot reach
that, document the achieved residual as the headline number and stop describing
the output as an eigenstate. A residual-vs-iteration plot belongs in every run.

### E4. ωT quantization is imposed, not solved as a BVP eigenvalue — `open`

EXPERIMENT.md: "`T` is the existence / quantization condition," to be found as a
BVP eigenvalue. Currently `n_t = 2` is hard-coded and `ω = 2π·n_t/n_slices` is
fixed kinematically by loop closure. The quantization is asserted by
construction, not discovered by the solver.

**Falsifier / fix.** Solve for the `(ω, T)` (equivalently `n_t`/period) that
admits a stationary bound worldtube, e.g. via `solver/breather.py`'s periodic
eigen-BVP, rather than fixing `n_t` a priori.

---

## Group II — Instrumentation: devices that report misleading verdicts

### E5. D3 winding device has inverted pass/fail and a self-contradictory verdict — `open`

For the `Y_1^1` seed the experiment *wants* `winding_z = m = 1`
(`winding_closure.json`: `winding_through_z_normal = 1.0`, `winding_ok = true`).
But `device_winding` sets `closure_ok = (max|winding| < 0.1)`, so it flags the
**desired** unit winding as a failure. The current run's `report.md` literally
prints:

> **D3 Winding**: net winding = 0 (closure **FAIL**); …  ← `max|winding| = 1.00`

— "net winding = 0" and "FAIL because winding = 1" in the same sentence. The note
("Net winding = 0 for a contractible vortex ring") is **stale**, left over from
the deleted smoke-ring seed; it is now actively wrong for `Y_1^1`. D3 is, as
written, non-discriminating and misleading.

**Falsifier / fix.** Rewrite D3's pass/fail to the actual target: `winding_z ≈ m`
(within tol) AND `winding_x, winding_y ≈ 0`, **constant across slices**. Delete
the contractible-ring note.

### E6. Condition-number claim is wrong by ~4 orders of magnitude — `open`

`DEPLOYMENT.md`, `orchestration/aws/RUNBOOK.md`, the experiment module's
docstrings/manifest, and several memory files all state PeriodicBC
`cond ~ 1e3–1e4`. The actual `condition_estimate` reported by the last eigensolve
is **`3.1e7`**. It is still far below Dirichlet's `1e14` (so PeriodicBC is the
right choice), but the specific quoted figure is false and has propagated
everywhere.

**Falsifier / fix.** Recompute `PeriodicBC.condition_estimate` across the
production grids and replace every `1e3–1e4` claim with the measured value (and
its grid/`n_slices` scaling).

### E7. `iterations` reporting is hard-coded; provenance text is over-optimistic — `open`

`solve_block`'s PeriodicBC path reports `"iterations": opts.max_iter`
unconditionally (`bvp.py:395`), regardless of what the solver did. The manifest
and RUNBOOK describe the relaxation as "residual drops ~100×" — the eigensolve's
actual drop was **44×** (78.2 → 1.78), to a residual still **6 orders of
magnitude** above `tol`. The provenance reads like a near-success; it was a stall.

**Falsifier / fix.** Report the true outer-iteration count and the achieved
`residual_final / tol` ratio; replace "drops ~100×" with the measured factor and
state `converged = False` prominently.

### E8. Berry/EM diagnostics on the seed are injector read-backs, not physics — `wip-expected`

> **WIP note.** Expected while diagnostics run on the seed/unconverged state
> (gated on E3). Not a defect — the read-backs correctly validate the injector
> round-trip. Only the `report.md` *wording* needs to wait for a solved state
> before claiming "emergent EM/Berry."

D4 (Berry phase) and D5 (EM `A_μ`,`F_μν`) are computed on `iter_0000`, where — by
the module's own docstring — "the carrier phase advance is prescribed kinematics,
not dynamics." So the accumulated Berry phase is exactly `2π·n_t` by construction,
the B-flux is the injected winding, and the E-field ≈ 0 (since `A_t ≈ 0` on the
seed). These are valid **round-trip checks of the injector**, but `report.md`'s
verdicts ("B field shows vortex topology") read as if they were emergent
measurements of the substrate.

**Falsifier / fix.** Label seed-state diagnostics as injector round-trip checks;
reserve "emergent EM/Berry" language for a converged solved state (gated on E3).

### E9. Winding preservation is asserted but unverified — and currently contradicted — `open`

The manifest/RUNBOOK claim "the carrier winding is preserved." The relaxed
`iter_0020` D3 reports `max|winding| = 2.00` (seed slice-0 was `winding_z = 1`).
Whether this is a genuine winding change, a different axis picking up circulation,
or a transient artifact is unknown — the suite does not certify preservation, and
the one number bearing on it shows a change buried under E5's inverted pass/fail.

**Falsifier / fix.** Add a per-slice winding-vs-iteration trace and a explicit
"winding conserved across the loop and across the solve" check; resolve whether
the `→2` is physical or a contour artifact.

### E10. Binding signals on the relaxed state are at the noise floor but reported as positive — `wip-expected`

> **WIP note.** Expected: the binding signals are measured on the non-converged
> relaxation (gated on E3), so noise-floor magnitudes are normal at this stage.
> Not a defect — the D8 device works; only the "binding ON" *verdict* should wait
> for a converged state and a stated per-probe noise threshold.

The eigensolve pretest declared "Part-1 time-link binding source turned ON" from
a relaxed `Δu_∥` extremum of **−9.9e-6** (seed: −1.7e-8), `sep_d = 3.3e-3`,
`R_kähler = 6e-6` (seed: 8.6e-18). These are round-off/transient-scale numbers on
a **non-converged** worldvolume (E3). Reading a binding verdict into them is
premature.

**Falsifier / fix.** Re-evaluate D8 only on a converged state with a stated
noise-floor threshold per probe; require signal ≫ floor before any "binding ON"
claim.

---

## Group III — Documentation / configuration drift (minor)

### E11. Production scale disagrees between docs — `open`

EXPERIMENT.md §"Scale ladder" targets production `96³ × 128`; RUNBOOK and the
module default to `n_slices = 32` and explicitly "prefer raising the spatial grid
over the slice count" (because the periodic-operator condition grows ~`P²`). Pick
one production target and make the docs agree.

### E12. `condition_estimate` is a linear modal proxy, mislabeled as the solve's conditioning — `open`

`PeriodicBC.condition_estimate` is a *linearized* modal ratio computed
independently of the actual nonlinear JFNK behavior. It is surfaced as the run's
"cond," but the solve stalled at `‖𝓡‖ = 1.78` for reasons the proxy does not
capture. Either compute an empirical conditioning estimate from the Krylov
iteration, or label this number as "linear modal estimate (not the achieved
nonlinear conditioning)."

---

*Generated from an audit of `vortex_seed_render.py`, `vortex_worldtube.py`,
`run_measurements.py`, `bvp.py`/`boundary.py`, the AWS runbook, and the
`runs/` artifacts (seed `vortex_seed_2026-06-07_091740`, eigensolve
`binding_eigensolve_2026-06-07_083140_local`).*