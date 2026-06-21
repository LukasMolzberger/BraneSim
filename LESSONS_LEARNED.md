# Lessons Learned

A durable record of mistakes not to repeat and results we trust, so the knowledge
survives without being expensively re-encoded in tests or frozen `test-runs/`.
(See `BACKBONE.md` for the theory, each paper's `derivations/<bridge>/status.md` *Open derivations* sections for open derivations.)

## Methodology — mistakes we kept making

1. **Clamp the pretension (periodic spatial BC).** The substrate has prestress
   `α<1`: bonds are stretched (rest length `αa < a`), so the vacuum is under
   tension. With **open/free boundaries** an unpinned tensioned lattice just
   relaxes — it contracts globally toward `αa`. The tell: an `A=0` control (no
   soliton) collapses identically. Use periodic BC for every nonlinear run; the
   tensioned square lattice is then a stable vacuum (`A=0 → E_excess ≡ 0`).
2. **A bound particle is an eigenstate — never seed-and-watch.** A hand-built lump
   marched/relaxed forward always disperses/unwinds/contracts because it is not a
   stationary solution. Solve for the self-consistent state; do not evolve a seed
   and hope it persists.
3. **Real diagnostics only.** The `compute_B_analytic` estimator hard-coded
   `F_inner = π`, so it reported `B≈1` for any vacuum-boundary field regardless of
   the interior — it never measured the winding. Always measure topological charge
   from the *actual* field (signed solid-angle degree), and measure energy as
   **vacuum-subtracted excess** against the stable periodic vacuum (raw `V*` is
   dominated by the prestress floor and is meaningless on its own).
4. **`ActionParams.r_t` defaults to 0 for a reason.** Defaulting it to the
   prestressed `α·β·dt` silently flips every default-constructed params onto the
   nonlinear spring path and hangs solvers (the breather eigen-solve spun forever).
   Prestressed is opt-in via an explicit `r_t`.
5. **Parameter sweeps rarely paid off.** Broad `(α, w, amplitude, …)` sweeps cost a
   lot and taught little; a derivation + one decisive, correctly-clamped test beat
   them every time. Prefer derive-then-test over scan.
6. **The discrete grid forbids Derrick collapse** (the core can't shrink below `a`)
   — but it does **not** prevent *unwinding* or *spreading*. "It can't collapse"
   was true and irrelevant to why solitons failed.

## Results we trust (machine-checked, periodic, small-amplitude or analytic)

- **Dispersion / light-cone:** `bvp_chiral` reproduces the analytic eigenmode to
  floating-point precision (residual/DOF ~1e-14). `c_T/c_L = √(1−α)`;
  `c_L² = k_s a²/m`, `c_T² = (1−α)k_s a²/m`. Cubic anisotropy `c_L([100])/c([111]) =
  1.074172`; `[110]` birefringence `1.06066`; `[111]` triplet degeneracy exact.
  Conditioning `κ` is N-independent (the chiral-BC guarantee).
- **Gauge layer:** the k-space Berry / Wilczek–Zee curvature is **identically zero
  ∀α** (real symmetric `D(k)` → trivial real eigenbundle). The linear envelope is a
  spin-1 vector. Spin-½ is therefore a *soliton-layer* `π₁(SO(3))=ℤ₂` effect, not
  linear. The only clean linear `α`-window is the color split (`∝α`). (P1/P2/P3.)
- **Geometric nonlinearity is exactly `∝α`:** the entire anharmonic sector is the
  norm term `−k_s αa|ΔR|`; `α=0` is exactly linear; quartic coeff `∝ k_s α/a`. StVK
  is the quadratic-order proxy of the spring (agree on wave speeds, differ at
  quartic). The substrate computes the exact spring, not StVK.
- **Kinematic color confinement (BACKBONE #24):** no linear direction is
  simultaneously coherent (needs branch degeneracy → only `[111]`) and color-active
  (needs `g(k̂)≠0` → only off `[111]`); colored content lives only in a nonlinear
  phase-locked soliton. This is an asymptotic-states statement, not a string-tension
  derivation.

## Dead ends (closed — do not redo without a new idea)

- **Breather (time-periodic Skyrme orbit):** converged 32³ orbits are robustly
  Floquet-unstable, `ρ≈3.5`, flat in `α`. The breathing/scale DOF carries the
  instability. Route closed.
- **Static B=1 hedgehog (FIRE minimization):** unwinds to vacuum. Two harness bugs
  inflated false positives (open-BC tension relaxation + the hard-coded B-estimator,
  lesson 1/3). All 2026-06-05 static/box-doubling/w-scan/saturation results were
  retracted.
- **Seed-and-watch worldtube (model-a / r_t=0):** the carrier was imposed by hand
  and radiated; with periodic BC it disperses. Nothing binds it along time without
  the temporal spring (`r_t>0`).
- **Worldtube with the temporal spring (`r_t=α·β·dt`):** binding is **untested**, not
  falsified — the test was confounded (no stable baseline soliton; a kinetic-kick
  seed). The real open question: does the substrate bind a B=1 object at all, on a
  clean periodic-clamped `r_t>0` footing (`papers/matter_mass/derivations/matter/status.md` C2).
