# BraneSim — Open Problems & Open Derivations (relocated)

**The open-derivation entries that used to live here have been distributed to the
per-paper "bridge" status files**, next to the HAVE / MISSING ledger for each
bridge. This file is now a thin index so that existing cross-references
(`OPEN_PROBLEMS.md A1`, `OPEN_PROBLEMS D6`, …) still resolve to the right place.

Each bridge's `status.md` now carries an `## Open derivations` section holding the
full entries (statement → why it matters → candidate approach → status), with the
original IDs preserved.

| ID(s) | Topic | Now lives in |
|-------|-------|--------------|
| **A1, A2, A3, A4, A4a** | Foundational solver — 4D block-variational formulation (saddle-point solve, two-time BVP well-posedness, `c²` Lorentz tuning, temporal-link form, single 4D prestress `α_t=α`) | `papers/core/derivations/status.md` |
| **B1, B2, B3** | Quantum-foundations / Bell stance (Tsirelson bound, no-signalling, baryon-to-photon ratio `η`) | `papers/bell/derivations/status.md` |
| **C1, C1a, C2, D1** | Constitutive law & soliton stability (StVK non-coercivity, StVK quartic α-scaling, rest baryon, EM-charge ⇄ spin-½ bridge) | `papers/matter_mass/derivations/matter/status.md` |
| **D2, D3, D5** | Gauge sector — EM/`U(1)` (prestress α from EM/colour ratio, dynamical Maxwell, massive-vs-massless screening) | `papers/gauge_color/derivations/gauge/status.md` |
| **D4, D6** | Gauge sector — `U(1)`↔`SU(3)` binding (vortex/texture binding & electron stabilization, what stops the proton splitting) | `papers/gauge_color/derivations/color/status.md` |

The single-dial co-variation note shared by **A4a** and **D6** (gravity strength
and binding strength co-vary in `α`) is also reflected in
`papers/lorentz_gravity/derivations/gravity/status.md`.

## Scoped non-claims (former group E — tracked elsewhere)

These are bounded-scope limits already stated honestly in each paper's
"Non-claims" / scope section; they are *scope statements*, not active derivation
tasks: derived value of `G` (Newtonian limit of the contraction channel), full
Maxwell / Yang–Mills dynamics in all regimes, a finished proton/neutron solution,
the full Standard-Model gauge group / generations / couplings, a closed
multi-particle scattering theory.

## Related documents

- `PRINCIPLES.md` — canonical non-negotiables
- `BACKBONE.md` — canonical project backbone
- `archive/VALIDATION_ROADMAP.md` — numerical validation sprints *(archived/historical)*
- `archive/BARYON_SIMULATION_ROADMAP.md` — soliton search program *(archived/historical)*