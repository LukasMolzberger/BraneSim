# Archived documents (historical — not current)

These are superseded planning documents, kept for the parts that are still
referenced (the VSH ansatz menu, the sprint structure) and for provenance. They
**predate** the 6-neighbor-only commitment, the 2026-06-05 retraction of the
nonlinear-soliton runs, and the 2026-06-06 U(1)-vortex/SU(3)-texture reframing,
so they describe the program in stale terms (legacy `components/` modules, the
breather solver vehicle, multi-shell `(k1,k2,k3)` couplings, an `L0–L8` layer
numbering). Read them as snapshots, not as the current plan.

For the current theory and program see, at the repo root:
`BACKBONE.md`, `PRINCIPLES.md`, `OPEN_PROBLEMS.md`, `EXPERIMENT.md`,
`LESSONS_LEARNED.md`.

- `VALIDATION_ROADMAP.md` — sprint-organized validation subtasks (L0→gravity).
- `BARYON_SIMULATION_ROADMAP.md` — soliton/baryon search program; the VSH ansatz
  menu (Candidates 1–5) is still referenced by `BACKBONE.md` #20.
- `run_pipeline.py` — the obsolete 4-separate-process pipeline (init → sim → viz →
  diag) that imported the deleted `components.*` modules. Superseded by the
  single entry point `branesim/run_experiment.py` (init+solve) +
  `branesim/diagnostics/run_measurements.py` (measure).
- `run_topological_breather_validation.py` — standalone driver for the
  topological-breather baryon candidate (open boundary). Part of the retracted
  2026-06-05 nonlinear-soliton line; the breather route is closed
  (`LESSONS_LEARNED.md`, `OPEN_PROBLEMS.md` C2).