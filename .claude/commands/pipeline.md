---
description: Run a BraneSim solver experiment (init+solve → worldvolume) from a JSON config. Defaults to orchestration/configs/branesim_chiral_smoke.json.
---

Run a BraneSim experiment via the package entry point `branesim.run_experiment`
(initialization + solve, writing `worldvolume.zip` + `summary.json`).

If `$ARGUMENTS` is empty, use `orchestration/configs/branesim_chiral_smoke.json`.
Otherwise, treat `$ARGUMENTS` as the path to a JSON config.

Execute:

```
python -m branesim.run_experiment --config <config> --output-dir ./run-out
```

Report `run-out/summary.json` (mode, `interior_residual_norm`, `solver_report`).

For the full measurement suite (energy / confinement / winding / Berry / EM /
per-colour SU(3) / spectra → CSV + paper-ready PNG + `report.md`) on an
experiment run folder, use:

```
python -m branesim.diagnostics.run_measurements <run_dir>
```