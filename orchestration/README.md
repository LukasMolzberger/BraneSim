# Orchestration

Orchestration runs outside the components and invokes each component as a separate process.

Run with config:

```bash
python orchestration/run_pipeline.py --config orchestration/configs/local_extensive.json
```

This writes one run directory with:
- `initial_state.npz`
- `trajectory.zip`
- `plots/*.mp4`
- `diagnostics/*`
- `pipeline_summary.json`
