# Baryon Pipeline (Modular Refactor)

This pipeline splits the baryon simulation workflow into independent components:

1. **Initialization** (`experiments/baryon_init_component.py`)
2. **Simulation** (`experiments/baryon_simulation_component.py`)
3. **Visualization** (`experiments/baryon_visualization_component.py`)
4. **Diagnostics** (`experiments/baryon_diagnostics_component.py`)

All handoff is file-based and compressed.

## File Flow

- Initialization output: compressed `.npz` initial-state package
- Simulation output: compressed `.zip` trajectory with frame checkpoints
- Visualization input: trajectory `.zip`
- Diagnostics input: trajectory `.zip`

## Single-Config Orchestration (Recommended)

Use one JSON config to run all components in sequence:

```bash
python experiments/baryon_pipeline_run.py \
  --config experiments/configs/baryon_pipeline.example.json
```

You can override output location:

```bash
python experiments/baryon_pipeline_run.py \
  --config experiments/configs/baryon_pipeline.example.json \
  --output-dir test-runs/my-baryon-run
```

## Example End-to-End (Local)

```bash
python experiments/baryon_init_component.py \
  --output test-runs/baryon/initial_state.npz \
  --nx 64 --ny 64 --nz 64 \
  --radius 10.0 --amplitude 0.25

python experiments/baryon_simulation_component.py \
  --input test-runs/baryon/initial_state.npz \
  --output test-runs/baryon/trajectory.zip \
  --num-steps 300 --checkpoint-interval 1

python experiments/baryon_visualization_component.py \
  --input test-runs/baryon/trajectory.zip \
  --mode volume \
  --output test-runs/baryon/volume.mp4

python experiments/baryon_visualization_component.py \
  --input test-runs/baryon/trajectory.zip \
  --mode slice --plane xy \
  --output test-runs/baryon/slice_xy.mp4

python experiments/baryon_diagnostics_component.py \
  --input test-runs/baryon/trajectory.zip \
  --output-dir test-runs/baryon/diagnostics
```

## Notes

- Initialization uses spherical coordinates + spherical harmonics as the center of the baryon seed.
- Diagnostics include Berry analysis and QCD-inspired metrics (triplet fractions, trace/traceless structure, mixing, confinement leakage).
- Visualization conversion is independent from simulation and can be rerun repeatedly from one stored trajectory.
