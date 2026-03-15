# 4 Independent Components

Each component is independently callable and communicates only via files.
A component must not import code from another component.

## Shared Scope

`components/shared/` is intentionally minimal and only contains:
- I/O contracts (`io.py`)
- Core brane state datastructure (`state.py`)
- Small common helpers (`utils.py`)

No simulation mechanics or diagnostics math live in shared.

## 1) Initialization

```bash
python -m components.initialization.run --help
```

Input: CLI parameters only.
Output: `initial_state.npz`.

## 2) Simulation

```bash
python -m components.simulation.run --help
```

Input: `initial_state.npz`.
Output: `trajectory.zip`.

## 3) Visualization

```bash
python -m components.visualization.run --help
```

Input: `trajectory.zip`.
Output: video file.

## 4) Diagnostics

```bash
python -m components.diagnostics.run --help
```

Input: `trajectory.zip`.
Output: diagnostics files.
