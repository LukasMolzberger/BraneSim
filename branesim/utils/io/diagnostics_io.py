"""
Export utilities for diagnostic results.

Provides systematic export of DiagnosticResult objects to:
- NPZ format (NumPy compressed archive)
- CSV format (flattened with coordinate columns)
- HDF5 format (optional, for large datasets)
"""

from __future__ import annotations
from pathlib import Path
import csv
import numpy as np
import torch

from branesim.diagnostics.types import DiagnosticResult, GridSpec


def save_result_npz(
    path: Path | str,
    result: DiagnosticResult,
    grid: GridSpec | None = None
) -> None:
    """
    Save a single DiagnosticResult to NPZ format.

    Parameters
    ----------
    path : Path | str
        Output file path (.npz)
    result : DiagnosticResult
        Diagnostic result to save
    grid : GridSpec | None
        Optional grid specification (saved as metadata)

    Examples
    --------
    >>> result = berry_connection_along_axis(...)
    >>> save_result_npz("output/berry_t10.npz", result, grid)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Convert all tensors to numpy
    arrays_to_save = {}

    # Add data fields
    for key, val in result.data.items():
        if isinstance(val, torch.Tensor):
            arrays_to_save[f"data_{key}"] = val.cpu().numpy()
        elif isinstance(val, np.ndarray):
            arrays_to_save[f"data_{key}"] = val
        else:
            # Store as metadata string
            arrays_to_save[f"meta_data_{key}"] = np.array([str(val)])

    # Add quality fields
    for key, val in result.quality.items():
        if isinstance(val, torch.Tensor):
            arrays_to_save[f"quality_{key}"] = val.cpu().numpy()
        elif isinstance(val, np.ndarray):
            arrays_to_save[f"quality_{key}"] = val

    # Add metadata as strings
    meta_dict = result.meta.copy()
    if grid is not None:
        meta_dict["grid_shape"] = grid.shape
        meta_dict["grid_spacing_sim"] = grid.spacing_sim
        meta_dict["grid_D"] = grid.D

    for key, val in meta_dict.items():
        arrays_to_save[f"meta_{key}"] = np.array([str(val)])

    # Add time information
    if result.t_sim is not None:
        arrays_to_save["t_sim"] = np.array([result.t_sim])
    if result.t_phys_s is not None:
        arrays_to_save["t_phys_s"] = np.array([result.t_phys_s])

    # Add result name
    arrays_to_save["name"] = np.array([result.name])

    # Save
    np.savez_compressed(path, **arrays_to_save)


def save_results_npz_timeseries(
    path: Path | str,
    results: list[DiagnosticResult],
    grid: GridSpec | None = None
) -> None:
    """
    Save a timeseries of DiagnosticResults to a single NPZ file.

    Each field becomes a 2D array with first dimension = time.

    Parameters
    ----------
    path : Path | str
        Output file path (.npz)
    results : list[DiagnosticResult]
        List of results (one per timestep)
    grid : GridSpec | None
        Optional grid specification

    Examples
    --------
    >>> results = pipeline.run(snapshots, grid)["berry_phase_axis0"]
    >>> save_results_npz_timeseries("output/berry_timeseries.npz", results, grid)
    """
    if not results:
        raise ValueError("Empty results list")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Collect all data/quality fields across timesteps
    arrays_to_save = {}
    n_times = len(results)

    # Time arrays
    t_sim_array = np.array([r.t_sim if r.t_sim is not None else np.nan for r in results])
    t_phys_array = np.array([r.t_phys_s if r.t_phys_s is not None else np.nan for r in results])
    arrays_to_save["t_sim"] = t_sim_array
    arrays_to_save["t_phys_s"] = t_phys_array

    # Data fields: stack over time
    data_keys = results[0].data.keys()
    for key in data_keys:
        # Collect arrays
        arrays = []
        for r in results:
            val = r.data[key]
            if isinstance(val, torch.Tensor):
                arrays.append(val.cpu().numpy())
            elif isinstance(val, np.ndarray):
                arrays.append(val)
            else:
                # Can't stack non-arrays
                continue

        if arrays:
            try:
                stacked = np.stack(arrays, axis=0)  # [n_times, ...]
                arrays_to_save[f"data_{key}"] = stacked
            except ValueError:
                # Arrays have inconsistent shapes, save separately
                for i, arr in enumerate(arrays):
                    arrays_to_save[f"data_{key}_t{i}"] = arr

    # Quality fields: stack over time
    quality_keys = results[0].quality.keys()
    for key in quality_keys:
        arrays = []
        for r in results:
            val = r.quality[key]
            if isinstance(val, torch.Tensor):
                arrays.append(val.cpu().numpy())
            elif isinstance(val, np.ndarray):
                arrays.append(val)
            else:
                continue

        if arrays:
            try:
                stacked = np.stack(arrays, axis=0)
                arrays_to_save[f"quality_{key}"] = stacked
            except ValueError:
                for i, arr in enumerate(arrays):
                    arrays_to_save[f"quality_{key}_t{i}"] = arr

    # Metadata (same for all timesteps)
    meta_dict = results[0].meta.copy()
    if grid is not None:
        meta_dict["grid_shape"] = grid.shape
        meta_dict["grid_spacing_sim"] = grid.spacing_sim
        meta_dict["grid_D"] = grid.D

    for key, val in meta_dict.items():
        arrays_to_save[f"meta_{key}"] = np.array([str(val)])

    # Name
    arrays_to_save["name"] = np.array([results[0].name])

    # Save
    np.savez_compressed(path, **arrays_to_save)


def save_result_csv_1d(
    path: Path | str,
    result: DiagnosticResult,
    grid: GridSpec,
    coord_name: str = "x"
) -> None:
    """
    Save 1D DiagnosticResult to CSV format.

    Creates a CSV with columns: [coord, data_field1, data_field2, ...].

    Parameters
    ----------
    path : Path | str
        Output file path (.csv)
    result : DiagnosticResult
        Diagnostic result (1D data)
    grid : GridSpec
        Grid specification (must be 1D)
    coord_name : str
        Name for coordinate column

    Examples
    --------
    >>> result = berry_phase_profile_along_x(...)
    >>> grid = GridSpec(shape=(100,), spacing_sim=1.0)
    >>> save_result_csv_1d("output/berry_profile.csv", result, grid)
    """
    if grid.D != 1:
        raise ValueError("save_result_csv_1d only supports 1D grids")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Generate coordinates if not provided
    if grid.coords_sim is not None:
        x = grid.coords_sim[0].cpu().numpy()
    else:
        x = np.arange(grid.shape[0]) * grid.spacing_sim

    # Build CSV data
    csv_data = {coord_name: x}

    # Add data fields
    for key, val in result.data.items():
        if isinstance(val, torch.Tensor):
            arr = val.cpu().numpy()
        elif isinstance(val, np.ndarray):
            arr = val
        else:
            continue

        # Check if field has same length as grid
        if arr.size == grid.num_points:
            csv_data[f"data_{key}"] = arr.flatten()
        elif arr.size == grid.num_points - 1:
            # Edge quantities (e.g., dphi)
            # Pad with NaN or skip
            csv_data[f"data_{key}"] = np.concatenate([arr.flatten(), [np.nan]])

    # Add quality fields
    for key, val in result.quality.items():
        if isinstance(val, torch.Tensor):
            arr = val.cpu().numpy()
        elif isinstance(val, np.ndarray):
            arr = val
        else:
            continue

        if arr.size == grid.num_points:
            csv_data[f"quality_{key}"] = arr.flatten()
        elif arr.size == grid.num_points - 1:
            csv_data[f"quality_{key}"] = np.concatenate([arr.flatten(), [np.nan]])

    # Write to CSV using built-in csv module
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(csv_data.keys())
        # Write rows
        for i in range(len(x)):
            row = [csv_data[key][i] for key in csv_data.keys()]
            writer.writerow(row)


def save_pipeline_results(
    output_dir: Path | str,
    results: dict[str, list[DiagnosticResult]],
    grid: GridSpec,
    format: str = "npz"
) -> None:
    """
    Save all results from a DiagnosticPipeline.

    Creates one file per diagnostic, containing timeseries data.

    Parameters
    ----------
    output_dir : Path | str
        Output directory
    results : dict[str, list[DiagnosticResult]]
        Results from DiagnosticPipeline.run()
    grid : GridSpec
        Grid specification
    format : str
        Output format: "npz" (default)

    Examples
    --------
    >>> pipeline = DiagnosticPipeline([...])
    >>> results = pipeline.run(snapshots, grid)
    >>> save_pipeline_results("output/diagnostics", results, grid)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for diag_name, diag_results in results.items():
        if format == "npz":
            filename = output_dir / f"{diag_name}.npz"
            save_results_npz_timeseries(filename, diag_results, grid)
        else:
            raise ValueError(f"Unknown format: {format}")

    print(f"Saved {len(results)} diagnostic results to {output_dir}")