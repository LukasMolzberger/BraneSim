"""Debug script to check lateralization measurement."""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import torch
import numpy as np

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.physics.forces import SpringForceComputer
from branesim.diagnostics.lateralization import (
    LateralizationMeasurement,
    LateralizationConfig,
)

def debug_lateralization_1d():
    """Debug lateralization measurement in 1D."""

    # Simple setup
    nx = 20
    h = 1.0
    device = torch.device('cpu')
    dtype = torch.float64

    # Create state and grid
    state = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(h)

    # Store reference
    initial_positions = state.positions.clone()

    grid = BraneGrid((nx,), Dimensionality.ONE_D, h, device)

    print("=== Setup ===")
    print(f"Grid shape: {nx}")
    print(f"Spacing: {h}")
    print(f"Device: {device}")
    print(f"Grid neighbors: {grid.neighbors}")
    print(f"Grid neighbors is None: {grid.neighbors is None}")

    if grid.neighbors is not None:
        print(f"Neighbors shape: {len(grid.neighbors)}")
        print(f"First few neighbors: {grid.neighbors[:3]}")

    # Add a simple wave
    x = torch.arange(nx, device=device, dtype=dtype) * h
    wavelength = 4 * h
    k = 2 * np.pi / wavelength
    amplitude = 0.1 * h
    center = nx * h / 2

    # Just positions, no velocities
    envelope = amplitude * torch.exp(-((x - center)**2) / (2 * wavelength**2))
    state.positions[:, 3] = envelope * torch.cos(k * (x - center))

    print(f"\n=== Wave Setup ===")
    print(f"Wavelength: {wavelength}")
    print(f"Amplitude: {amplitude}")
    print(f"Center: {center}")
    print(f"Max amplitude displacement: {state.positions[:, 3].abs().max().item():.6e}")

    # Create physics
    k_spring = 1.0
    rest_length = 0.0
    physics = SpringForceComputer(k_spring, rest_length)

    print(f"\n=== Physics ===")
    print(f"Spring constant: {k_spring}")
    print(f"Rest length: {rest_length}")

    # Create lateralization measurement
    lat_config = LateralizationConfig(
        amplitude_dim=3,
        lateral_dims=(0,),
    )

    lateralization = LateralizationMeasurement(
        config=lat_config,
        grid=grid,
        m_point=1.0,
        reference_positions=initial_positions,
    )

    print(f"\n=== Lateralization Config ===")
    print(f"Amplitude dim: {lat_config.amplitude_dim}")
    print(f"Lateral dims: {lat_config.lateral_dims}")
    print(f"Reference positions shape: {lateralization.reference_positions.shape}")
    print(f"Reference positions sample:")
    for i in range(min(3, nx)):
        print(f"  Point {i}: {lateralization.reference_positions[i].numpy()}")

    # Measure
    print(f"\n=== Measuring ===")
    R_lat_local, R_lat_global, diagnostics = lateralization.measure(state, physics)

    print(f"\nResults:")
    print(f"  Global R_lat: {R_lat_global:.6f}")
    print(f"  Local R_lat range: [{R_lat_local.min().item():.6f}, {R_lat_local.max().item():.6f}]")
    print(f"  Local R_lat mean: {R_lat_local.mean().item():.6f}")

    print(f"\nEnergy totals:")
    print(f"  E_amp_kin total: {diagnostics['E_amp_kin'].sum().item():.6e}")
    print(f"  E_lat_kin total: {diagnostics['E_lat_kin'].sum().item():.6e}")
    print(f"  E_amp_pot total: {diagnostics['E_amp_pot'].sum().item():.6e}")
    print(f"  E_lat_pot total: {diagnostics['E_lat_pot'].sum().item():.6e}")

    print(f"\nSample point energies (point 10):")
    i = 10
    print(f"  E_amp_kin[{i}]: {diagnostics['E_amp_kin'][i].item():.6e}")
    print(f"  E_lat_kin[{i}]: {diagnostics['E_lat_kin'][i].item():.6e}")
    print(f"  E_amp_pot[{i}]: {diagnostics['E_amp_pot'][i].item():.6e}")
    print(f"  E_lat_pot[{i}]: {diagnostics['E_lat_pot'][i].item():.6e}")
    print(f"  R_lat[{i}]: {R_lat_local[i].item():.6f}")

if __name__ == '__main__':
    debug_lateralization_1d()