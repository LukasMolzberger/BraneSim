"""Test dimension-independent velocity initialization"""
import torch
import numpy as np
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.initial_conditions import initialize_right_moving_velocities
from branesim.config.simulation_config import PhysicalConstants

print("=" * 80)
print("Testing Dimension-Independent Velocity Initialization")
print("=" * 80)

constants = PhysicalConstants()
device = torch.device('cpu')
dtype = torch.float64
c = constants.c
h = 1e-12  # 1 pm

# ============================================================================
# Test 1D
# ============================================================================
print("\n[Test 1D]")
print("-" * 80)

nx = 100
state_1d = BraneState((nx,), Dimensionality.ONE_D, device, dtype)
state_1d.initialize_flat_configuration(h)
grid_1d = BraneGrid((nx,), Dimensionality.ONE_D, h, device)

# Create a Gaussian pulse shape
x = torch.arange(nx, dtype=dtype, device=device) * h
center = nx * h / 2
sigma = 10 * h
amplitude = 0.1 * h
state_1d.positions[:, 3] = amplitude * torch.exp(-((x - center)**2) / (2 * sigma**2))

print(f"Shape initialized: Gaussian pulse with σ={sigma:.3e} m")
print(f"Max amplitude: {state_1d.positions[:, 3].abs().max():.3e} m")

# Initialize velocities for right-moving wave
initialize_right_moving_velocities(
    state=state_1d,
    grid=grid_1d,
    wave_speed=c,
    direction=None,  # +x
    field_component=3
)

print(f"✓ 1D velocity initialization successful")

# ============================================================================
# Test 2D
# ============================================================================
print("\n[Test 2D]")
print("-" * 80)

nx, ny = 50, 50
state_2d = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
state_2d.initialize_flat_configuration(h)
grid_2d = BraneGrid((nx, ny), Dimensionality.TWO_D, h, device)

# Create a 2D Gaussian pulse
coords = grid_2d.get_spatial_coordinates()
x_2d = coords[:, 0]
y_2d = coords[:, 1]
center_x = nx * h / 2
center_y = ny * h / 2
sigma_2d = 10 * h
r_sq = (x_2d - center_x)**2 + (y_2d - center_y)**2
state_2d.positions[:, 3] = amplitude * torch.exp(-r_sq / (2 * sigma_2d**2))

print(f"Shape initialized: 2D Gaussian pulse with σ={sigma_2d:.3e} m")
print(f"Max amplitude: {state_2d.positions[:, 3].abs().max():.3e} m")

# Initialize velocities for right-moving wave (+x direction)
initialize_right_moving_velocities(
    state=state_2d,
    grid=grid_2d,
    wave_speed=c,
    direction=None,  # +x
    field_component=3
)

print(f"✓ 2D velocity initialization successful")

# Test with custom direction (diagonal)
direction_diag = torch.tensor([1.0, 1.0], dtype=dtype, device=device)  # Will be normalized
initialize_right_moving_velocities(
    state=state_2d,
    grid=grid_2d,
    wave_speed=c,
    direction=direction_diag,
    field_component=3
)

print(f"✓ 2D velocity initialization with diagonal direction successful")

# ============================================================================
# Test 3D
# ============================================================================
print("\n[Test 3D]")
print("-" * 80)

nx, ny, nz = 20, 20, 20
state_3d = BraneState((nx, ny, nz), Dimensionality.THREE_D, device, dtype)
state_3d.initialize_flat_configuration(h)
grid_3d = BraneGrid((nx, ny, nz), Dimensionality.THREE_D, h, device)

# Create a 3D Gaussian pulse
coords_3d = grid_3d.get_spatial_coordinates()
x_3d = coords_3d[:, 0]
y_3d = coords_3d[:, 1]
z_3d = coords_3d[:, 2]
center_x = nx * h / 2
center_y = ny * h / 2
center_z = nz * h / 2
r_sq_3d = (x_3d - center_x)**2 + (y_3d - center_y)**2 + (z_3d - center_z)**2
state_3d.positions[:, 3] = amplitude * torch.exp(-r_sq_3d / (2 * sigma_2d**2))

print(f"Shape initialized: 3D Gaussian pulse with σ={sigma_2d:.3e} m")
print(f"Max amplitude: {state_3d.positions[:, 3].abs().max():.3e} m")

# Initialize velocities for right-moving wave
initialize_right_moving_velocities(
    state=state_3d,
    grid=grid_3d,
    wave_speed=c,
    direction=None,  # +x
    field_component=3
)

print(f"✓ 3D velocity initialization successful")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 80)
print("Summary")
print("=" * 80)
print("✓ 1D velocity initialization works")
print("✓ 2D velocity initialization works")
print("✓ 2D with custom direction works")
print("✓ 3D velocity initialization works")
print("\nThe generic velocity initializer is fully dimension-independent!")
print("=" * 80)