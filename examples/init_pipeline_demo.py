"""
Demonstration of the new layered initialization pipeline.

This example shows how to initialize a photon wave packet using the clean
initialization pipeline with full debug output.
"""

import torch
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.init import PhotonSpec, initialize_state_from_spec

# Physical parameters (example values)
spacing = 1e-15  # 1 fm
wave_speed = 3e8  # speed of light

# Create 1D grid
grid_shape = (256,)
state = BraneState(grid_shape, Dimensionality.ONE_D, dtype=torch.float64)
state.initialize_flat_configuration(spacing)

grid = BraneGrid(grid_shape, Dimensionality.ONE_D, spacing, state.device)

# Define photon spec
k_vec = torch.tensor([1e15], dtype=torch.float64)  # ~1/fm wave number
center = torch.tensor([128 * spacing], dtype=torch.float64)
sigma = 10 * spacing  # envelope width

spec = PhotonSpec(
    intrinsic_dim=1,
    center=center,
    sigma=sigma,
    amplitude=1e-16,  # 0.1 fm amplitude
    k_vector=k_vec,
    helicity="R",
    prefer_shear=True,
    velocity_init="directional_derivative",
)

print("\n" + "="*70)
print("PHOTON INITIALIZATION (New Pipeline)")
print("="*70)
print(f"Grid: {grid_shape} points, spacing={spacing:.3e} m")
print(f"Center: {center.item():.3e} m")
print(f"Sigma: {sigma:.3e} m")
print(f"k: {k_vec[0].item():.3e} rad/m")
print(f"λ: {2*torch.pi/k_vec[0].item():.3e} m")
print(f"Velocity method: {spec.velocity_init}")
print("="*70)

# Initialize state from spec
artifact = initialize_state_from_spec(
    state,
    grid,
    spec,
    wave_speed=wave_speed,
    debug_out_dir="./debug_output",
    debug_tag="photon_1d",
)

# Print summary
print("\nInitialization complete!")
print(f"Max displacement: {artifact.layer2.meta['max_displacement']:.3e} m")
print(f"Max velocity: {artifact.layer2.meta['max_velocity']:.3e} m/s")
print(f"\nDebug output saved to: ./debug_output/")
print("  - photon_1d_layer0_rest_geometry.png")
print("  - photon_1d_layer1_spec.png")
print("  - photon_1d_layer2_carrier.png")

# Validation checks
print("\n" + "="*70)
print("VALIDATION CHECKS")
print("="*70)

# Check 1: Rest geometry preserved
u_computed = state.positions - state.rest_positions
print(f"1. Rest geometry check:")
print(f"   max|positions - rest_positions - u0|: {torch.abs(u_computed - artifact.layer2.u0).max().item():.3e}")

# Check 2: Polarization basis correctness (for 1D, should be in embedding space)
from branesim.init.polarization import validate_polarization_basis
p1_valid = validate_polarization_basis(
    artifact.layer2.p1,
    artifact.layer2.p2,
    artifact.layer1.k_hat,
    artifact.layer0.intrinsic_dim,
)
print(f"2. Polarization basis check:")
print(f"   Orthonormal: {p1_valid['valid']}")
print(f"   |p1|: {p1_valid['p1_norm']:.6f}")
print(f"   |p2|: {p1_valid['p2_norm']:.6f}")
print(f"   p1·p2: {p1_valid['dot_p1_p2']:.6e}")

# Check 3: Energy distribution
E_total = torch.sum(state.positions ** 2).item()
E_lateral = torch.sum(state.positions[:, :1] ** 2).item()  # 1D: only x
E_amplitude = torch.sum(state.positions[:, 1:] ** 2).item()  # Y, Z, W
print(f"3. Energy distribution:")
print(f"   Total: {E_total:.3e}")
print(f"   Lateral (X): {E_lateral:.3e} ({100*E_lateral/E_total:.2f}%)")
print(f"   Amplitude (Y,Z,W): {E_amplitude:.3e} ({100*E_amplitude/E_total:.2f}%)")

print("\n" + "="*70)
print("Demo complete! Check debug_output/ for visualizations.")
print("="*70)