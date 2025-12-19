"""
Acceptance tests for the initialization pipeline.

These tests enforce the non-negotiable principles:
1. Never write carrier values into absolute coordinates
2. Momentum is part of the spec (top-down)
3. Diagnostics don't require hardcoded ω
4. Every layer produces artifacts
"""

import torch
from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.init import PhotonSpec, initialize_state_from_spec
from branesim.init.polarization import validate_polarization_basis


def simple_1d_setup():
    """Create a simple 1D setup for testing."""
    grid_shape = (64,)
    spacing = 1e-15
    device = torch.device('cpu')
    dtype = torch.float64

    state = BraneState(grid_shape, Dimensionality.ONE_D, device, dtype)
    state.initialize_flat_configuration(spacing)

    grid = BraneGrid(grid_shape, Dimensionality.ONE_D, spacing, device)

    return state, grid, spacing


def test_rest_geometry_preserved(simple_1d_setup):
    """
    Acceptance Check 1: Rest geometry preserved.

    After initialization, state.positions[:, :d] must equal
    state.rest_positions[:, :d] + u0[:, :d] exactly.
    """
    state, grid, spacing = simple_1d_setup

    k_vec = torch.tensor([1e15], dtype=state.dtype)
    center = torch.tensor([32 * spacing], dtype=state.dtype)

    spec = PhotonSpec(
        intrinsic_dim=1,
        center=center,
        sigma=5 * spacing,
        amplitude=1e-16,
        k_vector=k_vec,
        helicity="R",
        velocity_init="directional_derivative",
    )

    artifact = initialize_state_from_spec(
        state, grid, spec, wave_speed=3e8
    )

    # Check: positions = rest_positions + u0
    u_computed = state.positions - state.rest_positions
    error = torch.abs(u_computed - artifact.layer2.u0).max().item()

    assert error < 1e-14, f"Rest geometry not preserved: max error = {error}"

    # Check: lateral coordinates not directly overwritten for pure shear (1D)
    # In 1D with prefer_shear, polarization should be in embedding axes 1,2
    # So lateral (axis 0) should remain at rest
    lateral_displacement = state.positions[:, 0] - state.rest_positions[:, 0]
    assert torch.abs(lateral_displacement).max().item() < 1e-14, \
        "Lateral coordinates were modified (violates displacement-only principle)"


def test_polarization_basis_correctness(simple_1d_setup):
    """
    Acceptance Check 2: Polarization basis correctness.

    For photons:
    - p1, p2 must be orthonormal
    - In 3D with P1: p1·k̂ ≈ 0, p2·k̂ ≈ 0
    """
    state, grid, spacing = simple_1d_setup

    # Test 1D case
    k_vec = torch.tensor([1e15], dtype=state.dtype)
    center = torch.tensor([32 * spacing], dtype=state.dtype)

    spec = PhotonSpec(
        intrinsic_dim=1,
        center=center,
        sigma=5 * spacing,
        amplitude=1e-16,
        k_vector=k_vec,
        helicity="R",
        prefer_shear=True,
        velocity_init="directional_derivative",
    )

    artifact = initialize_state_from_spec(
        state, grid, spec, wave_speed=3e8
    )

    validation = validate_polarization_basis(
        artifact.layer2.p1,
        artifact.layer2.p2,
        artifact.layer1.k_hat,
        artifact.layer0.intrinsic_dim,
    )

    assert validation['valid'], f"Polarization basis invalid: {validation}"
    assert abs(validation['p1_norm'] - 1.0) < 1e-10
    assert abs(validation['p2_norm'] - 1.0) < 1e-10
    assert abs(validation['dot_p1_p2']) < 1e-10


def test_no_direct_position_writes():
    """
    Acceptance Check 3: No direct position writes.

    The initialization code must NEVER write directly to state.positions.
    It must only use state.set_kinematics(u, v).
    """
    # This is enforced by the API design:
    # - The compiler only produces u0, v0
    # - initialize_state_from_spec calls state.set_kinematics(u0, v0)
    # - set_kinematics internally does: positions = rest_positions + u

    # We can check that rest_positions is never modified
    grid_shape = (64,)
    spacing = 1e-15

    state = BraneState(grid_shape, Dimensionality.ONE_D, dtype=torch.float64)
    state.initialize_flat_configuration(spacing)

    # Save rest positions
    rest_copy = state.rest_positions.clone()

    grid = BraneGrid(grid_shape, Dimensionality.ONE_D, spacing, state.device)

    k_vec = torch.tensor([1e15], dtype=state.dtype)
    center = torch.tensor([32 * spacing], dtype=state.dtype)

    spec = PhotonSpec(
        intrinsic_dim=1,
        center=center,
        sigma=5 * spacing,
        amplitude=1e-16,
        k_vector=k_vec,
        helicity="R",
        velocity_init="directional_derivative",
    )

    initialize_state_from_spec(state, grid, spec, wave_speed=3e8)

    # Check rest_positions unchanged
    assert torch.all(state.rest_positions == rest_copy), \
        "rest_positions was modified (violates immutability principle)"


def test_artifacts_produced():
    """
    Acceptance Check 4: Every layer produces artifacts.

    The pipeline must return InitPipelineArtifact with all three layers.
    """
    grid_shape = (64,)
    spacing = 1e-15

    state = BraneState(grid_shape, Dimensionality.ONE_D, dtype=torch.float64)
    state.initialize_flat_configuration(spacing)

    grid = BraneGrid(grid_shape, Dimensionality.ONE_D, spacing, state.device)

    k_vec = torch.tensor([1e15], dtype=state.dtype)
    center = torch.tensor([32 * spacing], dtype=state.dtype)

    spec = PhotonSpec(
        intrinsic_dim=1,
        center=center,
        sigma=5 * spacing,
        amplitude=1e-16,
        k_vector=k_vec,
        helicity="R",
        velocity_init="directional_derivative",
    )

    artifact = initialize_state_from_spec(state, grid, spec, wave_speed=3e8)

    # Check all layers present
    assert artifact.layer0 is not None
    assert artifact.layer1 is not None
    assert artifact.layer2 is not None

    # Check Layer 0
    assert artifact.layer0.intrinsic_dim == 1
    assert artifact.layer0.rest_positions is not None
    assert artifact.layer0.coords.shape[0] == 64

    # Check Layer 1
    assert artifact.layer1.kind == "photon"
    assert artifact.layer1.k_mag > 0

    # Check Layer 2
    assert artifact.layer2.u0 is not None
    assert artifact.layer2.v0 is not None
    assert artifact.layer2.psi is not None
    assert 'max_displacement' in artifact.layer2.meta
    assert 'max_velocity' in artifact.layer2.meta


def test_momentum_integrated_not_separate_layer():
    """
    Acceptance Check 5: Momentum is part of Layer 2 compilation.

    There is no separate "Layer 3" for momentum. Velocities are computed
    as part of carrier compilation using the velocity_init method.
    """
    grid_shape = (64,)
    spacing = 1e-15

    state = BraneState(grid_shape, Dimensionality.ONE_D, dtype=torch.float64)
    state.initialize_flat_configuration(spacing)

    grid = BraneGrid(grid_shape, Dimensionality.ONE_D, spacing, state.device)

    k_vec = torch.tensor([1e15], dtype=state.dtype)
    center = torch.tensor([32 * spacing], dtype=state.dtype)

    spec = PhotonSpec(
        intrinsic_dim=1,
        center=center,
        sigma=5 * spacing,
        amplitude=1e-16,
        k_vector=k_vec,
        helicity="R",
        velocity_init="directional_derivative",
    )

    artifact = initialize_state_from_spec(state, grid, spec, wave_speed=3e8)

    # Check that v0 was computed and applied in one step
    assert torch.all(state.velocities == artifact.layer2.v0), \
        "Velocities not applied correctly from Layer 2"

    # Check that velocity metadata is in Layer 2
    assert 'velocity_method' in artifact.layer2.meta
    assert artifact.layer2.meta['velocity_method'] == 'directional_derivative'


if __name__ == "__main__":
    # Run tests manually
    print("Running acceptance tests...")

    setup = simple_1d_setup()

    print("\n1. Testing rest geometry preservation...")
    test_rest_geometry_preserved(setup)
    print("   ✓ PASSED")

    print("\n2. Testing polarization basis correctness...")
    test_polarization_basis_correctness(setup)
    print("   ✓ PASSED")

    print("\n3. Testing no direct position writes...")
    test_no_direct_position_writes()
    print("   ✓ PASSED")

    print("\n4. Testing artifacts produced...")
    test_artifacts_produced()
    print("   ✓ PASSED")

    print("\n5. Testing momentum integration...")
    test_momentum_integrated_not_separate_layer()
    print("   ✓ PASSED")

    print("\n" + "="*70)
    print("ALL ACCEPTANCE TESTS PASSED")
    print("="*70)