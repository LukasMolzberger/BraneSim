"""
Test suite for four-potential EM solver.

This module tests the plane wave propagation in vacuum to validate the
solver implementation, checking that:
    1. Plane waves propagate without dispersion
    2. E and B field magnitudes satisfy |E| = c|B|
    3. Energy is conserved (approximately)
    4. CFL condition is satisfied
"""

import torch
from branesim.core.state import Dimensionality
from branesim.core.grid import BraneGrid
from branesim.em.em_state import EMState
from branesim.em.potential_solver import FourPotentialVerletSolver
from branesim.em.initial_conditions_em import initialize_plane_wave_Ay


def test_plane_wave_propagates():
    """Test that a plane wave propagates correctly in vacuum."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dtype = torch.float64

    # Physical constants
    c = 299_792_458.0  # m/s
    mu0 = 4e-7 * 3.141592653589793  # H/m

    # Grid setup (1D)
    nx = 128
    h = 1e-3  # 1 mm spacing
    grid = BraneGrid((nx,), Dimensionality.ONE_D, spacing=h, device=device)

    # Create state
    state = EMState((nx,), Dimensionality.ONE_D, device=device, dtype=dtype)

    # Time step (CFL safe)
    dt = 0.4 * h / c
    solver = FourPotentialVerletSolver(
        dt=dt, c=c, mu0=mu0, grid=grid, bc="periodic", enforce_lorenz=False
    )

    # Initialize plane wave: A_y with wavelength = L/4
    wavelength = nx * h / 4
    amplitude = 1e-6
    initialize_plane_wave_Ay(state, grid, c=c, amplitude=amplitude, wavelength=wavelength)

    # Check CFL condition
    assert solver.cfl_ok(), "CFL condition violated"

    # Initialize accelerations
    solver.initialize_accelerations(state)

    # Compute initial fields
    E0, B0 = solver.compute_fields(state)
    E0_norm = torch.norm(E0).item()
    B0_norm = torch.norm(B0).item()

    # For plane wave: |E| = c|B|
    ratio0 = c * B0_norm / (E0_norm + 1e-30)
    print(f"Initial |E| = {E0_norm:.6e} V/m")
    print(f"Initial |B| = {B0_norm:.6e} T")
    print(f"Initial c|B|/|E| = {ratio0:.6f} (should be ~1)")

    # Evolve for several periods
    steps = 200
    for _ in range(steps):
        solver.step(state)

    # Compute final fields
    E1, B1 = solver.compute_fields(state)
    E1_norm = torch.norm(E1).item()
    B1_norm = torch.norm(B1).item()

    ratio1 = c * B1_norm / (E1_norm + 1e-30)
    print(f"\nAfter {steps} steps:")
    print(f"Final |E| = {E1_norm:.6e} V/m")
    print(f"Final |B| = {B1_norm:.6e} T")
    print(f"Final c|B|/|E| = {ratio1:.6f} (should be ~1)")

    # Check that field ratio is approximately 1 (vacuum plane wave)
    assert abs(ratio0 - 1.0) < 5e-2, f"Initial field ratio {ratio0} not close to 1"
    assert abs(ratio1 - 1.0) < 5e-2, f"Final field ratio {ratio1} not close to 1"

    # Check energy conservation (field magnitudes should be similar)
    energy_change = abs(E1_norm - E0_norm) / (E0_norm + 1e-30)
    print(f"Energy change: {energy_change:.2%}")
    assert energy_change < 0.1, f"Energy changed by {energy_change:.2%} (>10%)"


def test_em_state_basics():
    """Test basic EMState functionality."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Create 2D state
    state = EMState((10, 10), Dimensionality.TWO_D, device=device)

    assert state.num_points == 100
    assert state.potential.shape == (100, 4)
    assert state.velocity.shape == (100, 4)

    # Test view/flatten
    pot_grid = state.view_grid(state.potential)
    assert pot_grid.shape == (10, 10, 4)

    pot_flat = state.flatten_grid(pot_grid)
    assert pot_flat.shape == (100, 4)
    assert torch.allclose(pot_flat, state.potential)


def test_lorenz_gauge_enforcement():
    """Test Lorenz gauge constraint enforcement."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dtype = torch.float64

    c = 299_792_458.0
    mu0 = 4e-7 * 3.141592653589793

    # 1D grid
    nx = 64
    h = 1e-3
    grid = BraneGrid((nx,), Dimensionality.ONE_D, spacing=h, device=device)
    state = EMState((nx,), Dimensionality.ONE_D, device=device, dtype=dtype)

    # Initialize with non-gauge-satisfying random fields
    state.potential = torch.randn_like(state.potential) * 1e-6
    state.velocity = torch.randn_like(state.velocity) * 1e-3

    # Enforce gauge
    residual = state.enforce_lorenz_gauge(grid, c, bc="periodic")

    # Check that residual is small after enforcement
    residual_norm = torch.norm(residual).item()
    print(f"Lorenz gauge residual after enforcement: {residual_norm:.6e}")
    assert residual_norm < 1e-10, f"Gauge residual {residual_norm} too large"


def test_cfl_condition():
    """Test CFL condition checking."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    c = 299_792_458.0
    mu0 = 4e-7 * 3.141592653589793

    nx = 100
    h = 1e-3
    grid = BraneGrid((nx,), Dimensionality.ONE_D, spacing=h, device=device)

    # Safe time step
    dt_safe = 0.4 * h / c
    solver_safe = FourPotentialVerletSolver(dt=dt_safe, c=c, mu0=mu0, grid=grid)
    assert solver_safe.cfl_ok(), "Safe CFL should pass"

    # Unsafe time step
    dt_unsafe = 2.0 * h / c
    solver_unsafe = FourPotentialVerletSolver(dt=dt_unsafe, c=c, mu0=mu0, grid=grid)
    assert not solver_unsafe.cfl_ok(), "Unsafe CFL should fail"


if __name__ == "__main__":
    print("Running EM solver tests...\n")
    print("=" * 60)
    print("Test 1: Plane wave propagation")
    print("=" * 60)
    test_plane_wave_propagates()
    print("\n" + "=" * 60)
    print("Test 2: EMState basics")
    print("=" * 60)
    test_em_state_basics()
    print("Passed!\n")
    print("=" * 60)
    print("Test 3: Lorenz gauge enforcement")
    print("=" * 60)
    test_lorenz_gauge_enforcement()
    print("\n" + "=" * 60)
    print("Test 4: CFL condition")
    print("=" * 60)
    test_cfl_condition()
    print("\n" + "=" * 60)
    print("All tests passed!")
    print("=" * 60)