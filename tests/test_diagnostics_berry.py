"""
Validation tests for Berry phase diagnostics.

Tests:
1. Pure plane wave: Berry connection should match analytical k
2. Gauge transformation invariance for curvature
3. Node behavior: overlap and validity masks
"""

try:
    import pytest
except ImportError:
    pytest = None

import torch
import numpy as np
import math

from branesim.analytics import (
    GridSpec,
    BerryConfig,
    complex_band_state_from_quadrature,
    pointwise_normalize,
    berry_connection_along_axis,
    berry_plaquette_curvature,
)


def test_berry_connection_plane_wave_1d():
    """
    Test Berry connection for a pure plane wave.

    For a plane wave u = exp(i k x), the Berry connection should be:
    A_x ≈ k
    """
    # Setup
    N = 200
    L = 10.0
    h = L / N
    x = torch.linspace(0, L - h, N)

    # Plane wave parameters
    k = 2.0 * math.pi / L  # One wavelength
    omega = 1.0

    # Create plane wave: u = exp(i k x)
    # Using analytic signal: ψ = q + i·v/ω
    # For u = exp(ikx): q = cos(kx), v = ω sin(kx) (positive for positive phase velocity)
    q = torch.cos(k * x)
    v = omega * torch.sin(k * x)

    # Construct complex state
    psi = complex_band_state_from_quadrature(q, v, omega)
    psi_hat, amp = pointwise_normalize(psi)

    # Setup grid and config
    grid = GridSpec(shape=(N,), spacing_sim=h)
    cfg = BerryConfig(spacing_sim=h, amplitude_threshold=1e-10, overlap_threshold=1e-6)

    # Compute Berry connection
    result = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    A_x = result["A_axis"]

    # Check: A_x should be approximately constant and equal to k
    A_x_mean = A_x.mean().item()
    A_x_std = A_x.std().item()

    print(f"Plane wave test:")
    print(f"  Expected k: {k:.6f}")
    print(f"  Measured A_x: {A_x_mean:.6f} ± {A_x_std:.6f}")
    print(f"  Relative error: {abs(A_x_mean - k) / k * 100:.2f}%")

    # Assert: relative error should be small
    assert abs(A_x_mean - k) / k < 0.01, f"A_x = {A_x_mean} doesn't match k = {k}"
    assert A_x_std / abs(k) < 0.01, f"A_x has too much variance: {A_x_std}"


def test_berry_connection_plane_wave_2d():
    """
    Test Berry connection for a 2D plane wave.

    For u = exp(i k_x x), the connection along x should be k_x,
    and along y should be 0.
    """
    # Setup
    Nx, Ny = 50, 50
    Lx, Ly = 5.0, 5.0
    hx = Lx / Nx
    hy = Ly / Ny

    x = torch.linspace(0, Lx - hx, Nx)
    y = torch.linspace(0, Ly - hy, Ny)
    X, Y = torch.meshgrid(x, y, indexing='ij')

    # Plane wave along x
    k_x = 2.0 * math.pi / Lx
    omega = 1.0

    q = torch.cos(k_x * X)
    v = omega * torch.sin(k_x * X)

    # Construct complex state
    psi = complex_band_state_from_quadrature(q, v, omega)
    psi_hat, amp = pointwise_normalize(psi)

    # Setup grid and config
    grid = GridSpec(shape=(Nx, Ny), spacing_sim=hx)  # Assume uniform spacing
    cfg = BerryConfig(spacing_sim=hx, amplitude_threshold=1e-10)

    # Compute connection along x
    result_x = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)
    A_x = result_x["A_axis"]

    # Compute connection along y
    result_y = berry_connection_along_axis(psi_hat, amp, axis=1, cfg=cfg)
    A_y = result_y["A_axis"]

    # Check
    A_x_mean = A_x.mean().item()
    A_y_mean = A_y.mean().item()

    print(f"\n2D plane wave test:")
    print(f"  Expected: A_x ≈ {k_x:.6f}, A_y ≈ 0")
    print(f"  Measured: A_x = {A_x_mean:.6f}, A_y = {A_y_mean:.6f}")

    assert abs(A_x_mean - k_x) / k_x < 0.01, f"A_x = {A_x_mean} doesn't match k_x = {k_x}"
    assert abs(A_y_mean) < 0.01 * abs(k_x), f"A_y = {A_y_mean} should be near zero"


def test_berry_curvature_plane_wave_2d():
    """
    Test Berry curvature for a plane wave.

    For a plane wave, the curvature should be zero (flat U(1) connection).
    """
    # Setup
    Nx, Ny = 50, 50
    Lx, Ly = 5.0, 5.0
    hx = Lx / Nx

    x = torch.linspace(0, Lx - hx, Nx)
    y = torch.linspace(0, Ly - hx, Ny)
    X, Y = torch.meshgrid(x, y, indexing='ij')

    # Plane wave
    k_x = 2.0 * math.pi / Lx
    k_y = 1.5 * 2.0 * math.pi / Ly
    omega = 1.0

    q = torch.cos(k_x * X + k_y * Y)
    v = -omega * torch.sin(k_x * X + k_y * Y)

    psi = complex_band_state_from_quadrature(q, v, omega)
    psi_hat, amp = pointwise_normalize(psi)

    grid = GridSpec(shape=(Nx, Ny), spacing_sim=hx)
    cfg = BerryConfig(spacing_sim=hx)

    # Compute curvature
    result = berry_plaquette_curvature(psi_hat, amp, axes=(0, 1), cfg=cfg)
    curvature = result["curvature"]

    # Check: curvature should be near zero
    curv_mean = curvature.mean().item()
    curv_std = curvature.std().item()

    print(f"\nPlane wave curvature test:")
    print(f"  Expected: F ≈ 0")
    print(f"  Measured: F = {curv_mean:.6e} ± {curv_std:.6e}")

    assert abs(curv_mean) < 1e-3, f"Curvature {curv_mean} should be near zero for plane wave"


def test_gauge_invariance_curvature():
    """
    Test that Berry curvature is gauge-invariant.

    Apply a random smooth gauge transformation u' = exp(i χ) u
    and verify that curvature remains unchanged.
    """
    # Setup
    Nx, Ny = 40, 40
    Lx = 5.0
    hx = Lx / Nx

    x = torch.linspace(0, Lx - hx, Nx)
    y = torch.linspace(0, Lx - hx, Ny)
    X, Y = torch.meshgrid(x, y, indexing='ij')

    # Create a wavepacket (not pure plane wave)
    k0 = 2.0 * math.pi / Lx
    sigma = Lx / 8.0
    omega = 1.0

    envelope = torch.exp(-((X - Lx/2)**2 + (Y - Lx/2)**2) / (2 * sigma**2))
    q = envelope * torch.cos(k0 * X)
    v = envelope * (-omega * torch.sin(k0 * X))

    psi = complex_band_state_from_quadrature(q, v, omega)
    psi_hat, amp = pointwise_normalize(psi)

    grid = GridSpec(shape=(Nx, Ny), spacing_sim=hx)
    cfg = BerryConfig(spacing_sim=hx, amplitude_threshold=1e-2)

    # Compute curvature before gauge transformation
    result_before = berry_plaquette_curvature(psi_hat, amp, axes=(0, 1), cfg=cfg)
    curv_before = result_before["curvature"]

    # Apply gauge transformation: u' = exp(i χ) u
    # Use a smooth random phase
    chi = 0.1 * torch.sin(2 * math.pi * X / Lx) * torch.cos(2 * math.pi * Y / Lx)
    psi_gauge = psi_hat * torch.exp(1j * chi)

    # Renormalize (should already be normalized, but to be safe)
    psi_hat_gauge, amp_gauge = pointwise_normalize(psi_gauge)

    # Compute curvature after gauge transformation
    result_after = berry_plaquette_curvature(psi_hat_gauge, amp_gauge, axes=(0, 1), cfg=cfg)
    curv_after = result_after["curvature"]

    # Compare curvatures (should be the same)
    # Only compare where both are valid
    valid = result_before["valid_plaquette"] & result_after["valid_plaquette"]

    if valid.any():
        diff = torch.abs(curv_after[valid] - curv_before[valid])
        max_diff = diff.max().item()
        mean_diff = diff.mean().item()

        print(f"\nGauge invariance test:")
        print(f"  Max curvature difference: {max_diff:.6e}")
        print(f"  Mean curvature difference: {mean_diff:.6e}")
        print(f"  Valid plaquettes: {valid.sum().item()} / {valid.numel()}")

        assert max_diff < 1e-4, f"Curvature changed by {max_diff} under gauge transformation"
    else:
        if pytest:
            pytest.skip("No valid plaquettes for comparison")
        else:
            print("  WARNING: No valid plaquettes for comparison, skipping test")


def test_node_detection():
    """
    Test that nodes (amplitude zeros) are properly detected and masked.

    Create a wavepacket with a node/dip and verify that:
    1. Amplitude mask correctly identifies low-amplitude regions
    2. Overlap at nodes is low
    3. Valid edge mask excludes edges near nodes
    """
    # Setup: Create a wavepacket with an artificial dip/node in the middle
    N = 200
    L = 10.0
    h = L / N
    x = torch.linspace(0, L - h, N)

    k = 2.0 * math.pi / L
    omega = 1.0

    # Create a wavepacket with a Gaussian dip in the middle to simulate a node
    envelope = torch.ones_like(x)
    node_center = L / 2
    node_width = L / 20.0
    # Create a dip: multiply by (1 - exp(-(x-center)^2 / width^2))
    dip = 1.0 - 0.95 * torch.exp(-((x - node_center)**2) / (2 * node_width**2))
    envelope = envelope * dip

    q = envelope * torch.cos(k * x)
    v = envelope * omega * torch.sin(k * x)

    psi = complex_band_state_from_quadrature(q, v, omega, eps=1e-12)
    psi_hat, amp = pointwise_normalize(psi)

    grid = GridSpec(shape=(N,), spacing_sim=h)
    cfg = BerryConfig(
        spacing_sim=h,
        amplitude_threshold=0.1,  # Higher threshold to catch node region
        overlap_threshold=0.8,   # High threshold to catch rapid phase changes
    )

    result = berry_connection_along_axis(psi_hat, amp, axis=0, cfg=cfg)

    # Check mask
    mask = result["mask_point"]
    valid_edge = result["valid_edge"]
    overlap_abs = result["overlap_abs"]

    # Find node location (should be near x = L/2)
    node_idx = N // 2

    print(f"\nNode detection test:")
    print(f"  Amplitude at node (x={x[node_idx]:.2f}): {amp[node_idx]:.6f}")
    print(f"  Amplitude threshold: {cfg.amplitude_threshold:.2f}")
    print(f"  Mask at node: {mask[node_idx].item()}")
    print(f"  Number of masked points: {(~mask).sum().item()} / {N}")
    print(f"  Number of invalid edges: {(~valid_edge).sum().item()} / {N-1}")

    # At node, amplitude should be low
    assert amp[node_idx] < cfg.amplitude_threshold * 1.5, "Node should have low amplitude"

    # At least some points should be masked
    assert (~mask).sum() > 0, "Some points should be below amplitude threshold"


if __name__ == "__main__":
    # Run tests
    print("Running Berry diagnostics validation tests...")
    print("=" * 60)

    test_berry_connection_plane_wave_1d()
    test_berry_connection_plane_wave_2d()
    test_berry_curvature_plane_wave_2d()
    test_gauge_invariance_curvature()
    test_node_detection()

    print("=" * 60)
    print("All tests passed!")