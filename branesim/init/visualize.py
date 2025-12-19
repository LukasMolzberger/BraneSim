"""
Visualization functions for all layers of the initialization pipeline.

Each layer produces diagnostic plots and metadata JSON for debugging.
If you can't plot it, you can't debug it.
"""

import json
import torch
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional

from branesim.init.artifacts import (
    RestGeometryArtifact,
    SpecArtifact,
    CarrierArtifact,
)


def plot_rest_geometry(
    geom: RestGeometryArtifact,
    out_dir: str,
    tag: str,
) -> None:
    """
    Plot Layer 0: rest geometry diagnostics.

    Produces:
    - Scatter plot of spatial coordinates
    - Histogram of neighbor spacing (if applicable)
    - JSON metadata
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    d = geom.intrinsic_dim
    coords_np = geom.coords.cpu().numpy()

    fig = plt.figure(figsize=(12, 4))

    # Subplot 1: Spatial layout
    ax1 = fig.add_subplot(1, 2, 1)
    if d == 1:
        ax1.scatter(coords_np[:, 0], np.zeros_like(coords_np[:, 0]), s=1, alpha=0.5)
        ax1.set_xlabel('X')
        ax1.set_title(f'1D Rest Geometry (N={len(coords_np)})')
    elif d == 2:
        ax1.scatter(coords_np[:, 0], coords_np[:, 1], s=1, alpha=0.5)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title(f'2D Rest Geometry (N={len(coords_np)})')
        ax1.axis('equal')
    else:  # 3D
        from mpl_toolkits.mplot3d import Axes3D
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax1.scatter(coords_np[:, 0], coords_np[:, 1], coords_np[:, 2], s=1, alpha=0.5)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')
        ax1.set_title(f'3D Rest Geometry (N={len(coords_np)})')

    # Subplot 2: Spacing histogram
    ax2 = fig.add_subplot(1, 2, 2)
    # Compute neighbor distances along first axis
    if d >= 1:
        x = coords_np[:, 0]
        diffs = np.diff(np.sort(x))
        diffs = diffs[diffs > 0]  # remove zeros
        if len(diffs) > 0:
            ax2.hist(diffs, bins=50, alpha=0.7, edgecolor='black')
            ax2.set_xlabel('Spacing')
            ax2.set_ylabel('Count')
            ax2.set_title(f'Grid Spacing Histogram\n(mean={np.mean(diffs):.6e})')
            ax2.axvline(geom.spacing, color='r', linestyle='--', label=f'nominal={geom.spacing:.6e}')
            ax2.legend()

    plt.tight_layout()
    plt.savefig(out_path / f'{tag}_layer0_rest_geometry.png', dpi=150)
    plt.close()

    # Save metadata
    meta = {
        'intrinsic_dim': geom.intrinsic_dim,
        'embedding_dim': geom.embedding_dim,
        'grid_shape': geom.grid_shape,
        'spacing': geom.spacing,
        'num_points': coords_np.shape[0],
    }
    with open(out_path / f'{tag}_layer0_metadata.json', 'w') as f:
        json.dump(meta, f, indent=2)


def plot_spec(
    spec_art: SpecArtifact,
    out_dir: str,
    tag: str,
) -> None:
    """
    Plot Layer 1: specification artifact.

    Produces:
    - Text summary of spec parameters
    - Arrows/vectors for k_hat (if applicable)
    - JSON metadata
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('off')

    # Build summary text
    lines = [
        f"Specification: {spec_art.kind.upper()}",
        f"",
        f"Wave vector k:",
        f"  |k| = {spec_art.k_mag:.6e}",
        f"  k̂ = {spec_art.k_hat.cpu().numpy()}",
        f"",
    ]

    if spec_art.kind == "photon":
        spec = spec_art.spec
        lines.extend([
            f"Photon parameters:",
            f"  σ = {spec.sigma:.6e}",
            f"  A = {spec.amplitude:.6e}",
            f"  λ = {2*np.pi/spec_art.k_mag:.6e}",
            f"  helicity = {spec.helicity}",
            f"  velocity_init = {spec.velocity_init}",
        ])
    else:  # electron
        spec = spec_art.spec
        lines.extend([
            f"Electron parameters:",
            f"  R_major = {spec.torus_major_radius:.6e}",
            f"  R_minor = {spec.torus_minor_radius:.6e}",
            f"  (p, q) = ({spec.p}, {spec.q})",
            f"  k_long = {spec.longitudinal_k:.6e}",
            f"  helicity = {spec.helicity}",
        ])

    text = "\n".join(lines)
    ax.text(0.1, 0.5, text, fontsize=10, family='monospace',
            verticalalignment='center')

    plt.tight_layout()
    plt.savefig(out_path / f'{tag}_layer1_spec.png', dpi=150)
    plt.close()

    # Save metadata
    meta = {
        'kind': spec_art.kind,
        'k_mag': spec_art.k_mag,
        'k_hat': spec_art.k_hat.cpu().numpy().tolist(),
        'notes': spec_art.notes,
    }
    with open(out_path / f'{tag}_layer1_metadata.json', 'w') as f:
        json.dump(meta, f, indent=2)


def plot_carrier(
    geom: RestGeometryArtifact,
    spec_art: SpecArtifact,
    carrier: CarrierArtifact,
    out_dir: str,
    tag: str,
) -> None:
    """
    Plot Layer 2: compiled carrier diagnostics.

    Produces:
    - Envelope A(x)
    - Phase φ(x) along propagation axis
    - Polarization plane projection energy
    - Spectrum (FFT along prop axis)
    - JSON metadata
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    d = geom.intrinsic_dim
    coords = geom.coords.cpu().numpy()
    envelope = carrier.envelope.cpu().numpy()
    phase = carrier.phase.cpu().numpy()
    u0 = carrier.u0.cpu().numpy()
    v0 = carrier.v0.cpu().numpy()

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))

    # 1. Envelope plot
    ax1 = fig.add_subplot(2, 3, 1)
    if d == 1:
        ax1.plot(coords[:, 0], envelope, 'b-', linewidth=1)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Amplitude')
        ax1.set_title('Envelope A(x)')
        ax1.grid(True, alpha=0.3)
    elif d == 2:
        # 2D heatmap
        nx, ny = geom.grid_shape
        env_grid = envelope.reshape(nx, ny)
        im = ax1.imshow(env_grid.T, origin='lower', aspect='auto', cmap='viridis')
        plt.colorbar(im, ax=ax1, label='Amplitude')
        ax1.set_title('Envelope A(x, y)')
        ax1.set_xlabel('X index')
        ax1.set_ylabel('Y index')
    else:  # 3D
        # Central slice
        nx, ny, nz = geom.grid_shape
        env_grid = envelope.reshape(nx, ny, nz)
        slice_z = nz // 2
        im = ax1.imshow(env_grid[:, :, slice_z].T, origin='lower', aspect='auto', cmap='viridis')
        plt.colorbar(im, ax=ax1, label='Amplitude')
        ax1.set_title(f'Envelope A(x, y) at z={slice_z}')
        ax1.set_xlabel('X index')
        ax1.set_ylabel('Y index')

    # 2. Phase plot along propagation direction
    ax2 = fig.add_subplot(2, 3, 2)
    k_hat = spec_art.k_hat.cpu().numpy()
    # Project coordinates along k_hat
    x_proj = coords @ k_hat
    sort_idx = np.argsort(x_proj)
    x_sorted = x_proj[sort_idx]
    phase_sorted = phase[sort_idx]

    # Unwrap phase
    phase_unwrapped = np.unwrap(phase_sorted)

    ax2.plot(x_sorted, phase_unwrapped, 'r-', linewidth=0.5, alpha=0.7)
    ax2.set_xlabel('Position along k̂')
    ax2.set_ylabel('Phase φ (unwrapped)')
    ax2.set_title('Carrier Phase φ(x)')
    ax2.grid(True, alpha=0.3)

    # 3. Displacement component magnitudes
    ax3 = fig.add_subplot(2, 3, 3)
    u_mags = np.linalg.norm(u0, axis=1)
    if d == 1:
        ax3.plot(coords[:, 0], u_mags, 'g-', linewidth=1)
        ax3.set_xlabel('X')
    elif d == 2:
        u_grid = u_mags.reshape(geom.grid_shape)
        im = ax3.imshow(u_grid.T, origin='lower', aspect='auto', cmap='plasma')
        plt.colorbar(im, ax=ax3, label='|u|')
        ax3.set_xlabel('X index')
        ax3.set_ylabel('Y index')
    else:  # 3D slice
        u_grid = u_mags.reshape(geom.grid_shape)
        slice_z = geom.grid_shape[2] // 2
        im = ax3.imshow(u_grid[:, :, slice_z].T, origin='lower', aspect='auto', cmap='plasma')
        plt.colorbar(im, ax=ax3, label='|u|')
        ax3.set_xlabel('X index')
        ax3.set_ylabel('Y index')
    ax3.set_title('Displacement Magnitude |u|')

    # 4. Polarization plane projection energy
    ax4 = fig.add_subplot(2, 3, 4)
    p1 = carrier.p1.cpu().numpy()
    p2 = carrier.p2.cpu().numpy()

    # For photon, p1/p2 are [4]; for electron they might be per-node
    if len(p1.shape) == 1:
        # Project u0 onto polarization plane
        proj_p1 = (u0 @ p1).reshape(-1, 1) * p1.reshape(1, -1)
        proj_p2 = (u0 @ p2).reshape(-1, 1) * p2.reshape(1, -1)
        E_plane = np.linalg.norm(proj_p1 + proj_p2, axis=1) ** 2
    else:
        # Position-dependent polarization (electron)
        E_plane = np.array([np.dot(u0[i], p1[i])**2 + np.dot(u0[i], p2[i])**2
                           for i in range(len(u0))])

    E_total = np.linalg.norm(u0, axis=1) ** 2
    E_out = E_total - E_plane

    if d == 1:
        ax4.plot(coords[:, 0], E_plane, 'b-', label='In-plane', linewidth=1)
        ax4.plot(coords[:, 0], E_out, 'r-', label='Out-of-plane', linewidth=1)
        ax4.set_xlabel('X')
        ax4.legend()
    else:
        # Show ratio
        ratio = E_plane / (E_total + 1e-20)
        if d == 2:
            ratio_grid = ratio.reshape(geom.grid_shape)
            im = ax4.imshow(ratio_grid.T, origin='lower', aspect='auto', cmap='RdBu_r',
                           vmin=0, vmax=1)
            plt.colorbar(im, ax=ax4, label='Fraction in plane')
        else:
            ratio_grid = ratio.reshape(geom.grid_shape)
            slice_z = geom.grid_shape[2] // 2
            im = ax4.imshow(ratio_grid[:, :, slice_z].T, origin='lower', aspect='auto',
                           cmap='RdBu_r', vmin=0, vmax=1)
            plt.colorbar(im, ax=ax4, label='Fraction in plane')
    ax4.set_title('Polarization Plane Projection')
    ax4.set_ylabel('Energy')

    # 5. Spectrum (FFT along propagation axis)
    ax5 = fig.add_subplot(2, 3, 5)
    # Use first non-zero component
    comp_idx = np.argmax(np.abs(u0).max(axis=0))
    u_comp = u0[sort_idx, comp_idx]

    # FFT
    fft = np.fft.rfft(u_comp)
    fft_mag = np.abs(fft)
    fft_freq = np.fft.rfftfreq(len(u_comp), d=(x_sorted[1] - x_sorted[0] if len(x_sorted) > 1 else 1.0))

    ax5.semilogy(fft_freq, fft_mag, 'b-', linewidth=1)
    ax5.set_xlabel('Wave number k')
    ax5.set_ylabel('|FFT|')
    ax5.set_title(f'Spectrum (component {comp_idx})')
    ax5.grid(True, alpha=0.3)

    # Mark expected k
    if spec_art.k_mag > 0:
        ax5.axvline(spec_art.k_mag, color='r', linestyle='--',
                   label=f'k={spec_art.k_mag:.3e}')
        ax5.legend()

    # 6. Velocity magnitude
    ax6 = fig.add_subplot(2, 3, 6)
    v_mags = np.linalg.norm(v0, axis=1)
    if d == 1:
        ax6.plot(coords[:, 0], v_mags, 'm-', linewidth=1)
        ax6.set_xlabel('X')
        ax6.set_ylabel('|v|')
    elif d == 2:
        v_grid = v_mags.reshape(geom.grid_shape)
        im = ax6.imshow(v_grid.T, origin='lower', aspect='auto', cmap='magma')
        plt.colorbar(im, ax=ax6, label='|v|')
        ax6.set_xlabel('X index')
        ax6.set_ylabel('Y index')
    else:  # 3D slice
        v_grid = v_mags.reshape(geom.grid_shape)
        slice_z = geom.grid_shape[2] // 2
        im = ax6.imshow(v_grid[:, :, slice_z].T, origin='lower', aspect='auto', cmap='magma')
        plt.colorbar(im, ax=ax6, label='|v|')
        ax6.set_xlabel('X index')
        ax6.set_ylabel('Y index')
    ax6.set_title('Velocity Magnitude |v|')

    plt.tight_layout()
    plt.savefig(out_path / f'{tag}_layer2_carrier.png', dpi=150)
    plt.close()

    # Save metadata
    meta = carrier.meta.copy()
    meta['max_envelope'] = float(envelope.max())
    meta['energy_in_plane_fraction'] = float(E_plane.sum() / (E_total.sum() + 1e-20))

    with open(out_path / f'{tag}_layer2_metadata.json', 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"\nVisualization saved to {out_path}")
    print(f"  - Rest geometry: {tag}_layer0_rest_geometry.png")
    print(f"  - Spec: {tag}_layer1_spec.png")
    print(f"  - Carrier: {tag}_layer2_carrier.png")