"""
Carrier compilation layer (Layer 2): specs → substrate kinematics.

This module compiles high-level specs into concrete displacement and
velocity fields that can be applied to the substrate. It integrates:
- Envelope construction
- Polarization basis selection
- Phase patterns
- Velocity initialization (momentum is part of this layer, not separate)
"""

import torch
import numpy as np
from typing import Optional, Union
from pathlib import Path

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.init.artifacts import (
    RestGeometryArtifact,
    SpecArtifact,
    CarrierArtifact,
    InitPipelineArtifact,
)
from branesim.init.specs import PhotonSpec, ElectronSpec
from branesim.init.polarization import photon_polarization_basis
from branesim.init.envelopes import (
    gaussian_envelope,
    plane_wave_phase,
    tubular_envelope_electron,
)
from branesim.init import velocity
from branesim.initialization.geometry.tubular_electron_geometry import (
    sample_torus_knot_centerline,
    TorusKnotParameters,
)


def compile_photon(
    geom: RestGeometryArtifact,
    spec: PhotonSpec,
    physics=None,
    grid=None,
    m_point: Optional[float] = None,
    wave_speed: Optional[float] = None,
) -> CarrierArtifact:
    """
    Compile photon spec into displacement and velocity fields.

    This is Layer 2 compilation: takes the "what" (spec) and produces
    the "how" (u0, v0) using the substrate geometry.

    Momentum is integrated here via velocity initialization - there is
    no separate "Layer 3".

    Args:
        geom: Layer 0 rest geometry
        spec: PhotonSpec describing desired packet
        physics: Force computer (required for time_reversal_shift mode)
        grid: Grid topology (required for time_reversal_shift mode)
        m_point: Point mass (required for time_reversal_shift mode)
        wave_speed: Wave speed (required for velocity initialization)

    Returns:
        CarrierArtifact with psi, u0, v0, and metadata
    """
    device = geom.rest_positions.device
    dtype = geom.rest_positions.dtype
    coords = geom.coords  # [N, d]
    N = coords.shape[0]

    # Extract k vector parameters
    k_vec = spec.k_vector.to(device=device, dtype=dtype)
    k_mag = torch.linalg.norm(k_vec).item()
    k_hat = k_vec / k_mag

    # Build envelope
    center = spec.center.to(device=device, dtype=dtype)
    A = gaussian_envelope(coords, center, spec.sigma, spec.amplitude)  # [N]

    # Build phase
    phi = plane_wave_phase(coords, k_vec)  # [N]

    # Get polarization basis
    p1, p2 = photon_polarization_basis(
        spec.intrinsic_dim,
        k_hat,
        prefer_shear=spec.prefer_shear,
        device=device,
        dtype=dtype,
    )

    # Build complex carrier with helicity
    # Right/Left helicity: psi = A * exp(iφ) * (p1 + i*s*p2)
    # where s = +1 for R, s = -1 for L
    s = 1.0 if spec.helicity == "R" else -1.0

    # Expand A and phi for broadcasting
    A_exp = A.unsqueeze(1)  # [N, 1]
    phi_exp = phi.unsqueeze(1)  # [N, 1]

    # Build complex exponential
    exp_iphi = torch.complex(torch.cos(phi_exp), torch.sin(phi_exp))  # [N, 1]

    # Polarization: p1 + i*s*p2
    p_complex = torch.complex(p1, s * p2).unsqueeze(0)  # [1, 4]

    # Full carrier: A * exp(iφ) * p
    psi = A_exp * exp_iphi * p_complex  # [N, 4]

    # Displacement: u0 = Re(psi)
    u0 = psi.real

    # Velocity depends on method
    if spec.velocity_init == "time_reversal_shift":
        if physics is None or grid is None or m_point is None or wave_speed is None:
            raise ValueError(
                "time_reversal_shift requires physics, grid, m_point, and wave_speed"
            )
        v0 = velocity.velocities_time_reversal_shift(
            geom, u0, physics, grid, m_point, wave_speed,
            spec.shift_cells, k_hat, spec.periodic_shift
        )

    elif spec.velocity_init == "directional_derivative":
        if wave_speed is None:
            raise ValueError("directional_derivative requires wave_speed")
        v0 = velocity.velocities_directional_derivative(
            geom, u0, wave_speed, k_hat
        )

    elif spec.velocity_init == "complex_quadrature":
        # Estimate omega from k and wave_speed
        if wave_speed is None:
            raise ValueError("complex_quadrature requires wave_speed")
        omega = k_mag * wave_speed
        v0 = velocity.velocities_from_complex_quadrature(psi, omega)

    else:
        raise ValueError(f"Unknown velocity_init: {spec.velocity_init}")

    # Build metadata
    meta = {
        "k_mag": k_mag,
        "k_hat": k_hat.cpu().numpy(),
        "wavelength": 2 * np.pi / k_mag if k_mag > 0 else float('inf'),
        "helicity": spec.helicity,
        "velocity_method": spec.velocity_init,
        "max_displacement": torch.abs(u0).max().item(),
        "max_velocity": torch.abs(v0).max().item(),
    }

    return CarrierArtifact(
        envelope=A,
        phase=phi,
        p1=p1,
        p2=p2,
        psi=psi,
        u0=u0,
        v0=v0,
        meta=meta,
    )


def compile_electron(
    geom: RestGeometryArtifact,
    spec: ElectronSpec,
    wave_speed: Optional[float] = None,
) -> CarrierArtifact:
    """
    Compile electron spec into tubular displacement and velocity fields.

    This implements the double-loop tube with spinorial transport.

    Args:
        geom: Layer 0 rest geometry (must be 3D)
        spec: ElectronSpec describing tubular geometry
        wave_speed: Wave speed for velocity (optional)

    Returns:
        CarrierArtifact with tubular carrier
    """
    if geom.intrinsic_dim != 3:
        raise ValueError(f"Electron requires 3D geometry, got {geom.intrinsic_dim}D")

    device = geom.rest_positions.device
    dtype = geom.rest_positions.dtype
    coords = geom.coords  # [N, 3]
    N = coords.shape[0]

    # Sample centerline
    torus_params = TorusKnotParameters(
        major_radius=spec.torus_major_radius,
        minor_radius=spec.torus_minor_radius,
        core_windings=spec.p,
        tube_windings=spec.q,
    )
    centerline_np = sample_torus_knot_centerline(torus_params, spec.num_samples)
    centerline = torch.from_numpy(centerline_np).to(device=device, dtype=dtype)

    # Center the tubular structure
    if spec.center is not None:
        center = spec.center.to(device=device, dtype=dtype)
        centerline = centerline + center

    # Compute tubular envelope
    envelope, s_norm, alpha, alpha_half = tubular_envelope_electron(
        coords, centerline, spec.tube_sigma, device=device, dtype=dtype
    )

    # Scale envelope
    envelope = spec.amplitude * envelope  # [N]

    # Longitudinal phase
    k_long = spec.longitudinal_k
    phi = k_long * s_norm  # [N]

    # Polarization basis with half-angle rotation
    # p1(x) = (cos(α/2), sin(α/2), 0, 0) in local frame
    # p2(x) = (-sin(α/2), cos(α/2), 0, 0) in local frame
    # For simplicity, embed these directly
    cos_half = torch.cos(alpha_half)  # [N]
    sin_half = torch.sin(alpha_half)  # [N]

    # Build position-dependent polarization (not constant!)
    # p1[n] has shape [4], but here we need [N, 4]
    p1 = torch.zeros(N, 4, device=device, dtype=dtype)
    p1[:, 0] = cos_half
    p1[:, 1] = sin_half

    p2 = torch.zeros(N, 4, device=device, dtype=dtype)
    p2[:, 0] = -sin_half
    p2[:, 1] = cos_half

    # Helicity
    s = 1.0 if spec.helicity == "R" else -1.0

    # Complex carrier: A(x) * exp(iφ(x)) * (p1(x) + i*s*p2(x))
    A_exp = envelope.unsqueeze(1)  # [N, 1]
    phi_exp = phi.unsqueeze(1)  # [N, 1]
    exp_iphi = torch.complex(torch.cos(phi_exp), torch.sin(phi_exp))  # [N, 1]

    p_complex = torch.complex(p1, s * p2)  # [N, 4]
    psi = A_exp * exp_iphi * p_complex  # [N, 4]

    # Displacement
    u0 = psi.real

    # Velocity
    if spec.velocity_init == "complex_quadrature":
        # Estimate omega from longitudinal k and loop circumference
        loop_circ = 2 * np.pi * spec.torus_major_radius
        # Assume phase speed ~ wave_speed if provided
        if wave_speed is not None:
            phase_speed = wave_speed
        else:
            # Fallback: assume some reasonable value
            phase_speed = 1.0  # Document this assumption
            print(f"Warning: Using default phase_speed={phase_speed} for electron")

        # Frequency: ω = k * v_phase
        omega = k_long * phase_speed
        v0 = velocity.velocities_from_complex_quadrature(psi, omega)

    else:
        # For other methods, would need to implement loop-tangent derivative
        raise NotImplementedError(
            f"Electron velocity_init {spec.velocity_init} not yet implemented. "
            "Use 'complex_quadrature'."
        )

    # Metadata
    meta = {
        "torus_major_radius": spec.torus_major_radius,
        "torus_minor_radius": spec.torus_minor_radius,
        "p": spec.p,
        "q": spec.q,
        "longitudinal_k": k_long,
        "helicity": spec.helicity,
        "velocity_method": spec.velocity_init,
        "max_displacement": torch.abs(u0).max().item(),
        "max_velocity": torch.abs(v0).max().item(),
    }

    # For electron, polarization is position-dependent, so store first sample
    # (visualization will need to handle this differently)
    return CarrierArtifact(
        envelope=envelope,
        phase=phi,
        p1=p1[0],  # Store first sample as representative
        p2=p2[0],
        psi=psi,
        u0=u0,
        v0=v0,
        meta=meta,
    )


def build_rest_geometry_artifact(
    state: BraneState,
    grid: BraneGrid,
) -> RestGeometryArtifact:
    """
    Build Layer 0 artifact from existing state and grid.

    Args:
        state: BraneState with rest_positions set
        grid: Grid topology

    Returns:
        RestGeometryArtifact
    """
    if not hasattr(state, 'rest_positions') or state.rest_positions is None:
        raise ValueError(
            "BraneState must have rest_positions set. "
            "Call initialize_flat_configuration first."
        )

    # Extract intrinsic coordinates from rest positions
    d = len(grid.grid_shape)
    coords = state.rest_positions[:, :d]

    # Get dimension enum
    if d == 1:
        dim_val = 1
    elif d == 2:
        dim_val = 2
    else:
        dim_val = 3

    return RestGeometryArtifact(
        intrinsic_dim=dim_val,
        embedding_dim=4,
        grid_shape=grid.grid_shape,
        spacing=grid.spacing,
        rest_positions=state.rest_positions.clone(),
        coords=coords,
        fixed_mask=state.fixed_mask.clone() if state.fixed_mask is not None else None,
    )


def build_spec_artifact(
    spec: Union[PhotonSpec, ElectronSpec],
) -> SpecArtifact:
    """
    Build Layer 1 artifact from spec.

    Args:
        spec: PhotonSpec or ElectronSpec

    Returns:
        SpecArtifact
    """
    if isinstance(spec, PhotonSpec):
        kind = "photon"
        k_vec = spec.k_vector
        k_mag = torch.linalg.norm(k_vec).item()
        k_hat = k_vec / k_mag
        notes = {
            "helicity": spec.helicity,
            "sigma": spec.sigma,
            "amplitude": spec.amplitude,
        }
    elif isinstance(spec, ElectronSpec):
        kind = "electron"
        # For electron, there's no single k_hat (it's a loop)
        # Use a placeholder
        k_hat = torch.tensor([1.0, 0.0, 0.0])
        k_mag = spec.longitudinal_k
        notes = {
            "torus_major_radius": spec.torus_major_radius,
            "torus_minor_radius": spec.torus_minor_radius,
            "p": spec.p,
            "q": spec.q,
        }
    else:
        raise ValueError(f"Unknown spec type: {type(spec)}")

    return SpecArtifact(
        kind=kind,
        spec=spec,
        k_hat=k_hat,
        k_mag=k_mag,
        notes=notes,
    )


def initialize_state_from_spec(
    state: BraneState,
    grid: BraneGrid,
    spec: Union[PhotonSpec, ElectronSpec],
    physics=None,
    m_point: Optional[float] = None,
    wave_speed: Optional[float] = None,
    debug_out_dir: Optional[str] = None,
    debug_tag: str = "",
) -> InitPipelineArtifact:
    """
    Main entrypoint: compile spec and apply to state.

    This runs the full pipeline:
    1. Layer 0: Extract rest geometry from state
    2. Layer 1: Build spec artifact
    3. Layer 2: Compile carrier (displacement + velocity)
    4. Apply u0, v0 to state using set_kinematics
    5. Optionally generate debug visualizations

    Args:
        state: BraneState to initialize (must have rest_positions)
        grid: Grid topology
        spec: PhotonSpec or ElectronSpec
        physics: Force computer (for time_reversal_shift)
        m_point: Point mass (for time_reversal_shift)
        wave_speed: Wave speed (for velocity initialization)
        debug_out_dir: Optional directory for debug plots
        debug_tag: Tag for debug filenames

    Returns:
        InitPipelineArtifact with all layer outputs
    """
    # Layer 0: Rest geometry
    layer0 = build_rest_geometry_artifact(state, grid)

    # Layer 1: Spec
    layer1 = build_spec_artifact(spec)

    # Layer 2: Compile carrier
    if isinstance(spec, PhotonSpec):
        layer2 = compile_photon(
            layer0, spec, physics, grid, m_point, wave_speed
        )
    elif isinstance(spec, ElectronSpec):
        layer2 = compile_electron(layer0, spec, wave_speed)
    else:
        raise ValueError(f"Unknown spec type: {type(spec)}")

    # Apply to state
    state.set_kinematics(layer2.u0, layer2.v0)

    # Generate visualizations if requested
    if debug_out_dir is not None:
        from branesim.init import visualize
        out_path = Path(debug_out_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        visualize.plot_rest_geometry(layer0, str(out_path), debug_tag)
        visualize.plot_spec(layer1, str(out_path), debug_tag)
        visualize.plot_carrier(layer0, layer1, layer2, str(out_path), debug_tag)

    return InitPipelineArtifact(
        layer0=layer0,
        layer1=layer1,
        layer2=layer2,
    )