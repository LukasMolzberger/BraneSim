"""
Energy diagnostics for substrate evolution.

Provides tools for computing kinetic and potential energy from snapshots
or state objects. These diagnostics verify energy conservation and
detect numerical instabilities.

Note: Potential energy computation requires access to the force model,
which may not always be available from a snapshot alone.
"""

from __future__ import annotations
import torch
from .types import Snapshot, DiagnosticResult


def kinetic_energy_from_velocities(
    velocities: torch.Tensor,
    masses: torch.Tensor | float = 1.0
) -> torch.Tensor:
    """
    Compute kinetic energy from velocities.

    K = (1/2) Σ_p m_p |v_p|^2

    Parameters
    ----------
    velocities : torch.Tensor
        Velocity vectors, shape [N, D] where D is embedding dimension
    masses : torch.Tensor | float
        Particle masses, shape [N] or scalar (default: 1.0)

    Returns
    -------
    torch.Tensor
        Total kinetic energy (scalar)

    Examples
    --------
    >>> velocities = torch.randn(1000, 4)  # N=1000, 4D embedding
    >>> K = kinetic_energy_from_velocities(velocities)
    >>> K.shape
    torch.Size([])
    """
    # Compute squared velocity magnitude per particle
    v_squared = torch.sum(velocities**2, dim=-1)  # [N]

    # Kinetic energy per particle
    if isinstance(masses, float):
        K_per_particle = 0.5 * masses * v_squared
    else:
        K_per_particle = 0.5 * masses * v_squared

    # Total kinetic energy
    return K_per_particle.sum()


def kinetic_energy_from_snapshot(
    snapshot: Snapshot,
    masses: torch.Tensor | float = 1.0
) -> float:
    """
    Compute kinetic energy from a snapshot.

    Parameters
    ----------
    snapshot : Snapshot
        Snapshot containing velocity field
        Must have "vel" key in fields dict
    masses : torch.Tensor | float
        Particle masses

    Returns
    -------
    float
        Total kinetic energy

    Examples
    --------
    >>> snapshot = Snapshot(
    ...     t_sim=10.0,
    ...     fields={"vel": torch.randn(1000, 4)}
    ... )
    >>> K = kinetic_energy_from_snapshot(snapshot)
    """
    if "vel" not in snapshot.fields:
        raise ValueError("Snapshot must contain 'vel' field")

    velocities = snapshot.fields["vel"]
    K = kinetic_energy_from_velocities(velocities, masses)
    return K.item()


def potential_energy_placeholder(
    positions: torch.Tensor,
    force_model: object | None = None
) -> torch.Tensor:
    """
    Placeholder for potential energy computation.

    Actual implementation depends on the specific force model and
    energy functional used by the substrate solver.

    Parameters
    ----------
    positions : torch.Tensor
        Particle positions, shape [N, D]
    force_model : object | None
        Force model object that can compute potential energy

    Returns
    -------
    torch.Tensor
        Total potential energy (scalar)

    Notes
    -----
    This is a placeholder. Real implementation should:
    1. Access the solver's energy computation
    2. Or recompute energy from positions using the force model
    3. Or read cached energy values from the snapshot metadata
    """
    raise NotImplementedError(
        "Potential energy computation requires access to the force model. "
        "This should be implemented by accessing solver internals or "
        "reading pre-computed energy values from snapshot metadata."
    )


def total_energy_timeseries_from_snapshots(
    snapshots: list[Snapshot],
    masses: torch.Tensor | float = 1.0,
    force_model: object | None = None
) -> DiagnosticResult:
    """
    Compute energy timeseries from a list of snapshots.

    Parameters
    ----------
    snapshots : list[Snapshot]
        List of snapshots
    masses : torch.Tensor | float
        Particle masses
    force_model : object | None
        Force model for potential energy (optional)

    Returns
    -------
    DiagnosticResult
        Result containing:
        - data["kinetic"]: kinetic energy array
        - data["total"]: total energy array (if potential available)
        - data["t_sim"]: simulation time array
        - data["t_phys_s"]: physical time array (if available)
        - quality["energy_drift"]: relative energy drift from initial value

    Examples
    --------
    >>> snapshots = [...]  # List of snapshots
    >>> result = total_energy_timeseries_from_snapshots(snapshots)
    >>> K = result.data["kinetic"]
    >>> drift = result.quality["energy_drift"]
    """
    n_snapshots = len(snapshots)

    # Allocate arrays
    K_array = torch.zeros(n_snapshots)
    t_sim_array = torch.zeros(n_snapshots)
    t_phys_array = torch.zeros(n_snapshots)

    # Compute kinetic energy for each snapshot
    for i, snap in enumerate(snapshots):
        K_array[i] = kinetic_energy_from_snapshot(snap, masses)
        t_sim_array[i] = snap.t_sim if snap.t_sim is not None else float('nan')
        t_phys_array[i] = snap.t_phys_s if snap.t_phys_s is not None else float('nan')

    # Build result
    data = {
        "kinetic": K_array,
        "t_sim": t_sim_array,
    }

    # Add physical time if available
    if not torch.isnan(t_phys_array).all():
        data["t_phys_s"] = t_phys_array

    # Compute energy drift (relative to initial)
    if K_array[0] > 0:
        energy_drift = (K_array - K_array[0]) / K_array[0]
    else:
        energy_drift = K_array - K_array[0]

    quality = {
        "energy_drift": energy_drift,
    }

    meta = {
        "n_snapshots": n_snapshots,
        "masses": "uniform" if isinstance(masses, float) else "per-particle",
    }

    return DiagnosticResult(
        name="energy_timeseries",
        data=data,
        quality=quality,
        meta=meta,
    )