"""
Lateralization Ratio Measurement (Amplitude vs. Lateral dimensions)

For each brane point i we compute:

    E_amp(i)  = kinetic + potential energy in the amplitude dimension ξ
    E_lat(i)  = kinetic + potential energy in all lateral dimensions

and the local lateralization ratio

    R_lat(i) = E_lat(i) / (E_amp(i) + E_lat(i) + eps)

This is PURELY DIAGNOSTIC and never feeds back into the dynamics.
"""

from dataclasses import dataclass
from typing import Tuple, Dict, Any

import torch
import numpy as np

from branesim.core.state import BraneState
from branesim.core.grid import BraneGrid


@dataclass
class LateralizationConfig:
    """
    Split energy into amplitude vs lateral directions.

    amplitude_dim:
        Index of the amplitude / ξ coordinate in positions/velocities.
        In all photon examples this is 3.

    lateral_dims:
        Indices of lateral coordinates in positions/velocities:
          - 1D photon: (0,)       -> x
          - 2D photon: (0, 1)     -> x, y
          - 3D photon: (0, 1, 2)  -> x, y, z
    """
    amplitude_dim: int
    lateral_dims: Tuple[int, ...]

    include_kinetic: bool = True
    include_potential: bool = True
    eps: float = 1e-16


class LateralizationMeasurement:
    """
    Purely diagnostic: does NOT affect dynamics.

    Usage:
        meas = LateralizationMeasurement(config, grid, m_point, reference_positions)
        R_local, R_global, diag = meas.measure(state, physics)
    """

    def __init__(
        self,
        config: LateralizationConfig,
        grid: BraneGrid,
        m_point: float,
        reference_positions: torch.Tensor = None,
    ):
        self.config = config
        self.grid = grid
        self.m_point = float(m_point)

        # Reference configuration for "zero dynamic energy"
        # (flat, pre-stretched brane with ξ = 0)
        if reference_positions is not None:
            self.reference_positions = reference_positions.detach().clone()
        else:
            self.reference_positions = None

    # ---------- kinetic part --------------------------------------------
    def _measure_kinetic(self, state: BraneState) -> Tuple[torch.Tensor, torch.Tensor]:
        N = state.positions.shape[0]
        device = state.device
        dtype = state.dtype

        E_amp = torch.zeros(N, dtype=dtype, device=device)
        E_lat = torch.zeros(N, dtype=dtype, device=device)

        if not self.config.include_kinetic:
            return E_amp, E_lat

        v = state.velocities  # (N, D_embed)

        amp_idx = self.config.amplitude_dim
        lat_idx = torch.tensor(self.config.lateral_dims, dtype=torch.long, device=device)

        v_amp = v[:, amp_idx]                   # (N,)
        v_lat = v.index_select(1, lat_idx)      # (N, n_lat)
        v_lat_sq = torch.sum(v_lat**2, dim=1)   # (N,)

        m = self.m_point
        E_amp = 0.5 * m * v_amp**2
        E_lat = 0.5 * m * v_lat_sq
        return E_amp, E_lat

    # ---------- potential part ------------------------------------------
    def _measure_potential(self, state: BraneState, physics) -> Tuple[torch.Tensor, torch.Tensor]:
        N = state.positions.shape[0]
        device = state.device
        dtype = state.dtype

        E_amp = torch.zeros(N, dtype=dtype, device=device)
        E_lat = torch.zeros(N, dtype=dtype, device=device)

        if not self.config.include_potential or self.grid.neighbors is None:
            return E_amp, E_lat

        if self.reference_positions is None:
            # Without reference positions, we can't compute incremental distortions
            return E_amp, E_lat

        positions = state.positions
        neighbors = self.grid.neighbors

        amp_idx = self.config.amplitude_dim
        lat_idx = torch.tensor(self.config.lateral_dims, dtype=torch.long, device=device)

        k = float(physics.spring_constant)

        for i in range(N):
            for j in neighbors[i]:
                j_val = j.item() if isinstance(j, torch.Tensor) else int(j)
                if j_val < 0 or j_val <= i:   # avoid invalid / double counts
                    continue
                j = j_val

                # Current neighbor vector
                dX = positions[j] - positions[i]           # (D_embed,)

                # Baseline neighbor vector
                dX0 = self.reference_positions[j] - self.reference_positions[i]  # (D_embed,)

                # Incremental distortion (this is SMALL, no big-minus-big)
                delta_X = dX - dX0  # (D_embed,)

                # Split incremental distortion into amplitude and lateral components
                delta_X_amp = delta_X[amp_idx]  # scalar
                delta_X_lat = delta_X.index_select(0, lat_idx)  # (n_lat,)

                # Energy from incremental distortions
                E_spring_amp = 0.5 * k * delta_X_amp**2
                E_spring_lat = 0.5 * k * torch.sum(delta_X_lat**2)

                # Distribute half to each endpoint
                E_amp[i] += 0.5 * E_spring_amp
                E_amp[j] += 0.5 * E_spring_amp
                E_lat[i] += 0.5 * E_spring_lat
                E_lat[j] += 0.5 * E_spring_lat

        return E_amp, E_lat

    # ---------- public API ----------------------------------------------
    def measure(self, state: BraneState, physics) -> Tuple[torch.Tensor, float, Dict[str, Any]]:
        """
        Returns:
            R_lat_local : (N,) tensor per point
            R_lat_global: float
            diagnostics : dict of energy splits
        """
        E_amp_kin, E_lat_kin = self._measure_kinetic(state)
        E_amp_pot, E_lat_pot = self._measure_potential(state, physics)

        E_amp = E_amp_kin + E_amp_pot
        E_lat = E_lat_kin + E_lat_pot

        eps = self.config.eps
        E_tot = E_amp + E_lat + eps

        R_lat_local = E_lat / E_tot

        E_amp_global = torch.sum(E_amp)
        E_lat_global = torch.sum(E_lat)
        R_lat_global = float(E_lat_global / (E_amp_global + E_lat_global + eps))

        diagnostics = {
            "E_amp_kin": E_amp_kin,
            "E_lat_kin": E_lat_kin,
            "E_amp_pot": E_amp_pot,
            "E_lat_pot": E_lat_pot,
            "E_amp_total": E_amp,
            "E_lat_total": E_lat,
        }

        return R_lat_local, R_lat_global, diagnostics