"""
Electron Stability Measurement and Loss Computation

This module provides metrics for assessing the stability of the toroidal
electron soliton and computing loss functions for optimization.

The stability metrics include:
1. Energy leakage: Fraction of energy radiating away from tube region
2. Shape drift: Change in envelope structure over time
3. Mode purity: Dominance of Compton frequency in internal oscillation
4. Constraint penalties: Deviations from physical constraints (E, P, S, Q)

These are combined into a scalar loss function that can be minimized to
find stable electron configurations.
"""

import torch
import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

from branesim.core.state import BraneState
from branesim.config.simulation_config import PhysicalConstants


@dataclass
class StabilityMetrics:
    """
    Container for electron stability measurements.

    Attributes:
        energy_leakage: Fraction of energy outside tube [0, 1+]
        shape_drift: Mean squared change in envelope [0, ∞)
        mode_purity_loss: 1 - (power at ω_C / total power) [0, 1]
        energy_error: Relative error in total energy vs m_e c²
        momentum_error: Relative momentum magnitude
        spin_error: Relative error in spin magnitude vs ℏ/2
        charge_error: Relative error in charge vs -e (if available)
        E_initial_ratio: E_initial / (m_e c²) - Initial energy ratio
        S_initial_ratio: |S_initial| / (ℏ/2) - Initial spin ratio
        P_initial_abs: |P_initial| - Initial momentum magnitude [kg·m/s]
        spin_alignment: cos(θ) between S and spin_axis - Spin alignment [-1, 1]
    """
    energy_leakage: float
    shape_drift: float
    mode_purity_loss: float
    energy_error: float
    momentum_error: float
    spin_error: float
    charge_error: Optional[float] = None
    E_initial_ratio: Optional[float] = None
    S_initial_ratio: Optional[float] = None
    P_initial_abs: Optional[float] = None
    spin_alignment: Optional[float] = None

    def __repr__(self) -> str:
        lines = [
            "StabilityMetrics:",
            f"  Energy leakage: {self.energy_leakage:.4f}",
            f"  Shape drift: {self.shape_drift:.6e}",
            f"  Mode purity loss: {self.mode_purity_loss:.4f}",
            f"  Energy error: {self.energy_error:.4f}",
            f"  Momentum error: {self.momentum_error:.4f}",
            f"  Spin error: {self.spin_error:.4f}",
        ]
        if self.charge_error is not None:
            lines.append(f"  Charge error: {self.charge_error:.4f}")
        if self.E_initial_ratio is not None:
            lines.append(f"  E_initial / (m_e c²): {self.E_initial_ratio:.6f}")
        if self.S_initial_ratio is not None:
            lines.append(f"  |S_initial| / (ℏ/2): {self.S_initial_ratio:.6f}")
        if self.P_initial_abs is not None:
            lines.append(f"  |P_initial|: {self.P_initial_abs:.6e} kg·m/s")
        if self.spin_alignment is not None:
            lines.append(f"  Spin alignment cos(θ): {self.spin_alignment:.6f}")
        return "\n".join(lines)


def compute_energy_density(state: BraneState, m_point: float) -> torch.Tensor:
    """
    Compute energy density at each brane point (kinetic only for simplicity).

    For full implementation, should include elastic potential energy contribution.

    Args:
        state: BraneState
        m_point: Point mass [kg]

    Returns:
        Energy density [N] in Joules per point
    """
    # Kinetic energy: (1/2) m v²
    v_mag_sq = (state.velocities ** 2).sum(dim=1)  # [N]
    E_kin = 0.5 * m_point * v_mag_sq  # [N]

    # TODO: Add potential energy contribution from springs
    # For now, return kinetic only
    return E_kin


def compute_energy_leakage(
    states: List[BraneState],
    tube_mask: torch.Tensor,
    m_point: float
) -> float:
    """
    Compute fraction of energy that has leaked out of the electron tube region.

    Args:
        states: List of BraneState snapshots over time
        tube_mask: Boolean mask [N] for tube region
        m_point: Point mass [kg]

    Returns:
        Average leakage fraction over time [0, 1+]
    """
    E0_density = compute_energy_density(states[0], m_point)
    E0_total = E0_density.sum()

    leakages = []
    for state in states[1:]:
        E_density = compute_energy_density(state, m_point)
        E_tube = (E_density * tube_mask).sum()

        # Fraction leaked = (E0 - E_tube) / E0
        frac_leaked = (E0_total - E_tube) / (E0_total + 1e-30)
        leakages.append(frac_leaked)

    # Average over time
    L_E = torch.mean(torch.stack(leakages))
    return torch.clamp(L_E, 0.0, 10.0).item()  # Clamp to reasonable range


def compute_shape_drift(
    states: List[BraneState],
    tube_mask: torch.Tensor
) -> float:
    """
    Compute mean squared change in amplitude envelope over time.

    Measures how much the envelope |X⁴| changes from initial configuration.

    Args:
        states: List of BraneState snapshots
        tube_mask: Boolean mask [N] for tube region

    Returns:
        Mean normalized squared difference
    """
    # Extract amplitude field envelope (|X⁴|) at t=0
    X4_0 = states[0].positions[:, 3]
    F0 = torch.abs(X4_0) * tube_mask
    norm0 = torch.sqrt((F0 ** 2).sum() + 1e-30)

    drifts = []
    for state in states[1:]:
        X4_t = state.positions[:, 3]
        Ft = torch.abs(X4_t) * tube_mask

        # Normalized squared difference
        diff = (Ft - F0)
        D_t = (diff ** 2).sum() / (norm0 ** 2 + 1e-30)
        drifts.append(D_t)

    D_shape = torch.mean(torch.stack(drifts))
    return D_shape.item()


def compute_core_time_series(
    states: List[BraneState],
    core_mask: torch.Tensor
) -> torch.Tensor:
    """
    Compute time series of average X⁴ over core region.

    Args:
        states: List of BraneState snapshots
        core_mask: Boolean mask [N] for core region

    Returns:
        Time series [T] of average amplitude
    """
    series = []
    for state in states:
        X4 = state.positions[:, 3]
        # Average over core region
        val = (X4 * core_mask).sum() / (core_mask.sum() + 1e-30)
        series.append(val)
    return torch.stack(series)


def compute_mode_purity(
    states: List[BraneState],
    core_mask: torch.Tensor,
    dt: float,
    target_omega: float,
    omega_window_frac: float = 0.1
) -> float:
    """
    Compute mode purity: fraction of power spectrum near target_omega.

    Uses FFT of core time series to measure how much power is concentrated
    at the Compton frequency vs. spread across other frequencies.

    Args:
        states: List of BraneState snapshots
        core_mask: Boolean mask [N] for core region
        dt: Time step between states [s]
        target_omega: Target frequency (ω_C) [rad/s]
        omega_window_frac: Fractional window around target (default 0.1 = ±10%)

    Returns:
        Mode purity loss = 1 - (power near ω_C / total power) in [0, 1]
    """
    signal = compute_core_time_series(states, core_mask)  # [T]
    T = signal.shape[0]

    if T < 4:
        return 0.5  # Not enough data, return neutral value

    # Remove mean to suppress DC
    signal = signal - signal.mean()

    # Real FFT
    fft_vals = torch.fft.rfft(signal)
    freqs = torch.fft.rfftfreq(T, d=dt)  # Hz
    omegas = 2.0 * torch.pi * freqs  # rad/s

    # Power spectrum
    power = fft_vals.real ** 2 + fft_vals.imag ** 2

    # Window around target frequency
    window = omega_window_frac * target_omega
    mask = (omegas >= target_omega - window) & (omegas <= target_omega + window)

    P_main = power[mask].sum()
    P_total = power.sum() + 1e-30

    purity = P_main / P_total  # Fraction of power near ω_C
    M = 1.0 - purity  # Loss: 0 = pure, 1 = totally spread

    return M.item()


def compute_total_energy(state: BraneState, m_point: float) -> float:
    """
    Compute total energy (kinetic only for now).

    Args:
        state: BraneState
        m_point: Point mass [kg]

    Returns:
        Total energy [J]
    """
    E_density = compute_energy_density(state, m_point)
    return E_density.sum().item()


def compute_total_momentum(state: BraneState, m_point: float) -> torch.Tensor:
    """
    Compute total linear momentum.

    Args:
        state: BraneState
        m_point: Point mass [kg]

    Returns:
        Momentum vector [4] in [kg·m/s]
    """
    # P = Σ m_p v_p
    P = m_point * state.velocities.sum(dim=0)  # [4]
    return P


def compute_total_spin(
    state: BraneState,
    m_point: float,
    center: Tuple[float, float, float]
) -> torch.Tensor:
    """
    Compute total angular momentum about center.

    Only uses spatial components (X^0, X^1, X^2) and corresponding velocities.

    Args:
        state: BraneState
        m_point: Point mass [kg]
        center: Center point (x, y, z) [m]

    Returns:
        Spin vector [3] in [J·s]
    """
    center_t = torch.tensor(center, dtype=state.dtype, device=state.device)

    # Positions and velocities (spatial only)
    r = state.positions[:, :3] - center_t.unsqueeze(0)  # [N, 3]
    v = state.velocities[:, :3]  # [N, 3]

    # Angular momentum: L = Σ r × (m v)
    # Cross product for each point
    mv = m_point * v  # [N, 3]
    L = torch.cross(r, mv, dim=1)  # [N, 3]

    # Sum over all points
    L_total = L.sum(dim=0)  # [3]

    return L_total


def compute_effective_charge(state: BraneState) -> float:
    """
    Compute effective charge from time-averaged amplitude.

    This is a placeholder - actual implementation requires the charge-from-energy
    relation from the paper.

    Args:
        state: BraneState

    Returns:
        Effective charge [C] (currently returns 0.0 as placeholder)
    """
    # TODO: Implement charge measurement from X̄⁴
    # For now, return 0.0 as placeholder
    return 0.0


def compute_constraint_penalties(
    state: BraneState,
    constants: PhysicalConstants,
    m_point: float,
    center: Tuple[float, float, float]
) -> Tuple[float, float, float, float]:
    """
    Compute penalties for constraint violations.

    Constraints:
        - Total energy ≈ m_e c²
        - Total momentum ≈ 0 (rest frame)
        - Spin magnitude ≈ ℏ/2
        - Charge ≈ -e (not yet implemented)

    Args:
        state: BraneState (initial state)
        constants: PhysicalConstants
        m_point: Point mass [kg]
        center: Center for spin calculation

    Returns:
        (L_E, L_P, L_S, L_Q): Constraint penalties as floats
    """
    # Energy constraint
    E0 = compute_total_energy(state, m_point)
    E_target = constants.m_e * constants.c ** 2
    L_E = ((E0 - E_target) / (E_target + 1e-30)) ** 2

    # Momentum constraint
    P0 = compute_total_momentum(state, m_point)  # [4]
    P0_spatial = P0[:3]  # Only spatial components
    L_P = (P0_spatial ** 2).sum() / ((constants.m_e * constants.c) ** 2 + 1e-30)

    # Spin constraint
    S0 = compute_total_spin(state, m_point, center)  # [3]
    S_mag = S0.norm()
    S_target = 0.5 * constants.hbar
    L_S = ((S_mag - S_target) / (S_target + 1e-30)) ** 2

    # Charge constraint (placeholder)
    Q0 = compute_effective_charge(state)
    e = 1.602176634e-19  # Elementary charge [C]
    L_Q = ((Q0 + e) / (e + 1e-30)) ** 2

    # Convert to Python floats, handling both tensors and scalars
    def to_float(x):
        return x.item() if hasattr(x, 'item') else float(x)

    return to_float(L_E), to_float(L_P), to_float(L_S), to_float(L_Q)


def compute_electron_stability_loss(
    states: List[BraneState],
    tube_mask: torch.Tensor,
    core_mask: torch.Tensor,
    constants: PhysicalConstants,
    m_point: float,
    dt: float,
    target_omega: float,
    center: Tuple[float, float, float],
    spin_axis: Optional[Tuple[float, float, float]] = None,
    weights: Optional[Dict[str, float]] = None
) -> Tuple[float, StabilityMetrics]:
    """
    Compute combined stability loss for electron optimization.

    This is the main function for assessing electron stability. It combines
    multiple metrics into a single scalar loss that can be minimized.

    Args:
        states: List of BraneState snapshots from simulation
        tube_mask: Boolean mask [N] for electron tube region
        core_mask: Boolean mask [N] for core region (where double loop is)
        constants: PhysicalConstants
        m_point: Point mass [kg]
        dt: Time step between states [s]
        target_omega: Target Compton frequency [rad/s]
        center: Center of torus for spin calculation
        spin_axis: Expected spin axis direction (x, y, z) [unitless], defaults to (0, 0, 1)
        weights: Optional dict of weighting factors for loss components

    Returns:
        (total_loss, metrics): Total loss value and detailed metrics
    """
    if spin_axis is None:
        spin_axis = (0.0, 0.0, 1.0)  # Default to z-axis
    if weights is None:
        weights = {
            'w_leak': 1.0,
            'w_shape': 1.0,
            'w_mode': 0.5,
            'w_Ec': 1.0,
            'w_Pc': 0.5,
            'w_Sc': 0.5,
            'w_Qc': 0.0,  # Disabled until charge measurement implemented
        }

    # Stability metrics
    L_leak = compute_energy_leakage(states, tube_mask, m_point)
    L_shape = compute_shape_drift(states, tube_mask)
    L_mode = compute_mode_purity(states, core_mask, dt, target_omega)

    # Constraint penalties (use initial state)
    L_Ec, L_Pc, L_Sc, L_Qc = compute_constraint_penalties(
        states[0], constants, m_point, center
    )

    # Compute initial ratios for diagnostics
    E_initial = compute_total_energy(states[0], m_point)
    E_target = constants.m_e * constants.c ** 2
    E_initial_ratio = E_initial / E_target

    P_initial = compute_total_momentum(states[0], m_point)
    P_initial_abs = torch.norm(P_initial[:3]).item()

    S_initial = compute_total_spin(states[0], m_point, center)
    S_initial_mag = torch.norm(S_initial).item()
    S_target = 0.5 * constants.hbar
    S_initial_ratio = S_initial_mag / S_target

    # Compute spin alignment with expected axis
    spin_axis_t = torch.tensor(spin_axis, dtype=states[0].dtype, device=states[0].device)
    spin_axis_norm = torch.norm(spin_axis_t) + 1e-30
    spin_axis_unit = spin_axis_t / spin_axis_norm

    # cos(θ) = (S · axis) / (|S| |axis|)
    if S_initial_mag > 1e-30:
        cos_theta = torch.dot(S_initial, spin_axis_unit) / S_initial_mag
        spin_alignment = cos_theta.item()
    else:
        spin_alignment = 0.0  # No spin, undefined alignment

    # Combine into total loss
    loss = (
        weights['w_leak'] * L_leak
        + weights['w_shape'] * L_shape
        + weights['w_mode'] * L_mode
        + weights['w_Ec'] * L_Ec
        + weights['w_Pc'] * L_Pc
        + weights['w_Sc'] * L_Sc
        + weights['w_Qc'] * L_Qc
    )

    # Package metrics
    metrics = StabilityMetrics(
        energy_leakage=L_leak,
        shape_drift=L_shape,
        mode_purity_loss=L_mode,
        energy_error=np.sqrt(L_Ec),
        momentum_error=np.sqrt(L_Pc),
        spin_error=np.sqrt(L_Sc),
        charge_error=np.sqrt(L_Qc) if weights['w_Qc'] > 0 else None,
        E_initial_ratio=E_initial_ratio,
        S_initial_ratio=S_initial_ratio,
        P_initial_abs=P_initial_abs,
        spin_alignment=spin_alignment,
    )

    return loss, metrics


def compute_lateralization_ratio(
    state: BraneState,
    tube_mask: torch.Tensor,
    m_point: float,
) -> float:
    """
    Compute ratio of lateral kinetic energy to amplitude kinetic energy.

    This metric indicates how much of the electron's internal motion is
    going into lateral (X^0, X^1, X^2) vs amplitude (X^4) oscillations.

    Near the geometric nonlinearity threshold, we expect significant
    lateral motion (ratio ~ O(1)) but not dominant.

    Args:
        state: BraneState snapshot
        tube_mask: Boolean mask [N] for electron tube region
        m_point: Point mass [kg]

    Returns:
        lateralization_ratio: E_kin_lateral / E_kin_amplitude
    """
    v_lat = state.velocities[:, :3]  # [N, 3]
    v_amp = state.velocities[:, 3]   # [N]

    # Kinetic energies in tube region only
    E_lat = 0.5 * m_point * (v_lat[tube_mask]**2).sum()
    E_amp = 0.5 * m_point * (v_amp[tube_mask]**2).sum()

    ratio = E_lat / (E_amp + 1e-30)
    return ratio.item()


def build_electron_masks(
    positions: torch.Tensor,
    center: Tuple[float, float, float],
    R: float,
    tube_max_radius: float,
    core_radius_frac: float = 0.5
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Build tube and core masks for electron region.

    Args:
        positions: [N, 4] tensor of positions
        center: Center of torus (x, y, z)
        R: Major radius of torus
        tube_max_radius: Maximum transverse distance for tube
        core_radius_frac: Fraction of tube_max_radius for core (default 0.5)

    Returns:
        (tube_mask, core_mask): Boolean masks [N] for tube and core regions
    """
    center_t = torch.tensor(center, dtype=positions.dtype, device=positions.device)
    X_lat = positions[:, :3]  # [N, 3]

    v = X_lat - center_t
    vx, vy, vz = v[:, 0], v[:, 1], v[:, 2]
    r_xy = torch.sqrt(vx * vx + vy * vy + 1e-30)

    radial_offset = r_xy - R
    transverse = torch.sqrt(radial_offset * radial_offset + vz * vz)

    tube_mask = transverse <= tube_max_radius
    core_mask = transverse <= (core_radius_frac * tube_max_radius)

    return tube_mask, core_mask