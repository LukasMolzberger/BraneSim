"""
Material Calibration: Physics-Anchored Spring Parameter Derivation

This module implements the three-layer calibration approach for deriving
physically realistic spring network parameters from continuum physics targets.

The calibration process:
1. Compute continuum targets (ρ_m, T_target, κ_target) from fundamental constants
2. Derive analytic first guess for rest_length_frac from tension requirement
3. Refine using electron stability experiment to match geometric nonlinearity threshold

See CALIBRATION_NOTES.md for detailed physics background.
"""

import torch
import numpy as np
import warnings
from typing import Tuple, Dict, Optional
from dataclasses import dataclass

from branesim.config.simulation_config import PhysicalConstants
from branesim.physics.dimensional_mapping import DimensionalMapper


@dataclass
class ContinuumTargets:
    """
    Continuum elastic parameters derived from fundamental physics.

    Attributes:
        rho_m: Mass density [kg/m²] - from electron calibration
        T_target: Target tension [N/m] - for wave speed c
        kappa_target: Target bending stiffness [J·m] - for gravitational coupling
        c: Speed of light [m/s]
        G: Newton's constant [m³/(kg·s²)]
    """
    rho_m: float
    T_target: float
    kappa_target: float
    c: float
    G: float

    def __repr__(self) -> str:
        return (
            f"ContinuumTargets(\n"
            f"  ρ_m = {self.rho_m:.6e} kg/m²,\n"
            f"  T_target = {self.T_target:.6e} N/m,\n"
            f"  κ_target = {self.kappa_target:.6e} J·m,\n"
            f"  c = {self.c:.6e} m/s,\n"
            f"  G = {self.G:.6e} m³/(kg·s²)\n"
            f")"
        )


@dataclass
class MaterialCalibrationResult:
    """
    Result of material calibration procedure.

    Attributes:
        rest_length_frac: Calibrated dimensionless pre-strain factor
        rest_length_frac_analytic: Initial analytic guess from tension
        rest_length_phys: Physical rest length [m]
        kappa_phys: Physical bending stiffness [J·m]
        continuum_targets: ContinuumTargets used for calibration
        refinement_loss: Final loss from electron refinement (if performed)
    """
    rest_length_frac: float
    rest_length_frac_analytic: float
    rest_length_phys: float
    kappa_phys: float
    continuum_targets: ContinuumTargets
    refinement_loss: Optional[float] = None

    def __repr__(self) -> str:
        lines = [
            "MaterialCalibrationResult:",
            f"  rest_length_frac (analytic) = {self.rest_length_frac_analytic:.6f}",
            f"  rest_length_frac (final)    = {self.rest_length_frac:.6f}",
            f"  rest_length_phys = {self.rest_length_phys:.6e} m",
            f"  κ_phys = {self.kappa_phys:.6e} J·m",
        ]
        if self.refinement_loss is not None:
            lines.append(f"  refinement_loss = {self.refinement_loss:.6e}")
        lines.append(f"\n{self.continuum_targets}")
        return "\n".join(lines)


def compute_continuum_targets(
    rho_m: float,
    constants: PhysicalConstants,
) -> ContinuumTargets:
    """
    Compute continuum elastic parameter targets from fundamental physics.

    Args:
        rho_m: Mass density [kg/m²] - typically from electron calibration
        constants: PhysicalConstants with c, hbar, m_e, G

    Returns:
        ContinuumTargets with T_target and κ_target
    """
    c = constants.c
    G = constants.G

    # Tension required for wave speed = c
    # c² = T / ρ_m  ⟹  T = ρ_m c²
    T_target = rho_m * c**2

    # Bending stiffness for correct gravitational coupling
    # κ/2 ≈ c³/(16π G)  ⟹  κ = c³/(8π G)
    kappa_target = c**3 / (8.0 * np.pi * G)

    return ContinuumTargets(
        rho_m=rho_m,
        T_target=T_target,
        kappa_target=kappa_target,
        c=c,
        G=G,
    )


def solve_rest_length_frac_from_tension(
    T_target: float,
    k_phys: float,
    h_phys: float,
) -> float:
    """
    Solve for rest_length_frac such that discrete spring pre-strain
    reproduces the target continuum tension T_target.

    Derivation:
        For a cubic spring network with lattice spacing h, spring constant k,
        and rest length L_rest = rest_length_frac * h:

        Pre-stretch: ΔL = h - L_rest
        Pre-stretch force: F_0 = k ΔL
        Tension per area: T ≈ F_0 / h² = k(h - L_rest) / h²

        Setting T = T_target:
            T_target = k(h - L_rest) / h²
            h - L_rest = T_target h² / k
            1 - (L_rest / h) = T_target h² / k
            1 - rest_length_frac = T_target h² / k

        Therefore:
            rest_length_frac = 1 - (T_target h²) / k

    Args:
        T_target: Target continuum tension [N/m]
        k_phys: Spring constant in SI [N/m]
        h_phys: Grid spacing in SI [m]

    Returns:
        rest_length_frac: Dimensionless pre-strain factor [0, 1]
    """
    rest_length_frac = 1.0 - (T_target * h_phys**2) / k_phys
    return rest_length_frac


def validate_rest_length_frac(
    rest_length_frac: float,
    min_allowed: float = 0.8,
    max_allowed: float = 0.95,
) -> Tuple[float, bool]:
    """
    Validate and optionally clamp rest_length_frac to acceptable range.

    The range [0.8, 0.95] is empirically known to work well for the
    electron bottleneck geometry.

    Args:
        rest_length_frac: Candidate value
        min_allowed: Minimum acceptable value (default 0.8)
        max_allowed: Maximum acceptable value (default 0.95)

    Returns:
        (clamped_value, is_valid): Clamped value and whether original was valid
    """
    is_valid = (min_allowed <= rest_length_frac <= max_allowed)

    if not is_valid:
        warnings.warn(
            f"Derived rest_length_frac={rest_length_frac:.4f} outside "
            f"acceptable range [{min_allowed}, {max_allowed}]. "
            f"Consider adjusting k_phys or h_phys for better numerical properties."
        )

    clamped = np.clip(rest_length_frac, min_allowed, max_allowed)
    return clamped, is_valid


def combined_rest_length_loss(metrics: Dict[str, float]) -> float:
    """
    Combine stability metrics into scalar loss for rest_length optimization.

    This loss function balances:
    - Energy leakage (want minimal radiation)
    - Shape drift (want stable envelope)
    - Mode purity (want pure Compton frequency)
    - Lateralization ratio (want near-threshold, ~1.0)

    Args:
        metrics: Dictionary containing:
            - energy_leakage: [0, 1+] - fraction of energy leaked
            - shape_drift: [0, ∞) - envelope change
            - mode_purity_loss: [0, 1] - 1 - (power at ω_C / total)
            - lateralization_ratio: (optional) [0, ∞) - lateral/amplitude energy

    Returns:
        Combined loss (lower is better)
    """
    # Weights (tunable based on priorities)
    w_leak = 1.0   # Energy conservation is critical
    w_shape = 1.0  # Shape stability is critical
    w_mode = 1.0   # Mode purity is important
    w_lat = 0.5    # Lateralization is a soft constraint

    L = 0.0
    L += w_leak * metrics.get("energy_leakage", 0.0)
    L += w_shape * metrics.get("shape_drift", 0.0)
    L += w_mode * metrics.get("mode_purity_loss", 0.0)

    # Lateralization: target ~ 1.0 (near-threshold nonlinearity)
    # Too low (<0.5): nonlinearity too weak, electron disperses
    # Too high (>2.0): nonlinearity too strong, electron collapses
    if "lateralization_ratio" in metrics:
        lat = metrics["lateralization_ratio"]
        # Quadratic penalty centered at 1.0
        L += w_lat * (lat - 1.0)**2

    return L


def compute_lateralization_ratio(
    state,
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


def run_short_electron_stability_test(
    state,
    params,
    integrator,
    n_periods: float = 3.0,
) -> Dict[str, float]:
    """
    Run a short electron stability test for rest_length calibration.

    This is a lightweight version of the full electron experiment,
    optimized for scanning multiple rest_length candidates.

    Args:
        state: Initial BraneState with electron initialized
        params: ElectronInitParams with geometry
        integrator: Time integrator instance
        n_periods: Number of Compton periods to simulate (default 3.0)

    Returns:
        metrics: Dictionary with stability metrics
    """
    # This is a placeholder - actual implementation will depend on
    # how the full simulation is structured
    #
    # The full implementation should:
    # 1. Run simulation for n_periods of Compton oscillation
    # 2. Collect state snapshots at regular intervals
    # 3. Compute stability metrics using electron_stability.py functions
    # 4. Return metrics dict

    raise NotImplementedError(
        "run_short_electron_stability_test needs integration with "
        "the simulation framework. See material_calibration.py for "
        "implementation notes."
    )

    # Pseudocode for reference:
    # T_compton = 2 * np.pi / params.compton_omega
    # t_total = n_periods * T_compton
    # dt = integrator.dt
    # n_steps = int(t_total / dt)
    # n_snapshots = min(100, n_steps // 10)
    # snapshot_interval = n_steps // n_snapshots
    #
    # states = [state.clone()]
    # for i in range(n_steps):
    #     state = integrator.step(state)
    #     if i % snapshot_interval == 0:
    #         states.append(state.clone())
    #
    # # Build masks
    # tube_mask, core_mask = build_electron_masks(...)
    #
    # # Compute stability metrics
    # from branesim.diagnostics.electron_stability import (
    #     compute_energy_leakage, compute_shape_drift, compute_mode_purity
    # )
    #
    # metrics = {
    #     "energy_leakage": compute_energy_leakage(states, tube_mask, m_point),
    #     "shape_drift": compute_shape_drift(states, tube_mask),
    #     "mode_purity_loss": compute_mode_purity(
    #         states, core_mask, dt, params.compton_omega
    #     ),
    #     "lateralization_ratio": compute_lateralization_ratio(
    #         states[-1], tube_mask, m_point
    #     ),
    # }
    #
    # return metrics


def calibrate_rest_length_frac_with_electron(
    mapper: DimensionalMapper,
    k_phys: float,
    h_phys: float,
    rest_length_frac_initial: float,
    scan_radius: float = 0.05,
    n_scan: int = 5,
    run_refinement: bool = False,
) -> Tuple[float, Optional[float]]:
    """
    Refine rest_length_frac around analytic guess using electron stability.

    This function scans multiple rest_length_frac candidates around the
    analytic initial guess and selects the one with best electron stability.

    Process:
    1. Generate n_scan candidate values in [initial - radius, initial + radius]
    2. For each candidate:
       - Update material parameters
       - Run short electron stability test (~3 Compton periods)
       - Compute combined loss from stability metrics
    3. Return candidate with lowest loss

    Args:
        mapper: DimensionalMapper with physical constants
        k_phys: Spring constant in SI [N/m]
        h_phys: Grid spacing in SI [m]
        rest_length_frac_initial: Analytic guess from tension formula
        scan_radius: +/- range around initial guess (default 0.05)
        n_scan: Number of candidate values to test (default 5)
        run_refinement: If False, skip refinement and return initial guess
                       (default False - refinement requires full sim setup)

    Returns:
        (rest_length_frac_refined, best_loss): Best value and its loss
            If run_refinement=False, returns (initial, None)
    """
    if not run_refinement:
        warnings.warn(
            "Electron refinement disabled (run_refinement=False). "
            "Using analytic rest_length_frac without stability optimization."
        )
        return rest_length_frac_initial, None

    # Generate candidate values
    candidates = np.linspace(
        rest_length_frac_initial - scan_radius,
        rest_length_frac_initial + scan_radius,
        num=n_scan,
    )

    best_frac = rest_length_frac_initial
    best_loss = float("inf")

    print(f"\n=== Electron Refinement Scan ===")
    print(f"  Initial guess: {rest_length_frac_initial:.6f}")
    print(f"  Scan range: [{candidates[0]:.6f}, {candidates[-1]:.6f}]")
    print(f"  Testing {n_scan} candidates...\n")

    for i, frac in enumerate(candidates):
        # Enforce physical bounds [0.7, 0.98]
        if not (0.7 <= frac <= 0.98):
            print(f"  [{i+1}/{n_scan}] rest_length_frac={frac:.6f} - SKIPPED (out of bounds)")
            continue

        print(f"  [{i+1}/{n_scan}] Testing rest_length_frac={frac:.6f}...")

        # TODO: Update material parameters for this candidate
        # rest_length = frac * h_phys
        # ... update simulation config ...

        # Run short stability test
        try:
            metrics = run_short_electron_stability_test(...)
            loss = combined_rest_length_loss(metrics)

            print(f"    Loss: {loss:.6e}")
            print(f"      energy_leakage: {metrics['energy_leakage']:.4f}")
            print(f"      shape_drift: {metrics['shape_drift']:.6e}")
            print(f"      mode_purity_loss: {metrics['mode_purity_loss']:.4f}")

            if "lateralization_ratio" in metrics:
                print(f"      lateralization_ratio: {metrics['lateralization_ratio']:.4f}")

            if loss < best_loss:
                best_loss = loss
                best_frac = frac
                print(f"    *** New best! ***")

        except Exception as e:
            print(f"    ERROR: {e}")
            continue

    print(f"\n=== Refinement Complete ===")
    print(f"  Best rest_length_frac: {best_frac:.6f}")
    print(f"  Best loss: {best_loss:.6e}")
    print(f"  Δ from initial: {best_frac - rest_length_frac_initial:+.6f}\n")

    return best_frac, best_loss


def calibrate_material_from_physics(
    rho_m: float,
    constants: PhysicalConstants,
    k_phys: float,
    h_phys: float,
    mapper: DimensionalMapper,
    run_electron_refinement: bool = False,
) -> MaterialCalibrationResult:
    """
    Main entry point: calibrate all material parameters from physics.

    This function implements the complete three-layer calibration:
    1. Compute continuum targets (T, κ) from fundamental constants
    2. Derive analytic rest_length_frac from tension requirement
    3. Optionally refine using electron stability experiment

    Args:
        rho_m: Mass density [kg/m²] from electron calibration
        constants: PhysicalConstants with c, G, hbar, m_e
        k_phys: Spring constant [N/m] - chosen for numerical stability
        h_phys: Grid spacing [m] - chosen for resolution
        mapper: DimensionalMapper for unit conversions
        run_electron_refinement: Whether to run stability refinement (default False)

    Returns:
        MaterialCalibrationResult with all calibrated parameters
    """
    print("\n" + "="*70)
    print("MATERIAL CALIBRATION FROM PHYSICS")
    print("="*70)

    # Step 1: Compute continuum targets
    targets = compute_continuum_targets(rho_m, constants)
    print(f"\n{targets}")

    # Step 2: Analytic first guess from tension
    print("\n--- Analytic Rest Length Derivation ---")
    rest_length_frac_analytic = solve_rest_length_frac_from_tension(
        T_target=targets.T_target,
        k_phys=k_phys,
        h_phys=h_phys,
    )
    print(f"  Analytic rest_length_frac = {rest_length_frac_analytic:.6f}")

    # Validate and clamp if needed
    rest_length_frac_clamped, is_valid = validate_rest_length_frac(
        rest_length_frac_analytic
    )
    if not is_valid:
        print(f"  Clamped to: {rest_length_frac_clamped:.6f}")

    # Step 3: Optional electron refinement
    refinement_loss = None
    if run_electron_refinement:
        print("\n--- Electron Stability Refinement ---")
        rest_length_frac_final, refinement_loss = calibrate_rest_length_frac_with_electron(
            mapper=mapper,
            k_phys=k_phys,
            h_phys=h_phys,
            rest_length_frac_initial=rest_length_frac_clamped,
            run_refinement=True,
        )
    else:
        rest_length_frac_final = rest_length_frac_clamped
        print("\n--- Electron Refinement Skipped ---")
        print(f"  Using analytic value: {rest_length_frac_final:.6f}")

    # Compute physical rest length
    rest_length_phys = rest_length_frac_final * h_phys

    # Package results
    result = MaterialCalibrationResult(
        rest_length_frac=rest_length_frac_final,
        rest_length_frac_analytic=rest_length_frac_analytic,
        rest_length_phys=rest_length_phys,
        kappa_phys=targets.kappa_target,
        continuum_targets=targets,
        refinement_loss=refinement_loss,
    )

    print("\n" + "="*70)
    print("CALIBRATION COMPLETE")
    print("="*70)
    print(f"\n{result}\n")

    return result