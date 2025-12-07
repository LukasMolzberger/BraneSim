# Physical Calibration of `rest_length` Parameter

This document describes the three-layer calibration approach for deriving a physically realistic `rest_length` parameter for the discrete spring network, based on continuum physics targets.

## Overview

The goal is to replace the empirical "pick something between 0.8 and 0.95 and hope for the best" approach with a physics-anchored calibration that:

1. **Uses continuum physics targets** from fundamental constants (c, G, m_e)
2. **Derives an analytic first guess** from tension requirements
3. **Refines using the electron experiment** to ensure geometric nonlinearity threshold matches Compton-scale physics

## 1. Continuum Physics Targets

We derive the target continuum elastic parameters from fundamental physics requirements.

### 1.1 Mass Density (ρ_m)

The brane mass density is calibrated so that a Compton-scale localized mode (the electron soliton) has the correct rest energy:

```
E_electron = m_e c²
```

This is achieved through the existing electron calibration pipeline in `electron_initialization.py`. After `calibrate_electron_init_params`, the effective mass density ρ_m can be read from the `DimensionalMapper` and material config.

### 1.2 Target Tension (T_target) - For Speed of Light

Small-amplitude waves on a pre-tensioned brane should propagate at speed c. In the continuum limit, for an isotropic brane model:

```
c² = T / ρ_m  ⟹  T_target = ρ_m c²
```

This ensures that linear photon modes travel at the correct speed.

**Physical interpretation:** The brane is under uniform tension T (force per unit length), similar to a drumhead. Wave speed on the membrane is determined by the ratio of tension to mass density.

### 1.3 Target Bending Stiffness (κ_target) - For Gravity

From the paper's analysis of Newton's constant in the static weak-field limit:

```
κ/2 ≈ c³/(16π G)  ⟹  κ_target = c³/(8π G)
```

This bending rigidity κ appears in the brane Lagrangian as the coefficient of the curvature term. It sets the coupling strength between brane curvature and gravitational effects.

**Physical interpretation:** κ is analogous to the bending rigidity of a thin elastic sheet - it resists curvature. This resistance to bending is what generates the effective gravitational force in this emergent gravity model.

**Implementation note:** In the PyTorch simulation, κ is a separate material parameter from `rest_length`. It should be set directly to κ_target (via DimensionalMapper → dimensionless units). The `rest_length` parameter will be used to match the tension T.

### Summary of Continuum Targets

After the electron calibration determines ρ_m:

- **ρ_m**: Mass density [kg/m²] - from electron calibration
- **T_target = ρ_m c²**: Tension [N/m] - for correct wave speed
- **κ_target = c³/(8π G)**: Bending stiffness [J·m] - for correct gravitational coupling

## 2. Analytic First Guess for `rest_length`

Now we connect the continuum tension T to the discrete spring parameters.

### 2.1 Spring Network Conventions

From the codebase:

```python
rest_length = rest_length_frac * h  # h = grid spacing in physical units
```

where `rest_length_frac ∈ [0.8, 0.95]` is the dimensionless pre-strain knob.

### 2.2 Continuum Limit of Discrete Springs

For a cubic spring network with:
- Lattice spacing: h_phys
- Spring constant: k_phys [N/m]
- Rest length: L_rest = rest_length_frac × h_phys

Each spring is pre-stretched by:
```
ΔL_0 = h_phys - L_rest
```

The pre-stretch force in each spring:
```
F_0 = k_phys ΔL_0
```

Each spring spans an area ~h_phys², so the background tension per unit area:
```
T ≈ F_0 / h_phys² = k_phys(h_phys - L_rest) / h_phys²
```

### 2.3 Solving for rest_length_frac

Setting T = T_target:

```
T_target = k_phys(h_phys - L_rest) / h_phys²
```

Rearranging:
```
h_phys - L_rest = T_target h_phys² / k_phys
```

Dividing by h_phys and using L_rest = rest_length_frac × h_phys:
```
1 - rest_length_frac = T_target h_phys² / k_phys
```

**Final formula:**
```
rest_length_frac = 1 - (T_target h_phys²) / k_phys
```

### 2.4 Implementation: Analytic Solver

```python
def solve_rest_length_frac_from_tension(
    T_target: float,
    k_phys: float,
    h_phys: float,
) -> float:
    """
    Solve for rest_length_frac such that the discrete spring pre-strain
    reproduces the target continuum tension T_target.

    Uses: T ≈ k_phys * (h_phys - L_rest) / h_phys²
          L_rest = rest_length_frac * h_phys

    Args:
        T_target: Target continuum tension [N/m]
        k_phys: Spring constant in SI [N/m]
        h_phys: Grid spacing in SI [m]

    Returns:
        rest_length_frac: Dimensionless pre-strain factor [0, 1]
    """
    rest_length_frac = 1.0 - (T_target * h_phys**2) / k_phys
    return rest_length_frac
```

### 2.5 Validation and Constraints

The derived `rest_length_frac` should satisfy:
```
0.8 ≤ rest_length_frac ≤ 0.95
```

This range is empirically known to work well for the electron bottleneck geometry.

If the analytic value falls outside this range:
1. Warn the user
2. Adjust k_phys and/or h_phys, OR
3. Clamp rest_length_frac and flag that the continuum calibration is approximate

## 3. Refinement via Electron Stability Experiment

The analytic formula guarantees correct linear-wave speed (c) and gravitational coupling (via κ_target), but "realistic" also means the geometric nonlinearity threshold sits at the right scale.

### 3.1 The Electron as a Geometric Nonlinearity Probe

The electron torus is designed to be:
- A self-confining soliton via geometric nonlinearity
- Near the "bottleneck" threshold where lateral motion becomes significant
- Stable over multiple Compton periods

The pre-strain (rest_length_frac) affects:
- Where the geometric nonlinearity "kicks in"
- The balance between amplitude oscillations and lateral motion
- The stability of the toroidal electron configuration

### 3.2 Stability Metrics for rest_length Calibration

From `electron_stability.py`, we use:

1. **energy_leakage**: Fraction of energy radiating away (want < 0.1)
2. **shape_drift**: Change in envelope structure over time (want small)
3. **mode_purity_loss**: 1 - (power at ω_C / total power) (want small)
4. **lateralization_ratio** (new): Ratio of lateral kinetic energy to amplitude energy

The sweet spot has:
- Minimal energy leakage and shape drift
- High mode purity (dominated by Compton frequency)
- Lateralization ratio ~ O(1) (near-threshold, not collapsed or dispersed)

### 3.3 Combined Loss Function

```python
def combined_rest_length_loss(metrics: Dict[str, float]) -> float:
    """
    Combine stability metrics into scalar loss for rest_length optimization.

    Args:
        metrics: Dictionary containing:
            - energy_leakage: [0, 1+]
            - shape_drift: [0, ∞)
            - mode_purity_loss: [0, 1]
            - lateralization_ratio: (optional) [0, ∞)

    Returns:
        Combined loss (lower is better)
    """
    # Weights (tunable)
    w_leak = 1.0
    w_shape = 1.0
    w_mode = 1.0
    w_lat = 0.5

    L = 0.0
    L += w_leak * metrics["energy_leakage"]
    L += w_shape * metrics["shape_drift"]
    L += w_mode * metrics["mode_purity_loss"]

    # Lateralization: target ~ 1.0 (near-threshold)
    if "lateralization_ratio" in metrics:
        lat = metrics["lateralization_ratio"]
        L += w_lat * (lat - 1.0)**2

    return L
```

### 3.4 Refinement Procedure

```python
def calibrate_rest_length_frac_with_electron(
    mapper: DimensionalMapper,
    k_phys: float,
    h_phys: float,
    rest_length_frac_initial: float,
    scan_radius: float = 0.05,
    n_scan: int = 5,
) -> float:
    """
    Refine rest_length_frac around analytic guess using electron stability.

    1. Start from analytic initial guess
    2. Scan n_scan values in [initial - radius, initial + radius]
    3. For each candidate:
       - Update material parameters
       - Run short electron stability test (~3 Compton periods)
       - Compute combined loss
    4. Return rest_length_frac with lowest loss

    Args:
        mapper: DimensionalMapper with physical constants
        k_phys: Spring constant in SI [N/m]
        h_phys: Grid spacing in SI [m]
        rest_length_frac_initial: Analytic guess from tension formula
        scan_radius: +/- range around initial guess
        n_scan: Number of candidate values to test

    Returns:
        rest_length_frac_refined: Best value from stability test
    """
```

The refinement ensures that the pre-strain is consistent with both:
- **Linear wave physics** (via analytic T = ρ_m c²)
- **Nonlinear soliton physics** (via electron bottleneck stability)

## 4. Implementation Workflow

### Step 1: Material Calibration Setup

In the physical parameter layer (extend `DimensionalMapper` or create new `MaterialCalibration` class):

```python
# After electron calibration determines ρ_m:
rho_m_phys = compute_mass_density_from_electron(...)

# Compute continuum targets
c_phys = constants.c
G_phys = 6.67430e-11  # Newton's constant [m³/(kg·s²)]

T_target = rho_m_phys * c_phys**2  # Tension for wave speed
kappa_target = c_phys**3 / (8 * np.pi * G_phys)  # Bending for gravity
```

### Step 2: Analytic First Guess

```python
# Choose k_phys and h_phys (from grid setup)
rest_length_frac_initial = solve_rest_length_frac_from_tension(
    T_target=T_target,
    k_phys=k_phys,
    h_phys=h_phys,
)

# Validate range
if not (0.8 <= rest_length_frac_initial <= 0.95):
    warnings.warn(f"Derived rest_length_frac={rest_length_frac_initial:.3f} "
                  f"out of [0.8, 0.95]. Consider adjusting k_phys or h_phys.")
    rest_length_frac_initial = np.clip(rest_length_frac_initial, 0.8, 0.95)
```

### Step 3: Electron Refinement

```python
rest_length_frac_refined = calibrate_rest_length_frac_with_electron(
    mapper=mapper,
    k_phys=k_phys,
    h_phys=h_phys,
    rest_length_frac_initial=rest_length_frac_initial,
    scan_radius=0.05,
    n_scan=5,
)

logger.info(f"rest_length_frac: analytic={rest_length_frac_initial:.4f}, "
           f"refined={rest_length_frac_refined:.4f}")
```

### Step 4: Store in Material Config

```python
# Store final calibrated value
material_config.rest_length_frac = rest_length_frac_refined
material_config.rest_length = rest_length_frac_refined * h_phys

# Set bending stiffness separately
material_config.kappa_phys = kappa_target
material_config.kappa_sim = mapper.to_sim_bending_stiffness(kappa_target)
```

## 5. Physics Summary

This three-layer approach ensures:

1. **Photon sector consistency**: Linear waves propagate at c via T = ρ_m c²
2. **Gravitational consistency**: Curvature coupling gives correct G via κ = c³/(8π G)
3. **Electron sector consistency**: Geometric nonlinearity threshold matches Compton-scale soliton physics

The calibrated `rest_length` is physically anchored to (c, G, m_e) and the continuum brane model, rather than being an arbitrary tuning parameter.

## References

- Williamson & van der Mark toroidal electron model
- Paper section on "Connection between pre-tension and the gravitational constant"
- Existing code: `electron_initialization.py`, `electron_stability.py`, `dimensional_mapping.py`