# Compton-Cell Calibration Integration

This document describes the integration of the Compton-cell mass calibration procedure from the paper into the BraneSim Python codebase.

## Overview

The **amplitude scale calibration** (Route i from the paper, Section "Amplitude scale calibration" in `experimental-setting.tex`) has been implemented as a Python module that computes physically-motivated brane lattice parameters.

## Theory

The Compton-cell calibration ties the discrete brane lattice to fundamental physical constants:

### Continuum Assumptions

- **Reduced Compton wavelength**: λ_C = ℏ/(m_e c) ≈ 3.86×10⁻¹³ m
- **Brane volume mass density**: ρ = m_e / λ_C³ ≈ 1.58×10⁷ kg/m³
- **Brane bulk modulus**: K = ρ c² ≈ 1.42×10²⁴ Pa

### Discrete Mapping (cubic lattice with spacing h)

- **Point mass**: m_point = ρ h³
- **Axial spring constant**: k_spring = K h

This ensures the lattice reproduces the continuum wave speed: **c² = K/ρ**

## Implementation

### 1. Core Function: `brane_lattice_params_3d()`

**Location**: `branesim/physics/parameters.py`

```python
from branesim.physics import brane_lattice_params_3d

# Compute parameters for grid spacing h = 10 λ_C
params = brane_lattice_params_3d(
    grid_spacing_m=1e-12,
    use_compton_default=True  # Use Compton-cell assumption
)

# Returns:
# {
#     "lambda_C": 3.86e-13,      # Reduced Compton wavelength (m)
#     "rho": 1.58e7,             # Mass density (kg/m³)
#     "K": 1.42e24,              # Bulk modulus (Pa)
#     "m_point": 1.58e-29,       # Point mass (kg)
#     "k_spring": 1.42e12        # Spring constant (N/m)
# }
```

### 2. Convenience Method: `SimulationConfig.from_compton_calibration()`

**Location**: `branesim/config/simulation_config.py`

```python
from branesim.config.simulation_config import SimulationConfig
from branesim.core.state import Dimensionality

# Create a complete simulation configuration
config = SimulationConfig.from_compton_calibration(
    grid_shape=(32, 32, 32),           # 32³ grid
    dimension=Dimensionality.THREE_D,
    lambda_C_multiplier=10.0,          # h = 10 λ_C
    cfl_factor=0.4,                    # CFL = 0.4 for stability
    critical_strain=0.1,               # ε_cr for saturation (optional)
    device='cpu',
    dtype='float64'
)

# Config is ready to use with correct mass density, spring constant, etc.
```

### 3. Helper: `print_calibration_summary()`

```python
from branesim.physics import brane_lattice_params_3d, print_calibration_summary

params = brane_lattice_params_3d(1e-12)
print_calibration_summary(params, grid_spacing=1e-12)

# Prints formatted summary of all parameters
```

## Files Modified/Created

### New Files

1. **`branesim/physics/parameters.py`** (new)
   - Core calibration function `brane_lattice_params_3d()`
   - Summary printer `print_calibration_summary()`
   - Comprehensive documentation with paper references

2. **`examples/compton_calibration_demo.py`** (new)
   - Usage examples
   - Parameter scans
   - Custom density examples

3. **`test_calibration_simple.py`** (new)
   - Integration test
   - Verifies wave speed matches c to machine precision

### Modified Files

1. **`branesim/physics/__init__.py`**
   - Exports `brane_lattice_params_3d` and `print_calibration_summary`
   - Graceful handling of missing torch dependency

2. **`branesim/config/simulation_config.py`**
   - Added `SimulationConfig.from_compton_calibration()` class method
   - Integrated with existing configuration system
   - Automatic CFL checking and wave speed verification

## Usage Examples

### Example 1: Direct Calibration

```python
from branesim.physics import brane_lattice_params_3d

# Grid spacing at 10× Compton wavelength
lambda_C = 3.8616e-13  # meters
h = 10.0 * lambda_C

params = brane_lattice_params_3d(h)

# Use in your simulation setup
mass_per_point = params['m_point']      # 9.1094e-28 kg
spring_constant = params['k_spring']    # 5.4903e+12 N/m
```

### Example 2: Complete Simulation Config

```python
from branesim.config.simulation_config import SimulationConfig
from branesim.core.state import Dimensionality

config = SimulationConfig.from_compton_calibration(
    grid_shape=(64, 64, 64),
    dimension=Dimensionality.THREE_D,
    lambda_C_multiplier=20.0,  # h = 20 λ_C for better resolution
    critical_strain=0.05
)

# Verify wave speed
expected_c, computed_c, error = config.verify_wave_speed()
print(f"Wave speed error: {error*100:.6f}%")  # < 0.0001%
```

### Example 3: Parameter Scan

```python
from branesim.physics import brane_lattice_params_3d

lambda_C = 3.8616e-13

for multiplier in [5, 10, 20, 50, 100]:
    h = multiplier * lambda_C
    params = brane_lattice_params_3d(h)
    print(f"h = {multiplier:.0f} λ_C: m = {params['m_point']:.2e} kg, "
          f"k = {params['k_spring']:.2e} N/m")
```

## Testing

Run the integration test:

```bash
python test_calibration_simple.py
```

Expected output:
```
✓ PASSED: Wave speed matches c to machine precision!
✓ ALL TESTS PASSED!
```

## Physical Validation

The calibration has been verified to reproduce the speed of light **exactly** (to machine precision):

- **Input**: Compton-cell assumption ρ = m_e / λ_C³
- **Derived**: K = ρ c²
- **Verified**: c_computed = √(K/ρ) = c_expected (error < 10⁻¹⁰)

## Paper References

This implementation corresponds to:

- **Paper Section**: "Amplitude scale calibration" (experimental-setting.tex, lines 238-293)
- **Related Section**: "Related Works" addition on Compton-scale models (introduction.tex, lines 111-125)
- **Theoretical Basis**: Section "Amplitude scale and units of the fourth coordinate" (reconstructing-physics.tex)

## Next Steps

1. **Use in simulations**: Replace ad-hoc parameter choices with `from_compton_calibration()`
2. **Validate numerically**: Compare discrete wave propagation with continuum prediction
3. **Extend to Route (ii)**: Implement charge-based calibration (Section "Charge and internal energy density")

## Dependencies

- **Core function**: Python 3.8+ (no external dependencies)
- **SimulationConfig**: Requires numpy
- **Full simulation**: Requires torch (for force computation)

---

*This integration bridges the theoretical amplitude scale calibration from the paper with the practical discrete simulation framework.*