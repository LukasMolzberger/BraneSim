# Brane Simulation Calibration Tools

This directory contains standalone tools for calibrating physical parameters in the brane simulation.

## Rest-Length Calibration Tool

### Overview

The `calibrate_rest_length.py` script computes a physically plausible dimensionless spring `rest_length` parameter (L₀/a) by enforcing the continuum relations:

- **Wave speed constraint**: c² = T / ρ_m
- **Gravitational coupling**: κ ≈ c³ / (8π G)

The calibration uses the 1D spring-chain derivation to connect microscopic pre-strain to continuum tension, yielding an analytic solution for `rest_length` that satisfies the relativistic constraint.

### Physics Background

From the paper's 1D chain analysis with cross-section a²:

```
T₀ = k × Δ₀              (tension force in 1D chain)
T  ≈ T₀ / a² = k(a - L₀)/a²
   = (k/a) × (1 - rest_length)
```

where:
- `k` is the microscopic spring constant [N/m]
- `a` is the lattice spacing [m]
- `L₀` is the spring rest length [m]
- `rest_length := L₀/a` (dimensionless)
- `Δ₀ = a - L₀` is the pre-strain [m]

Given a target tension `T_target = ρ_m c²`, we can solve analytically:

```
rest_length = 1 - (T_target × a) / k
```

### Usage

#### Basic usage

```bash
python -m branesim.tools.calibrate_rest_length \
  --rho-m 1.0e6 \
  --spring-k 1.0e14 \
  --lattice-spacing 1.0e-10
```

#### With JSON output

```bash
python -m branesim.tools.calibrate_rest_length \
  --rho-m 1.0e6 \
  --spring-k 1.0e14 \
  --lattice-spacing 1.0e-10 \
  --output-json branesim/config/rest_length_calibration.json
```

#### Advanced parameters

```bash
python -m branesim.tools.calibrate_rest_length \
  --rho-m 1.0e6 \
  --spring-k 1.0e14 \
  --lattice-spacing 1.0e-10 \
  --mu-G 0.5 \
  --ell-star-factor 2.0 \
  --output-json config/calibration.json
```

### Parameters

#### Required

- `--rho-m`: Brane mass density ρ_m [kg/m³]
  - From Compton-scale calibration: ρ_m λ_C³ ≈ m_e
  - Typical range: 10⁶–10¹² kg/m³

- `--spring-k`: Microscopic spring constant k [N/m]
  - Must be large enough: k > ρ_m c² a
  - Typical range: 10¹⁰–10¹⁶ N/m

- `--lattice-spacing`: Lattice spacing a [m]
  - Should be sub-Compton: a < λ_C ≈ 3.86×10⁻¹³ m
  - Typical range: 10⁻¹³–10⁻⁹ m

#### Optional

- `--mu-G`: Bending-to-tension factor μ_G (default: 1.0)
  - Dimensionless geometric coefficient
  - Used in bending stiffness estimate: κ_micro ≈ μ_G T ℓ_*²

- `--ell-star-factor`: Coarse-graining length factor (default: 1.0)
  - Defines ℓ_* = ell_star_factor × a
  - Used in gravitational sector analysis

- `--output-json`: Path for JSON output (optional)
  - If provided, writes full calibration data to file

### Output

The script prints a detailed summary including:

1. **Input parameters**: ρ_m, k, a, μ_G, ℓ_* factor
2. **Target continuum constants**: T_target, κ_target from physical constants
3. **Calibration result**: `rest_length` (L₀/a)
4. **Achieved values**: T_micro, κ_micro
5. **Consistency checks**:
   - T_micro / T_target (should be 1.0)
   - κ_micro / κ_target ratio
   - Effective geometric prefactor ζ_G

Example output:

```
======================================================================
  REST-LENGTH CALIBRATION SUMMARY
======================================================================

Input Parameters:
  ρ_m (brane density)  : 1.000000e+06 kg/m³
  k (spring constant)  : 1.000000e+14 N/m
  a (lattice spacing)  : 1.000000e-10 m

Calibration Result:
  rest_length (L₀/a)       : 0.9101244821
  Pre-strain Δ₀/a          : 0.0898755179

Consistency Checks:
  ✓ Tension constraint satisfied (T = ρ_m c²)
  ⚠ ζ_G = 8.94e+50 suggests geometric prefactor outside O(1)
    Consider adjusting μ_G or ell_star_factor
======================================================================
```

### JSON Output Format

When `--output-json` is specified, the script writes:

```json
{
  "rest_length": 0.9101244821,
  "pre_strain_ratio": 0.0898755179,
  "rho_m": 1000000.0,
  "spring_k": 1e14,
  "lattice_spacing": 1e-10,
  "T_target": 8.987552e22,
  "T_micro": 8.987552e22,
  "kappa_target": 1.606263e34,
  "kappa_micro": 898.7551787,
  "zeta_G_effective": 8.936040e50,
  "calibration_method": "analytic_1D_chain",
  "physical_constants": {
    "c": 299792458.0,
    "G": 6.6743e-11
  }
}
```

## Loading Calibrated Parameters

Use the `calibration_loader` module to load calibrated values in your simulation code:

```python
from branesim.config.calibration_loader import load_rest_length_calibration

# Load calibration (returns default if file not found)
rest_length, cal_data = load_rest_length_calibration()

# Or just get the value
from branesim.config.calibration_loader import get_rest_length
rest_length = get_rest_length(default=0.95)

# Access full calibration data
if cal_data is not None:
    print(f"Rest length: {cal_data.rest_length}")
    print(f"Pre-strain: {cal_data.pre_strain_ratio}")
    print(f"Tension: {cal_data.T_micro} J/m³")
    print(f"Method: {cal_data.calibration_method}")
```

## Parameter Selection Guidelines

### For Compton-scale simulations

Starting from electron properties:
- λ_C ≈ 3.86×10⁻¹³ m (Compton wavelength)
- m_e ≈ 9.11×10⁻³¹ kg (electron mass)
- Calibration: ρ_m λ_C³ ≈ m_e → ρ_m ≈ 1.58×10¹³ kg/m³

Choose:
- `lattice_spacing`: a ≈ 0.1 λ_C = 3.86×10⁻¹⁴ m
- `rho_m`: 1.58×10¹³ kg/m³ (from Compton calibration)
- `spring_k`: Solve for desired rest_length

Example: For rest_length ≈ 0.9:
```
T_target = ρ_m c² ≈ 1.42×10³⁰ J/m³
rest_length = 1 - (T_target × a) / k
0.9 = 1 - (1.42×10³⁰ × 3.86×10⁻¹⁴) / k
k = 5.48×10¹⁶ / 0.1 ≈ 5.5×10¹⁷ N/m
```

### Troubleshooting

**Error: "rest_length is outside (0,1)"**

This means the parameters are incompatible. The spring constant is too weak to support the required tension. Solutions:

1. **Increase k** (stiffer springs)
2. **Decrease ρ_m** (lower mass density)
3. **Increase a** (coarser lattice)

The constraint is:
```
k > ρ_m c² a
```

For rest_length ≈ 0.9, you need:
```
k ≈ 10 × ρ_m c² a
```

## Physical Interpretation

### Rest length and pre-strain

- `rest_length = 0.9` → springs are stretched 10% beyond rest length
- `rest_length = 0.95` → springs are stretched 5% beyond rest length
- Smaller rest_length → larger pre-strain → higher tension

### Geometric prefactor ζ_G

The effective ζ_G tells you how the microscopic bending stiffness compares to the gravitational target:

- ζ_G ≈ 1: Microscopic model naturally reproduces correct G
- ζ_G >> 1: κ_micro too small; increase μ_G or ell_star_factor
- ζ_G << 1: κ_micro too large; decrease μ_G or ell_star_factor

This is diagnostic only - the script doesn't enforce ζ_G ≈ 1, since the simple coarse-graining ansatz is approximate.

## References

See paper sections:
- Parameter matching: Section 3.4
- 1D chain derivation: Experimental setting section
- Gravitational coupling: Appendix on emergent gravity