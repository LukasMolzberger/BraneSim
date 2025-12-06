# Electron Initialization Experiment

## Overview

This experiment implements the **Williamson & van der Mark (W&vdM) toroidal double-loop electron model** adapted to the BraneSim brane framework. It provides a complete pipeline for:

1. **Analytical initialization** of a toroidal electron soliton
2. **Stability measurement** during time evolution
3. **Optimization framework** for refining parameters

## Theoretical Background

### The W&vdM Toroidal Electron

The electron is modeled as a **topologically stable soliton** with:

- **Geometry**: Toroidal structure with major radius R ~ λ_C / 2π (Compton scale)
- **Double loop**: Two intertwined paths forming a torus knot (m=2, n=1)
- **Internal oscillation**: Compton frequency ω_C = m_e c² / ℏ
- **Self-confinement**: Via geometric nonlinearity ("bottleneck" mechanism)

### Key Physical Properties

From the rest frame:

- **Mass**: m_e = 9.109×10⁻³¹ kg (from total energy E = m_e c²)
- **Charge**: -e = -1.602×10⁻¹⁹ C (from time-averaged amplitude X̄⁴)
- **Spin**: ℏ/2 = 5.273×10⁻³⁵ J·s (from internal angular momentum)
- **Compton wavelength**: λ_C = 2.426×10⁻¹² m

### Tubular Coordinate System

The implementation uses a **change of coordinates** to simplify the toroidal structure:

```
Global brane coords (x, y, z, w) → Tubular coords (z, x, y, w)
```

Where:
- **z**: Arclength along the torus centerline [0, 2πR)
- **x**: Radial transverse coordinate (in-plane)
- **y**: Vertical transverse coordinate (binormal)
- **w**: Amplitude dimension (unchanged)

In these coordinates:
- The torus appears as a "straight tunnel"
- The double loop appears as two rotating peaks in the (x, y) cross-section
- Curvature effects enter through metric corrections g_ab(x, y)

### The Geometric Nonlinearity ("Bottleneck")

The key mechanism for confinement is the **amplitude-to-lateral energy transfer**:

1. High amplitude gradients (∂_z ξ) in the 4D embedding increase spring extensions
2. Beyond a critical threshold, it becomes energetically favorable to develop lateral distortions
3. This creates a self-consistent "tunnel" that confines the mode
4. The threshold is controlled by **rest_length** (pre-tension parameter)

Critical relation:
```
|∂_z ξ|_electron ~ |∂_z ξ|_crit(rest_length, k, h)
```

## Implementation Structure

### Core Modules

#### `branesim/physics/electron_initialization.py`

Main initialization module providing:

- **`ElectronInitParams`**: Dataclass for all initialization parameters
- **`compute_tubular_coords_*`**: Mapping to tubular coordinates
- **`double_loop_envelope`**: Cross-section profile f(x, y) with two lobes
- **`init_electron_amplitude`**: Initialize X⁴ field with Compton oscillation
- **`init_electron_state`**: Main entry point for electron setup
- **`calibrate_electron_init_params`**: Physical parameter calibration

#### `branesim/diagnostics/electron_stability.py`

Stability measurement and loss computation:

- **`compute_energy_leakage`**: Fraction of energy radiating away
- **`compute_shape_drift`**: Change in envelope over time
- **`compute_mode_purity`**: FFT-based frequency analysis
- **`compute_constraint_penalties`**: Energy, momentum, spin, charge errors
- **`compute_electron_stability_loss`**: Combined loss function for optimization

### Experiment Scripts

#### `experiments/electron_stability_test.py`

Basic stability test that:
1. Initializes a 3D grid at Compton scale
2. Creates an electron with calibrated parameters
3. Runs a short simulation (3 Compton periods)
4. Measures stability metrics
5. Produces visualizations

## Running the Experiment

### Basic Run

```bash
cd experiments
python electron_stability_test.py
```

This will:
- Create a 20×20×20 grid (small for testing)
- Initialize an electron at the center
- Run for 3 Compton periods
- Output stability metrics
- Generate `electron_evolution.png`

### Expected Output

```
=== Electron Tube Region ===
  Points in electron region: 234/8000 (2.93%)

=== Initialized Electron Amplitude Field ===
  Torus major radius R = 3.862e-13 m
  Cross-section radius rho0 = 1.159e-13 m
  Max |ξ| = 1.000e-13 m
  Max |v_ξ| = 2.398e+08 m/s

=== Stability Results ===
Total Loss: 0.52
  Energy leakage: 0.05
  Shape drift: 1.2e-03
  Mode purity loss: 0.42
```

### Interpreting Results

**Good electron**:
- Energy leakage < 0.1 (< 10% radiation)
- Shape drift < 0.01 (< 1% change)
- Mode purity loss < 0.2 (> 80% power at ω_C)

**Current limitations** (first-pass ansatz):
- Energy may not match m_e c² exactly
- Amplitude scale not yet calibrated to charge = -e
- No optimization yet - this is the raw analytical guess

## Parameter Calibration

### Geometric Parameters

From `calibrate_electron_init_params`:

```python
R = λ_C / (2π)           # Major radius ~ 386 pm
rho0 = 0.3 * R           # Lobe radius ~ 116 pm
sigma_r = 0.1 * R        # Lobe width ~ 39 pm
```

These are default choices. Optimization will refine them.

### Amplitude Scale

The amplitude A controls:
- Internal energy density
- Effective charge (via X̄⁴)

Initial guess: A ~ 0.1 pm

**Calibration procedure** (to be implemented):

1. Measure effective charge Q(A) for trial amplitude
2. Adjust A to match Q = -e
3. Check if energy E ≈ m_e c²
4. If not, adjust geometry (rho0, sigma_r) and iterate

### Rest Length

Critical parameter controlling the bottleneck:

```python
rest_length = rest_length_frac * h
```

Typical range: `rest_length_frac ∈ [0.8, 0.95]`

- **Too small** (< 0.8): Strong pre-tension, mode doesn't confine
- **Too large** (> 0.95): Weak pre-tension, floppy network
- **Sweet spot**: Mode operates near geometric nonlinearity threshold

## Optimization Framework

### Objective

Minimize stability loss:

```python
L = w_leak * L_leak
  + w_shape * L_shape
  + w_mode * L_mode
  + w_E * (E - m_e c²)²
  + w_P * |P|²
  + w_S * (|S| - ℏ/2)²
```

### Optimization Parameters (θ)

Suggested parameters to optimize:

1. **Geometry**: R, rho0, sigma_r (within bounds)
2. **Amplitude**: A (constrained by charge)
3. **Phase**: phase_offset
4. **Spin**: lateral_spin_scale (for angular momentum)

### Methods

#### Derivative-Free (Recommended First)

```python
from scipy.optimize import minimize
# or use CMA-ES for better global search

def loss_function(theta):
    params = decode_theta_to_params(theta)
    state = init_electron(params)
    states = run_simulation(state, n_steps=1000)
    loss, metrics = compute_electron_stability_loss(states, ...)
    return loss

result = minimize(loss_function, theta0, method='Nelder-Mead')
```

#### Gradient-Based (Advanced)

If the simulation is differentiable through time:

```python
import torch.optim as optim

theta = torch.tensor([R, rho0, sigma_r, A], requires_grad=True)
optimizer = optim.Adam([theta], lr=1e-14)

for epoch in range(100):
    optimizer.zero_grad()
    loss = compute_loss_differentiable(theta)
    loss.backward()
    optimizer.step()
```

## Grid Resolution Requirements

### Minimum Requirements

For a stable electron at Compton scale:

- **Spatial resolution**: h < λ_C / 5 (at least 5 points per wavelength)
- **Grid extent**: Domain > 3λ_C in each dimension (buffer for radiation)
- **Time step**: dt < 0.2 h / c (CFL condition)
- **Simulation time**: T > 10 T_C (multiple Compton periods)

### Production Run Suggestions

| Grid Size | Points | Memory | Time/step | Purpose |
|-----------|--------|---------|-----------|---------|
| 20³ | 8k | ~1 MB | Fast | Quick testing |
| 40³ | 64k | ~10 MB | Medium | Initial optimization |
| 80³ | 512k | ~80 MB | Slow | High resolution |
| 160³ | 4M | ~640 MB | Very slow | Production |

## Physical Scales

### SI Units (automatically handled by DimensionalMapper)

| Quantity | Value | Unit |
|----------|-------|------|
| λ_C | 2.426×10⁻¹² | m |
| T_C | 8.093×10⁻²¹ | s |
| m_e c² | 8.187×10⁻¹⁴ | J |
| ω_C | 7.763×10²⁰ | rad/s |
| c | 2.998×10⁸ | m/s |

### Typical Grid Values

For h = λ_C / 5:

- h = 4.85×10⁻¹³ m (0.485 pm)
- dt = 3.24×10⁻²² s
- R = 3.86×10⁻¹³ m (~ 0.8 h)

## Connection to Paper Sections

This implementation directly supports:

### Section 2: Conceptual Model
- Tubular coordinates → Equation (tubular metric)
- Geometric coupling → Curvature terms in Lagrangian
- Bottleneck mechanism → Critical strain threshold

### Section 3: W&vdM Toroidal Electron
- Double-loop structure → `double_loop_envelope`
- Compton frequency → `compton_omega` parameter
- Charge-from-energy → Calibration in diagnostic module

### Section 4: Numerical Methods
- Discrete implementation → Grid setup
- CFL condition → Time step selection
- Stability criteria → Loss function metrics

## Next Steps

### 1. Charge Calibration

Implement `compute_effective_charge` in `electron_stability.py`:

```python
def compute_effective_charge(state):
    # Time-average X⁴ field
    X4_bar = time_average(state.positions[:, 3], T_avg)
    # Apply charge-from-energy relation
    Q = κ_Φ * integrate(X4_bar, volume)
    return Q
```

Then refine amplitude A until Q ≈ -e.

### 2. Energy Refinement

Once charge is correct, check total energy:

```python
E_total = E_kinetic + E_elastic
```

Adjust geometry (rho0, sigma_r) to match E ≈ m_e c².

### 3. Optimization

Run CMA-ES or Nelder-Mead to minimize stability loss:

```python
from cma import CMAEvolutionStrategy

es = CMAEvolutionStrategy(theta0, sigma0)
while not es.stop():
    solutions = es.ask()
    losses = [evaluate_loss(theta) for theta in solutions]
    es.tell(solutions, losses)
```

### 4. Long-Term Stability

Once optimized, run for many Compton periods:

- T_sim = 100 T_C or more
- Check topological integrity (double loop survives)
- Measure spinor holonomy (4π test)
- Compute effective observables (g-factor, magnetic moment)

## Troubleshooting

### "Electron explodes immediately"

**Causes**:
- Amplitude A too large → reduce by factor of 10
- Time step too large → reduce dt
- Grid too coarse → increase resolution

**Fixes**:
```python
A = A * 0.1
dt = dt * 0.5
h = h * 0.5  # Double resolution
```

### "Energy leaks rapidly (> 50%)"

**Causes**:
- Tube too small → increase tube_max_radius
- Initial velocities inconsistent → check Compton mode setup
- Geometry not matched to nonlinearity threshold

**Fixes**:
```python
tube_max_radius = 5.0 * sigma_r  # Increase buffer
# Re-run with adjusted rest_length
```

### "Mode purity low (< 50%)"

**Causes**:
- Initial condition has multiple frequency components
- Simulation time too short for FFT resolution
- Numerical dispersion breaking Compton mode

**Fixes**:
```python
# Ensure single-mode initialization
phase_offset = 0.0  # Clean phase
# Run longer
N_periods = 10
```

## References

### Original W&vdM Paper
Williamson, J. G., & van der Mark, M. B. (1997). "Is the electron a photon with toroidal topology?" *Annales de la Fondation Louis de Broglie*, 22, 133.

### Related Work
- van der Mark, M. B., & 't Hooft, G. (2000). "Light is heavy." arXiv:gr-qc/9906084
- Close, F. (2025). "Spin density waves and elementary particles"

### Project Documentation
- `PROJECT_PRINCIPLES.md`: Core assumptions and constraints
- `docs/paper_draft.tex`: Full theoretical derivation

## Contact

For questions or issues with this implementation:
- Check PROJECT_PRINCIPLES.md for architectural constraints
- Review tubular coordinate derivation in paper draft
- Ensure substrate-only evolution principle is maintained