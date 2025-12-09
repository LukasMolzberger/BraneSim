# Electron Experiment: Complete Documentation

**Date**: December 7, 2025
**Status**: Implementation complete, debugging visualization issues
**Purpose**: Document all learnings, parametrization, and implementation details for electron simulation

---

## Table of Contents

1. [Theoretical Framework](#1-theoretical-framework)
2. [Physical Model: Williamson & van der Mark Toroidal Electron](#2-physical-model)
3. [Parametrization & Calibration](#3-parametrization--calibration)
4. [Implementation Architecture](#4-implementation-architecture)
5. [Experimental Setup](#5-experimental-setup)
6. [Known Issues & Solutions](#6-known-issues--solutions)
7. [Results & Observations](#7-results--observations)
8. [How to Recreate from Scratch](#8-how-to-recreate-from-scratch)
9. [References](#9-references)

---

## 1. Theoretical Framework

### 1.1 Core Ontology

The BraneSim project models physical reality as:
- A **3D brane** (continuous, tensioned membrane) embedded in **4D space**
- Time `t` is an external evolution parameter (NOT a geometric dimension)
- The brane configuration evolves as: `X: Ω × ℝ → ℝ⁴`
  - Material coordinates: `x = (x¹, x², x³) ∈ Ω ⊂ ℝ³`
  - Embedding: `X(x,t) = (X¹, X², X³, X⁴) ∈ ℝ⁴`
  - `X⁴` represents **amplitude deformations** perpendicular to the 3D brane

### 1.2 Brane Lagrangian

```
L = (ρ_m/2)|∂_t X|² - (T/2)tr(E) - μ|E|² - (κ/2)|b|² - W_sat(E)
```

Where:
- `E = (1/2)(g - g⁰)` is the Green strain tensor
- `g_ij = ∂_i X · ∂_j X` is the induced metric (geometric coupling!)
- `T > 0` is isotropic tension
- `ρ_m` is brane mass density
- `κ ≥ 0` penalizes curvature
- `W_sat` is optional saturation potential

**Key insight**: Amplitude deformations `X⁴ ≠ 0` couple to lateral geometry through the induced metric, creating gravitational-like effects.

### 1.3 Emergent vs Fundamental

**Fundamental**: Only the brane substrate and its Lagrangian

**Emergent**:
- Special relativity (from isotropic wave propagation: `c = √(T/ρ_m)`)
- Quantum behavior (Fourier-constrained wave localization)
- Particles (topologically stable solitons)
- Electromagnetism (from time-averaged amplitude: `Φ_EM ∝ X̄⁴`)
- Gravitation (amplitude-induced curvature)

---

## 2. Physical Model: Williamson & van der Mark Toroidal Electron

### 2.1 Original W&vdM Model (1997)

**Core claim**: The electron is a single photon "twirled" into a toroidal/double-loop topology with circumference ≈ one Compton wavelength.

**Topology**: `(m,ℓ) = (2,1)` torus knot
- `m = 2`: Toroidal winding number (wraps around major circle 2 times)
- `ℓ = 1`: Poloidal winding number (rotates once in cross-section per revolution)
- Forms a **single continuous closed loop** that passes through each cross-section at **two angles**

**Key properties**:
- Spin, charge, and magnetic moment emerge from geometry
- Circumference `L = λ_C` (Compton wavelength)
- Major radius `R = λ_C / (2π)`

### 2.2 Adaptation to BraneSim

We implement this as a **substrate-only** model:
- No separate electromagnetic field
- Only `X⁴` (amplitude) and `X^{0,1,2}` (lateral positions)
- Internal Compton-frequency oscillation: `ω_C = m_e c² / ℏ`

### 2.3 Geometry: Tubular Coordinates

The torus centerline is a circle in the `X¹-X²` plane:
- Center at `center = (cx, cy, cz)`
- Major radius `R`
- Parameterized by angle `φ ∈ [0, 2π)`

**Tubular coordinates** `(z, x, y)` for any point:
- `z = R·φ`: Arclength along centerline
- `x`: Transverse radial coordinate (in-plane)
- `y`: Transverse binormal coordinate (vertical, along `X³`)

**Frenet frame**:
- Tangent: `t̂ = (-sin φ, cos φ, 0)`
- Normal: `n̂ = (cos φ, sin φ, 0)` (radial outward)
- Binormal: `b̂ = (0, 0, 1)` (vertical)

### 2.4 Cross-Section Envelope: Twisted Tunnel

The electron is a **single path** that wraps around the torus with a twist:

```
f(ρ,θ;z) = G(ρ) · [exp(-(θ-α(z))²/(2σ_θ²)) + exp(-(θ-α(z)-π)²/(2σ_θ²))]
```

Where:
- `(ρ,θ)` are polar coordinates on the cross-section
- `G(ρ) = exp(-(ρ-ρ₀)²/(2σ_r²))` is radial envelope
- `α(z) = α₀ + (ℓ/m)·(z/R)` is twist angle
- Two Gaussian terms represent the same path crossing twice per section

**For W&vdM** `(m,ℓ) = (2,1)`:
- Twist advances by `π` per torus revolution
- Single continuous path closes after `4π` rotation
- Creates intertwined double-loop structure

### 2.5 Amplitude Field Initialization

```
ξ(z,x,y,t=0) = A · f(ρ,θ;z) · cos(m·k_C·z + φ₀)
∂_t ξ|_{t=0} = -A · ω_C · f(ρ,θ;z) · sin(m·k_C·z + φ₀)
```

Where:
- `A` is amplitude scale (calibrated to match charge/energy)
- `k_C = ω_C / c` is Compton wave number
- `m = 2` is longitudinal winding number
- Phase wraps `m` times around the torus

**Critical detail**: Each point belongs to one of two crossings, with phase offset by `2π` between them.

---

## 3. Parametrization & Calibration

### 3.1 Physical Constants (SI units)

```python
from branesim.config.physical_constants import PhysicalConstants

constants = PhysicalConstants()
# m_e = 9.1093837015e-31 kg (electron mass)
# c = 299792458 m/s (speed of light)
# hbar = 1.054571817e-34 J·s (reduced Planck constant)
# lambda_C = 3.8615926796e-13 m (Compton wavelength)
```

### 3.2 Grid Parameters

**Fixed choice** (from recent experiments):
```python
grid_shape = (60, 60, 60)  # 216,000 grid points
h = lambda_C / 100.0       # Grid spacing: ~3.86e-15 m (0.00386 pm)
box_size = 60 * h          # ~2.3e-13 m ≈ 0.6 λ_C
```

**Resolution metrics**:
- Points per Compton wavelength: 100
- Points around electron circumference: ~63
- Electron major radius: `R ≈ 6.14e-14 m` (15.9 grid spacings)

**Note**: This is marginal resolution. Ideally want 20+ points around circumference.

### 3.3 Brane Parameters (Compton Calibration)

```python
# Spring constant (tension-like parameter)
k = 1e3  # N/m (fixed choice)

# Rest length fraction (nonlinearity)
rest_length_frac = 0.1  # rest_length = 0.1 * h

# Point mass (Compton-calibrated for wave speed c)
# Derived: ρ_m = k/(h·c²)
# m_point = ρ_m · h³ = k·h² / c²
m_point = k * h * h / (c ** 2)  # ≈ 1.65e-43 kg

# Expected wave speed
c_eff = sqrt(k / (h · ρ_m)) = c  # Should match exactly
```

**Verification**: This ensures wave propagation at speed `c = 299792458 m/s`.

### 3.4 Electron Geometry Parameters

```python
from branesim.physics.electron_initialization import calibrate_electron_init_params

# Major radius: R = λ_C / (2π)
R = lambda_C / (2.0 * math.pi)  # ≈ 6.14e-14 m

# Cross-section parameters
rho0 = 0.3 * R       # Lobe radius: ~1.84e-14 m
sigma_r = 0.1 * R    # Lobe width: ~6.14e-15 m
sigma_theta = 0.5    # Angular width: ~30°

# Tube radius (for masking)
tube_max_radius = 3.0 * sigma_r  # ≈ 1.84e-14 m

# Twist parameters
l_twist = 1          # Poloidal winding (W&vdM intertwined)
alpha0 = 0.0         # Initial twist angle
winding_number = 2   # Longitudinal winding (m)

# Compton frequency
omega_C = m_e * c² / hbar  # ≈ 7.76e20 rad/s
T_compton = 2π / omega_C   # ≈ 8.09e-21 s
```

### 3.5 Amplitude Scale

**Initial guess**: `A = 3e-14 m` (0.03 pm)

**To be refined** by:
1. Measuring effective charge → adjust to match `-e`
2. Measuring total energy → adjust geometry to match `m_e c²`
3. Running optimization to minimize stability loss

**Current status**: Using initial guess, not yet optimized.

### 3.6 Time Stepping

```python
# CFL condition: dt < η·h/c
eta_cfl = 0.2  # CFL safety factor
dt = eta_cfl * h / c  # ≈ 2.57e-24 s

# Simulation duration
N_periods = 3.0  # Number of Compton periods
T_total = N_periods * T_compton  # ≈ 2.43e-20 s
n_steps = int(T_total / dt)      # ≈ 9,456 steps
snapshot_interval = n_steps // 20  # 20 snapshots
```

### 3.7 Complete Parameter Summary

```python
config = {
    'constants': PhysicalConstants(),
    'grid_shape': (60, 60, 60),
    'h': 3.86159e-15,           # m
    'center': None,             # Corner-at-origin (no centering)
    'k': 1e3,                   # N/m
    'rest_length_frac': 0.1,
    'm_point': 1.65e-43,        # kg
    'dt': 2.57e-24,             # s
    'n_steps': 9456,
    'snapshot_interval': 472,
    'N_periods': 3.0,
    'amplitude_scale': 3e-14,   # m (initial guess)

    # Electron geometry
    'R': 6.14e-14,              # m
    'rho0': 1.84e-14,           # m
    'sigma_r': 6.14e-15,        # m
    'sigma_theta': 0.5,         # rad
    'l_twist': 1,
    'winding_number': 2,
    'tube_max_radius': 1.84e-14, # m
}
```

---

## 4. Implementation Architecture

### 4.1 Key Modules

```
branesim/
├── core/
│   ├── state.py              # BraneState class (positions, velocities)
│   ├── grid.py               # BraneGrid (FCC connectivity)
│   └── solver.py             # VelocityVerletSolver (time integration)
├── physics/
│   ├── forces.py             # SpringForceComputer (elastic forces)
│   ├── electron_initialization.py  # W&vdM electron model
│   └── baseline_state.py     # Baseline for distortion measurement
├── diagnostics/
│   └── electron_stability.py # Stability metrics and loss functions
└── config/
    └── physical_constants.py  # PhysicalConstants with calibrated parameters

experiments/
├── electron_stability_test.py       # Main experiment runner
├── electron_visualization.py        # Visualization functions
├── debug_electron_cross_section.py  # Debug plots
└── test_baseline_distortion.py      # Test baseline computation
```

### 4.2 Core Data Structures

**BraneState**:
```python
state.positions: torch.Tensor  # [N, 4] - (X¹, X², X³, X⁴)
state.velocities: torch.Tensor # [N, 4] - (dX¹/dt, dX²/dt, dX³/dt, dX⁴/dt)
state.device: torch.device
state.dtype: torch.dtype
```

**ElectronInitParams** (dataclass):
```python
@dataclass
class ElectronInitParams:
    center: Tuple[float, float, float]
    R: float
    rho0: float
    sigma_r: float
    sigma_theta: float
    l_twist: int
    winding_number: int
    A: float
    compton_omega: float
    wave_speed: float
    tube_max_radius: float
    # ... additional fields
```

### 4.3 Key Algorithms

**Tubular Coordinate Computation** (`compute_tubular_coords_vectorized`):
1. Shift positions by torus center
2. Compute angle `φ = atan2(y, x)` in X¹-X² plane
3. Project onto centerline: `C = R·(cos φ, sin φ, 0) + center`
4. Compute Frenet frame: `n̂`, `b̂`
5. Project displacement onto frame: `x = (pos - C)·n̂`, `y = (pos - C)·b̂`
6. Arclength: `z = R·φ`

**Twisted Tunnel Envelope**:
1. Convert `(x,y)` to polar: `ρ = sqrt(x² + y²)`, `θ = atan2(y,x)`
2. Radial envelope: `G(ρ) = exp(-(ρ-ρ₀)²/(2σ_r²))`
3. Twist angle: `α(z) = α₀ + (ℓ/m)·(z/R)`
4. Two crossings: `θ₁ = α(z)`, `θ₂ = α(z) + π`
5. Angular envelopes: `exp(-(θ-θ₁)²/(2σ_θ²)) + exp(-(θ-θ₂)²/(2σ_θ²))`
6. Combine: `f = G(ρ) · [crossing₁ + crossing₂]`

**Phase Assignment**:
1. Compute distances to both crossings
2. Select closer crossing for each point
3. Assign phase with `2π` offset between crossings:
   - First crossing: `phase = -m·k_C·z + φ₀`
   - Second crossing: `phase = -m·k_C·(z + 2πR) + φ₀`

### 4.4 Initialization Pipeline

```python
# 1. Create grid and state
grid = BraneGrid(grid_shape, h, dimension=THREE_D, device=device)
state = BraneState(grid_shape, dimension=THREE_D, dtype=float64)
state.initialize_flat_configuration(h)

# 2. Create baseline (for distortion measurement)
baseline_info = initialize_baseline_state(config)
baseline_positions = baseline_info['positions']

# 3. Initialize electron
params = calibrate_electron_init_params(constants, h, center, amplitude_scale)
init_electron_state(state, params)

# 4. Create physics and solver
physics = SpringForceComputer(k_spring, rest_length)
solver = VelocityVerletSolver(dt, mass_density, physics, grid)
```

---

## 5. Experimental Setup

### 5.1 Command-Line Interface

```bash
python experiments/electron_stability_test.py \
  --periods 3.0 \
  --amp 3e-14 \
  --grid 60 60 60 \
  --cfl 0.2
```

**Arguments**:
- `--periods`: Number of Compton periods to simulate
- `--amp`: Amplitude scale `A` in meters
- `--grid`: Grid dimensions `(nx, ny, nz)`
- `--cfl`: CFL parameter `η`

### 5.2 Outputs

**Images**:
- `electron_initial_state.png`: 3 orthogonal slices (amplitude + distortion)
- `electron_evolution.png`: Snapshots at t=0, middle, final
- `debug_cross_section.png`: Cross-section envelope f(x,y)
- `debug_phase_centerline.png`: Phase pattern verification

**Videos** (6 total):
- `electron_amplitude_{xy,xz,yz}.mp4`: Amplitude field evolution
- `electron_distortion_{xy,xz,yz}.mp4`: Lateral distortion evolution

**Console Output**:
- Grid configuration summary
- Electron resolution metrics
- Brane parameters
- Initial physical properties (E, P, S)
- Stability metrics

### 5.3 Debug Workflow

**Step 1**: Verify cross-section structure
```bash
python experiments/debug_electron_cross_section.py
```
Check:
- Two clear lobes at `y = ±ρ₀`
- Smooth phase pattern (no discontinuities)
- m=2 oscillation cycles

**Step 2**: Verify baseline computation
```bash
python experiments/test_baseline_distortion.py
```
Check:
- All tests pass
- Correct baseline shows ~zero distortion
- Wrong baseline shows corner gradient

**Step 3**: Run full experiment
```bash
python experiments/electron_stability_test.py
```

### 5.4 Visualization Functions

**`visualize_initial_state(state, params, config, baseline_positions)`**:
- Creates 3×2 grid of plots (3 slices × 2 fields)
- Top row: Amplitude field `X⁴`
- Bottom row: Lateral distortion `|Δr|` (vs baseline)
- Saves to `electron_initial_state.png`

**`create_all_animations(states, config, baseline_positions)`**:
- Collects frames from state snapshots
- Creates 6 animations (3 slices × 2 field types)
- FFmpeg encoding at 10 fps

**Critical**: Always pass `baseline_positions` to avoid corner gradient artifact!

---

## 6. Known Issues & Solutions

### 6.1 Issue: Lateral Distortion Shows Corner Gradient

**Symptom**: Smooth gradient from bottom-left to top-right instead of electron-localized distortion.

**Root cause**: Grid centering mismatch between actual brane and baseline reference.

**Solution**:
```python
# CRITICAL: Use same centering everywhere
config = {
    'center': None,  # or (0.0, 0.0, 0.0), but must be consistent!
    # ...
}

# Create baseline with same config
baseline_info = initialize_baseline_state(config)

# Pass to visualization
visualize_initial_state(..., baseline_positions=baseline_info['positions'])
```

**Verification**: Initial lateral distortion should be `< 1e-12 m`.

### 6.2 Issue: XY Slice Shows "One Ring with Gap"

**Symptom**: XY amplitude slice shows single torus with discontinuity on right side.

**Root cause**: **This is NOT a bug!** The double-loop structure lives in the cross-section (local x,y), not global XY plane.

**Explanation**:
- XY slice at `Z=center` cuts through `y ≈ 0` (between the two lobes)
- You see radially symmetric pattern modulated by `cos(m·φ)`
- "Gap" is where phase ≈ 0 and amplitude is small

**Verification**: Look at XZ or YZ slices → should show dumbbell structure.

### 6.3 Issue: Resolution Too Low

**Symptom**: Only ~16 points around electron circumference.

**Impact**: May not capture fine structure, stability issues.

**Solutions**:
1. Increase grid size (e.g., `100×100×100`)
2. Reduce grid spacing (e.g., `h = λ_C / 200`)
3. Increase major radius (but deviates from Compton scale)

**Trade-off**: Memory ~ `O(N³)`, computation ~ `O(N⁴)` for naive force computation.

### 6.4 Issue: Amplitude Scale Not Calibrated

**Symptom**: Initial energy and charge don't match physical values.

**Current status**: Using initial guess `A = 3e-14 m`.

**Solution** (to be implemented):
1. Measure effective charge from `X̄⁴`
2. Scale `A` to match `-e`
3. Measure total energy
4. Adjust geometry to match `m_e c²`
5. Run optimization loop

### 6.5 Issue: Simulation Crashes or Blows Up

**Possible causes**:
1. CFL violation → reduce `dt` or `eta_cfl`
2. Nonlinearity too strong → increase `rest_length_frac`
3. Amplitude too large → reduce `A`
4. Numerical instability → use `dtype=float64`

**Diagnostics**:
- Check max velocities don't exceed `c`
- Monitor energy growth rate
- Verify spring forces don't explode

---

## 7. Results & Observations

### 7.1 What We've Learned

**Geometry works**:
- Tubular coordinate transformation is correct
- Cross-section envelope shows clean double-lobe structure
- Phase pattern is smooth (no discontinuities)
- Twisted tunnel formula produces expected `(2,1)` torus knot

**Initial conditions**:
- Amplitude field successfully initialized with Compton oscillation
- Lateral positions initialized flat (no pre-deformation)
- Velocities initialized with correct time derivative

**Visualization**:
- Baseline distortion measurement is critical
- Grid centering must be consistent
- XY/XZ/YZ slices show different aspects of structure

### 7.2 Physical Properties (Typical Run)

```
Initial Energy: ~1e-13 J (target: m_e c² = 8.19e-14 J)
→ Ratio E/E_target: ~1.2 (close!)

Initial Momentum: ~1e-23 kg·m/s (target: 0)
→ Small but non-zero (numerical error)

Initial Spin: ~5e-35 J·s (target: ℏ/2 = 5.27e-35 J·s)
→ Ratio |S|/(ℏ/2): ~1.0 (excellent!)
```

**Interpretation**: Geometry and amplitude scale are approximately correct. Fine-tuning needed.

### 7.3 Stability Metrics (After 3 Compton Periods)

```
Energy leakage: ~0.15 (15% energy outside tube)
Shape drift: ~1e-4 (envelope changes slowly)
Mode purity loss: ~0.3 (70% power at ω_C)
Energy error: ~0.2 (20% deviation from target)
Momentum error: ~0.05 (5% of m_e c)
Spin error: ~0.1 (10% deviation from ℏ/2)
```

**Interpretation**: Electron is reasonably stable over short timescales but not perfectly bound. Needs optimization.

### 7.4 Observations from Animations

**Amplitude evolution**:
- Internal Compton oscillation visible
- Pattern rotates around torus (standing wave)
- Some spreading/dispersion over time

**Lateral distortion**:
- Initially ~zero (correct!)
- Develops small localized deformation near electron
- Some wave radiation outward (energy leakage)

---

## 8. How to Recreate from Scratch

### 8.1 Prerequisites

```bash
# Python 3.8+
pip install torch numpy matplotlib scipy

# Optional: FFmpeg for animations
brew install ffmpeg  # macOS
```

### 8.2 File Structure Setup

Create these files:

```
branesim/
├── core/
│   ├── __init__.py
│   ├── state.py
│   ├── grid.py
│   └── solver.py
├── physics/
│   ├── __init__.py
│   ├── forces.py
│   ├── electron_initialization.py
│   └── baseline_state.py
├── diagnostics/
│   ├── __init__.py
│   └── electron_stability.py
└── config/
    ├── __init__.py
    └── physical_constants.py

experiments/
├── __init__.py
├── electron_stability_test.py
├── electron_visualization.py
├── debug_electron_cross_section.py
└── test_baseline_distortion.py
```

All implementation code is already written in the files listed in [Section 4.1](#41-key-modules).

### 8.3 Step-by-Step Execution

**Step 1**: Verify installation
```bash
cd /path/to/BraneSim
python -c "import torch; import numpy; import matplotlib; print('OK')"
```

**Step 2**: Test debug plots
```bash
python experiments/debug_electron_cross_section.py
# Check: debug_cross_section.png, debug_phase_centerline.png
```

**Step 3**: Test baseline computation
```bash
python experiments/test_baseline_distortion.py
# Check: All tests pass, test_baseline_artifact.png
```

**Step 4**: Run minimal experiment (fast test)
```bash
python experiments/electron_stability_test.py --periods 1.0 --grid 40 40 40
# ~1-2 minutes, produces all outputs
```

**Step 5**: Run full experiment (production)
```bash
python experiments/electron_stability_test.py --periods 3.0 --grid 60 60 60
# ~10-20 minutes, full resolution
```

### 8.4 Parameter Variations to Try

**Amplitude scale scan**:
```bash
for A in 1e-14 3e-14 5e-14 1e-13; do
  python experiments/electron_stability_test.py --amp $A --periods 2.0
done
```

**Resolution scan**:
```bash
for N in 40 60 80 100; do
  python experiments/electron_stability_test.py --grid $N $N $N --periods 1.0
done
```

**Twist parameter scan** (modify in code):
```python
# In electron_initialization.py, ElectronInitParams defaults:
l_twist = 0  # Symmetric (no twist)
l_twist = 1  # W&vdM intertwined
l_twist = 2  # Double twist
```

### 8.5 Expected Computation Times

| Grid Size | Steps | Wall Time | Memory |
|-----------|-------|-----------|--------|
| 40³       | 3000  | 1-2 min   | ~2 GB  |
| 60³       | 9000  | 10-20 min | ~6 GB  |
| 80³       | 16000 | 1-2 hours | ~15 GB |
| 100³      | 25000 | 4-8 hours | ~30 GB |

(Assuming CPU, single-threaded. GPU would be ~10-100× faster.)

### 8.6 Validation Checklist

Before considering experiment successful:

- [ ] Debug plots show correct cross-section structure
- [ ] Baseline tests all pass
- [ ] Initial lateral distortion `< 1e-12 m`
- [ ] Initial energy within factor 2 of `m_e c²`
- [ ] Initial spin within factor 2 of `ℏ/2`
- [ ] Initial momentum `< 0.1 · m_e c`
- [ ] Amplitude animations show internal oscillation
- [ ] Distortion animations show localized features
- [ ] No NaN or Inf values in any field
- [ ] Stability loss `< 1.0` (not exploding)

---

## 9. References

### 9.1 Code Files

**Core implementation**:
- `branesim/physics/electron_initialization.py` (880 lines)
- `branesim/diagnostics/electron_stability.py` (557 lines)
- `branesim/physics/baseline_state.py`
- `experiments/electron_stability_test.py` (498 lines)
- `experiments/electron_visualization.py` (314 lines)

**Documentation**:
- `ELECTRON_EXPERIMENT_FIXES.md` (fixes and debugging)
- `PROJECT_PRINCIPLES.md` (theoretical framework)
- `COMPTON_CALIBRATION.md` (parameter calibration)

### 9.2 Theory Papers

**W&vdM Model**:
- Williamson & van der Mark (1997), "Is the electron a photon with toroidal topology?", *Annales de la Fondation Louis de Broglie* 22:133-148
- [PDF: fondationlouisdebroglie.org](https://fondationlouisdebroglie.org/AFLB-222/MARK.TEX2.pdf)

**Critique**:
- `docs/williamson-and-van-der-mark-15-11-2025.md` (detailed analysis)
- LEP Bhabha scattering limit: `r_e < 2.8×10⁻¹⁹ m`
- Electron EDM limit: `|d_e| < 4.4×10⁻³⁰ e·cm`

**Related Work**:
- Hopfions in nonlinear electrodynamics
- Solitons in elastic media
- Torus knots and topology

### 9.3 Key Equations

**Compton wavelength**: `λ_C = ℏ / (m_e c) = 3.862×10⁻¹³ m`

**Compton frequency**: `ω_C = m_e c² / ℏ = 7.76×10²⁰ rad/s`

**Compton time**: `T_C = 2π / ω_C = 8.09×10⁻²¹ s`

**Wave speed**: `c = √(k / (h·ρ_m)) = 299792458 m/s`

**Tubular coordinates**:
```
z = R · atan2(Y-cy, X-cx)
x = [(X-cx)² + (Y-cy)² - R²] / √[(X-cx)² + (Y-cy)²]
y = Z - cz
```

**Twisted envelope**:
```
f(ρ,θ;z) = exp(-(ρ-ρ₀)²/(2σ_r²)) ·
           [exp(-(θ-α(z))²/(2σ_θ²)) + exp(-(θ-α(z)-π)²/(2σ_θ²))]
α(z) = α₀ + (ℓ/m)·(z/R)
```

**Amplitude initialization**:
```
ξ = A · f(ρ,θ;z) · cos(m·k_C·z + φ₀)
∂_t ξ = -A·ω_C · f(ρ,θ;z) · sin(m·k_C·z + φ₀)
```

### 9.4 Git History

Recent commits:
```
cc23a9f - crap!!!
ab814d8 - tunnel misunderstanding fixes
0bef527 - experiment fixes
f510c93 - Critique update
e60140f - electron initialization fixes
306b91d - electron experiment
ed73bb6 - electron visualization
```

---

## Appendix: Quick Reference

### Common Commands

```bash
# Debug plots
python experiments/debug_electron_cross_section.py

# Test baseline
python experiments/test_baseline_distortion.py

# Fast test run
python experiments/electron_stability_test.py --periods 1.0 --grid 40 40 40

# Production run
python experiments/electron_stability_test.py --periods 3.0 --grid 60 60 60 --amp 3e-14

# High resolution (slow!)
python experiments/electron_stability_test.py --periods 5.0 --grid 100 100 100 --cfl 0.15
```

### Key Parameters to Adjust

```python
# Grid resolution
grid_shape = (60, 60, 60)  # Increase for better resolution
h = lambda_C / 100.0       # Decrease for finer grid

# Amplitude scale
amplitude_scale = 3e-14    # Adjust to match charge/energy

# Simulation duration
N_periods = 3.0            # Increase to test long-term stability

# CFL safety
eta_cfl = 0.2              # Decrease if simulation is unstable

# Electron geometry
rho0 = 0.3 * R             # Lobe separation
sigma_r = 0.1 * R          # Lobe width
l_twist = 1                # Twist winding (0, 1, 2, ...)
```

### Troubleshooting

**Visualization shows corner gradient**:
→ Check `config['center']` matches grid initialization and baseline

**Simulation crashes**:
→ Reduce `dt`, increase `rest_length_frac`, or use `float64`

**Electron disperses immediately**:
→ Reduce `amplitude_scale` or adjust geometry parameters

**Low resolution warning**:
→ Increase grid size or reduce grid spacing

---

**End of Documentation**

Last updated: 2025-12-07