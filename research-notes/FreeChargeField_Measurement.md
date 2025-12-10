# Free Charge Field Measurement

## Overview

This measurement device computes the **free charge density field** ρ_free(x) from the brane state, implementing the free vs bound charge framework described in the paper (Section: *Charge as Amplitude Deformation* and *Discussion*).

## Conceptual Background

### Free vs Bound Charge in the Brane Picture

The brane's 4th-dimension displacement X^4(x,t) encodes the **total charge content** of the medium:

- **Bound charge**: Fast, oscillatory patterns (internal photon modes, wave packets)
  - Creates strong local X^4 oscillations
  - Dipole-like patterns that cancel on scales larger than Compton wavelength
  - Examples: internal electron photon, external photon wave packets

- **Free charge**: Slowly varying monopole component
  - Survives after bound oscillations are averaged out
  - Sources the long-range Coulomb field via Poisson equation
  - Examples: electron (-e), positron (+e), neutral for photons (≈0)

### Mathematical Relations

From the paper (tex/reconstructing-physics.tex):

```
Φ(x) = κ X̄^4(x)                    [Eq. phi-from-xi4]
∇²Φ(x) = -ρ_free(x) / ε_eff        [Eq. poisson-from-brane]
```

Combining these:
```
ρ_free(x) = -ε_eff κ ∇²X̄^4(x)
```

Where:
- X̄^4 = time/space averaged X^4 (smoothed field)
- ε_eff = effective permittivity of the brane medium
- κ = stiffness parameter coupling amplitude to tension

## Implementation in `sim.Measurements`

### Core Methods

#### 1. `computeSmoothedX4Field()`
Computes X̄^4 by spatial averaging over 3×3×3 neighborhoods:
```java
double[][][] smoothedX4 = measurements.computeSmoothedX4Field();
```
- Filters out rapid bound-charge oscillations
- Preserves slowly varying monopole component
- Returns 3D array [x][y][z]

#### 2. `computeLaplacianX4(smoothedX4)`
Computes ∇²X̄^4 using 7-point finite difference stencil:
```java
double[][][] laplacian = measurements.computeLaplacianX4(smoothedX4);
```
- Measures local curvature of amplitude deformation
- Uses second-order central differences
- Handles boundaries with Neumann conditions (∂X̄^4/∂n = 0)

#### 3. `computeFreeChargeDensity()`
Computes ρ_free = -ε_eff κ ∇²X̄^4:
```java
double[][][] rhoFree = measurements.computeFreeChargeDensity();
```
- Currently uses normalized units: ε_eff κ = 1.0
- Returns 3D charge density field
- Positive values = positron-like, negative = electron-like

#### 4. `computeTotalFreeCharge()`
Integrates ρ_free over entire brane volume:
```java
double totalCharge = measurements.computeTotalFreeCharge();
```
- ∫ρ_free d³x over the grid
- Should be ≈0 for neutral photons
- Should be ≈±e for electrons/positrons

#### 5. `computeFreeChargeStatistics()`
Comprehensive charge distribution statistics:
```java
Map<String, Double> stats = measurements.computeFreeChargeStatistics();
// Keys: "total_charge", "max_positive_density",
//       "max_negative_density", "rms_density"
```

#### 6. `printFreeChargeFieldDiagnostics()`
Pretty-printed diagnostic summary:
```java
measurements.printFreeChargeFieldDiagnostics();
```
Output example:
```
Free Charge Field Diagnostics at t=1000:
  Total free charge Q_free: -1.234567e-08
  Max positive ρ_free:      5.678901e-06
  Max negative ρ_free:      4.321098e-06
  RMS charge density:        2.345678e-07
  → Electrically neutral (photon-like)
```

## Usage Examples

### Example 1: Verify Photon Neutrality

Test that your `PhotonWavePacket` has zero net monopole charge:

```java
// In your experiment or runner

import core.Measurements;
import core.objects.PhotonWavePacket;

// Create photon and initialize brane
PhotonWavePacket photon = new PhotonWavePacket(...);
        brane.

        init(experiment);

// Let it propagate for a while
for(
        int step = 0;
        step< 1000;step++){
        brane.

        simulationUpdateStep();
}

        // Measure free charge
        Measurements m = new Measurements(brane);
        double qTotal = m.computeTotalFreeCharge();

System.out.

        printf("Photon total charge: %.3e (should be ≈0)\\n",qTotal);
m.

        printFreeChargeFieldDiagnostics();
```

Expected output:
```
Photon total charge: 2.345e-12 (should be ≈0)
Free Charge Field Diagnostics at t=1000:
  Total free charge Q_free: 2.345000e-12
  ...
  → Electrically neutral (photon-like)
```

### Example 2: Track Charge During Experiment

Add charge monitoring to your experiment runner:

```java
// In E8_PhotonPropagation or similar experiment
@Override
public void simulationUpdateStep(Matrix world) {
    super.simulationUpdateStep(world);

    // Every 100 steps, check charge conservation
    if (getBrane().getT() % 100 == 0) {
        Measurements m = new Measurements(getBrane());
        double q = m.computeTotalFreeCharge();
        System.out.printf("t=%d: Q_free = %.6e\\n", getBrane().getT(), q);
    }
}
```

### Example 3: Compare Internal vs External Photon

```java
// Test 1: External photon (should be neutral)
PhotonWavePacket externalPhoton = PhotonWavePacket.createElectronComptonPhoton(
    0.1,  // amplitude
    0.0   // polarization
);
brane1.init(experiment1);
Measurements m1 = new Measurements(brane1);
double q1 = m1.computeTotalFreeCharge();
System.out.printf("External photon Q_free: %.3e\\n", q1);

// Test 2: Toroidal soliton (should have net charge)
ToroidalSoliton electron = new ToroidalSoliton(...);
brane2.init(experiment2);
Measurements m2 = new Measurements(brane2);
double q2 = m2.computeTotalFreeCharge();
System.out.printf("Electron Q_free: %.3e\\n", q2);
```

Expected:
```
External photon Q_free: 1.234e-12  ← neutral
Electron Q_free: -1.602e-19        ← elementary charge
```

## Interpretation Guide

### What the Numbers Mean

| Total Charge | Physical Interpretation |
|--------------|------------------------|
| |Q| < 10⁻¹⁰  | Electrically neutral (photon-like, radiation mode) |
| Q ≈ -e       | Net negative free charge (electron-like soliton) |
| Q ≈ +e       | Net positive free charge (positron-like soliton) |

### Charge Density Maps

The 3D field ρ_free[x][y][z] shows:
- **Positive regions**: Peaks in X̄^4 (convex curvature)
- **Negative regions**: Troughs in X̄^4 (concave curvature)
- **Zero regions**: Flat or saddle points

For a **neutral photon**:
- ρ_free oscillates between positive and negative
- Spatial integral cancels: ∫ρ_free ≈ 0

For a **charged soliton**:
- ρ_free has net monopole component
- Spatial integral non-zero: ∫ρ_free ≈ ±e

## Connection to Paper

This implementation directly corresponds to the framework in:

1. **tex/reconstructing-physics.tex** (Section: Charge as Amplitude Deformation):
   - Eq. [phi-from-xi4]: Φ = κ X̄^4
   - Eq. [poisson-from-brane]: ∇²Φ = -ρ_free/ε_eff
   - Curvature relation: ρ_free = -ε_eff κ ∇²X̄^4

2. **tex/discussion.tex** (Paragraph: Free versus bound charge):
   - Internal photon → bound charge (no net monopole)
   - External photon → neutral (zero monopole)
   - Electron/positron → free charge (non-zero monopole)

## Future Extensions

### 1. Time-Averaging
Currently uses spatial averaging only. Could add:
```java
public double[][][] computeTimeAveragedX4Field(int windowSize) {
    // Average X^4 over last N timesteps
    // Better filtering of Compton-frequency oscillations
}
```

### 2. Calibrated Units
Use physical constants instead of normalized ε_eff κ = 1:
```java
public double[][][] computeFreeChargeDensity(double eps_eff, double kappa) {
    // ρ_free = -eps_eff * kappa * laplacian
}
```

### 3. Field Visualization
Export charge density for 3D visualization:
```java
public void exportFreeChargeField(String filename) {
    // Write ρ_free[x][y][z] to VTK/HDF5 for ParaView
}
```

### 4. Multipole Analysis
Decompose charge distribution into monopole, dipole, quadrupole:
```java
public Map<String, Double> computeChargeMultipoles() {
    // Q_monopole = ∫ρ d³x
    // p_dipole = ∫r ρ d³x
    // Q_quadrupole = ∫(3r_i r_j - r²δ_ij) ρ d³x
}
```

## References

- Paper section: `tex/reconstructing-physics.tex` (lines 508-641)
- Paper discussion: `tex/discussion.tex` (lines 31-56)
- Implementation: `src/main/java/sim/Measurements.java` (lines 439-668)
- Example usage: `src/main/java/sim/experiments/E8_PhotonPropagation.java`