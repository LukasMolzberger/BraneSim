# PhotonWavePacket Implementation

## Overview

This document describes the implementation of realistic photon wave packets in the brane simulation framework. The implementation models photons as localized Gaussian wave packets in the 4th dimension with physically motivated parameters.

## Core Components

### 1. PhotonWavePacket Class

**Location:** `src/main/java/sim/objects/PhotonWavePacket.java`

The `PhotonWavePacket` class implements the `PhysicalObject` interface and represents a photon as a Gaussian wave packet in the 4th brane dimension.

#### Key Features:

- **Wavelength**: Configurable, defaults to Compton wavelength (λ_C ≈ 2.43×10⁻¹² m)
- **Propagation direction**: Arbitrary unit vector (n̂)
- **Gaussian envelope**: Separate longitudinal (σ∥) and transverse (σ⊥) widths
- **Phase velocity**: Consistent with speed of light (c)
- **Polarization**: Stored parameter (currently conceptual for scalar field)

#### Mathematical Formulation:

The photon is initialized as:

```
X^4(x,t=0) = A₀ exp(-s²/(2σ∥²) - r⊥²/(2σ⊥²)) cos(k·s + φ₀)
```

where:
- `s = (x - x₀)·n̂` is the longitudinal coordinate
- `r⊥² = |x - x₀|² - s²` is the transverse distance squared
- `k = 2π/λ` is the wave number
- `A₀` is the amplitude

The 4th-dimension velocity is set to:

```
∂X^4/∂t|_{t=0} = A₀ ω exp(-...) sin(k·s + φ₀)
```

where `ω = c·k` ensures the wave packet propagates at speed of light.

#### Constructors:

1. **Full constructor**: Allows complete control over all parameters
   ```java
   PhotonWavePacket(amplitude, wavelength, centerX, centerY, centerZ,
                    dirX, dirY, dirZ, sigmaParallel, sigmaPerp,
                    polarizationAngle, initialPhase)
   ```

2. **Factory method**: Creates photon with Compton wavelength
   ```java
   PhotonWavePacket.createElectronComptonPhoton(amplitude, polarizationAngle)
   ```

### 2. E8_PhotonPropagation Experiment

**Location:** `src/main/java/sim/experiments/E8_PhotonPropagation.java`

Demonstration experiment that launches a photon through a non-cubic grid optimized for propagation studies.

#### Grid Configuration:

- **Dimensions**: 128×32×32 (long corridor along x-axis)
- **Grid spacing**: λ_C / 8 (8 cells per wavelength)
- **Time step**: 0.5 × (dx/c) for CFL stability
- **No nonlinearity**: Linear regime (amplitude = 0.01)
- **No damping**: Energy-conserving propagation

#### Photon Parameters:

- **Wavelength**: Compton wavelength (λ_C)
- **Longitudinal width**: 3 λ_C
- **Transverse width**: 0.5 λ_C
- **Initial position**: 25% along x-axis, centered in y,z
- **Propagation**: Along +x direction
- **Amplitude**: 0.01 (well below nonlinearity scale)

## Usage

### Interactive Mode (GUI)

1. Launch the GUI:
   ```bash
   java -cp target/classes sim.gui.Gui
   ```

2. Press '8' to activate the photon propagation experiment

3. Controls:
   - `r`: Reset to default
   - `s`: Pause/unpause
   - `g`: Toggle gravity visualization
   - `p`: Toggle plane highlighting
   - Mouse drag: Rotate view
   - Shift+drag: Additional rotation
   - Mouse wheel: Zoom

### Programmatic Usage

```java
// Create a photon with Compton wavelength
PhotonWavePacket photon = PhotonWavePacket.createElectronComptonPhoton(
    0.01,  // amplitude
    0.0    // polarization angle
);

// Or create with custom parameters
PhotonWavePacket customPhoton = new PhotonWavePacket(
    0.05,                              // amplitude
    1e-12,                             // wavelength (m)
    centerX, centerY, centerZ,         // position
    1.0, 0.0, 0.0,                     // direction (+x)
    3.0 * wavelength,                  // sigma parallel
    0.5 * wavelength,                  // sigma perpendicular
    0.0,                               // polarization angle
    0.0                                // initial phase
);

// Use in experiment
BraneConfig config = new BraneConfig(128, 32, 32);
Brane brane = new Brane(config);
photon.initializeBrane(brane);
```

## Physical Interpretation

### Mapping to Paper Theory

From the paper's 4D brane model:
- **Brane position**: X(x,t) = (X¹, X², X³, X⁴)
- **4th dimension**: X⁴ represents amplitude into extra dimension
- **Electric field**: E ~ ∇X⁴ (and ∂_t X⁴)
- **Photon**: Localized wave in X⁴ with λ = λ_C

### Compton Wavelength Choice

The Compton wavelength is used to match the "internal photon" geodesic length in the toroidal soliton (electron) model from Williamson & van der Mark:

```
λ_C = ℏ / (m_e c) ≈ 2.43×10⁻¹² m
```

This creates consistency between:
- Free photons (this implementation)
- Internal photon structure in electrons (E2_ToroidalSoliton)

### Speed of Light

The photon's propagation speed emerges from:
1. Initialization with 4th-component velocity: v₄ = ω sin(ks)
2. Angular frequency: ω = c·k
3. Brane dynamics that maintain this phase relationship

## Next Steps / Refinements

### Short Term:
1. Add diagnostics to track:
   - Center-of-energy position
   - Dispersion (width evolution)
   - Group velocity measurement
   - Energy conservation

2. Calibrate grid spacing/time step to minimize dispersion

3. Create test cases for multiple photons (interference)

### Medium Term:
1. Implement full vector polarization once E-field mapping is refined
2. Add photon-electron interaction experiments
3. Create photon scattering experiments (photon-photon, photon-soliton)

### Long Term:
1. Implement QED-inspired photon statistics
2. Model photon emission/absorption by solitons
3. Study emergent dispersion relations
4. Develop calibration to map simulation units to SI units precisely

## Implementation Notes

### Design Decisions:

1. **Additive initialization**: The `initializeBrane()` method adds to existing field values, allowing composition of multiple objects

2. **Scalar field**: Currently X⁴ is scalar; polarization is stored but not fully utilized (awaiting refined E-field mapping)

3. **4D velocity**: Both position and velocity in 4th dimension are initialized to ensure proper wave propagation

4. **Normalization**: Direction vector is automatically normalized in constructor

5. **Parameter validation**: Constructor checks for positive wavelength and widths

### Numerical Stability:

- **CFL condition**: dt < dx/c ensures stability for wave propagation
- **Grid resolution**: 8 cells per wavelength provides good balance between accuracy and performance
- **Amplitude**: Kept at 0.01 to stay in linear regime (well below nonlinearity scale)

## Files Modified/Created

### Created:
- `src/main/java/sim/objects/PhotonWavePacket.java`
- `src/main/java/sim/experiments/E8_PhotonPropagation.java`
- `docs/PhotonWavePacket_Implementation.md` (this file)

### Modified:
- `src/main/java/sim/gui/Gui.java` (added keyboard shortcut '8' for E8)

## References

- Paper Section 4.2.7: Experiments E1-E7
- PhysicsConstants.java: Physical constants (c, λ_C, etc.)
- Williamson & van der Mark: Toroidal electron model
- Existing physical objects: GaussianPulse, ToroidalSoliton, ProbeWavePacket