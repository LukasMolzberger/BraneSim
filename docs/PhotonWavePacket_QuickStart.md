# PhotonWavePacket - Quick Start Guide

## What Was Implemented

A realistic photon wave packet implementation based on 4D brane theory, modeling photons as Gaussian wave packets in the 4th dimension.

## Files Created

1. **`src/main/java/sim/objects/PhotonWavePacket.java`**
   - Core photon implementation
   - Implements `PhysicalObject` interface
   - Models photon as Gaussian wave packet with wavelength, direction, polarization
   - Defaults to Compton wavelength (λ_C ≈ 2.43×10⁻¹² m)

2. **`src/main/java/sim/experiments/E8_PhotonPropagation.java`**
   - Example experiment demonstrating photon propagation
   - Non-cubic grid (128×32×32) optimized for propagation along x-axis
   - Grid spacing: 8 cells per wavelength
   - Time step: CFL-stable (0.5 × dx/c)

3. **`docs/PhotonWavePacket_Implementation.md`**
   - Detailed implementation documentation
   - Physical interpretation and theory mapping
   - Usage examples and API reference

## Files Modified

1. **`src/main/java/sim/gui/Gui.java`**
   - Added keyboard shortcut '8' to launch E8_PhotonPropagation

2. **`CLAUDE.md`**
   - Updated with photon experiment and physical object documentation

## Quick Usage

### In GUI (Interactive)

```bash
# Compile
mvn compile

# Launch GUI
java -cp target/classes sim.gui.Gui

# Press '8' to activate photon propagation experiment
```

### Programmatically

```java
// Create photon with Compton wavelength
PhotonWavePacket photon = PhotonWavePacket.createElectronComptonPhoton(
    0.01,  // amplitude (linear regime)
    0.0    // polarization angle (rad)
);

// Or create with custom parameters
PhotonWavePacket customPhoton = new PhotonWavePacket(
    0.05,                      // amplitude
    1e-12,                     // wavelength (m)
    centerX, centerY, centerZ, // position
    1.0, 0.0, 0.0,             // direction (+x)
    3.0 * wavelength,          // sigma parallel
    0.5 * wavelength,          // sigma perpendicular
    0.0,                       // polarization angle
    0.0                        // initial phase
);

// Initialize brane
BraneConfig config = new BraneConfig(128, 32, 32);
Brane brane = new Brane(config);
photon.initializeBrane(brane);
```

## Key Features

### Physical Parameters
- **Wavelength**: λ (defaults to Compton wavelength)
- **Wave number**: k = 2π/λ
- **Angular frequency**: ω = c·k
- **Propagation direction**: Arbitrary unit vector n̂
- **Gaussian envelope**:
  - Longitudinal width: σ∥ (default: 3λ_C)
  - Transverse width: σ⊥ (default: 0.5λ_C)
- **Polarization**: Angle in plane ⊥ to n̂ (conceptual for scalar field)

### Initialization
Sets both position and velocity in 4th dimension:
- Position: X⁴ = A exp(-...) cos(k·s + φ)
- Velocity: ∂X⁴/∂t = A ω exp(-...) sin(k·s + φ)

This ensures the wave packet propagates at speed c.

## Theory Mapping

From paper's 4D brane model:
- **Brane**: X(x,t) = (X¹, X², X³, X⁴)
- **4th dimension**: X⁴ = amplitude into extra dimension
- **Electric field**: E ~ ∇X⁴ (and ∂_t X⁴)
- **Photon**: Localized wave in X⁴

The Compton wavelength choice matches the "internal photon" in the toroidal soliton (electron) model, ensuring consistency between free photons and bound photon structures.

## Grid Configuration for E8

The E8_PhotonPropagation experiment uses:
- **Grid**: 128×32×32 (long corridor along x)
- **Spacing**: λ_C / 8 (8 cells per wavelength)
- **Time step**: 0.5 × (dx/c) for stability
- **Spring constant**: 1.0 (default)
- **No nonlinearity**: Linear regime
- **No damping**: Energy-conserving

## Next Steps

1. Run the experiment and observe propagation
2. Add diagnostics (see E8_PhotonPropagation.onSimulationStep())
3. Measure group velocity and compare to c
4. Test with multiple photons for interference effects
5. Refine polarization implementation once E-field mapping is finalized

## References

- Full documentation: `docs/PhotonWavePacket_Implementation.md`
- Paper: Section 4.2.7 (Experiments E1-E7)
- Constants: `src/main/java/sim/PhysicsConstants.java`
- Interface: `src/main/java/sim/objects/PhysicalObject.java`