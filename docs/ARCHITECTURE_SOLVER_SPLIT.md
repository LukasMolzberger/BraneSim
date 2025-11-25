# Architecture: Physics/Solver Split

## Overview

The brane simulation has been refactored to cleanly separate the **physics model** (what forces act on the system) from the **numerical solver** (how to integrate the equations of motion).

This separation follows a fundamental principle in scientific computing: the **laws of physics** should be independent of the **numerical method** used to simulate them.

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                     Brane (Facade)                       │
│  - Grid management                                       │
│  - Initialization                                        │
│  - Public API (backward compatible)                      │
└───────────────┬─────────────────┬───────────────────────┘
                │                 │
      ┌─────────▼────────┐  ┌────▼──────────────┐
      │  BranePhysics    │  │   BraneSolver      │
      │   (Interface)    │  │   (Interface)      │
      └─────────┬────────┘  └────┬──────────────┘
                │                 │
      ┌─────────▼────────┐  ┌────▼──────────────┐
      │ SpringMeshPhysics│  │ VelocityVerlet     │
      │ (Implementation) │  │Solver              │
      │                  │  │ (Implementation)   │
      │ - computeAccel() │  │  - step()          │
      │ - calcPotential()│  │  - integrate()     │
      └──────────────────┘  └────────────────────┘
```

## Components

### 1. BranePhysics Interface (`sim.physics.BranePhysics`)

Defines **WHAT** forces act on the system.

**Responsibilities:**
- Compute accelerations from current positions (F = ma)
- Calculate potential energy
- Define the Lagrangian/Hamiltonian of the system

**Key methods:**
```java
void computeAccelerations(BranePoint[][][] points, BraneConfig config);
double calculateSpringPotential(double extension, BraneConfig config);
boolean isOutside(int nx, int ny, int nz, BraneConfig config);
boolean checkSymmetry(int nxRel, int nyRel, int nzRel);
```

**Properties:**
- Pure function: state → accelerations
- No knowledge of time integration
- No knowledge of timestep size (dt)
- Stateless (except for caching optimizations)

### 2. SpringMeshPhysics Implementation (`sim.physics.SpringMeshPhysics`)

Implements the 4D spring-mass mesh physics.

**Physics:**
- Each point has position R_p ∈ ℝ⁴
- Neighbors connected by springs with stiffness k
- Force: F = -∂U/∂R where U is elastic potential
- Linear regime: F = k * displacement (Hooke's law)
- Nonlinear regime: F = k * displacement / (1 + (ε/ε_cr)²)

**Key features:**
- Treats all 4 dimensions uniformly
- No special-casing for the 4th dimension
- All EM/GR interpretations are emergent (handled by measurements)

### 3. BraneSolver Interface (`sim.solvers.BraneSolver`)

Defines **HOW** to evolve the system forward in time.

**Responsibilities:**
- Integrate dx/dt = v, dv/dt = a(x,v,t)
- Maintain numerical stability
- Preserve conservation laws (energy, momentum)

**Key methods:**
```java
void step(BranePoint[][][] points, BranePhysics physics, BraneConfig config, Experiment experiment);
String getName();  // e.g., "Velocity Verlet"
int getOrder();    // Accuracy order (2 for Verlet, 4 for RK4)
```

**Properties:**
- Independent of physics model
- Only knows about time integration
- Can be swapped without changing physics

### 4. VelocityVerletSolver Implementation (`sim.solvers.VelocityVerletSolver`)

Implements the Velocity Verlet time integration scheme.

**Algorithm:**
1. x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
2. Compute a(t+dt) from x(t+dt) using physics model
3. v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt

**Properties:**
- 2nd order symplectic integrator
- Good energy conservation
- Time-reversible
- Explicit (no matrix inversion)

### 5. Brane Facade (`sim.Brane`)

Combines physics + solver + infrastructure.

**Responsibilities:**
- Grid point creation and management
- Initialization (delegates to experiments)
- Public API (backward compatible with existing code)
- Facade pattern: delegates to physics and solver

**Constructor options:**
```java
// Default: SpringMeshPhysics + VelocityVerletSolver
Brane brane = new Brane(N, M, K);

// Custom physics and/or solver
Brane brane = new Brane(config, customPhysics, customSolver);
```

## Benefits of This Architecture

### 1. **Separation of Concerns**

- Physics experts can modify forces without touching time integration
- Numerical analysts can implement new solvers without understanding physics
- Clear boundaries reduce coupling and complexity

### 2. **Extensibility**

**New physics models:**
- Alternative force laws (e.g., relativi stic corrections)
- Different nonlinearity models
- External fields

**New solvers:**
- Runge-Kutta 4 (4th order, more stable)
- Leapfrog (explicit symplectic)
- Backward Euler (implicit, unconditionally stable)
- Adaptive timestep solvers

### 3. **Testability**

- Test physics models independently (force correctness)
- Test solvers independently (convergence rates, stability)
- Test integration (end-to-end behavior)

### 4. **Comparison and Validation**

- Run same physics with different solvers → verify convergence
- Run different physics with same solver → isolate physics effects
- A/B testing for numerical accuracy vs. performance

### 5. **Backward Compatibility**

- Existing code continues to work unchanged
- Default behavior identical to before refactoring
- All tests pass without modification

## Usage Examples

### Default Usage (Unchanged)

```java
// Works exactly as before - no code changes needed
Brane brane = new Brane(35, 35, 35);
Experiment exp = new E1_ThresholdLocalization(brane);
```

### Custom Physics Model

```java
// Implement alternative physics
class RelativisticPhysics implements BranePhysics {
    @Override
    public void computeAccelerations(BranePoint[][][] points, BraneConfig config) {
        // Implement relativistic corrections to forces
    }
    // ...
}

Brane brane = new Brane(config, new RelativisticPhysics(), new VelocityVerletSolver());
```

### Custom Solver

```java
// Implement higher-order solver
class RK4Solver implements BraneSolver {
    @Override
    public void step(BranePoint[][][] points, BranePhysics physics,
                     BraneConfig config, Experiment experiment) {
        // Implement 4th order Runge-Kutta
    }
    // ...
}

Brane brane = new Brane(config, new SpringMeshPhysics(), new RK4Solver());
```

### Solver Comparison Study

```java
// Compare Verlet vs RK4 for accuracy and stability
Brane braneVerlet = new Brane(config, physics, new VelocityVerletSolver());
Brane braneRK4 = new Brane(config, physics, new RK4Solver());

// Run both and compare energy conservation, convergence, etc.
```

## Relation to PHILOSOPHY.md Layer Separation

This split is orthogonal to the layer separation in PHILOSOPHY.md:

```
Layer 5: Visualization (GUI)
  ↓ (read-only)
Layer 4: Measurements (MeasurementDevices)
  ↓ (read-only)
Layer 3: Experiments (orchestration)
  ↓ (initialization only)
Layer 2: PhysicalObjects (initial conditions)
  ↓
Layer 1: Brane Physics (THIS REFACTORING)
  ├─ BranePhysics (force computation)
  └─ BraneSolver (time integration)
  ↓
Layer 0: Data Structures (BranePoint, Vector, etc.)
```

**Key point:** Layer 1 (Brane Physics) is now split internally into:
- **Physics model** - defines the Lagrangian
- **Numerical solver** - integrates the equations of motion

This split maintains the existing layer boundaries while improving internal organization.

## Migration Notes

### For Developers

**No changes required for:**
- Existing experiments (all work unchanged)
- GUI code (all interfaces unchanged)
- Measurement devices (all work unchanged)
- Test code (all tests pass)

**Optional enhancements:**
- Implement custom physics models for research
- Implement higher-order solvers for better accuracy
- Create solver benchmarks for performance tuning

### For Researchers

**Benefits:**
- Easy to experiment with alternative physics models
- Easy to verify numerical convergence
- Clear separation makes code easier to understand and modify

**Example research directions:**
- Implement relativistic corrections to spring forces
- Test different nonlinearity models (exponential, power-law, etc.)
- Compare energy conservation across solvers
- Implement adaptive timestep for stiff problems

## Implementation Details

### File Organization

```
src/main/java/sim/
  ├── Brane.java                    (Facade)
  ├── physics/
  │   ├── BranePhysics.java         (Interface)
  │   └── SpringMeshPhysics.java    (Implementation)
  └── solvers/
      ├── BraneSolver.java          (Interface)
      └── VelocityVerletSolver.java (Implementation)
```

### Code Statistics

**Before refactoring:**
- Brane.java: ~400 lines (physics + solver mixed)

**After refactoring:**
- Brane.java: ~295 lines (facade only)
- SpringMeshPhysics.java: ~200 lines (pure physics)
- VelocityVerletSolver.java: ~100 lines (pure solver)
- **Total:** Slightly more code, but much better organized

### Performance

- **No performance impact** - same algorithms, just reorganized
- All tests pass with identical results
- Compiler optimizations should inline delegate calls

## Future Work

### Potential Enhancements

1. **Additional Solvers:**
   - Runge-Kutta 4 (4th order)
   - Adaptive timestep (RK45, Dormand-Prince)
   - Implicit solvers (Backward Euler, Crank-Nicolson)

2. **Advanced Physics:**
   - Relativistic corrections
   - Alternative elasticity models
   - External field interactions

3. **Analysis Tools:**
   - Solver convergence tests
   - Energy conservation benchmarks
   - Stability analysis utilities

4. **Optimization:**
   - SIMD vectorization of force computation
   - GPU acceleration for physics and/or solver
   - Multi-physics coupling

## References

- **Hairer, Lubich, Wanner:** "Geometric Numerical Integration" (2006)
  Standard reference for symplectic integrators

- **Press et al:** "Numerical Recipes" (3rd ed)
  Practical guide to numerical methods

- **Swope et al:** "A computer simulation method..." J. Chem. Phys. 76, 637 (1982)
  Original Velocity Verlet paper

## Summary

This refactoring achieves:
- ✅ Clean separation of physics model from numerical solver
- ✅ Improved code organization and maintainability
- ✅ Extensibility for research and experimentation
- ✅ Full backward compatibility
- ✅ All tests passing
- ✅ No performance degradation
- ✅ Well-documented architecture

The simulation is now more modular, testable, and ready for future enhancements!