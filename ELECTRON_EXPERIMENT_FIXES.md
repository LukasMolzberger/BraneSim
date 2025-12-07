# Electron Experiment Fixes - Action Items

This document explains the remaining issues with the electron visualization and provides concrete fixes.

## Problem Summary

### Issue 1: Lateral Distortion Shows "Corner Gradient" Artifact

**Symptom:** The lateral distortion plots show a smooth gradient from bottom-left (dark) to top-right (bright), instead of showing electron-localized distortion.

**Root Cause:** The experiment script needs to ensure that:
1. The baseline positions use the same centering as the actual brane grid
2. The `config['center']` is properly passed to the visualization functions

**Current Status:**
- ✓ `baseline_state.py` module works correctly (test with `test_baseline_distortion.py`)
- ✓ `electron_visualization.py` functions updated to accept `baseline_positions`
- ✗ **The experiment script calling these functions may not be passing the correct center**

### Issue 2: XY Amplitude Slice Shows Single Ring

**Symptom:** The XY slice shows a single torus with a "gap" or discontinuity on the right side.

**Root Cause:** This is actually **NOT a bug**! The double-loop structure lives in the *cross-section* (local x, y coordinates), not in the global XY plane. The XY slice at Z=center_z cuts through y≈0, which is between the two lobes, so you see a radially symmetric pattern.

**Solution:** Use the debug plots to verify the cross-section structure:
```bash
python experiments/debug_electron_cross_section.py
```

This will create:
- `debug_cross_section.png` - Shows the double-lobe dumbbell in the local (x,y) frame
- `debug_phase_centerline.png` - Shows the m-winding phase pattern is smooth

## Required Fixes

### Fix 1: Update Experiment Script Grid Initialization

In your main electron experiment script (e.g., `electron_stability_test.py` or similar), ensure that:

**BEFORE:**
```python
# Grid setup
grid_shape = (nx, ny, nz)
h = 1e-15

# Initialize state (positions at grid points)
positions = create_grid_positions(grid_shape, h)  # Might be centered or not?
```

**AFTER:**
```python
# Grid setup with explicit centering
grid_shape = (nx, ny, nz)
h = 1e-15
center = (0.0, 0.0, 0.0)  # Or whatever center you want

# Store in config for visualization
config = {
    'grid_shape': grid_shape,
    'h': h,
    'center': center,  # ← CRITICAL: Must match how positions are initialized
    # ... other config ...
}

# Initialize state with same centering
positions = create_grid_positions(grid_shape, h, center=center)
```

### Fix 2: Pass Baseline to Visualization Functions

**In the experiment script, when calling visualization:**

```python
from branesim.physics.baseline_state import initialize_baseline_state
from experiments.electron_visualization import visualize_initial_state, create_all_animations

# Create baseline ONCE at start
baseline_info = initialize_baseline_state(config)
baseline_positions = baseline_info['positions']

# Pass to visualization functions
visualize_initial_state(
    state=initial_state,
    params=electron_params,
    config=config,
    baseline_positions=baseline_positions,  # ← Pass the baseline!
)

# Later, for animations
create_all_animations(
    states=state_snapshots,
    config=config,
    baseline_positions=baseline_positions,  # ← Pass the same baseline!
)
```

### Fix 3: Verify Grid Centering Matches

The grid positions in your `BraneState` must be initialized consistently with the `config['center']`:

**Check your grid initialization code:**

```python
def create_flat_grid_positions(grid_shape, h, center=(0.0, 0.0, 0.0)):
    """
    Create flat grid positions centered at specified center.

    This must match the logic in baseline_state.compute_flat_baseline_positions!
    """
    nx, ny, nz = grid_shape
    cx, cy, cz = center

    # Create grid
    x = torch.arange(nx, dtype=torch.float32) * h
    y = torch.arange(ny, dtype=torch.float32) * h
    z = torch.arange(nz, dtype=torch.float32) * h

    # Center the grid
    x = x - (nx - 1) * h / 2 + cx
    y = y - (ny - 1) * h / 2 + cy
    z = z - (nz - 1) * h / 2 + cz

    # Create meshgrid
    X, Y, Z = torch.meshgrid(x, y, z, indexing='ij')

    # Flatten and stack
    positions = torch.stack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)], dim=1)

    return positions
```

**CRITICAL:** If your grid uses `center=None` or doesn't center at all (corner at origin), then you MUST pass `center=None` to both:
- `initialize_baseline_state(config)` in the experiment
- The grid initialization itself

The centering must be **identical** between baseline and actual grid, otherwise you get the corner gradient artifact.

## Testing Procedure

### Step 1: Verify Debug Plots

```bash
cd /Users/lukasmolzberger/PycharmProjects/BraneSim
python experiments/debug_electron_cross_section.py
```

**Expected output:**
- `debug_cross_section.png` shows:
  - Two clear Gaussian lobes at y = ±ρ₀
  - Lobes centered at x = 0
  - Dumbbell or figure-8 shape

- `debug_phase_centerline.png` shows:
  - m=2 complete oscillation cycles
  - Smooth phase everywhere (max |Δφ| < 0.1 rad)
  - Sinusoidal amplitude field

If these look correct → electron initialization is working properly!

### Step 2: Verify Baseline Computation

```bash
python experiments/test_baseline_distortion.py
```

**Expected output:**
- All tests pass (✓)
- `test_baseline_artifact.png` shows:
  - Left panel: uniform ~zero distortion (correct baseline)
  - Right panel: gradient from corner (wrong baseline)

This confirms the baseline functions work correctly.

### Step 3: Fix and Re-run Experiment

1. Update your experiment script with Fixes 1-3 above
2. Ensure `config['center']` matches grid initialization
3. Pass `baseline_positions` to visualization functions
4. Re-run the experiment

**Expected result after fix:**
- Lateral distortion at t=0 should be ~zero everywhere (or very small, < 1e-12 m)
- Lateral distortion later should show localized features near the electron
- No corner gradient artifact!

### Step 4: Interpret the XY Slice Correctly

The XY slice will still look like a single torus with angular modulation because:
- The double-loop is in the **cross-section** (local x, y)
- The XY slice at Z=center cuts through y≈0 (between the lobes)

To see the double-lobe structure:
- Look at **XZ** or **YZ** slices (should show dumbbell)
- Use the `debug_cross_section.png` plot

The "gap" on the right is NOT a discontinuity - it's just where cos(phase) ≈ 0 and the amplitude is small.

## Common Mistakes

### Mistake 1: Using Different Centers

```python
# WRONG - centering mismatch!
positions = create_grid(grid_shape, h, center=(5e-15, 5e-15, 0))
config = {'center': None}  # ← Doesn't match!
```

**Fix:** Use same center everywhere.

### Mistake 2: Not Passing Baseline

```python
# WRONG - baseline computed separately each time
visualize_initial_state(state, params, config)  # Uses its own baseline
create_all_animations(states, config)  # Uses different baseline!
```

**Fix:** Create baseline once, pass to all functions.

### Mistake 3: Expecting Two Rings in XY

```python
# MISUNDERSTANDING
# "The XY slice should show two separate rings at R±ρ₀"
# ↑ This is not how the W&vdM double-loop works!
```

**Reality:** The two lobes are in the **binormal direction** (Z), not radially. You see them in XZ/YZ slices, not XY.

## Quick Checklist

Before running the full experiment, verify:

- [ ] `debug_electron_cross_section.py` produces correct plots
- [ ] `test_baseline_distortion.py` all tests pass
- [ ] Experiment config has explicit `center` field
- [ ] Grid initialization uses same `center` as config
- [ ] Baseline created once with `initialize_baseline_state(config)`
- [ ] Baseline passed to `visualize_initial_state(..., baseline_positions=...)`
- [ ] Baseline passed to `create_all_animations(..., baseline_positions=...)`
- [ ] Initial lateral distortion is ~zero (< 1e-12 m)

If all checkboxes are ✓, the visualization should be correct!

## Summary

The core issue is **not** with the electron initialization or the baseline computation - those are working correctly. The issue is with:

1. **Grid centering consistency:** The experiment script must use the same centering for actual grid and baseline
2. **Passing baseline through:** The visualization functions need the pre-computed baseline

Once these are fixed, the lateral distortion will show only electron-induced deformation, not grid artifacts.

The XY slice looking like "one ring" is **expected behavior** - the double-loop structure is in the cross-section, not in the global XY projection.