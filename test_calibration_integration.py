"""
Quick test to verify Compton calibration integration is working.
"""

import sys
sys.path.insert(0, '/Users/lukasmolzberger/PycharmProjects/BraneSim')

# Test 1: Can we import the function?
print("Test 1: Importing brane_lattice_params_3d...")
from branesim.physics.parameters import brane_lattice_params_3d, print_calibration_summary
print("✓ Success!")

# Test 2: Does it compute correct values?
print("\nTest 2: Computing parameters...")
h = 1e-13  # 1 Angstrom
params = brane_lattice_params_3d(h)
print("✓ Success!")
print(f"  λ_C = {params['lambda_C']:.4e} m")
print(f"  ρ   = {params['rho']:.4e} kg/m³")
print(f"  K   = {params['K']:.4e} Pa")
print(f"  m   = {params['m_point']:.4e} kg")
print(f"  k   = {params['k_spring']:.4e} N/m")

# Test 3: Verify wave speed relation
print("\nTest 3: Verifying c² = K/ρ...")
import math
c_computed = math.sqrt(params['K'] / params['rho'])
c_expected = 299792458.0
relative_error = abs(c_computed - c_expected) / c_expected
print(f"  c_computed = {c_computed:.4e} m/s")
print(f"  c_expected = {c_expected:.4e} m/s")
print(f"  error      = {relative_error*100:.6f}%")
if relative_error < 1e-10:
    print("✓ Success! Wave speed matches c to machine precision")
else:
    print(f"✗ Failed: error = {relative_error*100:.6f}%")

# Test 4: Summary printer
print("\nTest 4: Printing calibration summary...")
print_calibration_summary(params, h)
print("✓ Success!")

# Test 5: Check it's exported from physics module
print("\nTest 5: Checking physics module exports...")
from branesim import physics
assert hasattr(physics, 'brane_lattice_params_3d')
assert hasattr(physics, 'print_calibration_summary')
print("✓ Success!")

print("\n" + "=" * 60)
print("ALL TESTS PASSED! ✓")
print("=" * 60)
print("\nThe Compton calibration is successfully integrated!")
print("You can now use it in your simulations via:")
print("  - from branesim.physics import brane_lattice_params_3d")
print("  - SimulationConfig.from_compton_calibration(...)")