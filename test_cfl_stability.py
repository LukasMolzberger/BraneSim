"""Test if smaller CFL factor provides stability"""
import sys
sys.path.append('.')

import torch
import numpy as np

from branesim.core.state import BraneState, Dimensionality
from branesim.core.grid import BraneGrid
from branesim.core.solver import VelocityVerletSolver
from branesim.physics.linear_tension_forces import LinearTensionForceComputer
from branesim.config.simulation_config import PhysicalConstants

constants = PhysicalConstants()

h = constants.lambda_C * 10.0
nx, ny = 400, 50
c = constants.c
tension = 1.0
sigma = tension / c**2

print("Testing different CFL factors...\n")

for cfl_factor in [0.001, 0.0001]:
    dt = cfl_factor * h / c
    print(f"CFL = {cfl_factor}: dt = {dt:.6e} s")
    
    device = torch.device('cpu')
    dtype = torch.float64
    
    state = BraneState((nx, ny), Dimensionality.TWO_D, device, dtype)
    state.initialize_flat_configuration(h)
    state.set_fixed_boundaries()
    
    grid = BraneGrid((nx, ny), Dimensionality.TWO_D, h, device)
    physics = LinearTensionForceComputer(tension, h)
    solver = VelocityVerletSolver(dt, sigma, physics, grid)
    
    # Small initial displacement
    center_idx = nx * ny // 2
    state.positions[center_idx, 3] = 1e-14  # 0.01 pm
    
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()
    
    # Run 1000 steps
    stable = True
    for step in range(1000):
        solver.step(state)
        max_xi = state.positions[:, 3].abs().max()
        
        if torch.isnan(max_xi) or max_xi > 1e-10:
            print(f"  UNSTABLE at step {step+1}: |ξ|={max_xi:.6e}")
            stable = False
            break
    
    if stable:
        print(f"  STABLE after 1000 steps: |ξ|={max_xi:.6e}")
    print()

print("Done")
