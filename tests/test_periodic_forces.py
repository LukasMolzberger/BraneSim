import unittest

import torch

from components.shared.state import BraneState3D
from components.simulation.forces import SpringForceComputer
from components.simulation.grid import BraneGrid3D


class PeriodicForceTests(unittest.TestCase):
    def test_flat_periodic_lattice_has_no_boundary_force(self):
        device = torch.device("cpu")
        dtype = torch.float64
        grid = BraneGrid3D((4, 4, 4), spacing=1.0, device=device, periodic_axes=(True, True, True))
        state = BraneState3D(grid.grid_shape, device=device, dtype=dtype)
        state.initialize_flat_configuration(grid.spacing)

        forces = SpringForceComputer(spring_constant=1.0, rest_length=0.2).compute_forces(state, grid)

        self.assertLess(float(torch.max(torch.abs(forces))), 1e-12)


if __name__ == "__main__":
    unittest.main()
