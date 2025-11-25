"""
SimulationConfig: Configuration management for brane simulations.

This module provides a dataclass-based configuration system with automatic
validation, including CFL stability checking.
"""

from dataclasses import dataclass
from typing import Tuple, Optional
import numpy as np

from branesim.core.state import Dimensionality


@dataclass
class PhysicalConstants:
    """Physical constants in SI units."""
    c: float = 299792458.0  # Speed of light [m/s]
    hbar: float = 1.054571817e-34  # Reduced Planck constant [J·s]
    m_e: float = 9.1093837015e-31  # Electron mass [kg]

    @property
    def lambda_C(self) -> float:
        """Compton wavelength [m]."""
        return self.hbar / (self.m_e * self.c)


@dataclass
class SimulationConfig:
    """
    Complete simulation configuration.

    Physical Parameters:
        tension: float (T), membrane tension [N/m]
        mass_density: float (ρ_m), mass per unit volume [kg/m³]
        grid_spacing: float (h), spatial resolution [m]
        time_step: float (dt) [s]

    Grid Parameters:
        grid_shape: tuple (nx,) or (nx, ny) or (nx, ny, nz)
        dimension: Dimensionality enum

    Spring Parameters:
        spring_constant: float (k) [N/m]
        rest_length: float (L_0) [m]
        critical_strain: float or None (ε_cr) for saturation

    Boundary Conditions:
        apply_boundary_tension: bool

    Computational:
        device: str ('cpu' or 'cuda')
        dtype: str ('float32' or 'float64')

    Output:
        save_interval: int (time steps between saves)
        output_dir: str
    """

    # Physical parameters
    tension: float
    mass_density: float
    grid_spacing: float
    time_step: float

    # Grid parameters
    grid_shape: Tuple[int, ...]
    dimension: Dimensionality

    # Spring parameters
    spring_constant: float
    rest_length: float
    critical_strain: Optional[float] = None

    # Boundary conditions
    apply_boundary_tension: bool = False

    # Computational
    device: str = 'cpu'
    dtype: str = 'float32'

    # Output
    save_interval: int = 10
    output_dir: str = './output'

    def __post_init__(self):
        """Validate configuration after initialization."""
        self._validate()

    def _validate(self):
        """
        Check parameter validity and stability.

        Raises:
            ValueError: If configuration is invalid or CFL condition violated
        """
        # Validate grid shape matches dimension
        if len(self.grid_shape) != self.dimension.value:
            raise ValueError(
                f"Grid shape {self.grid_shape} has {len(self.grid_shape)} dimensions, "
                f"but dimensionality is {self.dimension.value}D"
            )

        # Check positive values
        if self.tension <= 0:
            raise ValueError(f"Tension must be positive, got {self.tension}")
        if self.mass_density <= 0:
            raise ValueError(f"Mass density must be positive, got {self.mass_density}")
        if self.grid_spacing <= 0:
            raise ValueError(f"Grid spacing must be positive, got {self.grid_spacing}")
        if self.time_step <= 0:
            raise ValueError(f"Time step must be positive, got {self.time_step}")
        if self.spring_constant <= 0:
            raise ValueError(f"Spring constant must be positive, got {self.spring_constant}")

        # Check CFL stability condition
        wave_speed = self.compute_wave_speed()
        cfl_number = wave_speed * self.time_step / self.grid_spacing

        if cfl_number > 0.5:
            raise ValueError(
                f"CFL condition violated: c·dt/h = {cfl_number:.3f} > 0.5\n"
                f"Wave speed c = {wave_speed:.2e} m/s\n"
                f"Time step dt = {self.time_step:.2e} s\n"
                f"Grid spacing h = {self.grid_spacing:.2e} m\n"
                f"Reduce time step to dt < {0.5 * self.grid_spacing / wave_speed:.2e} s"
            )

        # Warn if spring constant doesn't match tension/rest_length
        expected_k = self.tension / self.rest_length
        relative_error = abs(self.spring_constant - expected_k) / expected_k
        if relative_error > 0.1:
            print(
                f"Warning: Spring constant k={self.spring_constant:.2e} differs "
                f"from T/L_0={expected_k:.2e} by {relative_error*100:.1f}%"
            )

        # Validate device
        if self.device not in ['cpu', 'cuda']:
            raise ValueError(f"Device must be 'cpu' or 'cuda', got '{self.device}'")

        # Validate dtype
        if self.dtype not in ['float32', 'float64']:
            raise ValueError(f"dtype must be 'float32' or 'float64', got '{self.dtype}'")

    def compute_wave_speed(self) -> float:
        """
        Compute wave speed based on dimensionality.

        1D: c = √(T₀/μ) where T₀ = k·(h - L₀), μ = mass_density [kg/m]
        2D/3D: c = √(T/ρ) where T is tension, ρ is mass density

        Returns:
            Wave speed in m/s
        """
        if self.dimension == Dimensionality.ONE_D:
            # For 1D: use pre-tension formula
            T_0 = self.spring_constant * (self.grid_spacing - self.rest_length)
            mu = self.mass_density  # Linear density [kg/m]
            return np.sqrt(T_0 / mu)
        else:
            # For 2D/3D: use standard formula
            return np.sqrt(self.tension / self.mass_density)

    def compute_cfl_number(self) -> float:
        """
        Compute CFL number: c·dt/h.

        Returns:
            CFL number (should be < 0.5 for stability)
        """
        c = self.compute_wave_speed()
        return c * self.time_step / self.grid_spacing

    @classmethod
    def from_physical_units(
        cls,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        lambda_C_multiplier: float = 10.0,
        cfl_factor: float = 0.4,
        **kwargs
    ) -> 'SimulationConfig':
        """
        Create config from physical constants.

        Automatically computes grid spacing, time step, and spring constant
        based on Compton wavelength and wave speed.

        Args:
            grid_shape: Grid dimensions
            dimension: Dimensionality enum
            lambda_C_multiplier: Grid spacing as multiple of Compton wavelength
            cfl_factor: CFL number (default 0.4 for stability margin)
            **kwargs: Override specific parameters

        Returns:
            SimulationConfig instance
        """
        constants = PhysicalConstants()

        # Derived parameters
        h = constants.lambda_C * lambda_C_multiplier  # Grid spacing
        c = constants.c  # Speed of light
        dt = cfl_factor * h / c  # Time step from CFL

        # Spring parameters
        L_0 = h  # Rest length = grid spacing
        T = 1.0  # Arbitrary tension (can be calibrated later)
        rho_m = T / c**2  # Mass density from wave speed c = √(T/ρ_m)
        k = T / L_0  # Spring constant

        # Build configuration dict
        config_dict = {
            'tension': T,
            'mass_density': rho_m,
            'grid_spacing': h,
            'time_step': dt,
            'grid_shape': grid_shape,
            'dimension': dimension,
            'spring_constant': k,
            'rest_length': L_0,
        }

        # Override with user-provided kwargs
        config_dict.update(kwargs)

        return cls(**config_dict)

    @classmethod
    def for_1d_test(
        cls,
        nx: int = 100,
        wave_speed: float = 1.0,
        spacing: float = 0.01,
        cfl_factor: float = 0.4,
        **kwargs
    ) -> 'SimulationConfig':
        """
        Create simple 1D test configuration.

        For a 1D pre-tensioned string, the wave speed is:
            c = √(T₀/μ) where μ = ρ_m·h is mass per unit length

        The pre-tension T₀ comes from pre-stretching the springs:
            T₀ = k·(h - L₀)

        This gives: c = √(k·(h - L₀)/(ρ_m·h))

        Args:
            nx: Number of grid points
            wave_speed: Desired wave speed [m/s]
            spacing: Grid spacing [m]
            cfl_factor: CFL number
            **kwargs: Override parameters

        Returns:
            SimulationConfig for 1D simulation
        """
        # Compute parameters
        dt = cfl_factor * spacing / wave_speed

        # For 1D: mass per unit length μ [kg/m]
        # Set μ = ρ_m where ρ_m is treated as linear mass density
        # (NOT volumetric density - units are kg/m, not kg/m³)
        mu = 1.0  # kg/m (mass per unit length)
        rho_m = mu  # For 1D, store as mass_density for consistency

        # Required tension for desired wave speed: T₀ = c²·μ
        T_0 = wave_speed**2 * mu

        # Target a reasonable pre-strain (e.g., 1%)
        target_prestrain = 0.01
        L_0_target = spacing / (1 + target_prestrain)

        # Required spring constant: k = T₀/(h - L₀)
        k = T_0 / (spacing - L_0_target)

        # Rest length from T₀ = k·(h - L₀)
        L_0 = spacing - T_0 / k

        # Validate that we have pre-stretch (L₀ < h)
        if L_0 >= spacing:
            raise ValueError(
                f"Cannot achieve wave speed {wave_speed} m/s with k={k} N/m. "
                f"Need larger k or smaller wave speed."
            )

        pre_strain = (spacing - L_0) / L_0

        # Tension parameter (not directly used in 1D, but keep for consistency)
        T = k * L_0

        config_dict = {
            'tension': T,
            'mass_density': rho_m,
            'grid_spacing': spacing,
            'time_step': dt,
            'grid_shape': (nx,),
            'dimension': Dimensionality.ONE_D,
            'spring_constant': k,
            'rest_length': L_0,
        }

        config_dict.update(kwargs)

        config = cls(**config_dict)

        # Print diagnostic info about pre-tension
        mass_per_point = mu * spacing
        print(f"\n1D Configuration with Pre-Tension:")
        print(f"  Grid spacing h = {spacing:.6f} m")
        print(f"  Rest length L₀ = {L_0:.6f} m")
        print(f"  Pre-strain ε₀ = {pre_strain:.6f} ({pre_strain*100:.2f}%)")
        print(f"  Spring constant k = {k:.2f} N/m")
        print(f"  Pre-tension T₀ = k·(h-L₀) = {T_0:.6f} N")
        print(f"  Linear mass density μ = {mu:.6f} kg/m")
        print(f"  Mass per point m = μ·h = {mass_per_point:.6f} kg")
        print(f"  Expected wave speed c = √(T₀/μ) = {np.sqrt(T_0/mu):.6f} m/s\n")

        return config

    def __repr__(self) -> str:
        """String representation."""
        c = self.compute_wave_speed()
        cfl = self.compute_cfl_number()
        return (
            f"SimulationConfig(\n"
            f"  {self.dimension.name}, grid_shape={self.grid_shape}\n"
            f"  c = {c:.2e} m/s, CFL = {cfl:.3f}\n"
            f"  h = {self.grid_spacing:.2e} m, dt = {self.time_step:.2e} s\n"
            f"  T = {self.tension:.2e} N/m, ρ_m = {self.mass_density:.2e} kg/m³\n"
            f"  k = {self.spring_constant:.2e} N/m, L_0 = {self.rest_length:.2e} m\n"
            f"  device={self.device}, dtype={self.dtype}\n"
            f")"
        )
