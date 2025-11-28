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

    CRITICAL CONSTRAINT: c = 3×10⁸ m/s is a PHYSICAL CONSTANT.
    The relation T = ρ_m · c² is ENFORCED, not computed.

    Physical Parameters:
        mass_density: float (ρ_m), mass per unit volume/area/length [kg/m³ or kg/m² or kg/m]
                      This is THE fundamental parameter. Tension is derived from this.
        grid_spacing: float (h), spatial resolution [m]
        time_step: float (dt) [s]

    Grid Parameters:
        grid_shape: tuple (nx,) or (nx, ny) or (nx, ny, nz)
        dimension: Dimensionality enum

    Spring Parameters:
        spring_constant: float (k) [N/m]
        rest_length: float (L_0) [m]
        critical_strain: float or None (ε_cr) for saturation

    Derived Parameters (computed automatically):
        tension: float (T), computed as T = ρ_m · c² [N/m]

    Computational:
        device: str ('cpu' or 'cuda')
        dtype: str ('float32' or 'float64')

    Output:
        save_interval: int (time steps between saves)
        output_dir: str
    """

    # Physical parameters - mass_density is fundamental
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

    # Computational
    device: str = 'cpu'
    dtype: str = 'float32'

    # Output
    save_interval: int = 10
    output_dir: str = './output'

    # Derived parameter (will be set in __post_init__)
    tension: float = None

    def __post_init__(self):
        """
        Derive tension from mass_density and validate configuration.

        ENFORCES: T = ρ_m · c²
        """
        constants = PhysicalConstants()

        # CRITICAL: Derive tension from mass density and physical constant c
        self.tension = self.mass_density * constants.c ** 2

        # Now validate
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
        wave_speed = self._get_effective_wave_speed()
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

    def verify_wave_speed(self) -> tuple[float, float, float]:
        """
        Verify that configured parameters reproduce the physical speed of light.

        CRITICAL: c = 3×10⁸ m/s is a PHYSICAL CONSTANT, not something we compute.
        This method checks if the discrete model parameters reproduce this value.

        Returns:
            Tuple of (expected_c, computed_c, relative_error)
        """
        constants = PhysicalConstants()
        expected_c = constants.c

        # Compute effective wave speed from current parameters
        if self.dimension == Dimensionality.ONE_D:
            # For 1D: use pre-tension formula
            T_0 = self.spring_constant * (self.grid_spacing - self.rest_length)
            mu = self.mass_density  # Linear density [kg/m]
            computed_c = np.sqrt(T_0 / mu)
        else:
            # For 2D/3D: use standard formula
            computed_c = np.sqrt(self.tension / self.mass_density)

        relative_error = abs(computed_c - expected_c) / expected_c

        return expected_c, computed_c, relative_error

    def _get_effective_wave_speed(self) -> float:
        """
        Internal helper to get effective wave speed for CFL checking.

        Returns:
            Computed wave speed in m/s
        """
        _, computed_c, _ = self.verify_wave_speed()
        return computed_c

    def compute_cfl_number(self) -> float:
        """
        Compute CFL number: c·dt/h.

        Returns:
            CFL number (should be < 0.5 for stability)
        """
        c = self._get_effective_wave_speed()
        return c * self.time_step / self.grid_spacing

    @classmethod
    def from_physical_units(
        cls,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        mass_density: float = None,
        lambda_C_multiplier: float = 10.0,
        cfl_factor: float = 0.4,
        **kwargs
    ) -> 'SimulationConfig':
        """
        Create config from physical constants.

        ENFORCES: c = 3×10⁸ m/s (physical constant)
        ENFORCES: T = ρ_m · c² (derived, not specified)

        Automatically computes grid spacing, time step, and spring constant
        based on Compton wavelength and the speed of light.

        Args:
            grid_shape: Grid dimensions
            dimension: Dimensionality enum
            mass_density: Mass density ρ_m [kg/m³ or kg/m² or kg/m]. If None, uses a default.
            lambda_C_multiplier: Grid spacing as multiple of Compton wavelength
            cfl_factor: CFL number (default 0.4 for stability margin)
            **kwargs: Override specific parameters (cannot override tension - it's derived!)

        Returns:
            SimulationConfig instance with physical wave speed c
        """
        constants = PhysicalConstants()

        # Derived parameters
        h = constants.lambda_C * lambda_C_multiplier  # Grid spacing
        c = constants.c  # Speed of light (PHYSICAL CONSTANT)
        dt = cfl_factor * h / c  # Time step from CFL

        # Mass density - use provided or default
        if mass_density is None:
            # Default: choose ρ_m such that T = 1.0 N
            # This gives ρ_m = T/c² = 1.0 / c²
            mass_density = 1.0 / c**2

        # Spring parameters
        L_0 = h  # Rest length = grid spacing
        # CRITICAL: Tension is DERIVED from mass density
        # (Will be computed in __post_init__ but we need k = T/L_0)
        T = mass_density * c**2
        k = T / L_0  # Spring constant

        # Build configuration dict
        config_dict = {
            'mass_density': mass_density,  # Fundamental parameter
            'grid_spacing': h,
            'time_step': dt,
            'grid_shape': grid_shape,
            'dimension': dimension,
            'spring_constant': k,
            'rest_length': L_0,
        }

        # Override with user-provided kwargs (except tension!)
        if 'tension' in kwargs:
            raise ValueError(
                "Cannot override 'tension' - it is derived from mass_density and c. "
                "Set 'mass_density' instead."
            )
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

    @classmethod
    def for_2d_test(
        cls,
        nx: int = 50,
        ny: int = 50,
        wave_speed: float = 1.0,
        spacing: float = 0.02,
        cfl_factor: float = 0.4,
        **kwargs
    ) -> 'SimulationConfig':
        """
        Create simple 2D test configuration.

        For a 2D membrane with springs, the wave speed is:
            c = √(T/ρ_s) where ρ_s [kg/m²] is surface mass density

        The effective 2D tension comes from the spring network:
            T = k·L₀ (effective surface tension)

        Args:
            nx: Number of grid points in x direction
            ny: Number of grid points in y direction
            wave_speed: Desired wave speed [m/s]
            spacing: Grid spacing [m]
            cfl_factor: CFL number
            **kwargs: Override parameters

        Returns:
            SimulationConfig for 2D simulation
        """
        # Compute parameters
        dt = cfl_factor * spacing / wave_speed

        # For 2D: surface mass density σ [kg/m²]
        sigma = 1.0  # kg/m² (mass per unit area)
        rho_s = sigma

        # Required tension for desired wave speed: T = c²·σ
        T = wave_speed**2 * sigma

        # For 2D, we set L₀ = h (no pre-tension in 2D springs)
        L_0 = spacing

        # Spring constant: k = T/L₀
        k = T / L_0

        config_dict = {
            'tension': T,
            'mass_density': rho_s,
            'grid_spacing': spacing,
            'time_step': dt,
            'grid_shape': (nx, ny),
            'dimension': Dimensionality.TWO_D,
            'spring_constant': k,
            'rest_length': L_0,
        }

        config_dict.update(kwargs)

        config = cls(**config_dict)

        # Print diagnostic info
        mass_per_point = sigma * spacing**2
        print(f"\n2D Configuration:")
        print(f"  Grid: {nx} × {ny} points")
        print(f"  Grid spacing h = {spacing:.6f} m")
        print(f"  Rest length L₀ = {L_0:.6f} m")
        print(f"  Spring constant k = {k:.2f} N/m")
        print(f"  Tension T = k·L₀ = {T:.6f} N")
        print(f"  Surface mass density σ = {sigma:.6f} kg/m²")
        print(f"  Mass per point m = σ·h² = {mass_per_point:.6f} kg")
        print(f"  Expected wave speed c = √(T/σ) = {np.sqrt(T/sigma):.6f} m/s")
        print(f"  Time step dt = {dt:.6f} s")
        print(f"  CFL = {cfl_factor:.3f}\n")

        return config

    @classmethod
    def from_compton_calibration(
        cls,
        grid_shape: Tuple[int, ...],
        dimension: Dimensionality,
        lambda_C_multiplier: float = 10.0,
        cfl_factor: float = 0.4,
        critical_strain: Optional[float] = None,
        **kwargs
    ) -> 'SimulationConfig':
        """
        Create 3D configuration using Compton-cell mass calibration (Route i from paper).

        This implements the amplitude scale calibration described in the paper:
        - Compton-cell assumption: ρ λ_C³ ≈ m_e
        - Bulk modulus: K = ρ c²
        - Point mass: m_point = ρ h³
        - Spring constant: k_spring = K h

        This ensures the lattice reproduces the continuum wave speed c² = K/ρ.

        Args:
            grid_shape: Grid dimensions (nx, ny, nz) for 3D
            dimension: Must be Dimensionality.THREE_D
            lambda_C_multiplier: Grid spacing as multiple of Compton wavelength
            cfl_factor: CFL number (default 0.4 for stability)
            critical_strain: Optional ε_cr for saturation nonlinearity
            **kwargs: Override specific parameters

        Returns:
            SimulationConfig with Compton-calibrated parameters

        Raises:
            ValueError: If dimension is not THREE_D

        Examples:
            >>> config = SimulationConfig.from_compton_calibration(
            ...     grid_shape=(32, 32, 32),
            ...     dimension=Dimensionality.THREE_D,
            ...     lambda_C_multiplier=10.0,
            ...     critical_strain=0.1
            ... )

        References:
            See paper Section "Amplitude scale calibration" (experimental-setting.tex)
        """
        if dimension != Dimensionality.THREE_D:
            raise ValueError(
                f"Compton calibration is only valid for 3D. Got {dimension}"
            )

        # Import here to avoid circular dependency
        from branesim.physics.parameters import compton_calibrated_brane_lattice_params

        constants = PhysicalConstants()

        # Grid spacing
        h = constants.lambda_C * lambda_C_multiplier

        # Compute calibrated parameters
        params = compton_calibrated_brane_lattice_params(
            grid_spacing_m=h,
            dimensionality=3,
            c=constants.c
        )

        # Time step from CFL condition
        dt = cfl_factor * h / constants.c

        # Mass density from Compton-cell calibration
        rho_m = params["rho_D"]

        # Spring constant from bulk modulus
        k = params["k_spring"]

        # Rest length = grid spacing (no pre-tension in 3D)
        L_0 = h

        # Build configuration
        config_dict = {
            'mass_density': rho_m,
            'grid_spacing': h,
            'time_step': dt,
            'grid_shape': grid_shape,
            'dimension': dimension,
            'spring_constant': k,
            'rest_length': L_0,
            'critical_strain': critical_strain,
        }

        # Override with user kwargs
        if 'tension' in kwargs:
            raise ValueError(
                "Cannot override 'tension' - it is derived from mass_density and c."
            )
        config_dict.update(kwargs)

        config = cls(**config_dict)

        # Print diagnostic info
        nx, ny, nz = grid_shape
        mass_per_point = params["m_point"]
        print(f"\n3D Compton-Calibrated Configuration:")
        print(f"  Grid: {nx} × {ny} × {nz} = {nx*ny*nz} points")
        print(f"  Reduced Compton wavelength λ_C = {constants.lambda_C:.4e} m")
        print(f"  Grid spacing h = {lambda_C_multiplier:.1f} λ_C = {h:.4e} m")
        print(f"  Mass density ρ = m_e/λ_C³ = {rho_m:.4e} kg/m³")
        print(f"  Bulk modulus K = ρ c² = {params['T_D']:.4e} Pa")
        print(f"  Spring constant k = K h = {k:.4e} N/m")
        print(f"  Point mass m = ρ h³ = {mass_per_point:.4e} kg")
        print(f"  Rest length L₀ = {L_0:.4e} m")
        if critical_strain is not None:
            print(f"  Critical strain ε_cr = {critical_strain:.4f}")
        print(f"  Time step dt = {dt:.4e} s")
        print(f"  CFL = {cfl_factor:.3f}\n")

        return config

    def __repr__(self) -> str:
        """String representation."""
        expected_c, computed_c, error = self.verify_wave_speed()
        cfl = self.compute_cfl_number()

        # Show warning if wave speed deviates significantly from physical constant
        if error > 0.01:  # > 1% error
            c_str = f"c_eff = {computed_c:.2e} m/s [WARNING: {error*100:.1f}% from c_physical]"
        else:
            c_str = f"c_eff = {computed_c:.2e} m/s"

        return (
            f"SimulationConfig(\n"
            f"  {self.dimension.name}, grid_shape={self.grid_shape}\n"
            f"  {c_str}, CFL = {cfl:.3f}\n"
            f"  h = {self.grid_spacing:.2e} m, dt = {self.time_step:.2e} s\n"
            f"  T = {self.tension:.2e} N/m (= ρ_m · c²), ρ_m = {self.mass_density:.2e} kg/m³\n"
            f"  k = {self.spring_constant:.2e} N/m, L_0 = {self.rest_length:.2e} m\n"
            f"  device={self.device}, dtype={self.dtype}\n"
            f")"
        )
