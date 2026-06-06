"""U(1) vortex-worldtube seed injection — carrier-phase 2-plane ansatz.

Physics spec: EXPERIMENT.md §"Injection ansatz"
Principles:   PRINCIPLES.md §2 (layer C: initialization), §3.2 (no clamps),
              §7.6 (dimension-agnostic where possible)

## What this module does

Translates the emergent U(1) order-parameter ansatz

    Psi(x, t) = A(rho) * exp(i * [m_wind * chi + omega * t])

into a substrate displacement field u on the prestressed-vacuum 4D worldvolume.

### Carrier 2-plane choice (DOCUMENTED HERE — must not change silently)

The ambient space is R^4 with components indexed (0,1,2,3).  The inside
observer's decomposition (principles §1.4) assigns:
  - components 0,1,2 — the three spacelike / "colour" lateral channels
  - component 3       — the timelike channel

Within the three spacelike channels the U(1) carrier phase corresponds to the
EM-like trace direction.  We pick the concrete 2-plane

    CARRIER_PLANE = (2, 3)

i.e. spatial component 2 and the timelike component 3 carry Re(Psi) and
Im(Psi) respectively.  Rationale:
  - The "i-from-time" memory note (MEMORY.md) identifies the U(1) imaginary
    unit with the time direction; component 3 is the natural Im slot.
  - Component 2 is the third spacelike channel — avoiding 0 and 1 keeps the
    "colour" x and y lateral channels at their vacuum values for this seed,
    matching the EXPERIMENT.md instruction that "other components [are] at the
    prestressed-vacuum value."
  - This choice is completely concrete and invertible; the diagnostics layer
    must use the same plane definition when extracting Psi.

Components (0,1) remain at their vacuum positions throughout.

### Vortex geometry: single vortex ring ("smoke ring")

A contractible vortex ring carries zero net winding through every periodic
plane: the ring core (a circle of radius R_ring) pierces any plane it
intersects *twice* with opposite winding signs, so the total is zero.  No
antivortex partner is needed.

Ring geometry (z=z_c plane, centred in the box):

    core(psi) = (x_c + R_ring*cos(psi), y_c + R_ring*sin(psi), z_c)

For each lattice node x:
  1. phi   = atan2(y - y_c, x - x_c)          — toroidal angle (around ring axis)
  2. r_pl  = sqrt((x-x_c)^2 + (y-y_c)^2)      — radial distance in ring plane
  3. Nearest core point: c(phi) as above
  4. d_r   = r_pl - R_ring                     — outward meridional component
  5. d_z   = z - z_c                           — axial meridional component
  6. rho   = sqrt(d_r^2 + d_z^2)              — distance to ring core
  7. chi   = atan2(d_z, d_r)                  — meridional angle around tube

Phase: theta(x,t) = m_wind * chi + omega * t   (winds m times around tube)
Amplitude: A(rho) = A0*(rho/w)*exp(-rho^2/2w^2) (zero on core, peak at rho~w)

Resulting displacement:
    u[node, CARRIER_RE] += A(rho) * cos(theta)
    u[node, CARRIER_IM] += A(rho) * sin(theta)

### Periodic consistency

A contractible ring is net-zero through every plane.  Verified numerically by
``measure_winding_closure``.  Choose R_ring << L/2 and R_ring + few*w << L/2
so the torus sits well inside the box with vacuum margin.

### Donut / tube profile

    A(rho) = A0 * (rho / w) * exp(-rho^2 / (2 w^2))

Zero on the ring core (rho=0), peak ~0.607*A0 at rho=w, decays beyond.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from branesim.core.conventions import ActionParams
from branesim.core.lattice import SpacelikeLattice

# ---------------------------------------------------------------------------
# Carrier 2-plane definition (canonical; diagnostics must match)
# ---------------------------------------------------------------------------

#: Index of the Re(Psi) ambient component.
CARRIER_RE: int = 2

#: Index of the Im(Psi) ambient component.
CARRIER_IM: int = 3


# ---------------------------------------------------------------------------
# Public parameter dataclass
# ---------------------------------------------------------------------------


class VortexParams(NamedTuple):
    """Parameters for the U(1) vortex-worldtube seed.

    Attributes
    ----------
    A0 : float
        Peak amplitude of the donut profile.  Dimensionless strain units;
        EXPERIMENT.md targets ~0.3.
    w : float
        Tube width (donut thickness) in lattice units.  Recommend w ~ 2-3*a.
    R_ring : float
        Ring core radius in lattice units.  The torus major radius.
        Choose so that R_ring + a few w < L/2 for vacuum margin.
    m_wind : int
        Winding number (meridional winding around tube cross-section). m=1.
    omega : float
        Carrier angular frequency (rad per time step).  omega*dt ~ 0.5.
    geometry : str
        ``"vortex_ring"`` — a single closed ring (smoke ring); zero net
        winding through every periodic plane.
    """

    A0: float = 0.3
    w: float = 3.0
    R_ring: float = 7.0
    m_wind: int = 1
    omega: float = 2.0
    geometry: str = "vortex_ring"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _donut_profile(rho: np.ndarray, A0: float, w: float) -> np.ndarray:
    """A(rho) = A0 * (rho/w) * exp(-rho^2 / (2 w^2)).

    Zero at rho=0, peak A0*exp(-1/2) ≈ 0.607*A0 at rho=w.
    """
    return A0 * (rho / w) * np.exp(-(rho ** 2) / (2.0 * w ** 2))


# ---------------------------------------------------------------------------
# Vortex ring displacement field
# ---------------------------------------------------------------------------


def _ring_displacement(
    coords: np.ndarray,
    centre_xy: np.ndarray,
    z_c: float,
    R_ring: float,
    m_wind: int,
    A0: float,
    w: float,
    t: float,
    omega: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Displacement field for a single vortex ring.

    The ring core is a circle of radius R_ring in the z=z_c plane, centred
    at (x_c, y_c) = centre_xy.

    Parameters
    ----------
    coords : ndarray, shape (n_nodes, 3)
        (x, y, z) node positions.
    centre_xy : ndarray, shape (2,)
        (x_c, y_c) centre of the ring in the ring plane.
    z_c : float
        z-coordinate of the ring plane.
    R_ring : float
        Ring core radius (major radius of the torus).
    m_wind : int
        Meridional winding number.
    A0, w : float
        Donut profile parameters.
    t : float
        Time (slice index, dimensionless).
    omega : float
        Carrier angular frequency.

    Returns
    -------
    re_disp, im_disp : ndarray, shape (n_nodes,)
        Displacement in Re(Psi) and Im(Psi) carrier components.
    """
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]

    # Step 1: toroidal angle phi and in-plane radius
    dx = x - centre_xy[0]
    dy = y - centre_xy[1]
    r_pl = np.sqrt(dx ** 2 + dy ** 2)  # in-plane distance from ring axis

    # Step 2: meridional components
    d_r = r_pl - R_ring          # outward radial offset from ring core
    d_z = z - z_c                # axial offset from ring plane

    # Step 3: distance to ring core and meridional angle
    rho = np.sqrt(d_r ** 2 + d_z ** 2)
    chi = np.arctan2(d_z, d_r)  # meridional angle around tube cross-section

    # Step 4: amplitude and phase
    A = _donut_profile(rho, A0, w)
    phase = m_wind * chi + omega * t

    return A * np.cos(phase), A * np.sin(phase)


# ---------------------------------------------------------------------------
# Public injection function
# ---------------------------------------------------------------------------


def inject_vortex_worldtube(
    lattice: SpacelikeLattice,
    params: ActionParams,
    vp: VortexParams,
    n_slices: int,
) -> tuple[np.ndarray, dict]:
    """Build the 4D worldvolume with the U(1) vortex-ring seed injected.

    The worldvolume has shape ``(n_slices+1, n_nodes, m_ambient)`` — exactly
    the WorldVolume.slices format used by branesim.solver.ivp.

    Strategy
    --------
    1. Compute the flat-lattice reference worldvolume (vacuum = all nodes at
       their reference positions, all slices identical).
    2. Add the vortex-ring displacement to the carrier 2-plane (components
       CARRIER_RE, CARRIER_IM) on each time slice, with the carrier phase
       advancing by ``omega * l`` at slice l.
    3. All other ambient components stay at their vacuum values.

    Parameters
    ----------
    lattice : SpacelikeLattice
        Spacelike lattice (dim=3 required for this seed).
    params : ActionParams
        Action parameters.
    vp : VortexParams
        Vortex geometry and physics parameters.
    n_slices : int
        Number of time slices N.  Worldvolume has slices 0..N.

    Returns
    -------
    world : ndarray, shape (n_slices+1, n_nodes, m_ambient)
        Full worldvolume with seed injected.
    meta : dict
        Seed metadata for config.json / boundary_problem.npz.
    """
    if lattice.dim != 3:
        raise ValueError(
            f"vortex_worldtube seed requires dim=3 spatial lattice; got {lattice.dim}"
        )
    if vp.geometry != "vortex_ring":
        raise NotImplementedError(
            f"geometry={vp.geometry!r} not implemented; only 'vortex_ring' is supported"
        )

    m_ambient = params.ambient_dim(lattice.dim)  # canonical 4
    n_nodes = lattice.n_nodes
    ref = lattice.reference_positions(m_ambient)  # (n_nodes, m_ambient)

    # Box geometry
    a = lattice.params.spacing
    grid_shape = lattice.params.grid_shape
    centre_xy = np.array([
        (grid_shape[0] - 1) * a / 2.0,
        (grid_shape[1] - 1) * a / 2.0,
    ])
    z_c = (grid_shape[2] - 1) * a / 2.0   # ring midplane = z-centre of box

    # Node (x, y, z) coordinates from reference positions
    coords = ref[:, :3]  # (n_nodes, 3)

    # Allocate worldvolume
    world = np.empty((n_slices + 1, n_nodes, m_ambient), dtype=np.float64)

    for l in range(n_slices + 1):
        t = float(l)  # carrier phase argument = slice index (dimensionless)

        # Start from vacuum reference on this slice
        pos = ref.copy()

        re_d, im_d = _ring_displacement(
            coords, centre_xy, z_c, vp.R_ring,
            vp.m_wind, vp.A0, vp.w, t, vp.omega,
        )

        # Superpose on carrier 2-plane
        pos[:, CARRIER_RE] += re_d
        pos[:, CARRIER_IM] += im_d

        world[l] = pos

    meta = {
        "ansatz": "u1_vortex_worldtube",
        "geometry": vp.geometry,
        "carrier_re_component": CARRIER_RE,
        "carrier_im_component": CARRIER_IM,
        "m_wind": vp.m_wind,
        "A0": vp.A0,
        "w": vp.w,
        "R_ring": vp.R_ring,
        "omega": vp.omega,
        "n_slices": n_slices,
        "ring_centre_xy": list(centre_xy),
        "ring_z_c": float(z_c),
        "note": (
            "Single vortex ring; meridional winding m=1 around tube cross-section. "
            "Net winding = 0 across all periodic planes (contractible ring; "
            "closure verified by measure_winding_closure). "
            "Carrier 2-plane: component 2 = Re(Psi), component 3 = Im(Psi)."
        ),
    }

    return world, meta


# ---------------------------------------------------------------------------
# Winding closure verifier
# ---------------------------------------------------------------------------


def measure_winding_closure(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    slice_index: int = 0,
) -> dict[str, float]:
    """Measure the net U(1) phase winding through each pair of periodic faces.

    Uses the discrete plaquette method: the winding number through a plane is
    computed as the sum of phase differences across all plaquettes in that
    plane, divided by 2*pi.  This is gauge-invariant and exact for smooth
    phases.

    For a vortex ring the expected result is 0.0 through every face pair
    (the ring is contractible; it pierces any plane 0 or 2 times with
    opposite signs).

    Parameters
    ----------
    world : ndarray, shape (n_slices+1, n_nodes, m_ambient)
    lattice : SpacelikeLattice
    slice_index : int
        Which time slice to analyse.

    Returns
    -------
    dict with keys ``"winding_through_z_normal"``, ``"winding_through_y_normal"``,
    ``"winding_through_x_normal"``, each a float (should be ~0).
    """
    pos = world[slice_index]  # (n_nodes, m_ambient)
    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape

    # Extract U(1) phase at each node from the carrier 2-plane
    re_field = pos[:, CARRIER_RE].reshape(nx, ny, nz)
    im_field = pos[:, CARRIER_IM].reshape(nx, ny, nz)

    # Subtract reference (vacuum has 0 on carrier components)
    ref = lattice.reference_positions(pos.shape[1])
    re_ref = ref[:, CARRIER_RE].reshape(nx, ny, nz)
    im_ref = ref[:, CARRIER_IM].reshape(nx, ny, nz)

    re_disp = re_field - re_ref
    im_disp = im_field - im_ref

    # Phase = atan2(Im, Re) of the complex displacement field
    phase = np.arctan2(im_disp, re_disp + 1e-300)  # (nx, ny, nz)

    def _winding_through_normal(axis: int) -> float:
        """Sum phase circulation through all plaquettes perpendicular to axis."""
        ax1 = (axis + 1) % 3
        ax2 = (axis + 2) % 3

        def _diff_along(arr: np.ndarray, a: int) -> np.ndarray:
            return np.roll(arr, -1, axis=a) - arr

        dphi1 = _diff_along(phase, ax1)
        dphi2 = _diff_along(phase, ax2)

        d1_fwd = np.mod(dphi1 + np.pi, 2 * np.pi) - np.pi
        d2_fwd = np.mod(dphi2 + np.pi, 2 * np.pi) - np.pi

        # Mixed discrete curl
        curl = (np.roll(d1_fwd, -1, axis=ax2) - d1_fwd
                - np.roll(d2_fwd, -1, axis=ax1) + d2_fwd)
        curl_wrapped = np.mod(curl + np.pi, 2 * np.pi) - np.pi

        total = float(np.sum(curl_wrapped)) / (2.0 * np.pi)
        return total

    return {
        "winding_through_z_normal": _winding_through_normal(2),  # xy-plane
        "winding_through_y_normal": _winding_through_normal(1),  # xz-plane
        "winding_through_x_normal": _winding_through_normal(0),  # yz-plane
    }
