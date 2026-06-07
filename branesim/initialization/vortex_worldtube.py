"""U(1) vortex-worldtube seed injection — spherical-harmonic carrier ansatz.

Physics spec: EXPERIMENT.md §"Injection ansatz"
Principles:   PRINCIPLES.md §2 (layer C: initialization), §3.2 (no clamps),
              §7.6 (dimension-agnostic where possible)
Memory:       spherical harmonics are the soliton-layer (L5) *description
              language* (project_soliton_layer_description_language); they seed
              the substrate field, they are not a structural commitment of it.

## What this module does

Translates the emergent U(1) order-parameter ansatz, written in the natural
soliton-layer basis (spherical harmonics about the object centre), into a
substrate displacement field u on the prestressed-vacuum 4D worldvolume.

    Psi(r, theta, phi, t) = A0 * Rhat(r) * Yhat_l^m(theta, phi) * exp(i*omega*t)

where (r, theta, phi) are spherical coordinates about the box centre, Y_l^m is
the (complex) spherical harmonic, Rhat is a localizing radial profile, and
exp(i*omega*t) is the closure-locked temporal carrier.

### Why spherical harmonics, not a hand-built torus (DESIGN, 2026-06-06)

The previous implementation built a literal "smoke ring" torus and wound the
phase *meridionally* around the tube cross-section.  That is wrong for a U(1)
line vortex: going around the ring azimuthally the phase was constant, so the
XY plan view showed no azimuthal winding.  The spherical-harmonic ansatz fixes
this at the root.  For the canonical EM/electron seed Y_1^1:

  - exp(i*phi) factor  -> the U(1) phase winds m times **azimuthally around the
    z-axis** (the vortex axis = the "donut hole").  On the XY midplane the phase
    advances 2*pi*m around the centre — the defining feature of a line vortex.
  - |Y_1^1| ~ sin(theta) -> amplitude zero on the axis (theta=0,pi), peak at the
    equator.  The energy density |Psi|^2 ~ sin^2(theta) is therefore a donut
    around the axis WITHOUT any literal torus being constructed.  The donut
    *emerges* from the angular harmonic.
  - the radial profile Rhat(r) only localizes the object; it sets no winding.

(l, m) are parameters.  l=1, m=1 is the canonical EM/U(1) vortex seed.

### Carrier 2-plane choice (DOCUMENTED HERE — must not change silently)

The ambient space is R^4 with components indexed (0,1,2,3).  The inside
observer's decomposition (principles §1.4) assigns:
  - components 0,1,2 — the three spacelike / "colour" lateral channels
  - component 3       — the timelike channel

The complex order parameter Psi populates the concrete 2-plane

    CARRIER_PLANE = (2, 3)

i.e. spatial component 2 carries Re(Psi) and timelike component 3 carries
Im(Psi).  Rationale (unchanged from the original seed):
  - The "i-from-time" memory note identifies the U(1) imaginary unit with the
    time direction; component 3 is the natural Im slot.
  - Component 2 is the third spacelike channel; components (0,1) remain at their
    vacuum positions, so the "colour" x/y lateral channels are untouched by the
    bare seed (the full U(3) field is then free to relax — SU(3) coexcitation is
    measured by the per-colour diagnostic).
  - The choice is concrete and invertible; the diagnostics layer must use the
    same plane definition when extracting Psi.

### Temporal carrier is QUANTIZED by loop closure (not a free knob)

EXPERIMENT.md: the phase must close around the time loop, omega*T + tumble =
2*pi*integer.  We therefore parameterize the carrier by an integer number of
full turns ``n_t`` over the whole ``n_slices`` loop; the per-slice increment is

    omega_per_slice = 2*pi*n_t / n_slices

so the worldtube is single-valued across the periodic time seam.  ``n_t=2`` is a
clean 720 deg carrier (the spin-1/2 double-cover target).

### Periodicity / closure

The seed is a *localized* excitation with a vacuum margin to the box faces, so
the field is continuous across the periodic faces (vacuum = vacuum).  The
azimuthal winding m is the vortex charge of the object, measured through the
central xy-plane by ``measure_winding_closure`` (expected ~m about the z-axis,
~0 about x and y).  In the full periodic 3-torus the localized object is
contractible (a *semilocal* vortex, EXPERIMENT.md §"What this answers": binding
is a dynamical, not topological, condition).
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np
from scipy.special import sph_harm_y

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
    """Parameters for the U(1) spherical-harmonic vortex seed.

    Attributes
    ----------
    A0 : float
        Peak amplitude of the order parameter (dimensionless strain units;
        EXPERIMENT.md targets ~0.3).  The angular and radial profiles are each
        normalized to unit peak, so ``A0`` is the true peak |Psi|.
    r0 : float
        Radial-shell peak (lattice units) — the radius at which the localizing
        profile peaks.  Sets the object's overall size, NOT a literal core
        circle.  Choose so r0 + a few*w < L/2 for a vacuum margin.
    w : float
        Radial-shell width (lattice units).
    l : int
        Spherical-harmonic degree.
    m : int
        Spherical-harmonic order = azimuthal U(1) winding around the z-axis.
        l=1, m=1 is the canonical EM/U(1) vortex seed (sin(theta) donut).
    n_t : int
        Number of full carrier turns over the whole time loop (temporal
        winding).  Quantized by loop closure; per-slice increment is
        2*pi*n_t/n_slices.  n_t=2 -> 720 deg (spin-1/2 double-cover target).
    geometry : str
        ``"spherical_harmonic"`` — the only supported geometry.
    """

    A0: float = 0.3
    r0: float = 6.0
    w: float = 2.5
    l: int = 1
    m: int = 1
    n_t: int = 2
    geometry: str = "spherical_harmonic"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _radial_profile(r: np.ndarray, r0: float, w: float) -> np.ndarray:
    """Unit-peak localizing radial shell  Rhat(r) = exp(-(r-r0)^2 / (2 w^2)).

    Peaks at 1 on the shell r=r0, decays with width w.  For r0>0 it is small
    at the origin, so combined with the sin(theta) angular zero the object is a
    localized donut with a vacuum core on the axis and at large r.
    """
    return np.exp(-((r - r0) ** 2) / (2.0 * w ** 2))


def _sph_harm_displacement(
    coords: np.ndarray,
    centre: np.ndarray,
    l: int,
    m: int,
    A0: float,
    r0: float,
    w: float,
    t: float,
    omega_per_slice: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Displacement field for the spherical-harmonic carrier seed.

    Parameters
    ----------
    coords : ndarray, shape (n_nodes, 3)
        (x, y, z) node positions.
    centre : ndarray, shape (3,)
        (x_c, y_c, z_c) object centre (box centre).
    l, m : int
        Spherical-harmonic degree and order.
    A0, r0, w : float
        Peak amplitude, radial-shell peak, radial-shell width.
    t : float
        Slice index (dimensionless); the carrier phase argument.
    omega_per_slice : float
        Closure-locked per-slice carrier increment (rad/slice).

    Returns
    -------
    re_disp, im_disp : ndarray, shape (n_nodes,)
        Re(Psi) and Im(Psi) displacement on the carrier 2-plane.
    """
    dx = coords[:, 0] - centre[0]
    dy = coords[:, 1] - centre[1]
    dz = coords[:, 2] - centre[2]

    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    r_safe = np.where(r < 1e-12, 1.0, r)

    # Spherical angles: theta = polar/colatitude in [0, pi], phi = azimuth.
    theta = np.arccos(np.clip(dz / r_safe, -1.0, 1.0))
    phi = np.arctan2(dy, dx)

    # Complex spherical harmonic Y_l^m(theta, phi) (scipy>=1.15 signature).
    Y = sph_harm_y(l, m, theta, phi)  # complex, shape (n_nodes,)

    # Normalize angular part to unit peak so A0 is the true peak amplitude.
    y_max = float(np.max(np.abs(Y)))
    if y_max < 1e-30:
        y_max = 1.0
    Yhat = Y / y_max

    R = _radial_profile(r, r0, w)  # unit-peak radial shell

    # Full complex order parameter, with the closure-locked carrier.
    carrier = np.exp(1j * omega_per_slice * t)
    Psi = A0 * R * Yhat * carrier

    # On the axis (r=0) sin(theta) -> 0 already; force exact vacuum there.
    Psi = np.where(r < 1e-12, 0.0 + 0.0j, Psi)

    return np.real(Psi), np.imag(Psi)


# ---------------------------------------------------------------------------
# Public injection function
# ---------------------------------------------------------------------------


def inject_vortex_worldtube(
    lattice: SpacelikeLattice,
    params: ActionParams,
    vp: VortexParams,
    n_slices: int,
) -> tuple[np.ndarray, dict]:
    """Build the 4D worldvolume with the U(1) spherical-harmonic seed injected.

    The worldvolume has shape ``(n_slices+1, n_nodes, m_ambient)`` — exactly
    the WorldVolume.slices format (branesim.solver.worldvolume).

    Strategy
    --------
    1. Flat-lattice reference worldvolume (vacuum = all nodes at reference,
       carrier components 0 on every slice).
    2. Add the spherical-harmonic displacement to the carrier 2-plane on each
       time slice, the carrier phase advancing by the closure-locked increment.
    3. All other ambient components stay at vacuum.

    Returns
    -------
    world : ndarray, shape (n_slices+1, n_nodes, m_ambient)
    meta : dict
        Seed metadata for config.json / boundary_problem.npz.
    """
    if lattice.dim != 3:
        raise ValueError(
            f"vortex_worldtube seed requires dim=3 spatial lattice; got {lattice.dim}"
        )
    if vp.geometry != "spherical_harmonic":
        raise NotImplementedError(
            f"geometry={vp.geometry!r} not implemented; "
            "only 'spherical_harmonic' is supported"
        )
    if abs(vp.m) > vp.l:
        raise ValueError(f"spherical harmonic requires |m|<=l; got l={vp.l}, m={vp.m}")

    m_ambient = params.ambient_dim(lattice.dim)  # canonical 4
    n_nodes = lattice.n_nodes
    ref = lattice.reference_positions(m_ambient)  # (n_nodes, m_ambient)

    # Box geometry — object centred in the box.
    a = lattice.params.spacing
    grid_shape = lattice.params.grid_shape
    centre = np.array([
        (grid_shape[0] - 1) * a / 2.0,
        (grid_shape[1] - 1) * a / 2.0,
        (grid_shape[2] - 1) * a / 2.0,
    ])

    coords = ref[:, :3]  # (n_nodes, 3)

    # Carrier rate quantized by time-loop closure (rad per slice).  Physical
    # angular frequency = omega_per_slice / dt.
    omega_per_slice = 2.0 * math.pi * vp.n_t / n_slices

    world = np.empty((n_slices + 1, n_nodes, m_ambient), dtype=np.float64)

    for l_slice in range(n_slices + 1):
        t = float(l_slice)  # carrier phase argument = slice index

        pos = ref.copy()
        re_d, im_d = _sph_harm_displacement(
            coords, centre, vp.l, vp.m, vp.A0, vp.r0, vp.w, t, omega_per_slice,
        )
        pos[:, CARRIER_RE] += re_d
        pos[:, CARRIER_IM] += im_d

        world[l_slice] = pos

    meta = {
        "ansatz": "u1_spherical_harmonic_vortex",
        "geometry": vp.geometry,
        "carrier_re_component": CARRIER_RE,
        "carrier_im_component": CARRIER_IM,
        "l": vp.l,
        "m": vp.m,
        "A0": vp.A0,
        "r0": vp.r0,
        "w": vp.w,
        "n_t": vp.n_t,
        "n_slices": n_slices,
        "omega_per_slice": omega_per_slice,
        "omega_phys": omega_per_slice / params.dt,
        "carrier_total_turns": vp.n_t,
        "carrier_total_deg": 360.0 * vp.n_t,
        "centre": list(centre),
        "note": (
            f"Spherical-harmonic U(1) vortex seed Y_{vp.l}^{vp.m}. "
            f"Azimuthal winding m={vp.m} around the z-axis (the vortex axis); "
            "energy density ~|Y|^2 forms a donut around that axis. "
            f"Carrier advances n_t={vp.n_t} full turns ({360 * vp.n_t} deg) over "
            "the time loop and closes exactly at the periodic seam "
            "(omega quantized by closure, not free). "
            "Carrier 2-plane: component 2 = Re(Psi), component 3 = Im(Psi)."
        ),
    }

    return world, meta


# ---------------------------------------------------------------------------
# Winding measurement
# ---------------------------------------------------------------------------


def measure_winding_closure(
    world: np.ndarray,
    lattice: SpacelikeLattice,
    slice_index: int = 0,
) -> dict[str, float]:
    """Measure the net U(1) azimuthal winding about the three central axes.

    The winding about an axis is the discrete phase circulation
    (1/2pi) * sum of wrapped phase increments around a square contour centred
    on the box centre, in the plane perpendicular to that axis.  A contour
    *enclosing* the vortex axis is required — summing vorticity over a whole
    periodic plane always gives 0 (total winding on a 2-torus cancels), so it
    cannot see the enclosed charge.  This is gauge-invariant and integer for a
    smooth phase with the contour in a single-valued region.

    For the spherical-harmonic seed Y_l^m the expected result is:
      - winding_through_z_normal (central xy-plane)  ~  m   (azimuthal winding
        around the z = vortex axis)
      - winding_through_y_normal (central xz-plane)  ~  0
      - winding_through_x_normal (central yz-plane)  ~  0

    Parameters
    ----------
    world : ndarray, shape (n_slices+1, n_nodes, m_ambient)
    lattice : SpacelikeLattice
    slice_index : int
        Which time slice to analyse.

    Returns
    -------
    dict with keys ``"winding_through_z_normal"``, ``"winding_through_y_normal"``,
    ``"winding_through_x_normal"`` (floats).
    """
    pos = world[slice_index]  # (n_nodes, m_ambient)
    grid_shape = lattice.params.grid_shape
    nx, ny, nz = grid_shape

    re_field = pos[:, CARRIER_RE].reshape(nx, ny, nz)
    im_field = pos[:, CARRIER_IM].reshape(nx, ny, nz)

    ref = lattice.reference_positions(pos.shape[1])
    re_ref = ref[:, CARRIER_RE].reshape(nx, ny, nz)
    im_ref = ref[:, CARRIER_IM].reshape(nx, ny, nz)

    re_disp = re_field - re_ref
    im_disp = im_field - im_ref

    phase = np.arctan2(im_disp, re_disp + 1e-300)  # (nx, ny, nz)

    def _wrap(d: np.ndarray | float) -> np.ndarray | float:
        return np.mod(d + np.pi, 2 * np.pi) - np.pi

    def _winding_about_axis(axis: int) -> float:
        """Phase circulation around a square contour enclosing the central axis."""
        ax1 = (axis + 1) % 3
        ax2 = (axis + 2) % 3
        centre_idx = (nx, ny, nz)[axis] // 2
        c1 = (nx, ny, nz)[ax1] // 2
        c2 = (nx, ny, nz)[ax2] // 2
        h = min((nx, ny, nz)[ax1], (nx, ny, nz)[ax2]) // 4  # contour half-width

        # Extract the central plane perpendicular to `axis`, as (ax1, ax2).
        sl = [slice(None), slice(None), slice(None)]
        sl[axis] = centre_idx
        plane = phase[tuple(sl)]
        if ax1 > ax2:
            plane = plane.T  # now indexed [ax1, ax2]

        # Walk the square contour [c1-h, c1+h] x [c2-h, c2+h] counter-clockwise.
        i0, i1 = c1 - h, c1 + h
        j0, j1 = c2 - h, c2 + h
        path = []
        for i in range(i0, i1):
            path.append((i, j0))
        for j in range(j0, j1):
            path.append((i1, j))
        for i in range(i1, i0, -1):
            path.append((i, j1))
        for j in range(j1, j0, -1):
            path.append((i0, j))

        total = 0.0
        for k in range(len(path)):
            a = plane[path[k]]
            b = plane[path[(k + 1) % len(path)]]
            total += float(_wrap(b - a))
        return total / (2.0 * np.pi)

    return {
        "winding_through_z_normal": _winding_about_axis(2),  # about z (xy contour)
        "winding_through_y_normal": _winding_about_axis(1),  # about y (xz contour)
        "winding_through_x_normal": _winding_about_axis(0),  # about x (yz contour)
    }