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

    CARRIER_PLANE = (trace direction of (0,1,2),  3)

i.e. Re(Psi) is written along the SYMMETRIC TRACE direction of the lateral
triplet — all three spacelike components equally, with unit weights
``CARRIER_RE_WEIGHTS = (1,1,1)/sqrt(3)`` on components
``CARRIER_RE_COMPONENTS = (0,1,2)`` — and timelike component 3 carries Im(Psi).
Rationale:
  - The "i-from-time" memory note identifies the U(1) imaginary unit with the
    time direction; component 3 is the natural Im slot.
  - **E1 fix (2026-06-07).** Writing Re(Psi) into a *single* lateral component
    (the old ``CARRIER_RE = 2``) is NOT a pure EM/U(1) vortex: under the U(3)
    projection (alpha_separability.projection_operators) a single component reads
    1/3 trace (U(1)/EM) + 2/3 traceless (SU(3)/colour), so the "EM/electron"
    object was mislabelled (D6 read ``u1_fraction = 0.333``).  Writing Re(Psi)
    along the trace direction (1,1,1)/sqrt(3) makes the bare seed *purely* U(1):
    P_SU3 @ (v * (1,1,1)/sqrt(3)) = 0, so D6 reads ``u1_fraction -> 1``.  The
    full U(3) field is then free to relax — SU(3) coexcitation, if any, is
    measured by the per-colour diagnostic, not injected by construction.
  - The choice is concrete and invertible: ``project_carrier_re`` recovers
    Re(Psi) by projecting the lateral displacement back onto the unit trace
    vector.  The diagnostics layer MUST use that same projection when extracting
    Psi (every consumer routes through ``project_carrier_re``).

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

#: Lateral components carrying Re(Psi) — the symmetric TRACE direction of the
#: spacelike triplet, so the bare seed is a pure U(1)/EM vortex (E1 fix).
CARRIER_RE_COMPONENTS: tuple[int, int, int] = (0, 1, 2)

#: Unit weights of the trace direction (1,1,1)/sqrt(3).  Re(Psi) is written as
#: ``re_d * CARRIER_RE_WEIGHTS`` across CARRIER_RE_COMPONENTS, and recovered by
#: projecting the lateral displacement back onto this unit vector.
_INV_SQRT3: float = 1.0 / math.sqrt(3.0)
CARRIER_RE_WEIGHTS: tuple[float, float, float] = (_INV_SQRT3, _INV_SQRT3, _INV_SQRT3)

#: Index of the Im(Psi) ambient component (the timelike channel; "i from time").
CARRIER_IM: int = 3

#: Trace-direction unit vector as an array (for projection / broadcasting).
_TRACE_WEIGHTS = np.asarray(CARRIER_RE_WEIGHTS)


def project_carrier_re(disp_lateral: np.ndarray) -> np.ndarray:
    """Re(Psi) from a lateral ``(..., 3)`` displacement, by projecting onto the
    unit trace direction ``(1,1,1)/sqrt(3)``.

    This is the exact inverse of the trace-direction write performed by the
    injector: a pure-trace displacement ``v * (1,1,1)/sqrt(3)`` projects back to
    exactly ``v``.  Every diagnostic that reads Re(Psi) MUST route through this
    helper so the injector and the measurement layer share one definition.
    """
    return np.asarray(disp_lateral) @ _TRACE_WEIGHTS


# ---------------------------------------------------------------------------
# Prestressed periodic-vacuum offset (derived 2026-06-07; physics-derivation)
# ---------------------------------------------------------------------------


def vacuum_offsets(n_slices: int, r_t: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-slice carrier-plane offsets of the prestressed closed-loop vacuum.

    The r_t>0 temporal central-force spring wants every consecutive pair of
    slices separated by its rest length r_t.  A naive seed (all slices = same
    reference) leaves the temporal bond = 0 in vacuum, so the spring's
    ``ΔR/|ΔR|`` direction term has a 1/|ΔR| Jacobian that blows up (the 64³ run
    failed this way).  The fix (derived): add a per-slice GLOBAL translation
    (the SAME offset on every node, so spatial bonds — within-slice differences —
    are untouched) tracing a regular N-gon in the carrier 2-plane (Re-direction,
    Im-direction):

        v_l = ρ ( cos(2π l/N) ê_Re + sin(2π l/N) ê₃ ),   ρ = r_t / (2 sin(π/N))

    where ``ê_Re`` is the UNIT trace direction (1,1,1)/√3 of components (0,1,2)
    (E1 fix) — so this ``off_re`` magnitude is distributed across comps 0,1,2 with
    ``CARRIER_RE_WEIGHTS``, preserving ‖v_l‖ = ρ — and ê₃ is the timelike Im slot.
    The N-gon edge (chord) is 2ρ·sin(π/N) = r_t, so |ΔR_temporal| = r_t uniformly
    on every node (including the vortex core, now the most-protected point), the
    loop closes (v_N = v_0), and the configuration stays codim-0 in ℝ⁴.  This is
    an n_vac=1 background turn in the carrier plane; the carrier (winding n_t)
    rides on top.  Diagnostics must subtract this offset to read the true carrier.

    For r_t <= 0 (the linear/Verlet limit, which has no 1/|ΔR| term) the offset
    is zero — no temporal stretch is needed.

    Returns
    -------
    off_re, off_im : ndarray, shape (n_slices+1,)
        Per-slice offset magnitudes for the Re channel (distributed along the
        trace direction (1,1,1)/sqrt(3) of components 0,1,2) and the Im channel
        (CARRIER_IM = 3), on each slice l = 0..n_slices (slice N duplicates
        slice 0 — closure).
    """
    l = np.arange(n_slices + 1)
    if r_t <= 0.0:
        z = np.zeros(n_slices + 1)
        return z, z
    rho = r_t / (2.0 * math.sin(math.pi / n_slices))
    ang = 2.0 * math.pi * l / n_slices
    return rho * np.cos(ang), rho * np.sin(ang)


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

    # Prestressed periodic-vacuum offset: a per-slice global N-gon translation in
    # the carrier plane so every temporal bond = r_t (no 1/|ΔR| singularity).
    off_re, off_im = vacuum_offsets(n_slices, params.r_t)
    rho = (params.r_t / (2.0 * math.sin(math.pi / n_slices))) if params.r_t > 0 else 0.0

    world = np.empty((n_slices + 1, n_nodes, m_ambient), dtype=np.float64)

    for l_slice in range(n_slices + 1):
        t = float(l_slice)  # carrier phase argument = slice index

        pos = ref.copy()
        # Vacuum N-gon offset (uniform across nodes) — establishes the r_t>0
        # prestressed timelike structure before the carrier rides on top.  The
        # Re offset is distributed along the trace direction (1,1,1)/sqrt(3) so
        # its norm stays off_re (the N-gon chord length r_t is preserved) and the
        # vacuum background remains pure-trace (no spurious SU(3) content).
        pos[:, 0:3] += off_re[l_slice] * _TRACE_WEIGHTS
        pos[:, CARRIER_IM] += off_im[l_slice]

        re_d, im_d = _sph_harm_displacement(
            coords, centre, vp.l, vp.m, vp.A0, vp.r0, vp.w, t, omega_per_slice,
        )
        # Re(Psi) along the trace direction -> pure U(1)/EM (E1 fix); Im(Psi) in
        # the timelike component.
        pos[:, 0:3] += re_d[:, np.newaxis] * _TRACE_WEIGHTS
        pos[:, CARRIER_IM] += im_d

        world[l_slice] = pos

    meta = {
        "ansatz": "u1_spherical_harmonic_vortex",
        "geometry": vp.geometry,
        "carrier_re_components": CARRIER_RE_COMPONENTS,
        "carrier_re_weights": list(CARRIER_RE_WEIGHTS),
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
        "vacuum_ngon_rho": rho,
        "r_t": params.r_t,
        "centre": list(centre),
        "note": (
            f"Spherical-harmonic U(1) vortex seed Y_{vp.l}^{vp.m}. "
            f"Azimuthal winding m={vp.m} around the z-axis (the vortex axis); "
            "energy density ~|Y|^2 forms a donut around that axis. "
            f"Carrier advances n_t={vp.n_t} full turns ({360 * vp.n_t} deg) over "
            "the time loop and closes exactly at the periodic seam "
            "(omega quantized by closure, not free). "
            "Carrier 2-plane: Re(Psi) along the trace direction (1,1,1)/sqrt(3) of "
            "components (0,1,2) -> pure U(1)/EM; component 3 = Im(Psi) (timelike)."
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
    r_t: float = 0.0,
) -> dict[str, float]:
    """Measure the net U(1) azimuthal winding about the three central axes.

    ``r_t`` (the temporal rest length, > 0 for the prestressed substrate) must
    be passed so the per-slice prestressed-vacuum N-gon offset is subtracted
    before reading the phase — otherwise the vacuum background (ρ ≈ r_t/(2 sin
    π/N), which exceeds the carrier amplitude) swamps the contour and the winding
    reads 0 instead of m.  ``r_t=0`` (default) subtracts nothing (linear limit).

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

    ref = lattice.reference_positions(pos.shape[1])

    # Re(Psi) = projection of the lateral displacement onto the unit trace
    # direction (E1 fix); Im(Psi) = the timelike component.
    disp_lat = pos[:, 0:3] - ref[:, 0:3]                     # (n_nodes, 3)
    re_node = project_carrier_re(disp_lat)                   # (n_nodes,)
    im_node = pos[:, CARRIER_IM] - ref[:, CARRIER_IM]        # (n_nodes,)

    # Subtract the prestressed-vacuum N-gon offset for this slice (else the
    # ρ-sized vacuum background swamps the carrier in the contour).  The Re
    # offset projects back to exactly off_re (the trace write preserves its norm).
    n_slices = world.shape[0] - 1
    off_re, off_im = vacuum_offsets(n_slices, r_t)

    re_disp = (re_node - off_re[slice_index]).reshape(nx, ny, nz)
    im_disp = (im_node - off_im[slice_index]).reshape(nx, ny, nz)

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