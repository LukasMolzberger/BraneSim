from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass
class TorusKnotParameters:
    """
    Parameters for the (p, q) torus-knot centerline C(z).

    For the electron model we use a (2, 1) torus knot, where
    `core_windings = 2` (around the torus core) and
    `tube_windings = 1` (around the minor circle).
    """
    major_radius: float = 1.0   # R
    minor_radius: float = 0.3   # r_0
    core_windings: int = 2      # p in (p, q)
    tube_windings: int = 1      # q in (p, q)


def sample_torus_knot_centerline(
    params: TorusKnotParameters,
    num_samples: int = 2000,
) -> np.ndarray:
    """
    Sample a closed (p, q) torus-knot centerline C(z) ⊂ R^3.

    This is the outer geometry of the toroidal electron in the
    brane's lateral coordinates (X^1, X^2, X^3). It corresponds
    to the centerline C(z) described in the tubular-coordinate
    appendix, arranged as a (2, 1) torus knot.

    Parameters
    ----------
    params : TorusKnotParameters
        Geometry of the torus and the (p, q) winding numbers.
    num_samples : int
        Number of samples along the closed curve.

    Returns
    -------
    points : (num_samples, 3) np.ndarray
        Discrete samples of C(z) in brane coordinates (X^1, X^2, X^3).
    """
    R = params.major_radius
    r0 = params.minor_radius
    p = params.core_windings
    q = params.tube_windings

    # Parameter u runs from 0 to 2π; the (p, q) torus knot closes once.
    u = np.linspace(0.0, 2.0 * np.pi, num_samples, endpoint=False)

    # φ: angle around the torus core (major circle)
    # θ: angle around the minor circle
    phi = p * u
    theta = q * u

    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Standard torus embedding:
    #   C(u) = ((R + r0 cos θ) cos φ,
    #           (R + r0 cos θ) sin φ,
    #           r0 sin θ)
    x = (R + r0 * cos_theta) * cos_phi
    y = (R + r0 * cos_theta) * sin_phi
    z = r0 * sin_theta

    points = np.stack([x, y, z], axis=-1)
    return points


def _normalize(v: np.ndarray, axis: int = -1, eps: float = 1e-9) -> np.ndarray:
    """Normalize vectors along a given axis."""
    norm = np.linalg.norm(v, axis=axis, keepdims=True)
    norm = np.maximum(norm, eps)
    return v / norm


def compute_frenet_frames(
    centerline: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute discrete Frenet–Serret frame (t, n, b) along a closed centerline.

    This is the discrete analogue of the continuum Frenet–Serret frame
    (t, n, b) with curvature κ and torsion τ used in the tubular
    embedding r(z, x, y) = C(z) + x n(z) + y b(z).

    Parameters
    ----------
    centerline : (N, 3) np.ndarray
        Discrete samples of C(z).

    Returns
    -------
    t : (N, 3) np.ndarray
        Unit tangent vectors.
    n : (N, 3) np.ndarray
        Unit normal vectors.
    b : (N, 3) np.ndarray
        Unit binormal vectors.
    """
    pts = np.asarray(centerline)
    N = pts.shape[0]

    # Periodic central differences for the tangent t
    forward = np.roll(pts, -1, axis=0)
    backward = np.roll(pts, 1, axis=0)
    tangents = _normalize(forward - backward, axis=-1)

    # Central difference of tangents for the normal n
    t_forward = np.roll(tangents, -1, axis=0)
    t_backward = np.roll(tangents, 1, axis=0)
    dn = t_forward - t_backward

    normals = _normalize(dn, axis=-1)

    # Binormal b = t × n
    binormals = _normalize(np.cross(tangents, normals), axis=-1)

    # Re-orthogonalize n to guarantee n ⟂ t numerically
    normals = _normalize(
        normals - np.sum(normals * tangents, axis=-1, keepdims=True) * tangents,
        axis=-1,
    )

    return tangents, normals, binormals


def construct_twisted_strip(
    centerline: np.ndarray,
    normals: np.ndarray,
    binormals: np.ndarray,
    strip_half_width: float = 0.1,
    num_width_samples: int = 16,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct a narrow twisted strip around a given centerline.

    This is the *physical strip* that supports the internal mode, as
    opposed to a full tube. It implements the embedding

        r(z, w) = C(z) + w * e_strip(z),

    which is the special case of the tubular map

        r(z, x, y) = C(z) + x n(z) + y b(z)

    with y = 0, x = w. Here we choose e_strip(z) = n(z), so the
    strip lies in the osculating plane of the curve.

    Parameters
    ----------
    centerline : (N, 3) np.ndarray
        Samples of C(z).
    normals, binormals : (N, 3) np.ndarray
        Frenet–Serret frame along the centerline.
        (binormals are kept for completeness / future extensions.)
    strip_half_width : float
        Half-width of the strip in the local w-coordinate.
    num_width_samples : int
        Number of samples across the strip width.

    Returns
    -------
    X, Y, Z : (N, num_width_samples) np.ndarray
        Grids of 3D coordinates describing the strip surface.
    """
    pts = np.asarray(centerline)
    n = np.asarray(normals)

    N = pts.shape[0]
    w = np.linspace(-strip_half_width, strip_half_width, num_width_samples)

    # Broadcast: (N, 1, 3) + (N, 1, 3) * (1, M, 1)
    pts_expanded = pts[:, np.newaxis, :]          # (N, 1, 3)
    n_expanded = n[:, np.newaxis, :]              # (N, 1, 3)
    w_expanded = w[np.newaxis, :, np.newaxis]     # (1, M, 1)

    strip_points = pts_expanded + w_expanded * n_expanded  # (N, M, 3)

    X = strip_points[..., 0]
    Y = strip_points[..., 1]
    Z = strip_points[..., 2]

    return X, Y, Z