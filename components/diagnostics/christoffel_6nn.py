"""Closed-form Christoffel matrix for the 6-neighbor axial-only lattice.

The canonical substrate (``paper/backbone.md`` #15) is the cubic lattice with
only the six axial bonds ``±êᵢ`` at each node. For this stencil the dynamical
matrix admits a closed form that is *diagonal in the Cartesian basis at every
``k`` and every ``α``*. This module exposes that closed form so other
diagnostics can use it without re-deriving.

Starting from ``paper/derivations/lattice_to_continuum.md`` §3.1 with the link
energy ``(s_δ²)⁽²⁾ = α (δ̂·Δξ)² + (1−α) |Δξ|²`` and restricting the bond set to
``δ ∈ {±êₓ, ±ê_y, ±ê_z}`` (all axial, ``|δ|=1``), the per-bond contribution to
``D_ac(k)`` reduces because ``δ̂_a δ̂_c = δ_{a,i}δ_{c,i}`` for ``δ̂ = ±êᵢ``. The
``+δ`` and ``−δ`` contributions are identical (``1−cos`` is even), giving

    D_ac(k) · ρ / k₀ = 2 Σᵢ hᵢ · [ α · δ_{a,i} δ_{c,i}  +  (1 − α) · δ_ac ]

with ``hᵢ ≡ 1 − cos(kᵢ a)``. The off-diagonal entries (``a ≠ c``) vanish because
``δ_{a,i} δ_{c,i}`` is zero for distinct ``a,c`` and the second piece is the
Kronecker delta. The three diagonal eigenvalues are

    ω²ₐ(k) = (2 k₀ / ρ) · [ α hₐ  +  (1 − α) (hₓ + h_y + h_z) ],   a ∈ {x,y,z}

with eigenvectors ``êₐ`` independent of ``k`` and ``α``.

Implications (used elsewhere):

- The eigenframe ``E(k) = I`` is constant in ``k``. Any k-space Berry / WZ
  connection built from this frame is identically zero (paper §5.6 caveat;
  backbone #16).
- At ``α = 0`` the bracket collapses to ``hₓ + h_y + h_z`` for every ``a``, so
  ``D(k) ∝ I`` — full 3-fold eigenvalue degeneracy at every ``k``. The
  Wilczek–Zee gauge group on the lateral triplet is ``U(3)``.
- At ``α = 1`` the bracket reduces to ``hₐ`` — pure Hookean longitudinal.
  Transverse modes (``a ≠ k̂``) are at zero frequency. ``U(1)³`` in the
  strongest sense.
- The current default ``α = 0.2`` sits close to the U(3) end: the relative
  eigenvalue spread at e.g. ``k`` along [100] is ``1 − √(1−α) ≈ 10.6 %``.
"""

from __future__ import annotations

import numpy as np


def christoffel_6nn(k, alpha, k0: float = 1.0, rho: float = 1.0, a: float = 1.0):
    """Return the 3x3 dynamical matrix ``D(k)`` for the 6-neighbor axial-only lattice.

    Parameters
    ----------
    k : array_like, shape (3,)
        Wavevector in inverse length units.
    alpha : float
        Prestretch ``α = rest_length / spacing`` (``α=1`` is no prestress).
    k0, rho, a : float
        Per-bond spring constant, mass density, lattice spacing.

    Returns
    -------
    D : ndarray, shape (3, 3)
        Real diagonal matrix with entries ``ω²ₐ`` given above.
    """
    k = np.asarray(k, dtype=float).reshape(3)
    h = 1.0 - np.cos(k * a)  # shape (3,)
    h_sum = float(h.sum())
    prefactor = 2.0 * k0 / rho
    diag = prefactor * (alpha * h + (1.0 - alpha) * h_sum)
    return np.diag(diag)


def eigvals_6nn(k, alpha, k0: float = 1.0, rho: float = 1.0, a: float = 1.0):
    """Return the three eigenvalues ``ω²ₐ(k)`` as a length-3 array in axis order (x,y,z).

    Eigenvectors are the Cartesian unit vectors at every ``k``; we therefore
    do not sort here, so the caller can identify which axis each eigenvalue
    belongs to.
    """
    k = np.asarray(k, dtype=float).reshape(3)
    h = 1.0 - np.cos(k * a)
    h_sum = float(h.sum())
    return (2.0 * k0 / rho) * (alpha * h + (1.0 - alpha) * h_sum)