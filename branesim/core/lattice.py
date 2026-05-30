"""Spacelike axial stencil and neighbor lookup — dimension-agnostic.

The canonical stencil (ARCHITECTURE.md §3.1, backbone #15) is the
``2*dim``-neighbor axial-only set:
    N_s = { +/-e_i  |  i in 0..dim-1 }
Diagonal bonds are intentionally absent; this keeps the spatial dynamical
matrix D(k) diagonal in the Cartesian basis at every k and alpha.

The ``l`` (timelike) index is managed separately by the solver; this module
handles only the spacelike topology.

No torch dependency; uses numpy only so the module is usable in any context.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from branesim.core.conventions import LatticeParams


# ---------------------------------------------------------------------------
# Axial offset enumeration (dimension-agnostic)
# ---------------------------------------------------------------------------


def axial_offsets(dim: int) -> np.ndarray:
    """Return all 2*dim unit axial offsets for a ``dim``-dimensional lattice.

    Result shape: (2*dim, dim).  Order: -e_0, +e_0, -e_1, +e_1, …

    Examples
    --------
    >>> axial_offsets(1)
    array([[-1], [ 1]])
    >>> axial_offsets(2)
    array([[-1,  0], [ 1,  0], [ 0, -1], [ 0,  1]])
    >>> axial_offsets(3)
    array([[-1, 0, 0], [1, 0, 0], [0,-1, 0], [0, 1, 0], [0, 0,-1], [0, 0, 1]])
    """
    offsets = []
    for axis in range(dim):
        for sign in (-1, +1):
            off = [0] * dim
            off[axis] = sign
            offsets.append(off)
    return np.array(offsets, dtype=np.int64)  # (2*dim, dim)


# ---------------------------------------------------------------------------
# Node-index helpers
# ---------------------------------------------------------------------------


def linear_index(multi_idx: np.ndarray, grid_shape: tuple[int, ...]) -> np.ndarray:
    """Convert multi-index array (..., dim) to flat index (...) using C-order.

    No bounds checking — callers are responsible for valid indices.
    """
    strides = np.ones(len(grid_shape), dtype=np.int64)
    for axis in range(len(grid_shape) - 2, -1, -1):
        strides[axis] = strides[axis + 1] * grid_shape[axis + 1]
    return (multi_idx * strides).sum(axis=-1)


def multi_index(flat: np.ndarray, grid_shape: tuple[int, ...]) -> np.ndarray:
    """Convert flat indices to multi-index array (n_nodes, dim)."""
    dim = len(grid_shape)
    result = np.empty((len(flat), dim), dtype=np.int64)
    remainder = np.asarray(flat, dtype=np.int64).copy()
    for axis in range(dim - 1, -1, -1):
        result[:, axis] = remainder % grid_shape[axis]
        remainder //= grid_shape[axis]
    return result


# ---------------------------------------------------------------------------
# Neighbor table
# ---------------------------------------------------------------------------


@dataclass
class SpacelikeLattice:
    """Precomputed neighbor table for a d-dimensional periodic/open lattice.

    Attributes
    ----------
    params : LatticeParams
    offsets : ndarray, shape (2*dim, dim)
        The 2*dim unit axial offset vectors.
    neighbors : ndarray, shape (n_nodes, 2*dim)
        Flat neighbor indices; -1 means "no neighbor" (open boundary).
    """

    params: LatticeParams
    offsets: np.ndarray = field(init=False)
    neighbors: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.offsets = axial_offsets(self.params.dim)  # (2*dim, dim)
        self.neighbors = self._build_neighbors()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def dim(self) -> int:
        return self.params.dim

    @property
    def n_nodes(self) -> int:
        return self.params.n_nodes

    @property
    def n_neighbors(self) -> int:
        return 2 * self.dim

    @cached_property
    def multi_indices(self) -> np.ndarray:
        """Multi-index for every node, shape (n_nodes, dim)."""
        flat = np.arange(self.n_nodes, dtype=np.int64)
        return multi_index(flat, self.params.grid_shape)

    # ------------------------------------------------------------------
    # Reference (unstressed) positions
    # ------------------------------------------------------------------

    def reference_positions(self, m_ambient: int | None = None) -> np.ndarray:
        """Flat-lattice reference position array, shape (n_nodes, m_ambient).

        The first ``dim`` components are ``a * (i0, i1, …, i_{d-1})``.
        The remaining ``m_ambient - dim`` components are zero (the "timelike"
        and any higher ambient directions are initialised to zero; temporal
        spacing is handled by the solver).

        Parameters
        ----------
        m_ambient : int, optional
            Number of ambient components.  Defaults to ``dim + 1``.
        """
        if m_ambient is None:
            m_ambient = self.dim + 1
        if m_ambient < self.dim:
            raise ValueError(
                f"m_ambient ({m_ambient}) must be >= dim ({self.dim})"
            )
        a = self.params.spacing
        ref = np.zeros((self.n_nodes, m_ambient), dtype=np.float64)
        mi = self.multi_indices  # (n_nodes, dim)
        ref[:, : self.dim] = mi * a
        return ref

    # ------------------------------------------------------------------
    # Neighbor-table construction (pure numpy, dimension-agnostic)
    # ------------------------------------------------------------------

    def _build_neighbors(self) -> np.ndarray:
        """Build the (n_nodes, 2*dim) neighbor index table.

        Uses vectorised numpy for efficiency — no Python loop over nodes.
        """
        grid_shape = self.params.grid_shape
        periodic = self.params.periodic_axes
        dim = self.dim
        n_nodes = self.n_nodes

        # Multi-index for all nodes: shape (n_nodes, dim)
        mi = multi_index(np.arange(n_nodes, dtype=np.int64), grid_shape)

        neighbors = -np.ones((n_nodes, 2 * dim), dtype=np.int64)

        for nb_idx, offset in enumerate(self.offsets):  # 2*dim iterations
            # Compute neighbour multi-index
            nb_mi = mi + offset[np.newaxis, :]  # (n_nodes, dim)

            # Validity mask per axis
            valid = np.ones(n_nodes, dtype=bool)
            for axis in range(dim):
                n_axis = grid_shape[axis]
                if periodic[axis]:
                    nb_mi[:, axis] = nb_mi[:, axis] % n_axis
                else:
                    in_bounds = (nb_mi[:, axis] >= 0) & (nb_mi[:, axis] < n_axis)
                    valid &= in_bounds

            # Convert valid entries to flat indices
            flat_nb = linear_index(nb_mi, grid_shape)  # (n_nodes,)
            neighbors[valid, nb_idx] = flat_nb[valid]

        return neighbors
