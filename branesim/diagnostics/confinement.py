"""Confinement diagnostics for the branesim world-volume format.

Metrics (ported from the validated legacy diagnostics, non-self-referential):

box_fill_radius
    RMS node radius about the geometric centre, uniform weight.
    Frame-invariant: depends only on the reference geometry, computed once.
    This is the value radius_rms approaches if the field disperses to fill the
    box uniformly, so it is the natural fixed reference for confinement ratios.

radius_rms
    Energy-weighted RMS radius about the energy-weighted centre (per slice).
    Energy proxy weight per node: |displacement_lateral|^2, where
    "lateral" = first ``dim`` ambient components (the spatial / gauge sector).
    Using the *energy-weighted* centre rather than the geometric centre makes
    the metric robust when the packet drifts off-centre.

spread_ratio  (per slice)
    radius_rms / box_fill_radius.
    -> 1.0 : field is dispersed to box-fill (deconfined).
    << 1.0 : field is localized/confined.

confined_fraction  (per slice)
    Fraction of total energy weight within confinement_radius_factor *
    box_fill_radius of the energy-weighted centre.  The radius threshold is
    a FIXED box scale (not the packet's own spread), so this metric detects
    dispersal even when the packet has expanded to fill the box.
    The distance is measured from the energy-weighted centre, keeping the
    metric robust to an off-centre packet.

radius_growth  (summary)
    radius_rms(last slice) / radius_rms(first slice).
    > 1 : the mode expanded; combined with spread_ratio -> 1 this signals
    dispersion to box-fill (the opposite of confinement).

Regression note
---------------
The legacy ``leakage_fraction`` (energy fraction OUTSIDE leakage_radius_factor
* radius_rms, i.e. outside a threshold proportional to the packet's own size)
was self-referential: it stays near 0 even when the packet fills the whole box,
because the threshold grows with the packet.  That metric is DEPRECATED and NOT
ported here.

Skyrme-aware confinement note
------------------------------
The default ``|lateral displacement|^2`` weight has a known deficiency for
Skyrme / topological solitons: the lateral component of the hedgehog
``sin(F(r)) * u0 * r_hat`` has algebraic (power-law) tails
``sin(pi/(1+(r/w)^2)) ~ pi*w^2/(r^2)`` which cause the energy-weighted centroid
to integrate over the whole box, yielding ``spread_ratio > 1`` even for a
perfectly localized topological object.

The fix is ``weight_mode="strain"``: weight each node by the sum of
``|displacement_p - displacement_q|^2`` over its incident bonds (the
displacement-direction variation).  For a Skyrme hedgehog the displacement
has constant magnitude ``u0`` everywhere (it lives on S^3), so this sum is
zero only where adjacent nodes have the same displacement vector (the flat
far-field) and is large at the soliton core where the hedgehog field rotates
rapidly.  Pass the lattice neighbor table (``SpacelikeLattice``) to
:func:`confinement_summary` when using this mode; ``k_s`` and ``alpha`` are
accepted for API convenience but are not used in the computation.

Principles compliance
---------------------
- Read-only: these functions take arrays as input and return dicts; they do not
  modify any solver state and carry no back-reaction.
- Dimension-agnostic: ``dim`` is passed explicitly; spatial coordinates are
  ``ref[:, :dim]`` and the energy proxy uses ``displacements[:, :dim]``.
  No hard-coded 3D logic.
- No clamps or artificial cutoffs: all thresholds are ratios of the reference
  geometry, not physics-imposed saturation rules.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from branesim.io.contracts import iter_slices, load_npy, load_manifest


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _box_fill_radius(ref_spatial: np.ndarray) -> float:
    """RMS node radius about the geometric centre (uniform weight).

    Parameters
    ----------
    ref_spatial : ndarray, shape (n_nodes, dim)
        Spatial (spacelike) coordinates of reference positions.

    Returns
    -------
    float
        The uniform-weight RMS radius.  Frame-invariant; depends only on the
        reference geometry.  This is the reference scale for confinement.
    """
    geom_centre = np.mean(ref_spatial, axis=0)          # (dim,)
    delta = ref_spatial - geom_centre[np.newaxis, :]    # (n_nodes, dim)
    r_sq = np.sum(delta * delta, axis=1)                # (n_nodes,)
    return float(np.sqrt(np.mean(r_sq)))


def _energy_weights(displacements: np.ndarray, dim: int) -> np.ndarray:
    """Per-node energy proxy weights: |displacement_lateral|^2.

    Parameters
    ----------
    displacements : ndarray, shape (n_nodes, m_ambient)
        Node displacements (positions - ref).
    dim : int
        Number of spatial (lateral) ambient components to use.

    Returns
    -------
    weights : ndarray, shape (n_nodes,)
        Non-negative energy proxy weights; a small floor avoids division by
        zero without introducing a physics-altering clamp (it is purely a
        numerical regulariser for the weighted-centroid formula when all
        displacements are exactly zero).
    """
    lateral = displacements[:, :dim]                   # (n_nodes, dim)
    return np.sum(lateral * lateral, axis=1) + 1e-40   # (n_nodes,)


def _direction_variation_weights(
    positions: np.ndarray,
    ref: np.ndarray,
    neighbor_table: np.ndarray,
) -> np.ndarray:
    """Per-node displacement-direction variation weight (Skyrme-aware confinement proxy).

    Computes the sum over incident bonds of ``|displacement_p - displacement_q|^2``,
    where ``displacement = position - ref``.  This measures how rapidly the
    displacement *direction* changes from node to node — it is large at the
    soliton core (where the hedgehog field rotates rapidly) and small far from
    the core (where the field is nearly uniform).

    This is the physically correct confinement weight for topological (Skyrme)
    solitons because:

    - The Skyrme hedgehog has ``|displacement|^2 = u0^2`` everywhere (constant
      magnitude on S^3), so ``|displacement|^2`` does not distinguish core from
      far-field.
    - The spring strain energy ``(|bond| - alpha*|ref_bond|)^2`` is nonzero
      everywhere for large-amplitude displacements, giving spread_ratio ~ 1 even
      for a perfectly localized soliton.
    - The displacement difference ``|disp_p - disp_q|^2`` is zero only when
      adjacent nodes have the same displacement vector, which happens only in the
      flat far-field; it peaks at the soliton core where the hedgehog field
      rotates most rapidly.

    For a scalar (small-amplitude) breather, ``|displacement_lateral|^2`` is the
    better proxy (available via ``weight_mode="lateral"``).

    Parameters
    ----------
    positions : ndarray, shape (n_nodes, m_ambient)
    ref : ndarray, shape (n_nodes, m_ambient)
    neighbor_table : ndarray, shape (n_nodes, n_neighbors), dtype int64
        Flat neighbor indices; -1 means "no neighbor" (open boundary).
        Typically ``lattice.neighbors``.

    Returns
    -------
    weights : ndarray, shape (n_nodes,)
        Per-node direction-variation weight; always >= 0; floor 1e-40 for
        numerical safety.
    """
    disp = positions - ref                             # (n_nodes, m_ambient)
    n_nodes, n_nb = neighbor_table.shape
    weights = np.zeros(n_nodes, dtype=np.float64)

    for nb_idx in range(n_nb):
        q_arr = neighbor_table[:, nb_idx]             # (n_nodes,)
        valid = q_arr >= 0
        p_idx = np.where(valid)[0]
        q_idx = q_arr[p_idx]

        diff = disp[p_idx] - disp[q_idx]             # (n_valid, m_ambient)
        weights[p_idx] += np.sum(diff * diff, axis=1) # (n_valid,)

    return weights + 1e-40


# ---------------------------------------------------------------------------
# Per-slice API
# ---------------------------------------------------------------------------


def confinement_metrics_per_slice(
    positions: np.ndarray,
    ref: np.ndarray,
    dim: int,
    confinement_radius_factor: float = 0.5,
    _box_fill_radius_cached: float | None = None,
    weight_mode: str = "lateral",
    _neighbor_table: np.ndarray | None = None,
    _k_s: float = 1.0,
    _alpha: float = 0.5,
) -> dict[str, float]:
    """Compute confinement metrics for a single spacelike slice.

    Parameters
    ----------
    positions : ndarray, shape (n_nodes, m_ambient)
        Node positions at this slice.
    ref : ndarray, shape (n_nodes, m_ambient)
        Reference (unstressed) node positions.
    dim : int
        Number of spatial (lateral / gauge) ambient components.
        Spatial coordinates are ``ref[:, :dim]``; energy proxy uses
        ``displacements[:, :dim]``.  Do NOT hard-code 3; pass ``lattice.dim``.
    confinement_radius_factor : float, optional
        Fraction of box_fill_radius defining the confinement sphere.
        Default 0.5.
    _box_fill_radius_cached : float or None, optional
        Pre-computed box_fill_radius (avoids redundant computation when
        called in a loop over slices).  If None, it is computed internally.
    weight_mode : str, optional
        ``"lateral"`` (default): weight by ``|displacement_lateral|^2``.
        Suitable for scalar (transverse) breathers and wavepackets.

        ``"strain"``: weight by per-node displacement-direction variation,
        ``sum_{bonds} |disp_p - disp_q|^2``.  Suitable for Skyrme /
        topological solitons whose displacement has constant magnitude on S^3
        (``|disp|^2 = u0^2`` everywhere) and whose lateral component has
        algebraic tails that inflate the ``"lateral"``-weighted spread_ratio
        above 1 even for a localized topological object.  Only the
        ``_neighbor_table`` argument is required; ``_k_s`` and ``_alpha``
        are accepted but ignored.
    _neighbor_table : ndarray or None
        Flat neighbor index table from ``lattice.neighbors``.  Required when
        ``weight_mode="strain"``.
    _k_s : float
        Accepted for API compatibility; not used in the direction-variation
        computation.
    _alpha : float
        Accepted for API compatibility; not used in the direction-variation
        computation.

    Returns
    -------
    dict with keys:
        radius_rms : float
            Energy-weighted RMS radius about the energy-weighted centre.
        box_fill_radius : float
            Uniform-weight RMS radius of the reference geometry.
        spread_ratio : float
            radius_rms / box_fill_radius.  -> 1 dispersed; << 1 localized.
        confined_fraction : float
            Fraction of energy weight within confinement_radius_factor *
            box_fill_radius of the energy-weighted centre.
    """
    ref_spatial = ref[:, :dim]                          # (n_nodes, dim)

    if _box_fill_radius_cached is None:
        bfr = _box_fill_radius(ref_spatial)
    else:
        bfr = float(_box_fill_radius_cached)

    if weight_mode == "strain":
        if _neighbor_table is None:
            raise ValueError(
                "weight_mode='strain' requires _neighbor_table (pass lattice.neighbors)"
            )
        weights = _direction_variation_weights(positions, ref, _neighbor_table)
    elif weight_mode == "lateral":
        displacements = positions - ref                 # (n_nodes, m_ambient)
        weights = _energy_weights(displacements, dim)   # (n_nodes,)
    else:
        raise ValueError(
            f"weight_mode must be 'lateral' or 'strain'; got {weight_mode!r}"
        )

    total_weight = float(np.sum(weights))

    # Energy-weighted centre (robust to off-centre packet)
    centre = np.sum(ref_spatial * weights[:, np.newaxis], axis=0) / total_weight

    # Radii from the energy-weighted centre
    delta = ref_spatial - centre[np.newaxis, :]        # (n_nodes, dim)
    r_sq = np.sum(delta * delta, axis=1)               # (n_nodes,)
    radii = np.sqrt(r_sq)                              # (n_nodes,)

    radius_rms = float(
        np.sqrt(np.sum(weights * r_sq) / total_weight)
    )

    # spread_ratio: how close to box-fill?
    if bfr > 1e-30:
        spread_ratio = radius_rms / bfr
    else:
        spread_ratio = float("nan")

    # confined_fraction: energy within a FIXED box-scale threshold
    confinement_radius = confinement_radius_factor * bfr
    mask = radii <= confinement_radius
    confined_fraction = float(np.sum(weights[mask]) / total_weight)

    return {
        "radius_rms": radius_rms,
        "box_fill_radius": bfr,
        "spread_ratio": spread_ratio,
        "confined_fraction": confined_fraction,
    }


# ---------------------------------------------------------------------------
# Whole world-volume API
# ---------------------------------------------------------------------------


def confinement_summary(
    slices: np.ndarray,
    ref: np.ndarray,
    dim: int,
    confinement_radius_factor: float = 0.5,
    weight_mode: str = "lateral",
    lattice: Any = None,
    k_s: float = 1.0,
    alpha: float = 0.5,
) -> dict[str, Any]:
    """Compute confinement metrics over all slices of a world-volume.

    Parameters
    ----------
    slices : ndarray, shape (n_slices, n_nodes, m_ambient)
        All spacelike slices (e.g. ``world_volume.slices`` from a block solve).
        First axis is the slice index.
    ref : ndarray, shape (n_nodes, m_ambient)
        Reference (unstressed) node positions.
    dim : int
        Number of spatial (lateral) ambient components.
    confinement_radius_factor : float, optional
        Fraction of box_fill_radius defining the confinement sphere.
    weight_mode : str, optional
        ``"lateral"`` (default): weight by ``|displacement_lateral|^2``.
        Suitable for scalar breathers and wavepackets.

        ``"strain"``: weight by per-node spacelike spring strain energy.
        Suitable for Skyrme / topological solitons whose lateral displacement
        has algebraic tails that mislead the ``"lateral"`` metric.
        Pass ``lattice``, ``k_s``, and ``alpha`` when using this mode.
    lattice : SpacelikeLattice or None
        Required when ``weight_mode="strain"``; provides ``lattice.neighbors``.
    k_s : float
        Spring constant (used when ``weight_mode="strain"``).
    alpha : float
        Rest-length fraction (used when ``weight_mode="strain"``).

    Returns
    -------
    dict with keys:
        box_fill_radius : float
            Frame-invariant reference scale (uniform-weight RMS of ref geometry).
        radius_rms : ndarray, shape (n_slices,)
            Energy-weighted RMS radius per slice.
        spread_ratio : ndarray, shape (n_slices,)
            spread_ratio per slice.
        confined_fraction : ndarray, shape (n_slices,)
            confined_fraction per slice.
        radius_growth : float
            radius_rms[-1] / radius_rms[0].  > 1 signals expansion.
        final : dict
            Metrics from the last slice.
        mean : dict
            Mean of per-slice metric arrays.
        weight_mode : str
            The weight mode used (echoed for traceability).
    """
    n_slices = slices.shape[0]
    ref_spatial = ref[:, :dim]
    bfr = _box_fill_radius(ref_spatial)

    neighbor_table = None
    if weight_mode == "strain":
        if lattice is None:
            raise ValueError(
                "weight_mode='strain' requires lattice (pass a SpacelikeLattice)"
            )
        neighbor_table = lattice.neighbors

    radius_rms_arr = np.empty(n_slices)
    spread_ratio_arr = np.empty(n_slices)
    confined_fraction_arr = np.empty(n_slices)

    for i in range(n_slices):
        m = confinement_metrics_per_slice(
            slices[i], ref, dim,
            confinement_radius_factor=confinement_radius_factor,
            _box_fill_radius_cached=bfr,
            weight_mode=weight_mode,
            _neighbor_table=neighbor_table,
            _k_s=k_s,
            _alpha=alpha,
        )
        radius_rms_arr[i] = m["radius_rms"]
        spread_ratio_arr[i] = m["spread_ratio"]
        confined_fraction_arr[i] = m["confined_fraction"]

    radius_growth = float(radius_rms_arr[-1] / radius_rms_arr[0]) if radius_rms_arr[0] > 1e-40 else float("nan")

    final = {
        "radius_rms": float(radius_rms_arr[-1]),
        "spread_ratio": float(spread_ratio_arr[-1]),
        "confined_fraction": float(confined_fraction_arr[-1]),
    }
    mean = {
        "radius_rms": float(np.mean(radius_rms_arr)),
        "spread_ratio": float(np.mean(spread_ratio_arr)),
        "confined_fraction": float(np.mean(confined_fraction_arr)),
    }

    return {
        "box_fill_radius": bfr,
        "radius_rms": radius_rms_arr,
        "spread_ratio": spread_ratio_arr,
        "confined_fraction": confined_fraction_arr,
        "radius_growth": radius_growth,
        "final": final,
        "mean": mean,
        "weight_mode": weight_mode,
    }


# ---------------------------------------------------------------------------
# World-volume zip wrapper
# ---------------------------------------------------------------------------


def confinement_from_worldvolume(
    path: str | Path,
    confinement_radius_factor: float = 0.5,
) -> dict[str, Any]:
    """Compute confinement metrics from a branesim worldvolume.zip file.

    Reads slices via ``branesim.io.contracts.iter_slices`` and the reference
    positions from ``aux/ref_positions.npy``.  The spatial dimension ``dim``
    is taken from the ``lattice.dim`` field in the manifest.

    Parameters
    ----------
    path : path-like
        Path to a ``worldvolume.zip`` written by ``WorldVolumeWriter``.
    confinement_radius_factor : float, optional
        Fraction of box_fill_radius for the confinement sphere threshold.

    Returns
    -------
    dict
        Same structure as :func:`confinement_summary`.
    """
    path = Path(path)
    manifest = load_manifest(path)
    dim = int(manifest["lattice"]["dim"])
    ref = load_npy(path, "aux/ref_positions.npy")

    all_positions = []
    for _index, _time, positions in iter_slices(path):
        all_positions.append(positions)

    if not all_positions:
        raise ValueError(f"No slices found in {path}")

    slices = np.stack(all_positions, axis=0)  # (n_slices, n_nodes, m_ambient)

    return confinement_summary(
        slices, ref, dim,
        confinement_radius_factor=confinement_radius_factor,
    )
