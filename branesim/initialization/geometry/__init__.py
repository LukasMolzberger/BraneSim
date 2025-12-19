"""
Geometry utilities for the brane simulation.

This module provides geometric constructions for analyzing brane configurations,
including tubular coordinates around toroidal solitons and photon modes.
"""

from .tubular_electron_geometry import (
    TorusKnotParameters,
    sample_torus_knot_centerline,
    compute_frenet_frames,
    construct_twisted_strip,
)
from .tubular_photon_mode import (
    PhotonModeParameters,
    compute_circular_polarization_EB,
    sample_gaussian_envelope,
)

__all__ = [
    "TorusKnotParameters",
    "sample_torus_knot_centerline",
    "compute_frenet_frames",
    "construct_twisted_strip",
    "PhotonModeParameters",
    "compute_circular_polarization_EB",
    "sample_gaussian_envelope",
]