"""
Geometry utilities for the brane simulation.

This module provides geometric constructions for analyzing brane configurations,
including tubular coordinates around toroidal solitons.
"""

from .tubular_electron_geometry import (
    TorusKnotParameters,
    sample_torus_knot_centerline,
    compute_frenet_frames,
    construct_twisted_strip,
)

__all__ = [
    "TorusKnotParameters",
    "sample_torus_knot_centerline",
    "compute_frenet_frames",
    "construct_twisted_strip",
]