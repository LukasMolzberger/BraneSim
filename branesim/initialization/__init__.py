"""
Initialization module for narrowband carrier packets.

This module provides preparation-first initialization routines for:
- Photon-like circular polarization packets (2D polarization subspace)
- Electron-like double-loop tube packets (spinorial transport)
"""

from branesim.initialization.carrier_packets import (
    make_photon_circular_packet,
    make_electron_double_loop_packet,
)

__all__ = [
    "make_photon_circular_packet",
    "make_electron_double_loop_packet",
]