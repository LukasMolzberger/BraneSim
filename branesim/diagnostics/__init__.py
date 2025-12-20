"""
Dimension-agnostic diagnostic tools for BraneSim.

This module provides:
- Data types: GridSpec, Snapshot, DiagnosticResult
- Time-domain Berry phase diagnostics (analytic signal + Berry connection/phase)
- Spectrum analysis
- Energy diagnostics
- Holonomy and degeneracy verification
- Symplectic band structure and Berry phase diagnostics (in bands subpackage)
"""

# Core data types
from branesim.diagnostics.types import (
    Axis,
    GridSpec,
    Snapshot,
    DiagnosticResult,
)




__all__ = [
    # Data types
    "Axis",
    "GridSpec",
    "Snapshot",
    "DiagnosticResult",
]