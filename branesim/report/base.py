"""Abstract report definitions for experiments."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Protocol, Union


Scalar = Union[str, int, float]


@dataclass
class FigureSpec:
    """Figure metadata for report inclusion."""

    path: str
    caption: str


@dataclass
class ReportData:
    """Structured data passed to report generators."""

    title: str
    experiment_name: str
    run_name: str
    summary: Optional[str] = None
    parameters: Dict[str, Scalar] = field(default_factory=dict)
    choices: List[str] = field(default_factory=list)
    assumptions: List[str] = field(default_factory=list)
    derived: Dict[str, Scalar] = field(default_factory=dict)
    dictionary: Dict[str, str] = field(default_factory=dict)
    paper_mapping: Dict[str, str] = field(default_factory=dict)
    figures: List[FigureSpec] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    metadata: Dict[str, Scalar] = field(default_factory=dict)


class ReportGenerator(Protocol):
    """Protocol for report generator implementations."""

    def generate(self, report: ReportData, output_path: str) -> None:
        """Generate a report from structured data."""
        raise NotImplementedError
