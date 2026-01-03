"""Experiment report generation utilities."""

from .base import FigureSpec, ReportData, ReportGenerator
from .latex_report import LatexReportGenerator

__all__ = [
    "FigureSpec",
    "ReportData",
    "ReportGenerator",
    "LatexReportGenerator",
]
