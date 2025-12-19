"""
Diagnostic pipeline for systematic analysis of simulation data.

Provides a framework for running multiple diagnostics on snapshots,
with caching support to avoid redundant computation (e.g., computing
complex band state once and reusing for multiple Berry diagnostics).

Key abstractions:
- Diagnostic: Protocol for a single diagnostic computation
- DiagnosticPipeline: Orchestrates multiple diagnostics with caching
"""

from __future__ import annotations
from typing import Protocol, Any
from dataclasses import dataclass
import torch

from .types import Snapshot, GridSpec, DiagnosticResult


class Diagnostic(Protocol):
    """
    Protocol for a diagnostic computation.

    Each diagnostic has a name and a compute method that takes a snapshot,
    grid specification, and shared cache, and returns a DiagnosticResult.

    The cache allows diagnostics to share intermediate computations:
    - One diagnostic can compute and cache expensive operations
    - Subsequent diagnostics can read from cache instead of recomputing
    """

    name: str

    def compute(
        self,
        snapshot: Snapshot,
        grid: GridSpec,
        cache: dict[str, Any]
    ) -> DiagnosticResult:
        """
        Compute diagnostic from snapshot.

        Parameters
        ----------
        snapshot : Snapshot
            Snapshot containing field data
        grid : GridSpec
            Grid specification
        cache : dict
            Shared cache for intermediate results
            Diagnostics can read/write to avoid redundant computation

        Returns
        -------
        DiagnosticResult
            Structured diagnostic output
        """
        ...


@dataclass
class AnalyticSignalDiagnostic:
    """
    Diagnostic that computes ω-free analytic signal and caches it.

    This is typically the first diagnostic in a Berry phase pipeline,
    as it computes psi, psi_hat, and amp that other diagnostics need.

    Unlike the old BandStateDiagnostic, this does NOT require omega or velocity.
    It constructs the complex state via positive-frequency projection along
    a chosen spatial axis.
    """
    name: str = "analytic_signal"
    axis: int = 0
    eps: float = 1e-12
    field_key: str = "xi"  # Fallback to "q" if not present

    def compute(
        self,
        snapshot: Snapshot,
        grid: GridSpec,
        cache: dict[str, Any]
    ) -> DiagnosticResult:
        """Compute and cache ω-free analytic signal."""
        from .analytic_signal import (
            analytic_signal_along_axis,
            pointwise_normalize_from_grid,
        )

        # Get position field (no velocity needed)
        q = snapshot.fields.get(self.field_key, None)
        if q is None:
            q = snapshot.fields.get("q", None)
        if q is None:
            raise ValueError(
                f"Snapshot must contain '{self.field_key}' or 'q' field for analytic signal"
            )

        # Compute analytic signal via positive-frequency projection
        psi = analytic_signal_along_axis(q, axis=self.axis, spatial_ndim=grid.D)
        psi_hat, amp = pointwise_normalize_from_grid(psi, grid, eps=self.eps)

        # Cache for other diagnostics
        cache["psi"] = psi
        cache["psi_hat"] = psi_hat
        cache["amp"] = amp

        # Return as diagnostic result
        return DiagnosticResult(
            name=f"{self.name}_axis{self.axis}",
            t_sim=snapshot.t_sim,
            t_phys_s=snapshot.t_phys_s,
            data={
                "psi": psi,
                "psi_hat": psi_hat,
                "amp": amp,
            },
            quality={
                "amp_mean": amp.mean().item(),
                "amp_std": amp.std().item(),
                "amp_max": amp.max().item(),
            },
            meta={
                "axis": self.axis,
                "eps": self.eps,
                "method": "analytic_signal_fft",
            },
        )


@dataclass
class BerryConnectionDiagnostic:
    """
    Diagnostic that computes Berry connection along a specified axis.
    """
    name: str = "berry_connection"
    axis: int = 0
    config: Any = None  # BerryConfig

    def compute(
        self,
        snapshot: Snapshot,
        grid: GridSpec,
        cache: dict[str, Any]
    ) -> DiagnosticResult:
        """Compute Berry connection from cached band state."""
        from .berry import berry_connection_along_axis, BerryConfig

        # Get cached band state (assumes AnalyticSignalDiagnostic ran first)
        if "psi_hat" not in cache or "amp" not in cache:
            raise ValueError(
                "Berry connection requires cached analytic signal. "
                "Run AnalyticSignalDiagnostic first."
            )

        psi_hat = cache["psi_hat"]
        amp = cache["amp"]

        # Use provided config or create default
        if self.config is None:
            cfg = BerryConfig(spacing_sim=grid.spacing_sim)
        else:
            cfg = self.config

        # Compute connection
        result = berry_connection_along_axis(psi_hat, amp, self.axis, cfg)

        # Return as diagnostic result
        return DiagnosticResult(
            name=f"{self.name}_axis{self.axis}",
            t_sim=snapshot.t_sim,
            t_phys_s=snapshot.t_phys_s,
            data={
                "dphi": result["dphi"],
                "A_axis": result["A_axis"],
            },
            quality={
                "overlap_abs": result["overlap_abs"],
                "mask_point": result["mask_point"],
                "valid_edge": result["valid_edge"],
            },
            meta={
                "axis": self.axis,
                "config": cfg,
                "spacing_sim": grid.spacing_sim,
            },
        )


@dataclass
class BerryPhaseDiagnostic:
    """
    Diagnostic that computes integrated Berry phase along an axis.

    This depends on BerryConnectionDiagnostic having run first,
    or it will compute the connection itself.
    """
    name: str = "berry_phase"
    axis: int = 0
    config: Any = None  # BerryConfig

    def compute(
        self,
        snapshot: Snapshot,
        grid: GridSpec,
        cache: dict[str, Any]
    ) -> DiagnosticResult:
        """Compute Berry phase profile."""
        from .berry import (
            berry_connection_along_axis,
            berry_phase_integrated_along_axis,
            BerryConfig,
        )

        # Check if connection was already computed
        conn_key = f"berry_connection_axis{self.axis}_dphi"
        if conn_key in cache:
            dphi = cache[conn_key]
        else:
            # Need to compute connection first
            if "psi_hat" not in cache or "amp" not in cache:
                raise ValueError(
                    "Berry phase requires cached analytic signal. "
                    "Run AnalyticSignalDiagnostic first."
                )

            psi_hat = cache["psi_hat"]
            amp = cache["amp"]

            if self.config is None:
                cfg = BerryConfig(spacing_sim=grid.spacing_sim)
            else:
                cfg = self.config

            conn_result = berry_connection_along_axis(psi_hat, amp, self.axis, cfg)
            dphi = conn_result["dphi"]
            # Cache for potential reuse
            cache[conn_key] = dphi
            cache[f"berry_connection_axis{self.axis}_A"] = conn_result["A_axis"]
            cache[f"berry_connection_axis{self.axis}_valid_edge"] = conn_result["valid_edge"]

        # Use provided config or create default
        if self.config is None:
            cfg = BerryConfig(spacing_sim=grid.spacing_sim)
        else:
            cfg = self.config

        # Integrate to get phase
        phase_result = berry_phase_integrated_along_axis(dphi, self.axis, cfg)

        return DiagnosticResult(
            name=f"{self.name}_axis{self.axis}",
            t_sim=snapshot.t_sim,
            t_phys_s=snapshot.t_phys_s,
            data={
                "gamma_unwrapped": phase_result["gamma_unwrapped"],
                "gamma_wrapped": phase_result["gamma_wrapped"],
            },
            quality={},
            meta={
                "axis": self.axis,
                "config": cfg,
            },
        )


class DiagnosticPipeline:
    """
    Pipeline for running multiple diagnostics with caching.

    Example
    -------
    >>> from branesim.diagnostics.pipeline import (
    ...     DiagnosticPipeline,
    ...     AnalyticSignalDiagnostic,
    ...     BerryConnectionDiagnostic,
    ...     BerryPhaseDiagnostic,
    ... )
    >>>
    >>> # Create pipeline
    >>> pipeline = DiagnosticPipeline([
    ...     AnalyticSignalDiagnostic(axis=0),
    ...     BerryConnectionDiagnostic(axis=0),
    ...     BerryPhaseDiagnostic(axis=0),
    ... ])
    >>>
    >>> # Run on snapshots
    >>> results = pipeline.run(snapshots, grid)
    >>>
    >>> # Access results
    >>> berry_results = results["berry_phase_axis0"]
    >>> gamma = berry_results[0].data["gamma_wrapped"]
    """

    def __init__(self, diagnostics: list[Diagnostic]):
        """
        Initialize pipeline with list of diagnostics.

        Parameters
        ----------
        diagnostics : list[Diagnostic]
            List of diagnostics to run (in order)
        """
        self.diagnostics = diagnostics

    def run(
        self,
        snapshots: list[Snapshot],
        grid: GridSpec,
        verbose: bool = False
    ) -> dict[str, list[DiagnosticResult]]:
        """
        Run all diagnostics on all snapshots.

        Parameters
        ----------
        snapshots : list[Snapshot]
            List of snapshots to analyze
        grid : GridSpec
            Grid specification
        verbose : bool
            If True, print progress

        Returns
        -------
        dict[str, list[DiagnosticResult]]
            Dictionary mapping diagnostic name to list of results
            (one result per snapshot)
        """
        results: dict[str, list[DiagnosticResult]] = {
            diag.name: [] for diag in self.diagnostics
        }

        for i, snapshot in enumerate(snapshots):
            if verbose:
                print(f"Processing snapshot {i+1}/{len(snapshots)}...")

            # Create cache for this snapshot (shared across diagnostics)
            cache: dict[str, Any] = {}

            # Run each diagnostic
            for diag in self.diagnostics:
                if verbose:
                    print(f"  Running {diag.name}...")

                result = diag.compute(snapshot, grid, cache)
                results[diag.name].append(result)

        return results