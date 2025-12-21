"""
Test Run Manager

Handles creation and management of test run directories with timestamps.
Each test/experiment gets its own timestamped folder for outputs.
"""

from datetime import datetime
from pathlib import Path
from typing import Optional


class TestRunManager:
    """
    Manages test run directories and file organization.

    Creates directory structure:
        test-runs/
            experiment_name_YYYYMMDD_HHMMSS/
                plots/
                data/
                logs/
    """

    def __init__(self, base_dir: str = "test-runs", experiment_name: Optional[str] = None):
        """
        Initialize test run manager.

        Args:
            base_dir: Base directory for all test runs (default: "test-runs")
            experiment_name: Name of the experiment/test (default: auto-detect from script)
        """
        self.base_dir = Path(base_dir)

        # Auto-detect experiment name from calling script if not provided
        if experiment_name is None:
            import inspect
            frame = inspect.currentframe()
            caller_frame = frame.f_back
            caller_file = caller_frame.f_globals.get('__file__', 'unknown')
            experiment_name = Path(caller_file).stem

        self.experiment_name = experiment_name

        # Create timestamped run directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_name = f"{experiment_name}_{timestamp}"
        self.run_dir = self.base_dir / self.run_name

        # Subdirectories
        self.plots_dir = self.run_dir / "plots"
        self.data_dir = self.run_dir / "data"
        self.logs_dir = self.run_dir / "logs"

        # Create directories
        self._create_directories()

    def _create_directories(self):
        """Create all necessary directories."""
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(exist_ok=True)
        self.data_dir.mkdir(exist_ok=True)
        self.logs_dir.mkdir(exist_ok=True)

    def get_plot_path(self, filename: str) -> str:
        """
        Get full path for a plot file.

        Args:
            filename: Name of the plot file (e.g., "amplitude_xy.png")

        Returns:
            Full path to save the plot
        """
        return str(self.plots_dir / filename)

    def get_data_path(self, filename: str) -> str:
        """
        Get full path for a data file.

        Args:
            filename: Name of the data file (e.g., "results.npz")

        Returns:
            Full path to save the data
        """
        return str(self.data_dir / filename)

    def get_log_path(self, filename: str) -> str:
        """
        Get full path for a log file.

        Args:
            filename: Name of the log file (e.g., "experiment.log")

        Returns:
            Full path to save the log
        """
        return str(self.logs_dir / filename)

    def save_config(self, config_dict: dict, filename: str = "config.txt"):
        """
        Save configuration/parameters to a text file.

        Args:
            config_dict: Dictionary of configuration parameters
            filename: Name of config file (default: "config.txt")
        """
        config_path = self.run_dir / filename
        with open(config_path, 'w') as f:
            f.write(f"Test Run: {self.run_name}\n")
            f.write(f"{'=' * 70}\n\n")
            for key, value in config_dict.items():
                f.write(f"{key}: {value}\n")

    def get_summary(self) -> str:
        """Get summary of the test run directory structure."""
        summary = f"\nTest Run Directory Structure:\n"
        summary += f"{'=' * 70}\n"
        summary += f"Run name: {self.run_name}\n"
        summary += f"Base directory: {self.run_dir}\n"
        summary += f"  - Plots: {self.plots_dir}\n"
        summary += f"  - Data: {self.data_dir}\n"
        summary += f"  - Logs: {self.logs_dir}\n"
        return summary

    def __str__(self) -> str:
        """String representation."""
        return f"TestRunManager(run_name='{self.run_name}', run_dir='{self.run_dir}')"

    def __repr__(self) -> str:
        """Representation."""
        return self.__str__()


def get_latest_run(base_dir: str = "test-runs", experiment_name: Optional[str] = None) -> Optional[Path]:
    """
    Get the path to the most recent test run directory.

    Args:
        base_dir: Base directory for all test runs
        experiment_name: Filter by experiment name (optional)

    Returns:
        Path to most recent run directory, or None if not found
    """
    base_path = Path(base_dir)
    if not base_path.exists():
        return None

    # Get all subdirectories
    run_dirs = [d for d in base_path.iterdir() if d.is_dir()]

    # Filter by experiment name if provided
    if experiment_name:
        run_dirs = [d for d in run_dirs if d.name.startswith(experiment_name)]

    if not run_dirs:
        return None

    # Sort by modification time (most recent first)
    run_dirs.sort(key=lambda x: x.stat().st_mtime, reverse=True)

    return run_dirs[0]


def list_runs(base_dir: str = "test-runs", experiment_name: Optional[str] = None) -> list[Path]:
    """
    List all test run directories.

    Args:
        base_dir: Base directory for all test runs
        experiment_name: Filter by experiment name (optional)

    Returns:
        List of run directory paths, sorted by date (most recent first)
    """
    base_path = Path(base_dir)
    if not base_path.exists():
        return []

    # Get all subdirectories
    run_dirs = [d for d in base_path.iterdir() if d.is_dir()]

    # Filter by experiment name if provided
    if experiment_name:
        run_dirs = [d for d in run_dirs if d.name.startswith(experiment_name)]

    # Sort by modification time (most recent first)
    run_dirs.sort(key=lambda x: x.stat().st_mtime, reverse=True)

    return run_dirs