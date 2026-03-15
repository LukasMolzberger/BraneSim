"""Common shared helpers."""

from __future__ import annotations

import torch


def choose_device(device: str) -> torch.device:
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    return torch.device(device)


def choose_dtype(dtype_name: str, device: torch.device) -> torch.dtype:
    if dtype_name not in {"float32", "float64"}:
        raise ValueError("dtype must be one of: float32, float64")
    dtype = torch.float32 if dtype_name == "float32" else torch.float64
    if device.type == "mps" and dtype == torch.float64:
        return torch.float32
    return dtype


def parse_bool_triple(text: str) -> tuple[bool, bool, bool]:
    values = tuple(v.strip().lower() == "true" for v in text.split(","))
    if len(values) != 3:
        raise ValueError("Expected three comma-separated booleans")
    return values
