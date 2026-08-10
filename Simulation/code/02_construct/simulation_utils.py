"""Shared estimators and utilities for the cluster-quota RD simulations."""

from __future__ import annotations

import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "simulation.yaml"


def load_config() -> dict[str, Any]:
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def ensure_output_dirs() -> None:
    for relative in [
        "data/processed",
        "output/tables",
        "output/figures",
        "output/logs",
    ]:
        (PROJECT_ROOT / relative).mkdir(parents=True, exist_ok=True)


def triangular_kernel(u: np.ndarray) -> np.ndarray:
    """Triangular kernel K(u)=(1-|u|)1{|u|<1}."""
    return np.clip(1.0 - np.abs(u), 0.0, None)


def midpoint_quota_cutoff(x: np.ndarray, treated_count: int) -> float:
    """Midpoint between the last treated and first untreated order statistics."""
    m = x.size
    if not 1 <= treated_count <= m - 1:
        raise ValueError("treated_count must be interior")
    ordered = np.sort(x)
    first_treated = ordered[m - treated_count]
    first_untreated = ordered[m - treated_count - 1]
    return float((first_treated + first_untreated) / 2.0)


def boundary_outcome(
    x: np.ndarray,
    y0: np.ndarray,
    y1: np.ndarray,
    boundary: float,
) -> np.ndarray:
    """Sharp-assignment outcome generated at a supplied boundary."""
    return np.where(x > boundary, y1, y0)


def _local_linear_side(
    x: np.ndarray,
    y: np.ndarray,
    boundary: float,
    bandwidth: float,
    side: int,
) -> tuple[float, float, int]:
    """One-sided local-linear intercept and HC3 conditional variance."""
    u = (x - boundary) / bandwidth
    if side == 1:
        mask = (u > 0.0) & (u < 1.0)
    elif side == -1:
        mask = (u < 0.0) & (u > -1.0)
    else:
        raise ValueError("side must be +1 or -1")

    u_side = u[mask]
    y_side = y[mask]
    weights = triangular_kernel(u_side)
    positive = weights > 0.0
    u_side = u_side[positive]
    y_side = y_side[positive]
    weights = weights[positive]
    n_side = int(u_side.size)

    if n_side < 4:
        return math.nan, math.nan, n_side

    design = np.column_stack((np.ones(n_side), u_side))
    gram = design.T @ (weights[:, None] * design)
    if not np.all(np.isfinite(gram)) or np.linalg.cond(gram) > 1.0e12:
        return math.nan, math.nan, n_side

    gram_inv = np.linalg.inv(gram)
    beta = gram_inv @ (design.T @ (weights * y_side))
    residual = y_side - design @ beta
    leverage = weights * np.einsum(
        "ij,jk,ik->i", design, gram_inv, design
    )
    leverage = np.clip(leverage, 0.0, 1.0 - 1.0e-8)
    adjusted_squared_residual = (
        residual / (1.0 - leverage)
    ) ** 2
    meat = design.T @ (
        (weights**2 * adjusted_squared_residual)[:, None] * design
    )
    vcov = gram_inv @ meat @ gram_inv
    return float(beta[0]), float(vcov[0, 0]), n_side


def local_linear_rd(
    x: np.ndarray,
    y: np.ndarray,
    boundary: float,
    bandwidth: float,
) -> dict[str, float]:
    """Local-linear RD estimate with side-specific HC3 variance."""
    right, var_right, n_right = _local_linear_side(
        x, y, boundary, bandwidth, side=1
    )
    left, var_left, n_left = _local_linear_side(
        x, y, boundary, bandwidth, side=-1
    )
    estimate = right - left
    variance = var_right + var_left
    se = math.sqrt(variance) if variance >= 0.0 else math.nan
    return {
        "estimate": estimate,
        "se": se,
        "n_right": n_right,
        "n_left": n_left,
    }


def local_window_overlap(
    x: np.ndarray,
    first_boundary: float,
    second_boundary: float,
    bandwidth: float,
) -> float:
    """Jaccard overlap of the two compact local windows."""
    first = np.abs(x - first_boundary) < bandwidth
    second = np.abs(x - second_boundary) < bandwidth
    union = np.count_nonzero(first | second)
    if union == 0:
        return math.nan
    return float(np.count_nonzero(first & second) / union)


def normal_density_at_zero() -> float:
    return 1.0 / math.sqrt(2.0 * math.pi)


def triangular_local_linear_variance(
    error_sd_right: float = 1.0,
    error_sd_left: float = 1.0,
) -> float:
    """Asymptotic variance V(c) for N(0,1) scores and triangular kernel.

    For one-sided local-linear estimation with a triangular kernel,
    the squared equivalent-kernel integral is 24/5 on each side.
    """
    equivalent_kernel_l2 = 24.0 / 5.0
    return (
        equivalent_kernel_l2
        * (error_sd_right**2 + error_sd_left**2)
        / normal_density_at_zero()
    )


def write_dataframe(df: pd.DataFrame, relative_path: str) -> Path:
    path = PROJECT_ROOT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def write_json(payload: dict[str, Any], relative_path: str) -> Path:
    path = PROJECT_ROOT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path


def software_versions() -> dict[str, str]:
    return {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "pyyaml": yaml.__version__,
    }
