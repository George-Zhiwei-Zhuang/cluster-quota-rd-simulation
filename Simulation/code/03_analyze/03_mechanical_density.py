"""Simulation of mechanical density features after quota recentering."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "code" / "02_construct"))

from simulation_utils import (  # noqa: E402
    ensure_output_dirs,
    load_config,
    write_dataframe,
    write_json,
)


def main() -> None:
    started = time.time()
    ensure_output_dirs()
    config = load_config()
    settings = config["mechanical_density"]
    master_seed = int(config["project"]["master_seed"]) + 20000
    rng = np.random.default_rng(master_seed)

    clusters = int(settings["clusters"])
    m = int(settings["cluster_size"])
    k = int(settings["treated_count"])
    fixed_cutoff = float(settings["fixed_cutoff"])

    scores = rng.uniform(size=(clusters, m))
    ordered = np.sort(scores, axis=1)
    first_treated = ordered[:, m - k]
    first_untreated = ordered[:, m - k - 1]
    midpoint = (first_treated + first_untreated) / 2.0
    recentered = scores - midpoint[:, None]
    fixed_centered = scores - fixed_cutoff

    bin_min = float(settings["histogram_min"])
    bin_max = float(settings["histogram_max"])
    bin_width = float(settings["histogram_bin_width"])
    bin_edges = np.arange(
        bin_min, bin_max + bin_width * 1.0001, bin_width
    )
    fixed_counts, _ = np.histogram(
        fixed_centered.ravel(), bins=bin_edges
    )
    recentered_counts, _ = np.histogram(
        recentered.ravel(), bins=bin_edges
    )
    denominator = clusters * m * bin_width
    histogram = pd.DataFrame(
        {
            "bin_left": bin_edges[:-1],
            "bin_right": bin_edges[1:],
            "bin_midpoint": (bin_edges[:-1] + bin_edges[1:]) / 2.0,
            "fixed_cutoff_density": fixed_counts / denominator,
            "recentered_density": recentered_counts / denominator,
        }
    )
    write_dataframe(
        histogram,
        "data/processed/density_histogram.csv",
    )

    local_rows: list[dict[str, float]] = []
    for t in settings["local_windows"]:
        t = float(t)
        local_rows.append(
            {
                "window_t": t,
                "fixed_right_density": np.mean(
                    (fixed_centered > 0.0) & (fixed_centered <= t)
                )
                / t,
                "fixed_left_density": np.mean(
                    (fixed_centered >= -t) & (fixed_centered < 0.0)
                )
                / t,
                "recentered_right_density": np.mean(
                    (recentered > 0.0) & (recentered <= t)
                )
                / t,
                "recentered_left_density": np.mean(
                    (recentered >= -t) & (recentered < 0.0)
                )
                / t,
                "primitive_theoretical_density": 1.0,
                "recentered_theoretical_limit": 2.0,
            }
        )
    write_dataframe(
        pd.DataFrame.from_records(local_rows),
        "data/processed/density_local_limits.csv",
    )

    rank_grid = pd.DataFrame(
        {
            "rank": np.arange(1, m + 1),
            "normalized_rank": np.arange(1, m + 1) / m,
            "within_cluster_mass": np.repeat(1.0 / m, m),
        }
    )
    write_dataframe(
        rank_grid,
        "data/processed/normalized_rank_grid.csv",
    )

    realized_gap = first_treated - first_untreated
    write_json(
        {
            "runtime_seconds": time.time() - started,
            "master_seed": master_seed,
            "clusters": clusters,
            "cluster_size": m,
            "mean_boundary_spacing": float(realized_gap.mean()),
            "analytic_mean_boundary_spacing": 1.0 / (m + 1.0),
            "minimum_absolute_recentered_score": float(
                np.min(np.abs(recentered))
            ),
        },
        "output/logs/mechanical_density_run.json",
    )
    print(
        f"mechanical-density: clusters={clusters}, m={m}, "
        f"mean gap={realized_gap.mean():.6f}"
    )


if __name__ == "__main__":
    main()
