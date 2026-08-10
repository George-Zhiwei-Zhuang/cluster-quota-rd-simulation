"""Monte Carlo showing that precision does not eliminate the rank gap."""

from __future__ import annotations

import math
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
    settings = config["many_small"]
    master_seed = int(config["project"]["master_seed"]) + 10000
    replications = int(settings["replications"])
    m = int(settings["cluster_size"])
    k = int(settings["treated_count"])
    tau = float(settings["treatment_effect"])
    error_sd = float(settings["error_sd"])

    combinations = [
        (float(beta), int(g_count))
        for beta in settings["baseline_slopes"]
        for g_count in settings["cluster_counts"]
    ]
    child_seeds = np.random.SeedSequence(master_seed).spawn(
        len(combinations)
    )
    records: list[dict[str, float | int]] = []

    for (beta, g_count), child_seed in zip(combinations, child_seeds):
        rng = np.random.default_rng(child_seed)
        baseline_gap = beta / (m + 1.0)
        delta = tau + baseline_gap

        for replication in range(1, replications + 1):
            adjacent_spacing = rng.beta(1.0, float(m), size=g_count)
            residual_difference = rng.normal(
                scale=math.sqrt(2.0) * error_sd,
                size=g_count,
            )
            cluster_contrast = (
                tau
                + beta * adjacent_spacing
                + residual_difference
            )
            estimate = float(cluster_contrast.mean())
            se = float(cluster_contrast.std(ddof=1) / math.sqrt(g_count))
            lower = estimate - 1.959963984540054 * se
            upper = estimate + 1.959963984540054 * se

            records.append(
                {
                    "baseline_slope": beta,
                    "m": m,
                    "k": k,
                    "G": g_count,
                    "replication": replication,
                    "tau": tau,
                    "baseline_gap_B": baseline_gap,
                    "assignment_contrast_Delta": delta,
                    "estimate": estimate,
                    "estimated_se": se,
                    "coverage_Delta_95": int(lower <= delta <= upper),
                    "coverage_tau_95": int(lower <= tau <= upper),
                }
            )

        print(
            f"many-small: beta={beta:.2f}, G={g_count}, "
            f"tau={tau:.3f}, B={baseline_gap:.3f}, "
            f"Delta={delta:.3f}"
        )

    replications_df = pd.DataFrame.from_records(records)
    write_dataframe(
        replications_df,
        "data/processed/many_small_replications.csv",
    )

    summary_rows: list[dict[str, float | int]] = []
    for (beta, g_count), group in replications_df.groupby(
        ["baseline_slope", "G"], sort=False
    ):
        delta = group["assignment_contrast_Delta"].iloc[0]
        estimate = group["estimate"]
        summary_rows.append(
            {
                "baseline_slope": beta,
                "G": int(g_count),
                "m": m,
                "tau": tau,
                "baseline_gap_B": group["baseline_gap_B"].iloc[0],
                "assignment_contrast_Delta": delta,
                "mean_estimate": estimate.mean(),
                "bias_relative_Delta": estimate.mean() - delta,
                "rmse_relative_Delta": math.sqrt(
                    np.mean((estimate - delta) ** 2)
                ),
                "mc_sd": estimate.std(ddof=1),
                "mean_estimated_se": group["estimated_se"].mean(),
                "coverage_Delta_95": group[
                    "coverage_Delta_95"
                ].mean(),
                "coverage_tau_95": group["coverage_tau_95"].mean(),
                "estimate_q025": estimate.quantile(0.025),
                "estimate_q975": estimate.quantile(0.975),
            }
        )

    summary = pd.DataFrame.from_records(summary_rows)
    summary = summary.sort_values(["baseline_slope", "G"])
    write_dataframe(
        summary,
        "data/processed/many_small_summary.csv",
    )
    write_json(
        {
            "runtime_seconds": time.time() - started,
            "master_seed": master_seed,
            "replications": replications,
            "analytic_spacing_mean": 1.0 / (m + 1.0),
        },
        "output/logs/many_small_run.json",
    )


if __name__ == "__main__":
    main()
