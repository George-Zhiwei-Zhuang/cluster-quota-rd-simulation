"""Monte Carlo for marginal validity versus pathwise coupling."""

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
    boundary_outcome,
    ensure_output_dirs,
    load_config,
    local_linear_rd,
    local_window_overlap,
    midpoint_quota_cutoff,
    triangular_local_linear_variance,
    write_dataframe,
    write_json,
)


def main() -> None:
    started = time.time()
    ensure_output_dirs()
    config = load_config()
    settings = config["large_cluster"]
    master_seed = int(config["project"]["master_seed"])
    replications = int(settings["replications"])
    tau = float(settings["treatment_effect"])
    error_sd = float(settings["error_sd"])
    alpha = float(settings["treatment_share"])
    theoretical_v = triangular_local_linear_variance(error_sd, error_sd)

    seed_sequence = np.random.SeedSequence(master_seed)
    combinations = [
        (regime, int(m))
        for regime in settings["regimes"]
        for m in settings["cluster_sizes"]
    ]
    child_seeds = seed_sequence.spawn(len(combinations))
    records: list[dict[str, float | int | str]] = []

    for (regime, m), child_seed in zip(combinations, child_seeds):
        regime_config = settings["regimes"][regime]
        exponent = float(regime_config["bandwidth_exponent"])
        constant = float(regime_config["bandwidth_constant"])
        bandwidth = constant * m ** (-exponent)
        treated_count = int(round(alpha * m))
        rng = np.random.default_rng(child_seed)

        for replication in range(1, replications + 1):
            x = rng.normal(size=m)
            epsilon0 = rng.normal(scale=error_sd, size=m)
            epsilon1 = rng.normal(scale=error_sd, size=m)
            baseline = 0.30 * x + 0.15 * x**2
            y0 = baseline + epsilon0
            y1 = baseline + tau + epsilon1

            empirical_cutoff = midpoint_quota_cutoff(x, treated_count)
            population_cutoff = 0.0
            feasible_y = boundary_outcome(
                x, y0, y1, empirical_cutoff
            )
            reference_y = boundary_outcome(
                x, y0, y1, population_cutoff
            )

            feasible = local_linear_rd(
                x, feasible_y, empirical_cutoff, bandwidth
            )
            reference = local_linear_rd(
                x, reference_y, population_cutoff, bandwidth
            )

            success = int(
                np.isfinite(feasible["estimate"])
                and np.isfinite(feasible["se"])
                and feasible["se"] > 0.0
                and np.isfinite(reference["estimate"])
            )
            if success:
                scaled_difference = math.sqrt(m * bandwidth) * (
                    feasible["estimate"] - reference["estimate"]
                )
                z_feasible = (
                    feasible["estimate"] - tau
                ) / feasible["se"]
                coverage = int(abs(z_feasible) <= 1.959963984540054)
            else:
                scaled_difference = math.nan
                z_feasible = math.nan
                coverage = 0

            records.append(
                {
                    "regime": regime,
                    "m": m,
                    "replication": replication,
                    "bandwidth": bandwidth,
                    "kappa_sqrt_m_h": math.sqrt(m) * bandwidth,
                    "local_rate_mh": m * bandwidth,
                    "empirical_cutoff": empirical_cutoff,
                    "cutoff_displacement_bw": (
                        empirical_cutoff / bandwidth
                    ),
                    "feasible_estimate": feasible["estimate"],
                    "feasible_se": feasible["se"],
                    "reference_estimate": reference["estimate"],
                    "scaled_difference": scaled_difference,
                    "z_feasible": z_feasible,
                    "coverage_95": coverage,
                    "window_overlap": local_window_overlap(
                        x,
                        empirical_cutoff,
                        population_cutoff,
                        bandwidth,
                    ),
                    "n_feasible_right": feasible["n_right"],
                    "n_feasible_left": feasible["n_left"],
                    "success": success,
                }
            )

        print(
            f"large-cluster: regime={regime}, m={m}, "
            f"h={bandwidth:.6f}, kappa={math.sqrt(m)*bandwidth:.3f}"
        )

    replications_df = pd.DataFrame.from_records(records)
    write_dataframe(
        replications_df,
        "data/processed/large_cluster_replications.csv",
    )

    successful = replications_df.loc[
        replications_df["success"].eq(1)
    ].copy()
    summary_rows: list[dict[str, float | int | str]] = []
    for (regime, m), group in successful.groupby(
        ["regime", "m"], sort=False
    ):
        feasible = group["feasible_estimate"]
        reference = group["reference_estimate"]
        scaled_difference = group["scaled_difference"]
        summary_rows.append(
            {
                "regime": regime,
                "m": int(m),
                "replications_requested": replications,
                "replications_successful": int(group.shape[0]),
                "success_rate": group.shape[0] / replications,
                "bandwidth": group["bandwidth"].iloc[0],
                "kappa_sqrt_m_h": group["kappa_sqrt_m_h"].iloc[0],
                "local_rate_mh": group["local_rate_mh"].iloc[0],
                "feasible_bias": feasible.mean() - tau,
                "feasible_rmse": math.sqrt(
                    np.mean((feasible - tau) ** 2)
                ),
                "feasible_mc_sd": feasible.std(ddof=1),
                "mean_estimated_se": group["feasible_se"].mean(),
                "coverage_95": group["coverage_95"].mean(),
                "z_mean": group["z_feasible"].mean(),
                "z_sd": group["z_feasible"].std(ddof=1),
                "reference_bias": reference.mean() - tau,
                "reference_mc_sd": reference.std(ddof=1),
                "corr_feasible_reference": feasible.corr(reference),
                "mean_window_overlap": group["window_overlap"].mean(),
                "var_scaled_difference": scaled_difference.var(ddof=1),
                "var_difference_ratio_to_2V": (
                    scaled_difference.var(ddof=1)
                    / (2.0 * theoretical_v)
                ),
                "theoretical_V": theoretical_v,
            }
        )

    summary = pd.DataFrame.from_records(summary_rows)
    summary = summary.sort_values(["regime", "m"])
    write_dataframe(
        summary,
        "data/processed/large_cluster_summary.csv",
    )
    write_json(
        {
            "runtime_seconds": time.time() - started,
            "master_seed": master_seed,
            "theoretical_V": theoretical_v,
            "replications": replications,
        },
        "output/logs/large_cluster_run.json",
    )


if __name__ == "__main__":
    main()
