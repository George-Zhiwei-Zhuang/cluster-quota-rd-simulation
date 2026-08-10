"""Automated integrity and substantive checks for the replication package."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "code" / "02_construct"))

from simulation_utils import load_config, software_versions  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    config = load_config()
    processed = PROJECT_ROOT / "data" / "processed"
    output = PROJECT_ROOT / "output"
    large = pd.read_csv(processed / "large_cluster_summary.csv")
    small = pd.read_csv(processed / "many_small_summary.csv")
    density = pd.read_csv(processed / "density_local_limits.csv")

    checks: dict[str, bool] = {}
    checks["large_all_cells_present"] = (
        large.shape[0]
        == len(config["large_cluster"]["cluster_sizes"])
        * len(config["large_cluster"]["regimes"])
    )
    checks["large_success_rate_above_95pct"] = bool(
        (large["success_rate"] >= 0.95).all()
    )
    checks["large_coverage_reasonable"] = bool(
        large["coverage_95"].between(0.85, 0.99).all()
    )

    largest_m = int(large["m"].max())
    largest = large.loc[large["m"].eq(largest_m)].set_index("regime")
    checks["coupling_correlation_orders_regimes"] = bool(
        largest.loc[
            "strong_overlap", "corr_feasible_reference"
        ]
        > largest.loc[
            "partial_overlap", "corr_feasible_reference"
        ]
        > largest.loc[
            "vanishing_overlap", "corr_feasible_reference"
        ]
    )
    checks["window_overlap_orders_regimes"] = bool(
        largest.loc["strong_overlap", "mean_window_overlap"]
        > largest.loc["partial_overlap", "mean_window_overlap"]
        > largest.loc["vanishing_overlap", "mean_window_overlap"]
    )

    positive_gap = small.loc[small["baseline_gap_B"] > 0.0]
    no_gap = small.loc[small["baseline_gap_B"].eq(0.0)]
    max_g = int(small["G"].max())
    checks["many_small_Delta_coverage_reasonable"] = bool(
        small["coverage_Delta_95"].between(0.90, 0.99).all()
    )
    checks["many_small_tau_coverage_when_B_zero"] = bool(
        no_gap["coverage_tau_95"].between(0.90, 0.99).all()
    )
    checks["many_small_tau_coverage_collapses_when_B_positive"] = bool(
        positive_gap.loc[
            positive_gap["G"].eq(max_g), "coverage_tau_95"
        ].max()
        < 0.20
    )

    smallest_t = density.sort_values("window_t").iloc[0]
    checks["fixed_density_near_one"] = bool(
        abs(smallest_t["fixed_right_density"] - 1.0) < 0.10
        and abs(smallest_t["fixed_left_density"] - 1.0) < 0.10
    )
    checks["recentered_density_near_two"] = bool(
        abs(smallest_t["recentered_right_density"] - 2.0) < 0.15
        and abs(smallest_t["recentered_left_density"] - 2.0) < 0.15
    )

    required_outputs = [
        output / "tables" / "table_1_large_cluster.csv",
        output / "tables" / "table_2_many_small.csv",
        output / "tables" / "table_3_mechanical_density.csv",
        output / "figures" / "figure_1_large_cluster_phase_transition.png",
        output / "figures" / "figure_1b_scaled_coupling_difference.png",
        output / "figures" / "figure_2_many_small_precision.png",
        output / "figures" / "figure_3_mechanical_density.png",
    ]
    checks["required_outputs_exist"] = all(
        path.exists() and path.stat().st_size > 0
        for path in required_outputs
    )

    report = {
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "software_versions": software_versions(),
    }
    report_path = output / "validation_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    checksum_targets = sorted(
        [
            path
            for folder in [
                processed,
                output / "tables",
                output / "figures",
            ]
            for path in folder.glob("*")
            if path.is_file()
        ]
    )
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}"
        for path in checksum_targets
    ]
    (output / "checksums.sha256").write_text(
        "\n".join(checksum_lines) + "\n",
        encoding="utf-8",
    )
    (PROJECT_ROOT / "environment" / "software_versions.json").write_text(
        json.dumps(software_versions(), indent=2, sort_keys=True),
        encoding="utf-8",
    )

    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["all_checks_pass"]:
        raise SystemExit("One or more validation checks failed.")


if __name__ == "__main__":
    main()
