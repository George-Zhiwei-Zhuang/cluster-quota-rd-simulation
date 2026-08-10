"""Create publication-ready CSV and LaTeX tables from simulation results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT_ROOT / "data" / "processed"
TABLES = PROJECT_ROOT / "output" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)


def save_table(
    frame: pd.DataFrame,
    stem: str,
    caption: str,
    label: str,
) -> None:
    frame.to_csv(TABLES / f"{stem}.csv", index=False)
    alignment = "l" + "r" * (frame.shape[1] - 1)
    header = " & ".join(str(column) for column in frame.columns)
    rows = []
    for row in frame.itertuples(index=False, name=None):
        formatted = []
        for value in row:
            if isinstance(value, float):
                formatted.append(f"{value:.3f}")
            else:
                formatted.append(str(value))
        rows.append(" & ".join(formatted) + r" \\")
    latex = "\n".join(
        [
            r"\begin{table}[!htbp]",
            r"\centering",
            rf"\caption{{{caption}}}",
            rf"\label{{{label}}}",
            rf"\begin{{tabular}}{{{alignment}}}",
            r"\toprule",
            header + r" \\",
            r"\midrule",
            *rows,
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    (TABLES / f"{stem}.tex").write_text(latex, encoding="utf-8")


def main() -> None:
    large = pd.read_csv(PROCESSED / "large_cluster_summary.csv")
    regime_labels = {
        "strong_overlap": "Strong overlap",
        "partial_overlap": "Partial overlap",
        "vanishing_overlap": "Vanishing overlap",
    }
    table_large = large[
        [
            "regime",
            "m",
            "kappa_sqrt_m_h",
            "feasible_bias",
            "feasible_rmse",
            "coverage_95",
            "corr_feasible_reference",
            "mean_window_overlap",
            "var_difference_ratio_to_2V",
        ]
    ].copy()
    table_large["regime"] = table_large["regime"].map(regime_labels)
    table_large = table_large.rename(
        columns={
            "regime": "Regime",
            "m": "m",
            "kappa_sqrt_m_h": "kappa",
            "feasible_bias": "Bias(F)",
            "feasible_rmse": "RMSE(F)",
            "coverage_95": "Coverage(F)",
            "corr_feasible_reference": "Corr(F,R)",
            "mean_window_overlap": "Window overlap",
            "var_difference_ratio_to_2V": "Var(diff)/(2V)",
        }
    )
    save_table(
        table_large,
        "table_1_large_cluster",
        (
            "Marginal inference and pathwise coupling across "
            "bandwidth regimes."
        ),
        "tab:large-cluster-simulation",
    )

    small = pd.read_csv(PROCESSED / "many_small_summary.csv")
    table_small = small[
        [
            "baseline_slope",
            "G",
            "baseline_gap_B",
            "assignment_contrast_Delta",
            "mean_estimate",
            "bias_relative_Delta",
            "mean_estimated_se",
            "coverage_Delta_95",
            "coverage_tau_95",
        ]
    ].copy()
    table_small = table_small.rename(
        columns={
            "baseline_slope": "Baseline slope",
            "G": "G",
            "baseline_gap_B": "B",
            "assignment_contrast_Delta": "Delta",
            "mean_estimate": "Mean estimate",
            "bias_relative_Delta": "Bias wrt Delta",
            "mean_estimated_se": "Mean SE",
            "coverage_Delta_95": "Coverage(Delta)",
            "coverage_tau_95": "Coverage(tau)",
        }
    )
    save_table(
        table_small,
        "table_2_many_small",
        (
            "Increasing the number of small clusters improves "
            "precision for the assignment contrast but does not "
            "eliminate the adjacent-rank baseline gap."
        ),
        "tab:many-small-simulation",
    )

    density = pd.read_csv(PROCESSED / "density_local_limits.csv")
    table_density = density.rename(
        columns={
            "window_t": "t",
            "fixed_right_density": "Fixed right",
            "fixed_left_density": "Fixed left",
            "recentered_right_density": "Recentered right",
            "recentered_left_density": "Recentered left",
            "primitive_theoretical_density": "Fixed limit",
            "recentered_theoretical_limit": "Recentered limit",
        }
    )
    save_table(
        table_density,
        "table_3_mechanical_density",
        (
            "Local density estimates under a fixed cutoff and a "
            "same-sample midpoint cutoff for iid uniform scores."
        ),
        "tab:mechanical-density-simulation",
    )


if __name__ == "__main__":
    main()
