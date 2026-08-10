"""Build all simulation figures from processed CSV files."""

from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[2] / ".mplconfig"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROCESSED = PROJECT_ROOT / "data" / "processed"
FIGURES = PROJECT_ROOT / "output" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

COLORS = {
    "strong_overlap": "#1b9e77",
    "partial_overlap": "#d95f02",
    "vanishing_overlap": "#7570b3",
}
LABELS = {
    "strong_overlap": "Strong overlap",
    "partial_overlap": "Partial overlap",
    "vanishing_overlap": "Vanishing overlap",
}


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(
        FIGURES / f"{stem}.pdf",
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


def figure_large_cluster_summary() -> None:
    summary = pd.read_csv(PROCESSED / "large_cluster_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.3), sharex=True)
    panels = [
        ("coverage_95", "95% CI coverage", 0.95),
        ("mean_window_overlap", "Mean local-window overlap", None),
        (
            "corr_feasible_reference",
            "Correlation: feasible vs. reference",
            None,
        ),
        (
            "var_difference_ratio_to_2V",
            r"$\mathrm{Var}\{\sqrt{mh}(\hat\tau_F-\hat\tau_R)\}/(2V)$",
            1.0,
        ),
    ]

    for axis, (column, ylabel, reference) in zip(axes.ravel(), panels):
        for regime in LABELS:
            group = summary.loc[summary["regime"].eq(regime)].sort_values(
                "m"
            )
            axis.plot(
                group["m"],
                group[column],
                marker="o",
                linewidth=2,
                color=COLORS[regime],
                label=LABELS[regime],
            )
        if reference is not None:
            axis.axhline(
                reference,
                color="black",
                linestyle="--",
                linewidth=1,
                alpha=0.7,
            )
        axis.set_xscale("log")
        axis.set_xlabel("Cluster size $m$")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.2)

    axes[0, 0].set_ylim(0.80, 1.00)
    axes[0, 1].set_ylim(-0.02, 1.02)
    axes[1, 0].set_ylim(-0.15, 1.05)
    axes[1, 1].set_ylim(-0.05, 1.20)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.suptitle(
        "Marginal inference remains stable as pathwise coupling changes",
        y=1.055,
        fontsize=13,
    )
    fig.tight_layout()
    save_figure(fig, "figure_1_large_cluster_phase_transition")


def figure_scaled_differences() -> None:
    replications = pd.read_csv(
        PROCESSED / "large_cluster_replications.csv"
    )
    max_m = int(replications["m"].max())
    selected = replications.loc[
        replications["m"].eq(max_m)
        & replications["success"].eq(1)
    ]
    theoretical_v = float(
        pd.read_csv(PROCESSED / "large_cluster_summary.csv")[
            "theoretical_V"
        ].iloc[0]
    )
    x_grid = np.linspace(
        -4.5 * math.sqrt(2.0 * theoretical_v),
        4.5 * math.sqrt(2.0 * theoretical_v),
        500,
    )
    normal_density = (
        np.exp(-(x_grid**2) / (4.0 * theoretical_v))
        / math.sqrt(4.0 * math.pi * theoretical_v)
    )

    fig, axes = plt.subplots(1, 3, figsize=(11.3, 3.4), sharey=True)
    for axis, regime in zip(axes, LABELS):
        values = selected.loc[
            selected["regime"].eq(regime), "scaled_difference"
        ].dropna()
        axis.hist(
            values,
            bins=35,
            density=True,
            color=COLORS[regime],
            alpha=0.75,
            edgecolor="white",
        )
        if regime == "vanishing_overlap":
            axis.plot(
                x_grid,
                normal_density,
                color="black",
                linestyle="--",
                linewidth=1.5,
                label=r"$N(0,2V)$",
            )
            axis.legend(frameon=False)
        axis.axvline(0.0, color="black", linewidth=1)
        axis.set_title(LABELS[regime])
        axis.set_xlabel(
            r"$\sqrt{mh}(\hat\tau_F-\hat\tau_R)$"
        )
        axis.grid(alpha=0.15)
    axes[0].set_ylabel("Monte Carlo density")
    fig.suptitle(
        f"Coupling discrepancy at the largest simulated cluster size "
        f"($m={max_m:,}$)",
        y=1.02,
        fontsize=12,
    )
    fig.tight_layout()
    save_figure(fig, "figure_1b_scaled_coupling_difference")


def figure_many_small() -> None:
    summary = pd.read_csv(PROCESSED / "many_small_summary.csv")
    slopes = sorted(summary["baseline_slope"].unique())
    fig, axes = plt.subplots(1, len(slopes), figsize=(10.5, 4.0), sharey=True)
    if len(slopes) == 1:
        axes = [axes]

    for axis, beta in zip(axes, slopes):
        group = summary.loc[
            summary["baseline_slope"].eq(beta)
        ].sort_values("G")
        axis.fill_between(
            group["G"],
            group["estimate_q025"],
            group["estimate_q975"],
            color="#80b1d3",
            alpha=0.35,
            label="Monte Carlo 95% range",
        )
        axis.plot(
            group["G"],
            group["mean_estimate"],
            color="#1f78b4",
            marker="o",
            linewidth=2,
            label=r"Mean $\hat\Delta$",
        )
        tau = group["tau"].iloc[0]
        delta = group["assignment_contrast_Delta"].iloc[0]
        axis.axhline(
            tau,
            color="#e31a1c",
            linestyle=":",
            linewidth=2,
            label=r"Own-treatment effect $\tau$",
        )
        axis.axhline(
            delta,
            color="black",
            linestyle="--",
            linewidth=1.5,
            label=r"Assignment contrast $\Delta$",
        )
        axis.set_xscale("log")
        axis.set_xlabel("Number of clusters $G$")
        axis.set_title(
            rf"Baseline slope $\beta={beta:.1f}$; "
            rf"$B={group['baseline_gap_B'].iloc[0]:.2f}$"
        )
        axis.grid(alpha=0.2)

    axes[0].set_ylabel("Estimate")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.06),
    )
    fig.suptitle(
        "More small clusters tighten inference around "
        r"$\Delta=\tau+B$, not around $\tau$",
        y=1.16,
        fontsize=13,
    )
    fig.tight_layout()
    save_figure(fig, "figure_2_many_small_precision")


def figure_mechanical_density() -> None:
    histogram = pd.read_csv(PROCESSED / "density_histogram.csv")
    rank_grid = pd.read_csv(PROCESSED / "normalized_rank_grid.csv")
    local = histogram["bin_midpoint"].between(-0.18, 0.18)

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7))
    axes[0].plot(
        histogram.loc[local, "bin_midpoint"],
        histogram.loc[local, "fixed_cutoff_density"],
        color="#1b9e77",
        linewidth=2,
    )
    axes[0].axhline(
        1.0, color="black", linestyle="--", linewidth=1.2
    )
    axes[0].set_title("Original scores, fixed cutoff")
    axes[0].set_xlabel(r"$X-c$")
    axes[0].set_ylabel("Density")

    axes[1].plot(
        histogram.loc[local, "bin_midpoint"],
        histogram.loc[local, "recentered_density"],
        color="#d95f02",
        linewidth=2,
    )
    axes[1].axhline(
        2.0, color="black", linestyle="--", linewidth=1.2
    )
    axes[1].set_title("Scores, same-sample midpoint")
    axes[1].set_xlabel(r"$X-\hat c$")
    axes[1].set_ylabel("Density")

    markerline, stemlines, baseline = axes[2].stem(
        rank_grid["normalized_rank"],
        rank_grid["within_cluster_mass"],
        linefmt="-",
        markerfmt="o",
        basefmt=" ",
    )
    plt.setp(markerline, color="#7570b3", markersize=4)
    plt.setp(stemlines, color="#7570b3", linewidth=1.5)
    axes[2].set_title("Normalized ranks")
    axes[2].set_xlabel(r"$R/m$")
    axes[2].set_ylabel("Mass at each grid point")
    axes[2].set_ylim(0.0, 0.065)

    for axis in axes:
        axis.grid(alpha=0.18)
    fig.suptitle(
        "Same-sample quota recentering changes the diagnostic target",
        y=1.03,
        fontsize=13,
    )
    fig.tight_layout()
    save_figure(fig, "figure_3_mechanical_density")


def main() -> None:
    figure_large_cluster_summary()
    figure_scaled_differences()
    figure_many_small()
    figure_mechanical_density()


if __name__ == "__main__":
    main()
