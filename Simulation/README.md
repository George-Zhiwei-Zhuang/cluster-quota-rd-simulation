# Cluster-Quota RD Simulation Replication Package

This package reproduces the Monte Carlo evidence for *Regression
Discontinuity with Cluster-Level Quotas: Random Score Boundaries,
Finite-Rank Margins, and an Illustration from Formula One Qualifying*.

## Research questions

The simulations isolate three claims.

1. A feasible local-linear estimator centered at a same-sample
   order-statistic boundary can retain conventional marginal coverage
   even when it is not pathwise coupled to a population-boundary
   reference estimator.
2. With many bounded clusters, increasing the number of clusters
   estimates the finite-rank assignment contrast
   \(\Delta=\tau+B\) increasingly precisely but does not eliminate the
   adjacent-rank baseline gap \(B\).
3. Recentring iid uniform scores at a same-sample midpoint cutoff changes
   the local density target mechanically, while normalized ranks lie on
   a deterministic grid.

## Data

All inputs are generated synthetically from the documented
data-generating processes. No external or restricted data are required.
The raw-data directories are intentionally empty.

## Estimands and estimators

The large-cluster experiment uses iid standard-normal scores, a
fixed 50 percent quota, constant treatment effect \(\tau=1\), and a
smooth quadratic baseline regression. It computes:

- the feasible local-linear estimator at the realized midpoint cutoff;
- a coupled population-reference estimator at the population median
  using the same primitive scores and potential outcomes;
- an HC3 conditional sandwich standard error for the feasible estimator;
- local-window overlap and the scaled pathwise difference.

The many-small experiment uses \(m=20\), \(k=10\), iid uniform scores,
and the exact uniform-spacing result
\(E[S_k]=1/(m+1)\). With baseline slope \(\beta\),
\(B=\beta/(m+1)\) and \(\Delta=\tau+B\).

The density experiment uses 100,000 independent clusters of 20 iid
uniform scores and compares a fixed cutoff at 0.5 with the midpoint
between ranks 10 and 11.

## Reproduction

Python 3.12 or later is recommended. Install the pinned dependencies:

```bash
python3 -m pip install -r environment/requirements.txt
```

Run the entire package from the project root:

```bash
PROJECT_PYTHON=python3 ./run_all.sh
```

In the Codex primary runtime used to produce the archived outputs:

```bash
PROJECT_PYTHON="$CODEX_PRIMARY_RUNTIME_PYTHON" ./run_all.sh
```

The master seed is `20260727`. Component scripts use deterministic child
seeds recorded in their JSON run logs. A full run takes roughly one to
three minutes on a contemporary laptop, depending on CPU speed.

## Dependency map

```text
config/simulation.yaml
    -> code/03_analyze/01_large_cluster_phase_transition.py
    -> code/03_analyze/02_many_small_clusters.py
    -> code/03_analyze/03_mechanical_density.py
    -> data/processed/*.csv
    -> code/04_tables/01_build_tables.py
    -> code/05_figures/01_build_figures.py
    -> output/tables/*
    -> output/figures/*
    -> code/06_diagnostics/01_validate_outputs.py
    -> output/validation_report.json
    -> output/checksums.sha256
```

## Primary outputs

- `output/tables/table_1_large_cluster.*`
- `output/tables/table_2_many_small.*`
- `output/tables/table_3_mechanical_density.*`
- `output/figures/figure_1_large_cluster_phase_transition.*`
- `output/figures/figure_1b_scaled_coupling_difference.*`
- `output/figures/figure_2_many_small_precision.*`
- `output/figures/figure_3_mechanical_density.*`

Both CSV/LaTeX tables and PNG/PDF figures are generated. Processed
replication-level draws are stored under `data/processed/`.

## Validation

`code/06_diagnostics/01_validate_outputs.py` checks output existence,
Monte Carlo cell counts, estimator success, coverage, ordering of
coupling measures across regimes, the collapse of coverage for
\(\tau\) when \(B\ne0\), and the analytical density limits. It writes:

- `output/validation_report.json`;
- `output/checksums.sha256`;
- `environment/software_versions.json`.

## Scope and limitations

The simulations illustrate finite-sample implications of the paper's
theorems; they do not replace the proofs. The large-cluster DGP imposes
the conditional independence and smoothness conditions under which the
random-boundary expansion applies. The many-small experiment is designed
to isolate estimand mismatch, not interference. The density experiment
demonstrates a mechanical change in the local target; it does not claim
that every conventional density-discontinuity statistic must reject.
