# Variable dictionary

## `large_cluster_replications.csv`

| Variable | Meaning |
|---|---|
| `regime` | Bandwidth-rate regime |
| `m` | Cluster size |
| `bandwidth` | Local-linear bandwidth \(h_m\) |
| `kappa_sqrt_m_h` | Phase index \(\sqrt m h_m\) |
| `local_rate_mh` | Expected local-sample rate \(m h_m\) |
| `empirical_cutoff` | Midpoint sample-quantile boundary |
| `cutoff_displacement_bw` | \((\widehat c-c)/h\), with \(c=0\) |
| `feasible_estimate` | RD estimate centered at \(\widehat c\) |
| `reference_estimate` | Coupled RD estimate centered at \(c=0\) |
| `scaled_difference` | \(\sqrt{mh}(\widehat\tau_F-\widehat\tau_R)\) |
| `feasible_se` | Side-specific HC3 standard error |
| `coverage_95` | Indicator that the feasible CI covers \(\tau=1\) |
| `window_overlap` | Jaccard overlap of the two local windows |
| `success` | Indicator for a nondegenerate fit on both sides |

## `many_small_replications.csv`

| Variable | Meaning |
|---|---|
| `baseline_slope` | \(\beta\) in the untreated score gradient |
| `G` | Number of independent small clusters |
| `tau` | Own-treatment effect |
| `baseline_gap_B` | \(B=\beta/(m+1)\) |
| `assignment_contrast_Delta` | \(\Delta=\tau+B\) |
| `estimate` | Equal-cluster average of boundary contrasts |
| `estimated_se` | Cluster-level standard error |
| `coverage_Delta_95` | CI coverage for \(\Delta\) |
| `coverage_tau_95` | The same CI's coverage for \(\tau\) |

## Density outputs

`density_histogram.csv` contains common-bin density estimates for
fixed-cutoff and same-sample-recentered scores.
`density_local_limits.csv` reports one-sided local density estimates for
several windows \(t\). `normalized_rank_grid.csv` records the exact
rank support and its within-cluster point mass.
