# Data-source ledger

| Input | Status | Unit | Access | Provenance |
|---|---|---|---|---|
| Large-cluster Monte Carlo draws | Synthetic | Unit within simulated cluster | Generated locally | `01_large_cluster_phase_transition.py` |
| Many-small-cluster contrasts | Synthetic | Simulated cluster | Generated locally | `02_many_small_clusters.py` |
| Mechanical-density draws | Synthetic | Unit within simulated cluster | Generated locally | `03_mechanical_density.py` |

No external data, manual downloads, account access, or restricted inputs
are used. All random inputs are generated from the distributions and
seeds recorded in `config/simulation.yaml` and the JSON run logs.
