#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_PYTHON="${PROJECT_PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-${PROJECT_ROOT}/.mplconfig}"
mkdir -p "${MPLCONFIGDIR}"

"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/03_analyze/01_large_cluster_phase_transition.py"
"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/03_analyze/02_many_small_clusters.py"
"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/03_analyze/03_mechanical_density.py"
"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/04_tables/01_build_tables.py"
"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/05_figures/01_build_figures.py"
"${PROJECT_PYTHON}" "${PROJECT_ROOT}/code/06_diagnostics/01_validate_outputs.py"
