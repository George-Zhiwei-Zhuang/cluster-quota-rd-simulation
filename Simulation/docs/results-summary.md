# Simulation results

## Large-cluster random boundary

The feasible estimator has negligible Monte Carlo bias in every cell.
At the largest simulated cluster size, its 95 percent coverage is 0.948,
0.949, and 0.946 in the strong-, partial-, and vanishing-overlap
sequences, respectively.

Marginal stability coexists with sharply different pathwise behavior.
At \(m=16{,}000\):

- strong overlap has mean window overlap 0.949, feasible-reference
  correlation 0.811, and scaled-difference variance equal to 0.197 of
  the separated-window benchmark \(2V\);
- partial overlap has mean window overlap 0.414, correlation 0.081, and
  variance ratio 0.933;
- vanishing overlap has mean window overlap 0.273, correlation 0.046,
  and variance ratio 1.077.

The finite-sample patterns therefore reproduce the paper's distinction:
coverage of the feasible estimator approaches its nominal level in all
three sequences, while coupling measures move toward different limits.
The strong-overlap convergence of the pathwise discrepancy is visibly
slower than the convergence of marginal coverage.

## Many small clusters

When the baseline slope is zero, \(B=0\), the estimator is centered at
\(\tau=\Delta=1\) and coverage remains approximately 0.95 as \(G\)
increases.

When the baseline slope is 4.2,

\[
B=\frac{4.2}{21}=0.2,\qquad \Delta=1.2.
\]

The mean estimate is 1.204, 1.201, 1.200, and 1.200 for
\(G=50,200,800,3200\). Coverage for \(\Delta\) remains between 0.946
and 0.949. Coverage for the own-treatment effect \(\tau=1\), using the
same increasingly narrow confidence interval, falls from 0.819 to
0.484, 0.024, and finally 0.000. Thus increasing \(G\) removes sampling
uncertainty around the wrong causal target rather than removing the
baseline adjacent-rank gap.

## Mechanical density

With iid uniform primitive scores, the fixed-cutoff one-sided density
estimates in the smallest window \(t=0.005\) are 1.000 and 1.010,
close to the primitive density one. The corresponding midpoint-
recentered estimates are 1.940 and 1.932, close to the analytical limit
two. The simulation therefore displays a symmetric mechanical pileup,
not evidence of manipulation and not necessarily a density jump.

## Recommended main-text use

- Use `figure_1_large_cluster_phase_transition` as the principal
  simulation figure.
- Use `figure_2_many_small_precision` as the estimand-mismatch figure.
- Use `figure_3_mechanical_density` either at the end of the simulation
  section or alongside Section 9.
- Place `figure_1b_scaled_coupling_difference` and the full tables in the
  simulation appendix if the main text is constrained.
