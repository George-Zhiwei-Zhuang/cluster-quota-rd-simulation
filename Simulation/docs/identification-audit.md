# Simulation estimand and interpretation audit

## Large-cluster experiment

- **Primitive population:** iid \(X\sim N(0,1)\) with smooth conditional
  potential-outcome means and independent mean-zero residuals.
- **Assignment:** the highest \(k=m/2\) scores are treated; the feasible
  cutoff is the midpoint of order statistics \(k\) and \(k+1\).
- **Target:** \(\tau(c)=1\) at the population median \(c=0\).
- **Finite-sample object:** local-linear RD centered at the realized
  cutoff, compared pathwise with a coupled population-reference
  estimator using the same primitive draw.
- **What varies:** bandwidth exponent, and hence
  \(\kappa_m=\sqrt m h_m\), while \(h_m\to0\) and \(mh_m\to\infty\).
- **Diagnostic implication:** marginal coverage need not move in tandem
  with window overlap or estimator correlation.

## Many-small-cluster experiment

- **Primitive population:** iid uniform scores within clusters of
  \(m=20\).
- **Assignment:** the top \(k=10\) ranks are treated.
- **Own-treatment effect:** \(\tau=1\).
- **Observed boundary contrast:** \(\Delta=\tau+B\), where
  \(B=\beta E[S_k]=\beta/(m+1)\).
- **Weakest identifying restriction for interpreting \(\Delta\) as
  \(\tau\):** \(B=0\). Increasing \(G\) does not supply this restriction.

## Mechanical-density experiment

- **Primitive density:** \(f(x)=1\) on \([0,1]\).
- **Fixed-cutoff target:** the density of \(X-0.5\), equal to one near
  zero.
- **Recentered target:** the induced distribution of
  \(X-\widehat c\), where \(\widehat c\) is the same-sample midpoint.
- **Interpretation:** the experiment evaluates the induced diagnostic
  distribution, not manipulation or causal identification.

## Claims deliberately not made

- Simulation does not prove the asymptotic theorems.
- Good Monte Carlo coverage under this DGP does not establish robustness
  to irregular score densities, outcome-selected cutoffs, or
  interference.
- The mechanical-density exercise does not imply that every density
  discontinuity statistic must overreject.
