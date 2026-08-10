# Deviations from blueprint

## 2026-07-27

- The main simulation package omits a separate mixed small-large pooling
  experiment. The pooled probability limits are algebraic, whereas the
  finite-sample questions most directly requiring simulation concern
  marginal coverage, pathwise coupling, and estimand mismatch.
- The mechanical-density experiment reports induced local density levels
  and histograms. It does not report rejection frequencies for a
  conventional McCrary statistic because such frequencies would depend
  on a particular implementation and variance calibration and are not
  required for Proposition 4.
- Large-cluster inference uses an HC3 conditional sandwich variance.
  The DGP gives both potential-outcome regressions the same curvature,
  so the leading jump bias cancels; the partial- and vanishing-overlap
  sequences additionally undersmooth. Robust bias correction is outside
  this minimal package and can be added as a supplementary experiment.
