# Simulation design

## Experiment 1: marginal validity and coupling

Scores satisfy \(X_i\sim N(0,1)\). The quota treats half the cluster.
Potential outcomes are

\[
Y_i(0)=0.30X_i+0.15X_i^2+\varepsilon_{i0},
\qquad
Y_i(1)=Y_i(0)-\varepsilon_{i0}+1+\varepsilon_{i1},
\]

where \(\varepsilon_{i0},\varepsilon_{i1}\) are independent standard
normal variables. Hence the conditional mean treatment effect equals
one at every score.

For \(m\in\{1000,4000,16000\}\), the bandwidth is \(h_m=m^{-a}\), with

| Regime | \(a\) | \(\sqrt m h_m\) | \(mh_m\) |
|---|---:|---:|---:|
| Strong overlap | 0.20 | \(m^{0.30}\to\infty\) | \(m^{0.80}\to\infty\) |
| Partial overlap | 0.50 | \(1\) | \(m^{0.50}\to\infty\) |
| Vanishing overlap | 0.60 | \(1.5m^{-0.10}\to0\) | \(1.5m^{0.40}\to\infty\) |

The experiment uses 1,000 replications per cell. It reports feasible
coverage, local-window overlap, feasible-reference correlation, and the
variance of the scaled pathwise difference relative to its separated
limit \(2V(c)\).

## Experiment 2: many small clusters

Each cluster has \(m=20\) and \(k=10\). Uniform adjacent spacings satisfy

\[
S_k\sim\operatorname{Beta}(1,m),
\qquad E[S_k]=\frac{1}{m+1}.
\]

The simulated cluster contrast is

\[
A_g=\tau+\beta S_{k,g}+\eta_g,
\qquad
\eta_g\sim N(0,2),
\]

so

\[
E[A_g]=\Delta=\tau+\frac{\beta}{m+1}.
\]

The experiment uses \(\beta\in\{0,4.2\}\),
\(G\in\{50,200,800,3200\}\), and 5,000 replications. When
\(\beta=4.2\), \(B=0.2\) exactly.

## Experiment 3: diagnostic target

The experiment draws 100,000 clusters of 20 iid
\(\operatorname{Unif}(0,1)\) scores. It compares:

1. \(X_i-0.5\);
2. \(X_i-\widehat c\), with \(\widehat c\) the rank-10/rank-11 midpoint;
3. the deterministic normalized-rank grid.

For a fixed cutoff, the one-sided density limit equals one. For midpoint
recentring, the selected adjacent spacing has
\(q_S(0)=m\), giving the induced one-sided limit

\[
\frac{2q_S(0)}{m}=2.
\]

The simulation records common-bin histograms and direct local-window
estimates of these limits.
